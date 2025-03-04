#include <Rcpp.h>
using namespace Rcpp;

// Declare that the function is defined in another file
List estimate_alpha_beta_cpp(NumericVector X, bool known_alpha, double alpha);

// [[Rcpp::export]]
List monte_carlo_loop_down_cpp(double alpha, double beta, double alpha0_est, double beta0_est,
                               double beta1, int n_I, bool known_alpha, double H_delta,
                               double H_minus, double k_minus, double K_l, int delay, int tau, int max_iter = 1e4) {

  int RL = 0;
  int previous_RL = 0;
  bool IC = true;
  double Cminus = 0.0;
  double Cminus_l = 0.0;
  int lastupdate = 0;
  NumericVector muestra_nl(max_iter);
  List results_list(max_iter);

  // Initial Phase I sample generation
  NumericVector X = rgamma(n_I, alpha, beta);
  List params = estimate_alpha_beta_cpp(X, known_alpha, alpha);
  beta0_est = params["beta0.est"];

  double H_delta_n = H_delta;
  double H_minus_c = H_minus - H_delta_n; // Detection threshold for downward deviations

  while (IC && RL < max_iter) {
    double current_beta = (RL >= tau) ? beta1 : beta;
    double x = R::rgamma(alpha, current_beta);

    Cminus = std::min(0.0, Cminus + (x / beta0_est) - k_minus);
    Cminus_l = std::min(0.0, Cminus_l + (x / beta0_est) - K_l);

    muestra_nl[RL] = x;
    lastupdate++;

    // Dynamic update of beta0_est
    if ((Cminus_l == 0) && (lastupdate >= delay)) {
      lastupdate = 0;
      NumericVector combined(RL + X.size());
      std::copy(muestra_nl.begin(), muestra_nl.begin() + RL, combined.begin());
      std::copy(X.begin(), X.end(), combined.begin() + RL);
      double X1 = mean(combined);
      beta0_est = X1 / alpha0_est;
      H_delta_n = H_delta * pow(sqrt(n_I / static_cast<double>(n_I + RL)), 1.7);
      H_minus_c = H_minus - H_delta_n;
    }

    results_list[RL] = List::create(_["iteration"] = RL, _["Cminus"] = Cminus,
                                    _["Cminus_l"] = Cminus_l, _["H_minus_c"] = H_minus_c,
                                    _["beta0.est"] = beta0_est);

    RL++;

    if (Cminus < H_minus_c) {
      IC = false;
    }
  }

  // Adjust RL considering the tau condition
  if (RL < tau) {
    RL = previous_RL;
  } else {
    previous_RL = RL;
    RL = RL - tau + 1;
  }

  return List::create(_["RL"] = RL, _["results_df"] = results_list);
}
