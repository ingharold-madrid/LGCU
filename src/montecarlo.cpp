#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List estimate_alpha_beta_cpp(NumericVector X, bool known_alpha, double alpha) {

  double m = mean(X);
  double alpha0_est, beta0_est;

  if (known_alpha) {
    alpha0_est = alpha;
  } else {
    double v = var(X);
    alpha0_est = pow(m, 2) / v;
  }
  beta0_est = m / alpha0_est;
  return List::create(_["alpha0.est"] = alpha0_est, _["beta0.est"] = beta0_est);
}

// [[Rcpp::export]]
List monte_carlo_loop_cpp(double alpha, double beta, double alpha0_est, double beta0_est,
                          double beta1, int n_I, bool known_alpha, double H_delta,
                          double H_plus, double k_plus, double K_l, int delay, int tau, int max_iter = 1e4) {

  int RL = 0;
  int previous_RL = 0;
  bool IC = true;
  double Cplus = 0.0;
  double Cplus_l = 0.0;
  int lastupdate = 0;
  NumericVector muestra_nl(max_iter);
  List results_list(max_iter);

  // Initial Phase I sample generation
  NumericVector X = rgamma(n_I, alpha, beta);
  List params = estimate_alpha_beta_cpp(X, known_alpha, alpha);
  //alpha0_est = params["alpha0.est"];
  beta0_est = params["beta0.est"];

  double H_delta_n = H_delta;
  double H_plus_c = H_plus + H_delta_n;

  while (IC && RL < max_iter) {
    double current_beta = (RL >= tau) ? beta1 : beta;
    double x = R::rgamma(alpha, current_beta);

    Cplus = std::max(0.0, Cplus + (x / beta0_est) - k_plus);
    Cplus_l = std::max(0.0, Cplus_l + (x / beta0_est) - K_l);
    muestra_nl[RL] = x;
    lastupdate++;

    // Dynamic update of beta0_est
    if ((Cplus_l == 0) && (lastupdate >= delay)) {
      lastupdate = 0;
      NumericVector combined(RL + X.size());
      std::copy(muestra_nl.begin(), muestra_nl.begin() + RL, combined.begin());
      std::copy(X.begin(), X.end(), combined.begin() + RL);
      double X1 = mean(combined);
      beta0_est = X1 / alpha0_est;
      H_delta_n = H_delta * pow(sqrt(n_I / static_cast<double>(n_I + RL)), 1.7);
      H_plus_c = H_plus + H_delta_n;
    }

    results_list[RL] = List::create(_["iteration"] = RL, _["Cplus"] = Cplus,
                                    _["Cplus_l"] = Cplus_l, _["H_plus_c"] = H_plus_c,
                                    _["beta0.est"] = beta0_est);

    RL++;

    if (Cplus > H_plus_c) {
      IC = false;
    }
  }

  // Adjust RL considering the tau condition
  if (RL < tau) {
    RL = previous_RL; // Use the previous RL value if RL is less than tau
  } else {
    previous_RL = RL; // Update previous_RL for the next iteration
    RL = RL - tau + 1; // Adjust RL to consider only after the change
  }

  return List::create(_["RL"] = RL, _["results_df"] = results_list);
}

