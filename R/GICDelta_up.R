#' @title Estimation of the optimal `H_delta` value to guarantee performance in the upward CUSUM control chart
#'
#' @description This function calculates the optimal `H_delta` value that ensures specific performance in the Gamma CUSUM
#' control chart for upward detection. It relies on Monte Carlo simulations and an iterative adjustment process
#' to determine the appropriate value.
#'
#' Following the methodology proposed by Madrid-Alvarez, Garcia-Diaz, and Tercero-Gomez (2024),
#' this function allows adjusting `H_delta` for different sample size scenarios, ensuring that the control chart
#' maintains the expected performance in terms of ARL.
#'
#' ## Features:
#' - Implements Monte Carlo simulations to estimate `H_delta`.
#' - Based on parameter estimates obtained in Phase I.
#' - Iteratively adjusts `H_delta` until the specified ARL is reached.
#' - Displays total execution time using `tictoc`.
#'
#' ## **Recommendations**
#' - This function is useful for estimating `H_delta` values in scenarios where the sample size differs from the values reported in the reference article:
#'
#'   **Madrid-Alvarez, H. M., Garcia-Diaz, J. C., & Tercero-Gomez, V. G. (2024).**
#'   *A CUSUM control chart for gamma distribution with guaranteed performance.*
#'   Quality and Reliability Engineering International, 40(3), 1279-1301.
#'
#' - **The adjustment process is iterative and computationally demanding**, as its execution time depends on the number of iterations (`N_init + N_final`) and the sample size (`n_I`).
#' - It is recommended to establish an appropriate convergence criterion to optimize execution time without compromising the accuracy of `H_delta` estimation.
#' - For selecting values of `H_plus`, `a`, and `b`, refer to the reference article, which presents specific strategies and recommendations for calibration.
#'
#' @param n_I Sample size in Phase I.
#' @param alpha Shape parameter of the Gamma distribution.
#' @param beta Scale parameter of the Gamma distribution.
#' @param beta_ratio Ratio between beta and its estimate.
#' @param H_plus Initial upper limit of the CUSUM chart.
#' @param a Tolerance level for the expected ARL (0 <= a < 1).
#' @param b Tolerance level for the expected ARL (0 < b < 1).
#' @param ARL_esp Desired expected ARL value.
#' @param m Number of states in the Markov matrix.
#' @param N_init Number of initial iterations.
#' @param N_final Number of final iterations.
#' @param known_alpha Indicates whether `alpha` is known (`TRUE`) or needs to be estimated (`FALSE`).
#'
#' @return A numeric value corresponding to the optimal `H_delta` for the upward CUSUM control chart, ensuring the expected performance.
#' @export
#' @importFrom stats rgamma
#' @importFrom MASS fitdistr
#'
#' @examples
#' \donttest{
#' getDeltaH_up(n_I = 100, alpha = 1, beta = 1, beta_ratio = 2, H_plus = 6.8313,
#'              a = 0.1, b = 0.05, ARL_esp = 370, m = 100,
#'              N_init = 10, N_final = 1000, known_alpha = TRUE)
#'              }
#'
getDeltaH_up <- function(n_I, alpha, beta, beta_ratio, H_plus,
                         a, b, ARL_esp, m, N_init, N_final, known_alpha) {

  # Verify that the ARL function is available
  if (!exists("GICARL_CUSUM_up")) {
    stop("The function 'GICARL_CUSUM_up()' is not defined in the environment. Ensure it is included in the package.")
  }

  # Verify that the MASS library is installed
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("The 'MASS' library is required but not installed. Use install.packages('MASS') to install it.")
  }

  # Verify that the tictoc library is installed
  if (!requireNamespace("tictoc", quietly = TRUE)) {
    stop("The 'tictoc' library is required but not installed. Use install.packages('tictoc') to install it.")
  }

  # Start the timer
  tictoc::tic("Total execution time")

  # Initialize variables
  Dvector <- rep(NA, (N_init + N_final))
  D0 <- H_plus / 10
  gamma <- H_plus / 100
  D <- D0  # Starting point in the search for H_delta
  verbose = TRUE

  # Main iteration
  for (i in 1:(N_init + N_final)) {

    # Generate sample in Phase I
    X <- rgamma(n_I, shape = alpha, scale = beta)
    estimador <- MASS::fitdistr(x = X, 'gamma', method = "Nelder-Mead")

    # Estimate parameters
    if (known_alpha) {
      alpha_est <- alpha  # Known alpha
    } else {
      alpha_est <- as.numeric(estimador$estimate[1])
    }
    beta_est <- mean(X) / alpha_est  # Beta estimation

    # Calculate ARL using GICARL_CUSUM_up()
    ARL <- GICARL_CUSUM_up(alpha = alpha, beta = beta,
                           alpha_est = alpha_est, beta_est = beta_est,
                           beta_ratio = beta_ratio,
                           H_plus = H_plus, H_delta = D, m = m)

    # Adjust H_delta based on the calculated ARL
    y <- b
    if (ARL < ((1 - a) * ARL_esp)) {
      y <- b - 1
    }
    D <- max(0, (D - gamma * y))
    Dvector[i] <- D

    # Display progress every 50 iterations if verbose is TRUE
    if (verbose && i %% 50 == 0) {
      message(sprintf("Iteration: %d | H_delta: %.4f", i, D))
    }
  }

  # Compute the average H_delta over the last iterations
  H_delta <- mean(Dvector[(N_init + 1):(N_init + N_final)])

  # Stop the timer and display the total execution time
  tictoc::toc()

  return(H_delta)
}


