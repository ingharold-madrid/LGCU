#' @title Estimation of the Optimal `H_delta` Value to Guarantee Performance in the Downward CUSUM Control Chart
#'
#' @description This function calculates the optimal value of `H_delta` that guarantees a specific performance in the Gamma CUSUM control chart
#' for downward detection. It employs a Monte Carlo simulation approach and an iterative adjustment process to determine the appropriate value.
#'
#' Following the methodology presented by Madrid‐Alvarez, García‐Díaz, and Tercero‐Gómez (2024),
#' this function allows adjusting `H_delta` for different sample size configurations, ensuring that the control chart
#' maintains the desired performance in terms of expected ARL.
#'
#' ## Features:
#' - Implements Monte Carlo simulations to estimate `H_delta`.
#' - Based on parameter estimates obtained in Phase I.
#' - Iteratively adjusts `H_delta` until the specified ARL is achieved.
#' - Displays the total execution time using `tictoc`.
#'
#' ## **Recommendations**
#' - This function is useful for estimating `H_delta` values in scenarios where the sample size differs from
#'   the values reported in the reference paper:
#'
#'   **Madrid‐Alvarez, H. M., García‐Díaz, J. C., & Tercero‐Gómez, V. G. (2024).**
#'   *A CUSUM control chart for gamma distribution with guaranteed performance.*
#'   Quality and Reliability Engineering International, 40(3), 1279-1301.
#'
#' - **The adjustment process is iterative and computationally demanding**, as its execution time depends on the number of iterations (`N_init + N_final`) and the sample size (`n_I`).
#' - It is recommended to establish an appropriate convergence criterion to optimize execution time without compromising the accuracy of `H_delta` estimation.
#' - For selecting values of `H_minus`, `a`, and `b`, it is advisable to consult the reference paper, which provides specific calibration strategies and recommendations.
#'
#' @param n_I Sample size in Phase I.
#' @param alpha Shape parameter of the Gamma distribution.
#' @param beta Scale parameter of the Gamma distribution.
#' @param beta_ratio Ratio between beta and its estimate.
#' @param H_minus Initial lower limit of the CUSUM chart.
#' @param a Tolerance level for the expected ARL (0 <= a < 1).
#' @param b Tolerance level for the expected ARL (0 < b < 1).
#' @param ARL_esp Desired expected ARL value.
#' @param m Number of states in the Markov matrix.
#' @param N_init Number of initial iterations.
#' @param N_final Number of final iterations.
#' @param known_alpha Indicates whether `alpha` is known (`TRUE`) or should be estimated (`FALSE`).
#'
#' @return A numerical value corresponding to the optimal `H_delta` for the downward CUSUM control chart, ensuring the expected performance.
#' @export
#' @importFrom stats rgamma
#' @importFrom MASS fitdistr
#'
#' @examples
#' \donttest{
#' getDeltaH_down(n_I = 100, alpha = 1, beta = 1, beta_ratio = 1/2,
#'                H_minus = -4.1497, a = 0.1, b = 0.05, ARL_esp = 370,
#'                m = 100, N_init = 10, N_final = 1000, known_alpha = TRUE)
#'                }
#'
getDeltaH_down <- function(n_I, alpha, beta, beta_ratio, H_minus,
                           a, b, ARL_esp, m, N_init, N_final, known_alpha) {

  # Verify that the ARL function is available
  if (!exists("GICARL_CUSUM_down")) {
    stop("The function 'GICARL_CUSUM_down()' is not defined in the environment. Ensure it is included in the package.")
  }

  # Check that the MASS library is installed
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("The 'MASS' package is required but not installed. Use install.packages('MASS') to install it.")
  }

  # Check that the tictoc library is installed
  if (!requireNamespace("tictoc", quietly = TRUE)) {
    stop("The 'tictoc' package is required but not installed. Use install.packages('tictoc') to install it.")
  }

  # Start the timer
  tictoc::tic("Total execution time")

  # Initialize variables
  Dvector <- rep(NA, (N_init + N_final))
  D0 <- abs(H_minus) / 10  # Use abs() to avoid negative initial values
  gamma <- abs(H_minus) / 100
  D <- D0  # Starting point in the H_delta search
  verbose =TRUE

  # Main iteration
  for (i in 1:(N_init + N_final)) {

    # Generate sample in Phase I
    X <- rgamma(n_I, shape = alpha, scale = beta)
    estimador <- MASS::fitdistr(x = X, 'gamma', method = "Nelder-Mead")  # Use fitdistr()

    # Estimate parameters
    if (known_alpha) {
      alpha_est <- alpha  # Known alpha
    } else {
      alpha_est <- as.numeric(estimador$estimate[1])
    }
    beta_est <- mean(X) / alpha_est  # Beta estimation

    # Compute ARL using GICARL_CUSUM_down()
    ARL <- GICARL_CUSUM_down(alpha = alpha, beta = beta,
                             alpha_est = alpha_est, beta_est = beta_est,
                             beta_ratio = beta_ratio,
                             H_minus = H_minus, H_delta = D, m = m)

    # Adjust H_delta based on the computed ARL
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

  # Compute the average of H_delta over the last iterations
  H_delta <- mean(Dvector[(N_init + 1):(N_init + N_final)])

  # Stop the timer and display the total execution time
  tictoc::toc()

  return(H_delta)
}

