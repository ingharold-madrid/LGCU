#' @title Estimation of the `H_delta` parameter with learning for upward detection in Gamma CUSUM control charts
#'
#' @description This function calculates the optimal value of `H_delta` using a dynamic learning scheme
#' based on the `ARL_Clplus` function, iteratively adjusting `H_delta` to achieve an **expected ARL**
#' with higher accuracy and adaptability.
#'
#' Based on the methodology proposed by Madrid-Alvarez, García-Díaz, and Tercero-Gómez (2024),
#' this function allows adjusting `H_delta` in different sample size scenarios, ensuring
#' that the control chart progressively adapts to changes in the Gamma distribution.
#'
#' ## Features:
#' - Implements Monte Carlo simulations to estimate `H_delta`.
#' - Relies on parameter estimates obtained in Phase I.
#' - Iteratively adjusts `H_delta` until the specified ARL is reached.
#' - Incorporates a cautious learning mechanism to improve adjustment precision.
#' - Displays the total execution time using `tictoc`.
#'
#' ## **Recommendations**
#' - This function is useful for estimating `H_delta` values when the sample size differs from the values reported in the reference article:
#'
#'   **Madrid-Alvarez, H. M., García-Díaz, J. C., & Tercero-Gómez, V. G. (2024).**
#'   *A CUSUM control chart for the Gamma distribution with cautious parameter learning.*
#'   Quality Engineering, 1-23.
#'
#' - **The adjustment process is iterative and computationally intensive**, as execution time depends on the number of iterations (`N_init + N_final`)
#'   and the sample size (`n_I`).
#' - It is recommended to define a proper convergence criterion to optimize execution time without compromising the accuracy of `H_delta` estimation.
#' - For selecting values of `a`, `b`, `k_l`, `delay`, `tau`, and `H_plus`, consulting the reference article is recommended, as it provides specific strategies
#'   for their calibration in different scenarios.
#'
#' @param n_I Sample size in Phase I.
#' @param alpha Shape parameter of the Gamma distribution.
#' @param beta Scale parameter of the Gamma distribution.
#' @param beta_ratio Ratio between `beta` and its posterior estimate.
#' @param H_plus Initial limit of the CUSUM chart.
#' @param a Tolerance level for the expected ARL. (0 <= a < 1).
#' @param b Tolerance level for the expected ARL. (0 < b < 1)
#' @param ARL_esp Desired expected ARL value.
#' @param replicates Number of replications in the Monte Carlo simulation.
#' @param N_init Number of initial iterations for adjustment.
#' @param N_final Number of final iterations for averaging `H_delta`.
#' @param known_alpha `TRUE` if `alpha` is fixed, `FALSE` if it should be estimated.
#' @param K_l Secondary control threshold for parameter updating.
#' @param delay Number of observations before updating `beta0_est`.
#' @param tau Point in time where `beta` changes.
#'
#' @return A numeric value corresponding to the optimal `H_delta` estimated with learning for the upward CUSUM control chart.
#'
#' @export
#' @importFrom stats rgamma
#' @importFrom MASS fitdistr
#'
#' @examples
#' \donttest{
#' getDeltaHL_up(
#'            n_I = 200, alpha = 1, beta = 1, beta_ratio = 2,
#'             H_plus = 6.8313, a = 0.1, b = 0.05, ARL_esp = 370,
#'              replicates = 100, N_init = 100, N_final = 500,
#'              known_alpha = TRUE, K_l = 2, delay = 25, tau = 1
#'             )
#'             }
#'
#'
getDeltaHL_up <- function(
                          n_I, alpha, beta, beta_ratio, H_plus,
                          a, b, ARL_esp, replicates, N_init, N_final, known_alpha,
                          K_l, delay, tau
                          ) {

  # Verify that the ARL function is available
  if (!exists("ARL_Clplus")) {
    stop("The function 'ARL_Clplus()' is not defined in the environment. Ensure it is included in the package.")
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


  # Main iteration loop
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

    # Calculate ARL using the ARL_Clplus function
    ARL <- ARL_Clplus(alpha = alpha, beta = beta,
                      alpha0_est = alpha_est, beta0_est = beta_est,
                      beta_ratio = beta_ratio,
                      H_plus = H_plus, H_delta = D, n_I = n_I,
                      replicates = replicates, K_l = K_l, delay = delay, tau = tau,
                      known_alpha = known_alpha)

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

  # Calculate the average H_delta over the last iterations
  H_delta <- mean(Dvector[(N_init + 1):(N_init + N_final)])

  # Stop the timer and display the total execution time
  tictoc::toc()

  return(H_delta)
}
