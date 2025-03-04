#' @title ARL Estimation in CUSUM Control Charts with Gamma Distribution and Cautious Learning
#' for downward detection
#'
#' @description This function calculates the **Average Run Length (ARL)** of a CUSUM control chart based on the Gamma distribution,
#' incorporating a cautious learning scheme for the dynamic update of parameters.
#'
#' The function allows the evaluation of the CUSUM chart’s performance under different parameterization scenarios,
#' ensuring efficient detection of process changes.
#'
#' Based on the methodology presented in the work of Madrid-Alvarez, García-Díaz, and Tercero-Gómez (2024),
#' this implementation uses Monte Carlo simulations optimized in C++ for efficient execution
#' and progressive adjustment of the control chart parameters.
#'
#' The values for `H_minus`, `H_delta`, `K_l`, `delay`, and `tau` can be referenced in the tables from the article:
#'
#' **Madrid-Alvarez, H. M., García-Díaz, J. C., & Tercero-Gómez, V. G. (2024).**
#' *A CUSUM control chart for the Gamma distribution with cautious parameter learning.*
#' Quality Engineering, 1-23.
#'
#' ## Usage Scenarios:
#'
#' **Scenario 1: Known `alpha` and estimated `beta`**
#' - The `alpha` parameter is assumed to be fixed and known in advance.
#' - `beta` is estimated from a dataset or provided by the user.
#' - The user must specify `alpha` and an initial estimate of `beta` (`beta0_est`).
#'
#' **Scenario 2: Both `alpha` and `beta` are estimated**
#' - Both `alpha` and `beta` are estimated from an external dataset.
#' - The user must calculate `alpha0_est` and `beta0_est` before calling the function.
#' - `beta0_est` is dynamically updated during the simulation when a predefined condition is met.
#'
#' ## Features:
#' - Implements Monte Carlo simulations for ARL estimation.
#' - Allows dynamic updating of `beta0_est` to improve model adaptation.
#' - Uses C++ optimization for efficient and precise execution.
#' - Compatible with scenarios where `alpha` is either known or estimated.
#' - Recommended values for `H_minus`, `H_delta`, `K_l`, `delay`, and `tau` can be found in the reference article.
#'
#' This function is ideal for quality control studies where reliable detection of process changes
#' modeled with Gamma distributions is required.
#'
#' @param alpha Shape parameter of the Gamma distribution.
#' @param beta Scale parameter of the Gamma distribution.
#' @param alpha0_est Initial estimate of the shape parameter `alpha`. If `known_alpha` is `TRUE`, this value will be equal to `alpha`.
#' @param beta0_est Initial estimate of the scale parameter `beta`. This value is updated dynamically during the simulation.
#' @param known_alpha `TRUE` if `alpha0_est` is fixed, `FALSE` if it must be estimated.
#' @param beta_ratio Ratio between `beta` and its posterior estimate.
#' @param H_delta Increment of the lower control limit in the CUSUM chart.
#' @param H_minus Initial control limit of the CUSUM chart for downward detection.
#' @param n_I Sample size in Phase I.
#' @param replicates Number of Monte Carlo simulations.
#' @param K_l Secondary control threshold for parameter updating.
#' @param delay Number of observations before updating `beta0_est`.
#' @param tau Time point at which `beta` changes.
#'
#' @return A numeric value corresponding to the **ARL** estimate for the downward CUSUM control chart with cautious learning.
#'
#' @examples
#' # Option 1: Provide parameters directly
#' ARL_Clminus(
#'    alpha = 1,
#'    beta = 1,
#'    alpha0_est = 1.067,  # alpha = known_alpha
#'    beta0_est = 0.2760,   # Estimated Beta
#'    known_alpha = TRUE,
#'    beta_ratio = 1/2,
#'    H_delta = 0.6946,
#'    H_minus = -4.8272,
#'    n_I = 500,
#'    replicates = 1000,
#'    K_l = 0.5,
#'    delay = 25,
#'    tau = 1
#' )
#'
#' # Option 2: Use generated data
#' set.seed(123)
#' datos_faseI <- rgamma(n = 500, shape = 1, scale = 1)
#' alpha0_est <- mean(datos_faseI)^2 / var(datos_faseI)  # Alpha estimation
#' beta0_est <- mean(datos_faseI) / alpha0_est  # Beta estimation
#'
#' ARL_Clminus(
#'    alpha = 1,
#'    beta = 1,
#'    alpha0_est = 1.067,  # alpha = known_alpha
#'    beta0_est = 0.2760,   # Estimated Beta
#'    known_alpha = FALSE,
#'    beta_ratio = 1/2,
#'    H_delta = 0.6946,
#'    H_minus = -4.8272,
#'    n_I = 500,
#'    replicates = 1000,
#'    K_l = 0.5,
#'    delay = 25,
#'    tau = 1
#' )
#' @export
#'
#' @useDynLib LGCU, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rgamma
#' @importFrom MASS fitdistr
#'
ARL_Clminus <- function(alpha, beta, alpha0_est, beta0_est, known_alpha, beta_ratio, H_delta, H_minus, n_I,
                        replicates, K_l, delay, tau) {

  # Ensure alpha0_est is equal to alpha when alpha is known
  if (known_alpha) {
    alpha0_est <- alpha
  }

  beta1 <- beta
  k_minus <- -((alpha0_est * beta * beta_ratio * log(beta / beta_ratio)) / (beta_ratio - beta))

  RL_results <- numeric(replicates)
  for (i in 1:replicates) {
    resultado <- monte_carlo_loop_down_cpp(alpha, beta, alpha0_est, beta0_est, beta1, n_I,
                                           known_alpha, H_delta, H_minus, k_minus, K_l, delay, tau)
    RL_results[i] <- resultado$RL
  }

  return(mean(RL_results))
}

