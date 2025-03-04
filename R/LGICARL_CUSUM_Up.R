#' @title ARL Estimation in CUSUM Control Charts with Gamma Distribution and Cautious Learning
#' for upward detection
#'
#' @description This function calculates the **Average Run Length (ARL)** of a CUSUM control chart based on the Gamma distribution,
#' incorporating a cautious learning scheme for the progressive update of parameters and optimization of performance in upward detection.
#'
#' The function allows for the evaluation of the CUSUM chart’s behavior under different parameterization scenarios,
#' ensuring efficient detection of process changes.
#'
#' Following the methodology presented in the work of Madrid-Alvarez, García-Díaz, and Tercero-Gómez (2024),
#' this implementation utilizes Monte Carlo simulations in C++ for efficient execution,
#' ensuring a dynamic adjustment of parameters based on the evolution of the process.
#'
#' The values of `H_plus`, `H_delta`, `K_l`, `delay`, and `tau` can be referenced in the tables from the article:
#'
#' **Madrid-Alvarez, H. M., García-Díaz, J. C., & Tercero-Gómez, V. G. (2024).**
#' *A CUSUM control chart for the Gamma distribution with cautious parameter learning.*
#' Quality Engineering, 1-23.
#'
#' ## Usage Scenarios:
#'
#' **Scenario 1: Known `alpha` and estimated `beta`**
#' - The `alpha` parameter is assumed to be fixed and known in advance.
#' - `beta` is estimated from a dataset or defined by the user.
#' - The user must provide `alpha` and an initial estimate of `beta` (`beta0_est`).
#'
#' **Scenario 2: Both `alpha` and `beta` are estimated**
#' - Both `alpha` and `beta` are estimated from a dataset or external data source.
#' - The user must calculate `alpha0_est` and `beta0_est` before running the function.
#' - `beta0_est` is dynamically updated during the simulation when a predefined condition is met.
#'
#' ## Features:
#' - Implements Monte Carlo simulations optimized in C++ for ARL estimation.
#' - Allows dynamic updating of `beta0_est` to improve the model's adaptability.
#' - Compatible with scenarios where `alpha` is known or estimated.
#' - Ensures stable and reliable performance in detecting changes in processes modeled with Gamma distributions.
#' - Recommended values for `H_plus`, `H_delta`, `K_l`, `delay`, and `tau` can be found in the reference article.
#'
#' @param alpha Shape parameter of the Gamma distribution.
#' @param beta Scale parameter of the Gamma distribution.
#' @param alpha0_est Initial estimate of the shape parameter `alpha`. If `known_alpha` is `TRUE`, this value will be equal to `alpha`.
#' @param beta0_est Initial estimate of the scale parameter `beta`. This value is updated dynamically during the simulation.
#' @param known_alpha `TRUE` if `alpha0_est` is fixed, `FALSE` if it must be estimated.
#' @param beta_ratio Ratio between `beta` and its posterior estimate.
#' @param H_delta Increment of the upper control limit in the CUSUM chart.
#' @param H_plus Initial control limit of the CUSUM chart.
#' @param n_I Sample size in Phase I.
#' @param replicates Number of Monte Carlo simulations.
#' @param K_l Secondary control threshold for parameter updating.
#' @param delay Number of observations before updating `beta0_est`.
#' @param tau Time point at which `beta` changes. A value of 1 is recommended for IC states.
#'
#' @return A numeric value corresponding to the **ARL** estimate for the upward CUSUM control chart with cautious learning.
#'
#' @examples
#' # Option 1: Provide parameters directly
#' ARL_Clplus(
#'   alpha = 1,
#'   beta = 1,
#'   alpha0_est = 1,  # alpha = known_alpha
#'   beta0_est = 1.1,   # Estimated Beta
#'   known_alpha = TRUE,
#'   beta_ratio = 2,
#'   H_delta = 4.2433,
#'   H_plus = 8.7434,
#'   n_I = 200,
#'   replicates = 100,
#'   K_l = 2,
#'   delay = 25,
#'   tau = 1
#' )
#'
#' # Option 2: Use generated data
#' set.seed(123)
#' datos_faseI <- rgamma(n = 200, shape = 1, scale = 1)
#' alpha0_est <- mean(datos_faseI)^2 / var(datos_faseI)  # Alpha estimation
#' beta0_est <- mean(datos_faseI) / alpha0_est  # Beta estimation
#'
#' ARL_Clplus(
#'   alpha = 1,
#'   beta = 1,
#'   alpha0_est = alpha0_est,
#'   beta0_est = beta0_est,
#'   known_alpha = FALSE,
#'   beta_ratio = 2,
#'   H_delta = 4.2433,
#'   H_plus = 8.7434,
#'   n_I = 200,
#'   replicates = 1000,
#'   K_l = 2,
#'   delay = 25,
#'   tau = 1
#' )
#' @export
#'
#' @useDynLib LGCU, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rgamma
#' @importFrom MASS fitdistr
#'
ARL_Clplus <- function(alpha, beta, alpha0_est, beta0_est, known_alpha, beta_ratio, H_delta, H_plus, n_I,
                       replicates, K_l, delay, tau) {

  # Ensure alpha0_est is equal to alpha when alpha is known
  if (known_alpha) {
    alpha0_est <- alpha
  }

  beta1 <- beta
  k_plus <- -((alpha0_est * beta * beta_ratio * log(beta / beta_ratio)) / (beta_ratio - beta))

  RL_results <- numeric(replicates)
  for (i in 1:replicates) {
    resultado <- monte_carlo_loop_cpp(alpha, beta, alpha0_est, beta0_est, beta1, n_I,
                                      known_alpha, H_delta, H_plus, k_plus, K_l, delay, tau)
    RL_results[i] <- resultado$RL
  }

  return(mean(RL_results))
}
