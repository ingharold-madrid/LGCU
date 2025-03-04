#' @title ARL Calculation for an Upward CUSUM Control Chart with a Gamma Distribution
#'
#' @description This function estimates the average run length (**ARL**) of an upward CUSUM control chart applied to a Gamma distribution with guaranteed performance (GIC), considering both known and estimated parameters.
#'
#' This approach follows the methodology described in the work of Madrid‐Alvarez, García‐Díaz, and Tercero‐Gómez (2024),
#' which provides a detailed analysis of the performance of CUSUM control charts for Gamma distributions with guaranteed efficiency.
#' Specifically, the method implemented in this function enables the precise evaluation of ARL under different parameter settings,
#' ensuring appropriate calibration and monitoring of controlled processes.
#'
#' ### **Recommendations**
#'
#' For further consultation and to review values of `H_delta` and `H_plus`, it is recommended to refer to the following article:
#' Madrid‐Alvarez, H. M., García‐Díaz, J. C., & Tercero‐Gómez, V. G. (2024).
#' **A CUSUM control chart for gamma distribution with guaranteed performance**. Quality and Reliability Engineering.
#'
#' ## Key Considerations:
#' - The control chart is calibrated with `beta = 1`.
#' - When `alpha` and `beta` are known, it is recommended to use the same values for `alpha_est = alpha` and `beta_est = 1`.
#' - Higher values of `m` increase the accuracy of the results (`m = {10, 100}`).
#'
#' @param alpha Shape parameter of the Gamma distribution.
#' @param beta Scale parameter of the Gamma distribution.
#' @param alpha_est Estimated shape parameter.
#' @param beta_est Estimated scale parameter.
#' @param beta_ratio Ratio between `beta` and its estimation.
#' @param H_plus Upper control limit of the upward CUSUM chart.
#' @param H_delta Increment of the GIC threshold.
#' @param m Number of divisions for the probability matrix.
#'
#' @return A numeric value representing the average run length (**ARL**) of the upward CUSUM control chart.
#' @export
#' @importFrom stats pgamma
#'
#' @examples
#' # Example with known parameters
#' GICARL_CUSUM_up(alpha = 0.9, beta = 2.136, alpha_est = 0.9, beta_est = 1,
#'                 beta_ratio = 2.67, H_plus = 25.1592, H_delta = 0, m = 100)
#'
#' # Example with estimated parameters
#' GICARL_CUSUM_up(alpha = 1, beta = 1, alpha_est = 1.2, beta_est = 0.8,
#'                 beta_ratio = 2, H_plus = 6.5081, H_delta = 2.9693, m = 100)
#'
GICARL_CUSUM_up <- function(alpha, beta, alpha_est, beta_est, beta_ratio, H_plus, H_delta, m) {

  k_plus <- -((alpha_est * beta * beta_ratio * log(beta / beta_ratio)) / (beta_ratio - beta))
  H_plus_c <- H_plus + H_delta
  w <- (2 * H_plus_c) / (2 * m - 1)

  P <- matrix(nrow = m + 1, ncol = m + 1)

  for (i in 1:m) {
    P[i, 1] <- pgamma(w / 2 - (i - 1) * w + k_plus, shape = alpha, scale = beta / beta_est, lower.tail = TRUE)
  }

  for (i in 1:m) {
    for (j in 2:m) {
      P[i, j] <- pgamma(((j - 1) - (i - 1)) * w + (1 / 2) * w + k_plus, shape = alpha, scale = beta / beta_est, lower.tail = TRUE) -
        pgamma(((j - 1) - (i - 1)) * w - (1 / 2) * w + k_plus, shape = alpha, scale = beta / beta_est, lower.tail = TRUE)
    }
  }

  for (i in 1:m) {
    P[i, m + 1] <- pgamma(H_plus - (i - 1) * w + k_plus, shape = alpha, scale = beta / beta_est, lower.tail = FALSE)
  }

  for (j in 1:m) {
    P[m + 1, j] <- 0
  }

  P[m + 1, m + 1] <- 1

  # Constructing probability matrix
  R <- P[1:m, 1:m]
  I <- diag(1, ncol(R))
  Rest <- I - R
  Inv <- solve(Rest)

  ones <- matrix(rep(1, m), ncol = 1)
  d <- rep(0, m)
  d[1] <- 1
  d <- matrix(d, ncol = 1)

  ARL <- t(d) %*% Inv %*% ones

  return(as.numeric(ARL))
}
