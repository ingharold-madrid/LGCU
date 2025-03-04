#' @title ARL Calculation for a Downward CUSUM Control Chart with a Gamma Distribution
#'
#' @description This function estimates the average run length (**ARL**) of a downward CUSUM control chart applied to a Gamma distribution with guaranteed performance (GIC), considering both known and estimated parameters.
#'
#' This approach follows the methodology described in the work of Madrid‐Alvarez, García-díaz, and Tercero‐Gómez(2024),
#' which provides a detailed analysis of the performance of CUSUM control charts for Gamma distributions with guaranteed efficiency.
#' Specifically, the method implemented in this function enables the precise evaluation of ARL under different parameter settings,
#' ensuring appropriate calibration and monitoring of controlled processes.
#'
#' ### **Recommendations**
#'
#' For further consultation and to review values of `H_delta` and `H_minus`, it is recommended to refer to the following article:
#' Madrid‐Alvarez, H. M., García‐Díaz, J. C., & Tercero‐Gómez, V. G. (2024).
#' **A CUSUM control chart for gamma distribution with guaranteed performance**. Quality and Reliability Engineering.
#'
#' ## Considerations:
#' - The control chart is calibrated with `beta = 1`.
#' - When `alpha` and `beta` are known, it is recommended to use the same values for `alpha_est = alpha` and `beta_est = 1`.
#' - Higher values of `m` increase the accuracy of the results (`m = {10, 100}`).
#'
#' @param alpha Shape parameter of the Gamma distribution.
#' @param beta Scale parameter of the Gamma distribution.
#' @param alpha_est Estimated shape parameter.
#' @param beta_est Estimated scale parameter.
#' @param beta_ratio Ratio between `beta` and its estimation.
#' @param H_minus Lower control limit of the downward CUSUM chart.
#' @param H_delta Increment of the GIC threshold.
#' @param m Number of divisions for the probability matrix.
#'
#' @return A numeric value representing the average run length (**ARL**) of the downward CUSUM control chart.
#' @export
#' @importFrom stats pgamma

#'
#' @examples
#' # Example with known parameters
#' GICARL_CUSUM_down(alpha = 1, beta = 1, alpha_est = 1, beta_est = 1,
#'                   beta_ratio = 1/2.5, H_minus = -2.792, H_delta = 0, m = 100)
#'
#' # Example with estimated parameters
#' GICARL_CUSUM_down(alpha = 1, beta = 1, alpha_est = 1, beta_est = 1.1,
#'                   beta_ratio = 1/2, H_minus = -4.1497, H_delta = 1.5167,
#'                   m = 100)
#'
GICARL_CUSUM_down <- function(alpha, beta, alpha_est, beta_est, beta_ratio, H_minus, H_delta, m) {

  k_minus <- -((alpha_est * beta * beta_ratio * log(beta / beta_ratio)) / (beta_ratio - beta))
  H_minus_c <- H_minus - H_delta
  w <- (-2 * H_minus_c) / (2 * m - 1)

  P <- matrix(nrow = m + 1, ncol = m + 1)

  for (i in 1:m) {
    P[i, 1] <- pgamma(-w / 2 + (i - 1) * w + k_minus, shape = alpha, scale = beta / beta_est, lower.tail = FALSE)
  }

  for (i in 1:m) {
    for (j in 2:m) {
      P[i, j] <- pgamma(((i - 1) - (j - 1)) * w - (1 / 2) * w + k_minus, shape = alpha, scale = beta / beta_est, lower.tail = FALSE) -
        pgamma(((i - 1) - (j - 1)) * w + (1 / 2) * w + k_minus, shape = alpha, scale = beta / beta_est, lower.tail = FALSE)
    }
  }

  for (i in 1:m) {
    P[i, m + 1] <- pgamma(H_minus + (1/2) * w + k_minus, shape = alpha, scale = beta / beta_est, lower.tail = TRUE)
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
