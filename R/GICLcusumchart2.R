#' @title CUSUM Control Chart with Cautious Learning and Guaranteed Performance
#'
#' @description This function generates a bidirectional (upward and downward) CUSUM control chart for a Gamma distribution,
#' incorporating a cautious parameter update mechanism with guaranteed performance. Its purpose is to enhance sensitivity
#' and precision in detecting changes in dynamic processes.
#'
#' Based on the methodology presented by Madrid-Alvarez, García-Díaz, and Tercero-Gómez (2024), this implementation allows
#' control limits to adapt according to the evolution of the process, ensuring early detection of variations while minimizing
#' the risk of false alarms.
#'
#' ## Features:
#' - If the user does not provide Phase I and Phase II data, the function automatically generates them.
#' - If `known_alpha = TRUE`, `alpha` is fixed and not estimated.
#' - If `known_alpha = FALSE`, `alpha` is estimated from Phase I data.
#' - Includes dynamic control limits and a summary table of parameters.
#' - Enables the detection of both upward and downward deviations, progressively adjusting the control limits.
#'
#' ## **Recommendations**
#' - The parameters `k_l`, `delay`, and `tau` are crucial for the learning process in the control chart.
#'   They regulate the progressive update of control limits, allowing the dynamic update of `beta0_est`, `H_plus_c`, and `H_minus_c`, ensuring that the control chart
#'   gradually adjusts to changes in the process. It is recommended to use reference values presented in:
#'
#'   **Madrid-Alvarez, H. M., García-Díaz, J. C., & Tercero-Gómez, V. G. (2024).**
#'   *A CUSUM control chart for the Gamma distribution with cautious parameter learning.*
#'   Quality Engineering, 1-23.
#'
#' - Similar to the parameters above, for proper selection of `H_plus`, `H_minus`, `H_delta_plus`, and `H_delta_minus` values,
#'   it is recommended to review the reference article, where detailed calibration strategies for different scenarios are presented.
#'
#' @param H_delta_plus Increment of the upper control limit.
#' @param H_delta_minus Increment of the lower control limit.
#' @param beta_ratio_plus Ratio between `beta` and its estimate for upward detection.
#' @param beta_ratio_minus Ratio between `beta` and its estimate for downward detection.
#' @param alpha Shape parameter of the Gamma distribution (if `known_alpha = TRUE`).
#' @param beta Scale parameter of the Gamma distribution.
#' @param H_plus Initial upper limit of the CUSUM chart.
#' @param H_minus Initial lower limit of the CUSUM chart.
#' @param known_alpha Indicates whether `alpha` is known (`TRUE`) or should be estimated (`FALSE`).
#' @param k_l Secondary control threshold used in the learning logic.
#' @param delay Number of observations before updating `beta0_est`, `H_plus_c`, and `H_minus_c`.
#' @param tau Time point at which the `beta` parameter changes.
#' @param n_I Sample size in Phase I (if `faseI` is not provided).
#' @param n_II Sample size in Phase II (if `faseII` is not provided).
#' @param faseI Data sample from Phase I (numeric vector). If `NULL`, it is generated internally.
#' @param faseII Data sample from Phase II (numeric vector). If `NULL`, it is generated internally.
#'
#' @return A plot showing the evolution of the CUSUM statistic with cautious learning, including:
#' - Dynamically adjusted accumulated values of the CUSUM statistic.
#' - Progressively updated control limits with guaranteed performance.
#' - A summary of the parameters used in the control chart.
#' @export
#' @importFrom stats rgamma
#' @importFrom graphics lines
#' @importFrom graphics abline layout legend par rect text
#' @importFrom utils install.packages
#' @importFrom MASS fitdistr
#'
#' @examples
#' # Option 1: Automatically generated data
#' plot_GICCL_chart2(alpha = 1, beta = 1,
#'                  beta_ratio_plus = 2, beta_ratio_minus = 0.5,
#'                  H_delta_plus = 3.0, H_plus = 6.5,
#'                  H_delta_minus = 2.0, H_minus = -5.0,
#'                  known_alpha = TRUE, k_l = 2, delay = 25, tau = 1,
#'                  n_I = 200, n_II = 700,
#'                  faseI = NULL, faseII = NULL)
#'
#' # Option 2: User-provided data
#' datos_faseI <- rgamma(n = 200, shape = 1, scale = 1)
#' datos_faseII <- rgamma(n = 700, shape = 1, scale = 1)
#' plot_GICCL_chart2(alpha = 1, beta = 1,
#'                  beta_ratio_plus = 2, beta_ratio_minus = 0.5,
#'                  H_delta_plus = 3.0, H_plus = 6.5,
#'                  H_delta_minus = 2.0, H_minus = -5.0,
#'                  known_alpha = FALSE, k_l = 2, delay = 25, tau = 1,
#'                  n_I = 200, n_II = 700,
#'                  faseI = datos_faseI, faseII = datos_faseII)
#'

plot_GICCL_chart2 <- function(alpha, beta, beta_ratio_plus, beta_ratio_minus, H_delta_plus, H_plus,
                              H_delta_minus, H_minus, known_alpha, k_l, delay, tau,
                              n_I, n_II, faseI = NULL, faseII = NULL) {

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("The 'MASS' package is required but not installed. Use install.packages('MASS') to install it.")
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if (is.null(faseI)) {
    faseI <- rgamma(n = n_I, shape = alpha, scale = beta)
  } else {
    n_I <- length(faseI)
  }

  if (is.null(faseII)) {
    faseII <- rgamma(n = n_II, shape = alpha, scale = beta)
  } else {
    n_II <- length(faseII)
  }

  estimador <- MASS::fitdistr(x = faseI, 'gamma', method = "Nelder-Mead")
  alpha0_est <- ifelse(known_alpha, alpha, as.numeric(estimador$estimate[1]))
  beta0_est <- mean(faseI) / alpha0_est

  k_plus <- -((alpha0_est * beta * beta_ratio_plus * log(beta / beta_ratio_plus)) / (beta_ratio_plus - beta))
  k_minus <- -((alpha0_est * beta * beta_ratio_minus * log(beta / beta_ratio_minus)) / (beta_ratio_minus - beta))

  H_plus_c <- H_plus + H_delta_plus
  H_minus_c <- H_minus - H_delta_minus

  Cplus <- numeric(n_II)
  Cminus <- numeric(n_II)
  H_plus_c_values <- numeric(n_II)
  H_minus_c_values <- numeric(n_II)

  Cplus[1] <- max(0, (faseII[1] / beta0_est) - k_plus)
  Cminus[1] <- min(0, (faseII[1] / beta0_est) - k_minus)
  H_plus_c_values[1] <- H_plus_c
  H_minus_c_values[1] <- H_minus_c

  for (i in 2:n_II) {
    Cplus[i] <- max(0, Cplus[i - 1] + (faseII[i] / beta0_est) - k_plus)
    Cminus[i] <- min(0, Cminus[i - 1] + (faseII[i] / beta0_est) - k_minus)

    if (i > delay && i %% delay == 0) {
      beta0_est <- mean(faseII[1:i]) / alpha0_est
      H_plus_c <- H_plus + H_delta_plus * (sqrt(n_I / (n_I + i)))^1.7
      H_minus_c <- H_minus - H_delta_minus * (sqrt(n_I / (n_I + i)))^1.7
    }
    H_plus_c_values[i] <- H_plus_c
    H_minus_c_values[i] <- H_minus_c
  }

  H_plus_c_final <- H_plus_c
  H_minus_c_final <- H_minus_c

  par(mfrow = c(2,1))

  # Combined chart
  par(mar = c(4, 4, 2, 1))
  ylim_range <- range(c(Cplus, Cminus, H_plus_c_values, H_minus_c_values))

  plot(Cplus, type = "l", col = "blue", ylim = ylim_range,
       main = "CUSUM Control Chart with Cautious Learning",
       xlab = expression(bold("Observations (Phase II)")), ylab = expression(bold("CUSUM Statistic")),
       cex.main = 1.2, cex.lab = 0.9, cex.axis = 1.2,font.lab = 2, las = 1)

  lines(Cminus, col = "darkgreen")
  lines(H_plus_c_values, col = "red", lty = 2)
  lines(H_minus_c_values, col = "red", lty = 2)

  legend("topright", legend = c("C+", "C-", "Lower and Upper Limit"),
         col = c("blue", "darkgreen", "red", "red"),
         lwd = c(2, 2, 2, 2), lty = c(1, 1, 1, 1), cex = 0.9, bg = "white")

  # Summary box with gray background and border
  par(mar = c(1, 2, 1, 2))
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  rect(0.75, 0.6, 1.3, 1.1, col = "lightgray", border = "black", lwd = 1)  # Adjust rectangle size

  text(1, 1.05, "Control Chart Summary", cex = 1, font = 2, adj = 0.5)  # Adjusted title position
  text(0.76, 0.95, sprintf("Phase I Sample Size: %d", n_I), cex = 1, adj = 0)
  text(0.76, 0.85, sprintf("Estimated Alpha: %.2f", alpha0_est), cex = 1, adj = 0)
  text(0.76, 0.75, sprintf("Estimated Beta: %.2f", beta0_est), cex = 1, adj = 0)

  text(1.01, 0.95, sprintf("Updated C.L. (Final): %.2f", H_plus_c_final), cex = 1, adj = 0)
  text(1.01, 0.85, sprintf("Updated C.L. (Final): %.2f", H_minus_c_final), cex = 1, adj = 0)
  text(1.01, 0.75, sprintf("k_plus Value: %.4f", k_plus), cex = 1, adj = 0)
  text(1.01, 0.65, sprintf("k_minus Value: %.4f", k_minus), cex = 1, adj = 0)

  message("Execution completed successfully.")
}
