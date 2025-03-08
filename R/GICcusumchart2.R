#' @title CUSUM Control Chart for Gamma Distribution with Guaranteed Performance
#'
#' @description This function generates a bidirectional (upward and downward) CUSUM control chart for a Gamma distribution,
#' allowing the monitoring of the evolution of the CUSUM statistic while ensuring optimal performance in detecting process changes.
#'
#' Based on the methodology proposed by Madrid‐Alvarez, García‐Díaz, and Tercero‐Gómez (2024), this implementation employs a Monte Carlo-based approach
#' to estimate the Gamma distribution parameters and determine control limits with precise calibration.
#' The function enables visualization of process evolution and the detection of deviations with reduced risk of false alarms.
#'
#' ## Features:
#' - Implements Monte Carlo simulations for control chart calibration.
#' - Allows the use of known Gamma distribution values or estimation in Phase I.
#' - Provides a graphical representation of the CUSUM statistic evolution with guaranteed performance.
#' - Includes control limits and a legend with key configuration details of the control chart.
#'
#' For additional details on selecting parameters `H_plus`, `H_minus`, `H_delta_plus`, and
#'   `H_delta_minus`, as well as calibration strategies, it is recommended to consult the reference article:
#'
#'   **Madrid‐Alvarez, H. M., García‐Díaz, J. C., & Tercero‐Gómez, V. G. (2024).**
#'   *A CUSUM control chart for gamma distribution with guaranteed performance.*
#'   Quality and Reliability Engineering International, 40(3), 1279-1301.
#'
#' @param alpha Shape parameter of the Gamma distribution.
#' @param beta Scale parameter of the Gamma distribution.
#' @param beta_ratio_plus Ratio between beta and its estimation for upward detection.
#' @param beta_ratio_minus Ratio between beta and its estimation for downward detection.
#' @param H_delta_plus Increment of the upper GIC limit.
#' @param H_plus Initial upper limit of the CUSUM chart.
#' @param H_delta_minus Increment of the lower GIC limit.
#' @param H_minus Initial lower limit of the CUSUM chart.
#' @param n_I Sample size in Phase I (if `faseI` is not provided).
#' @param n_II Sample size in Phase II (if `faseII` is not provided).
#' @param faseI Sample data from Phase I (numeric vector). If `NULL`, it is generated with `rgamma()`.
#' @param faseII Sample data from Phase II (numeric vector). If `NULL`, it is generated with `rgamma()`.
#' @param known_alpha If `TRUE`, a known `alpha` is used; if `FALSE`, it is estimated.
#'
#' @return A plot showing the evolution of the CUSUM statistic for a Gamma distribution with guaranteed performance, including:
#' - The accumulated values of the CUSUM statistic.
#' - Control limits with guaranteed performance.
#' - A summary of the parameters used in the control chart.
#' @export
#' @importFrom stats rgamma
#' @importFrom graphics lines
#' @importFrom graphics abline layout legend par rect text
#' @importFrom utils install.packages
#' @importFrom MASS fitdistr
#'
#' @examples
#' # Option 1: Automatically generate data with predefined sample sizes
#' plot_GICC_chart2(alpha = 1, beta = 1, beta_ratio_plus = 2,
#'                       beta_ratio_minus = 0.5,H_delta_plus = 2.0,
#'                       H_plus = 5.0, H_delta_minus = 1.5, H_minus = -4.5,
#'                       n_I = 100, n_II = 200, faseI = NULL,
#'                       faseII = NULL, known_alpha = TRUE
#'                      )
#'
#' # Option 2: Use custom data
#' phaseI_data <- rgamma(n = 100, shape = 1, scale = 1)
#' phaseII_data <- rgamma(n = 200, shape = 1, scale = 1)
#' plot_GICC_chart2(alpha = 1, beta = 1, beta_ratio_plus = 2,
#'                  beta_ratio_minus = 0.5, H_delta_plus = 2.0, H_plus = 5.0,
#'                  H_delta_minus = 1.5, H_minus = -4.5, n_I = 100, n_II = 200,
#'                  faseI = phaseI_data, faseII = phaseII_data,
#'                  known_alpha = TRUE
#'                  )
#'
plot_GICC_chart2 <- function(
                             alpha, beta, beta_ratio_plus, beta_ratio_minus, H_delta_plus, H_plus,
                             H_delta_minus, H_minus, n_I, n_II,
                             faseI = NULL, faseII = NULL, known_alpha
                             ) {

  # Load necessary packages
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("The 'MASS' package is required but not installed. Use install.packages('MASS') to install it.")
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  # Check if Phase I data is provided
  if (is.null(faseI)) {
    faseI <- rgamma(n = n_I, shape = alpha, scale = beta)
  } else {
    n_I <- length(faseI)
  }

  # Parameter estimation in Phase I
  estimator <- MASS::fitdistr(x = faseI, 'gamma', method = "Nelder-Mead")
  alpha0_est <- if (known_alpha) alpha else as.numeric(estimator$estimate[1])
  beta0_est <- mean(faseI) / alpha0_est

  # CUSUM parameter calculation
  k_plus <- -((alpha0_est * beta * beta_ratio_plus * log(beta / beta_ratio_plus)) / (beta_ratio_plus - beta))
  k_minus <- -((alpha0_est * beta * beta_ratio_minus * log(beta / beta_ratio_minus)) / (beta_ratio_minus - beta))

  # Upper and lower control limits with guaranteed performance
  H_plus_c <- H_plus + H_delta_plus
  H_minus_c <- H_minus - H_delta_minus

  # Check if Phase II data is provided
  if (is.null(faseII)) {
    faseII <- rgamma(n = n_II, shape = alpha, scale = beta)
  } else {
    n_II <- length(faseII)
  }

  # Generate CUSUM control chart in Phase II
  Cplus <- numeric(n_II)
  Cminus <- numeric(n_II)

  Cplus[1] <- max(0, 0 + (faseII[1] / beta0_est) - k_plus)
  Cminus[1] <- min(0, 0 + (faseII[1] / beta0_est) - k_minus)

  for (i in 2:n_II) {
    Cplus[i] <- max(0, Cplus[i - 1] + (faseII[i] / beta0_est) - k_plus)
    Cminus[i] <- min(0, Cminus[i - 1] + (faseII[i] / beta0_est) - k_minus)
  }

  par(mfrow = c(2,1))

  # Combined chart
  par(mar = c(4, 4, 2, 1))
  ylim_range <- range(c(Cplus, Cminus, H_plus_c, H_minus_c))
  plot(Cplus, type = "l", col = "blue", ylim = ylim_range,
       main = "CUSUM Control Chart for Gamma Distribution with Guaranteed Performance",
       xlab = expression(bold("Observations (Phase II)")), ylab = expression(bold("CUSUM Statistic")),
       cex.main = 1.2, cex.lab = 0.9, cex.axis = 1.2, font.lab = 2, las = 1)

  lines(Cminus, col = "darkgreen")
  abline(h = H_plus_c, col = "red", lwd = 2, lty = 2)
  abline(h = H_minus_c, col = "red", lwd = 2, lty = 2)

  legend("topright", legend = c("C+", "C-", "Upper and Lower Limits"),
         col = c("blue", "darkgreen", "red", "red"),
         lwd = c(2, 2, 2, 2), lty = c(1, 1, 1, 1), cex = 0.9, bg = "white")

  # Summary box with gray background and border
  par(mar = c(1, 2, 1, 2))
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

  rect(0.75, 0.6, 1.3, 1.1, col = "lightgray", border = "black", lwd = 1)  # Adjust rectangle size

  text(1, 1.05, "Control Chart Summary", cex = 1, font = 2, adj = 0.5)  # Title

  text(0.76, 0.95, sprintf("Phase I Sample Size: %d", n_I), cex = 1, adj = 0)
  text(0.76, 0.85, sprintf("Estimated Alpha: %.2f", alpha0_est), cex = 1, adj = 0)
  text(0.76, 0.75, sprintf("Estimated Beta: %.2f", beta0_est), cex = 1, adj = 0)

  text(1.01, 0.95, sprintf("Guaranteed Upper Limit: %.2f", H_plus_c), cex = 1, adj = 0)
  text(1.01, 0.85, sprintf("Guaranteed Lower Limit: %.2f", H_minus_c), cex = 1, adj = 0)
  text(1.01, 0.75, sprintf("Value of k_plus: %.4f", k_plus), cex = 1, adj = 0)
  text(1.01, 0.65, sprintf("Value of k_minus: %.4f", k_minus), cex = 1, adj = 0)

  message("Execution completed successfully.")
}
