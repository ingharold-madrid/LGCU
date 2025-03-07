#' @title Upward CUSUM Control Chart for Gamma Distribution with Guaranteed Performance
#'
#' @description This function generates an upward CUSUM control chart for a Gamma distribution, displaying
#' the evolution of the CUSUM statistic, control limits, and a summary of the parameters.
#'
#' Based on the approach presented by Madrid‐Alvarez, García‐Díaz, and Tercero‐Gómez (2024), this implementation enables the evaluation and visualization
#' of the monitored process using a CUSUM chart adapted to Gamma distributions with guaranteed performance.
#'
#' Specifically, the library incorporates a Monte Carlo model for simulating the control chart behavior,
#' allowing the Gamma distribution to be estimated in Phase I or using predefined values. Additionally,
#' it provides a clear graphical representation of the CUSUM statistic's evolution, ensuring proper
#' calibration and process control.
#'
#' ### **Recommendations**
#'
#' To check specific values for `H_delta` and `H_plus`, it is recommended to review the reference article:
#' Madrid‐Alvarez, H. M., García‐Díaz, J. C., & Tercero‐Gómez, V. G. (2024).
#' **A CUSUM control chart for gamma distribution with guaranteed performance**. Quality and Reliability Engineering International, 40(3), 1279-1301.
#'
#' ## Features:
#' - Based on a Monte Carlo model.
#' - Estimates the Gamma distribution in Phase I or uses predefined values.
#' - Plots the accumulated values of the CUSUM statistic with guaranteed performance.
#' - Includes control limits and a summary table.
#'
#' @param alpha Shape parameter of the Gamma distribution.
#' @param beta Scale parameter of the Gamma distribution.
#' @param beta_ratio Ratio between beta and its estimation.
#' @param H_delta Increment of the upper GIC limit.
#' @param H_plus Initial upper limit of the CUSUM chart.
#' @param n_I Sample size in Phase I (if `faseI` is not provided).
#' @param n_II Sample size in Phase II (if `faseII` is not provided).
#' @param faseI Sample data from Phase I (numeric vector). If `NULL`, it is generated with `rgamma()`.
#' @param faseII Sample data from Phase II (numeric vector). If `NULL`, it is generated with `rgamma()`.
#' @param known_alpha If `TRUE`, a known `alpha` is used; if `FALSE`, it is estimated.
#'
#' @return A plot displaying the evolution of the upward CUSUM statistic, including:
#' - The accumulated values of the CUSUM statistic.
#' - Control limits with guaranteed performance.
#' - A summary of the parameters used in the control chart.
#' @export
#' @importFrom stats rgamma
#' @importFrom graphics abline layout legend par rect text
#' @importFrom utils install.packages
#' @importFrom MASS fitdistr
#'
#' @examples
#' # Option 1: Automatically generate data with defined sample sizes
#' plot_GICCup_chart(
#'                   alpha = 1, beta = 1, beta_ratio = 2, H_delta = 0,
#'                   H_plus = 5.16, n_I = 100, n_II = 200, faseI = NULL,
#'                   faseII = NULL, known_alpha = TRUE
#'                   )
#'
#' # Option 2: Use custom data
#' phaseI_data <- rgamma(n = 100, shape = 1, scale = 1)
#' phaseII_data <- rgamma(n = 200, shape = 1, scale = 1)
#' plot_GICCup_chart(
#'                   alpha = 1, beta = 1, beta_ratio = 2, H_delta = 2.9693,
#'                   H_plus = 6.5081, n_I = 100, n_II = 200,
#'                   faseI = phaseI_data, faseII = phaseII_data,
#'                   known_alpha = TRUE
#'                   )
plot_GICCup_chart <- function(
                              alpha, beta, beta_ratio, H_delta, H_plus,
                              n_I, n_II, faseI = NULL, faseII = NULL, known_alpha
                              ) {

  # Load required packages
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("The 'MASS' package is required but not installed. Use install.packages('MASS') to install it.")
  }

  # Check if Phase I data is provided
  if (is.null(faseI)) {
    faseI <- rgamma(n = n_I, shape = alpha, scale = beta)  # Generate if no data is provided
  } else {
    n_I <- length(faseI)  # Adjust muestra_n if user provides data
  }

  # Parameter estimation in Phase I
  estimator <- MASS::fitdistr(x = faseI, 'gamma', method = "Nelder-Mead")

  if (known_alpha) {
    alpha0_est <- alpha  # Known alpha
  } else {
    alpha0_est <- as.numeric(estimator$estimate[1])
  }

  beta0_est <- mean(faseI) / alpha0_est

  # Compute CUSUM parameters
  k_plus <- -((alpha0_est * beta * beta_ratio * log(beta / beta_ratio)) / (beta_ratio - beta))

  # Upper control limit with guaranteed performance
  H_plus_c <- H_plus + H_delta

  # Check if Phase II data is provided
  if (is.null(faseII)) {
    faseII <- rgamma(n = n_II, shape = alpha, scale = beta)  # Generate if no data is provided
  } else {
    n_II <- length(faseII)  # Adjust n_II if user provides data
  }

  # Generate the CUSUM control chart in Phase II
  Cplus <- numeric(n_II)
  Cplus[1] <- max(0, 0 + (faseII[1] / beta0_est) - k_plus)

  for (i in 2:n_II) {
    Cplus[i] <- max(0, Cplus[i - 1] + (faseII[i] / beta0_est) - k_plus)
  }

  # Adjust layout for better visualization
  layout(matrix(c(1,2), nrow = 2), heights = c(2, 1))

  # Plot the control chart
  par(mar = c(5, 5, 4, 2))
  ylim <- c(0, H_plus_c + 5)

  plot(Cplus, ylim = ylim, type = "l", col = "blue",
       main = "Gamma CUSUM Control Chart with Guaranteed Performance",
       xlab = "Observations (Phase II)", ylab = "CUSUM Statistic",
       cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2)

  abline(h = H_plus_c, col = "red", lwd = 2, lty = 2)  # Upper control limit
  legend("topright", legend = c("CUSUM", "Control Limit"),
         col = c("blue", "red"), lwd = c(2, 2), lty = c(1, 2), cex = 1)

  # Summary box with gray background and border
  par(mar = c(2, 2, 2, 2))
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

  # Draw a gray rectangle for better readability
  rect(0.3, 0.4, 1.7, 1.4, col = "lightgray", border = "black", lwd = 1)

  # Summary title
  text(1.0, 1.35, "Control Chart Summary", cex = 1, font = 2, adj = 0.5)

  # First column (Phase I values)
  text(0.6, 1.15, sprintf("Phase I Sample Size: %d", n_I), cex = 1, adj = 0)
  text(0.6, 0.95, sprintf("Estimated Alpha: %.2f", alpha0_est), cex = 1, adj = 0)
  text(0.6, 0.75, sprintf("Estimated Beta: %.2f", beta0_est), cex = 1, adj = 0)

  # Second column (CUSUM control values)
  text(1.0, 1.15, sprintf("Guaranteed Control Limit: %.2f", H_plus_c), cex = 1, adj = 0)
  text(1.0, 0.95, sprintf("Value of k_plus: %.4f", k_plus), cex = 1, adj = 0)

  message("Execution completed successfully.")
}


