#' @title Downward CUSUM Control Chart with Cautious Learning and Guaranteed Performance
#'
#' @description This function generates a downward CUSUM control chart for a Gamma distribution, incorporating a
#' cautious parameter updating mechanism based on guaranteed performance.
#'
#' It enables dynamic process monitoring, ensuring progressive adaptation to changes in the distribution.
#'
#' This approach follows the methodology presented in the work of Madrid-Alvarez, García-Díaz, and Tercero-Gómez (2024), where
#' a cautious learning scheme for parameter updating in CUSUM control charts applied to Gamma distributions is proposed.
#'
#' The implementation captures changes in the distribution and adjusts the
#' control limits to enhance the detection of process variations.
#'
#' ## Features:
#' - If the user does not provide Phase I and Phase II data, the function generates them automatically.
#' - If `known_alpha = TRUE`, `alpha` is fixed and not estimated.
#' - If `known_alpha = FALSE`, `alpha` is estimated from Phase I data.
#' - Dynamic control limits and a summary of parameters are included.
#' - Integrates a cautious learning scheme using the parameters `k_l`, `delay`, and `tau`.
#'
#' ## **Recommendations**
#' The parameters `k_l`, `delay`, and `tau` are part of the cautious learning mechanism of the CUSUM chart. These values enable the dynamic
#' updating of `beta0_est` and `H_minus`, ensuring that the control chart progressively adapts to process changes, improving sensitivity in detecting deviations.
#'
#' For proper implementation, it is recommended to reference the values proposed in:
#'
#' **Madrid-Alvarez, H. M., García-Díaz, J. C., & Tercero-Gómez, V. G. (2024).**
#' *A CUSUM control chart for the Gamma distribution with cautious parameter learning.*
#' Quality Engineering, 1-23.
#'
#' While these parameters have been tested and validated in the referenced article, users can adjust them based on the specific characteristics
#' of their process, considering factors such as system variability and desired update frequency.
#'
#' Additionally, if detailed guidance on selecting values for `H_delta` and `H_minus` is needed, it is recommended to review the referenced article,
#' which presents calibration and adjustment strategies for these limits to ensure optimal control chart performance.
#'
#' @param H_delta Increment of the lower control limit.
#' @param beta_ratio Ratio between `beta` and its estimation.
#' @param alpha Shape parameter of the Gamma distribution (if `alpha_conocido = TRUE`).
#' @param beta Scale parameter of the Gamma distribution.
#' @param H_minus Initial lower control limit of the CUSUM chart.
#' @param known_alpha Indicates whether `alpha` is known (`TRUE`) or should be estimated (`FALSE`).
#' @param k_l Secondary control threshold used in the learning logic.
#' @param delay Number of observations before updating `beta0_est` and `H_minus_c`.
#' @param tau Time point at which the `beta` parameter changes.
#' @param n_I Sample size in Phase I (if `faseI` is not provided).
#' @param n_II Sample size in Phase II (if `faseII` is not provided).
#' @param faseI Data sample from Phase I (numeric vector). If `NULL`, it is generated internally.
#' @param faseII Data sample from Phase II (numeric vector). If `NULL`, it is generated internally.
#'
#' @return A plot showing the evolution of the downward CUSUM statistic with cautious learning, including:
#' - The dynamically adjusted accumulated values of the CUSUM statistic.
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
#' # Option 1: Providing Phase I and Phase II data
#' phaseI_data <- rgamma(n = 200, shape = 1, scale = 1)
#' phaseII_data <- rgamma(n = 710, shape = 1, scale = 1)
#' plot_GICCLdown_Chart(alpha = 1, beta = 1, beta_ratio = 1/2, H_delta = 4.2433,
#'                      H_minus= -4.8257, known_alpha = FALSE, k_l = 0.739588,
#'                      delay = 25, tau = 1, n_I = 200, n_II = 700,
#'                      faseI = phaseI_data, faseII = phaseII_data)
#'
#' # Option 2: Without providing data, the function automatically generates them
#' plot_GICCLdown_Chart(alpha = 1, beta = 1, beta_ratio = 1/2, H_delta = 1.6763,
#'                      H_minus = -4.8257, known_alpha = FALSE, k_l = 0.739588,
#'                      delay = 25, tau = 1, n_I = 200,
#'                      n_II = 710, faseI = NULL, faseII = NULL)
#'
#'
#'
plot_GICCLdown_Chart <- function(alpha, beta, beta_ratio, H_delta, H_minus,
                                 known_alpha, k_l, delay, tau,
                                 n_I, n_II, faseI = NULL, faseII = NULL) {

  # Verify that the MASS package is installed
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("The 'MASS' package is required but not installed. Use install.packages('MASS') to install it.")
  }

  # Generate Phase I data if not provided
  if (is.null(faseI)) {
    faseI <- rgamma(n = n_I, shape = alpha, scale = beta)
  } else {
    n_I <- length(faseI)
  }

  # Generate Phase II data if not provided
  if (is.null(faseII)) {
    faseII <- rgamma(n = n_II, shape = alpha, scale = beta)
  } else {
    n_II <- length(faseII)
  }

  # Parameter estimation in Phase I
  estimator <- MASS::fitdistr(x = faseI, 'gamma', method = "Nelder-Mead")
  alpha0_est <- ifelse(known_alpha, alpha, as.numeric(estimator$estimate[1]))
  beta0_est <- mean(faseI) / alpha0_est

  # Downward CUSUM parameter calculation
  k_minus <- -((alpha0_est * beta * beta_ratio * log(beta / beta_ratio)) / (beta_ratio - beta))
  H_minus_c <- H_minus - H_delta
  H_minus_c_initial <- H_minus_c

  # Second Phase: Monte Carlo Learning Logic
  N <- length(faseII)
  Cminus <- rep(NA, N)
  Cminus_l <- rep(NA, N)
  H_minus_c_values <- rep(NA, N)
  beta0_est_values <- rep(beta0_est, N)

  Cminus[1] <- min(0, (faseII[1] / beta0_est) - k_minus)
  Cminus_l[1] <- min(0, (faseII[1] / beta0_est) - k_l)
  H_minus_c_values[1] <- H_minus_c

  for (i in 2:N) {
    Cminus[i] <- min(0, Cminus[i - 1] + (faseII[i] / beta0_est) - k_minus)
    Cminus_l[i] <- min(0, Cminus_l[i - 1] + (faseII[i] / beta0_est) - k_l)

    if (Cminus_l[i] == 0 && i > delay && i %% delay == 0) {
      beta0_est <- mean(faseII[1:i]) / alpha0_est
      H_delta_n <- H_delta * (sqrt(n_I / (n_I + i)))^1.7
      H_minus_c <- H_minus - H_delta_n
    }

    H_minus_c_values[i] <- H_minus_c
    beta0_est_values[i] <- beta0_est

    if (i == tau) {
      beta <- beta_ratio
    }
  }

  H_minus_c_final <- H_minus_c

  # Adjust layout for graph and summary box
  layout(matrix(c(1,2), nrow = 2), heights = c(2, 1))

  # CUSUM chart with H_minus_c
  ylim <- c(min(c(Cminus, H_minus_c_values), na.rm = TRUE) - 2, 0)
  plot(Cminus, ylim = ylim, type = "l", col = "blue",
       main = "CUSUM Control Chart with Cautious Learning - Downward Detection",
       xlab = "Observations (Phase II)", ylab = expression(C^"-"),
       cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.5)
  lines(1:N, H_minus_c_values, col = "red", type = "l")

  # Summary box
  par(mar = c(2, 2, 2, 2))
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  rect(0.3, 0.5, 1.7, 1.4, col = "lightgray", border = "black", lwd = 1)

  text(1, 1.35, "Control Chart Summary", cex = 1, font = 2, adj = 0.5)
  text(0.6, 1.15, sprintf("Initial Control Limit: %.2f", H_minus_c_initial), cex = 1, adj = 0)
  text(0.6, 0.95, sprintf("Updated Control Limit (Final): %.2f", H_minus_c_final), cex = 1, adj = 0)
  text(0.6, 0.75, sprintf("Phase I Sample Size: %d", n_I), cex = 1, adj = 0)

  text(1.2, 1.15, sprintf("Estimated Alpha: %.2f", alpha0_est), cex = 1, adj = 0)
  text(1.2, 0.95, sprintf("Estimated Beta: %.2f", beta0_est), cex = 1, adj = 0)
  text(1.2, 0.75, sprintf("Value of k_minus: %.4f", k_minus), cex = 1, adj = 0)

  message("Execution completed successfully.")
}

