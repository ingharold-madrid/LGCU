#' @title Upward CUSUM Control Chart with Cautious Learning and Guaranteed Performance
#'
#' @description This function generates an upward CUSUM control chart for a Gamma distribution, incorporating a
#' cautious parameter update mechanism based on guaranteed performance.
#'
#' It enables dynamic process monitoring, ensuring progressive adaptation to distribution changes.
#' This approach follows the methodology presented in the work of Madrid-Alvarez, García-Díaz, and Tercero-Gómez (2024),
#' where a cautious learning scheme for parameter updates in CUSUM control charts applied to Gamma distributions is proposed.
#'
#' The implementation captures distribution changes and adjusts the control limits to improve process variation detection.
#'
#' ## Features:
#' - If the user does not provide Phase I and Phase II data, the function automatically generates them.
#' - If `known_alpha = TRUE`, `alpha` is fixed and not estimated.
#' - If `known_alpha = FALSE`, `alpha` is estimated from Phase I data.
#' - Includes dynamic control limits and a summary of parameters.
#' - Integrates a cautious learning scheme using the parameters `k_l`, `delay`, and `tau`.
#'
#' ## **Recommendations**
#' The parameters `k_l`, `delay`, and `tau` are part of the cautious learning mechanism of the CUSUM chart. These values enable the dynamic
#' updating of `beta0_est` and `H_plus`, ensuring that the control chart progressively adapts to process changes, thus improving sensitivity in detecting deviations.
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
#' Additionally, if detailed guidance on selecting values for `H_delta` and `H_plus` is required, it is recommended to review the referenced article,
#' which presents calibration and adjustment strategies for these limits to ensure optimal control chart performance.
#'
#' @param H_delta Increment of the upper control limit.
#' @param beta_ratio Ratio between `beta` and its estimation.
#' @param alpha Shape parameter of the Gamma distribution (if `known_alpha = TRUE`).
#' @param beta Scale parameter of the Gamma distribution.
#' @param H_plus Initial upper control limit of the CUSUM chart.
#' @param known_alpha Indicates whether `alpha` is known (`TRUE`) or should be estimated (`FALSE`).
#' @param k_l Secondary control threshold used in the learning logic.
#' @param delay Number of observations before updating `beta0_est` and `H_plus_c`.
#' @param tau Time point at which the `beta` parameter changes.
#' @param n_I Sample size in Phase I (if `faseI` is not provided).
#' @param n_II Sample size in Phase II (if `faseII` is not provided).
#' @param faseI Data sample from Phase I (numeric vector). If `NULL`, it is generated internally.
#' @param faseII Data sample from Phase II (numeric vector). If `NULL`, it is generated internally.
#'
#' @return A plot showing the evolution of the upward CUSUM statistic with cautious learning, including:
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
#' plot_GICCLup_Chart(alpha = 1, beta = 1, beta_ratio = 2, H_delta = 4.2433,
#'                    H_plus = 8.9345, known_alpha = FALSE, k_l = 2, delay = 25,
#'                    tau = 1,n_I = 200, n_II = 700, faseI = phaseI_data,
#'                    faseII = phaseII_data)
#'
#' # Option 2: Without providing data, the function automatically generates them
#' plot_GICCLup_Chart(alpha = 1, beta = 1, beta_ratio = 2, H_delta = 2.9819,
#'                    H_plus = 6.5081, known_alpha = TRUE, k_l = 2, delay = 25,
#'                    tau = 1, n_I = 200, n_II = 710, faseI = NULL,
#'                    faseII = NULL)
#'
plot_GICCLup_Chart <- function(alpha, beta, beta_ratio, H_delta, H_plus,
                               known_alpha, k_l, delay, tau,
                               n_I, n_II, faseI = NULL, faseII = NULL) {

  # Verify that the MASS package is installed
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("The 'MASS' package is required but not installed. Use install.packages('MASS') to install it.")
  }


  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

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

  if (known_alpha) {
    alpha0_est <- alpha  # If known, do not estimate
  } else {
    alpha0_est <- as.numeric(estimator$estimate[1])  # Estimate from Phase I
  }

  beta0_est <- mean(faseI) / alpha0_est

  # Compute CUSUM parameters
  k_plus <- -((alpha0_est * beta * beta_ratio * log(beta / beta_ratio)) / (beta_ratio - beta))
  H_plus_c <- H_plus + H_delta
  H_plus_c_initial <- H_plus_c

  # Second Phase: Monte Carlo Learning Logic
  N <- length(faseII)
  Cplus <- rep(NA, N)
  Cplus_l <- rep(NA, N)
  H_plus_c_values <- rep(NA, N)
  beta0_est_values <- rep(beta0_est, N)  # Store estimated beta value at each iteration

  Cplus[1] <- max(0, (faseII[1] / beta0_est) - k_plus)
  Cplus_l[1] <- max(0, (faseII[1] / beta0_est) - k_l)
  H_plus_c_values[1] <- H_plus_c

  for (i in 2:N) {
    Cplus[i] <- max(0, Cplus[i - 1] + (faseII[i] / beta0_est) - k_plus)
    Cplus_l[i] <- max(0, Cplus_l[i - 1] + (faseII[i] / beta0_est) - k_l)

    if (Cplus_l[i] == 0 && i > delay && i %% delay == 0) {
      beta0_est <- mean(faseII[1:i]) / alpha0_est
      H_delta_n <- H_delta * (sqrt(n_I / (n_I + i)))^1.7
      H_plus_c <- H_plus + H_delta_n
    }

    H_plus_c_values[i] <- H_plus_c
    beta0_est_values[i] <- beta0_est  # Store the new estimated beta value

    if (i == tau) {
      beta <- beta_ratio
    }
  }

  H_plus_c_final <- H_plus_c  # Store final value

  par(mfrow = c(2,1))

  # CUSUM chart with H_plus_c
  par(mar = c(4, 4, 2, 1))
  ylim <- c(0, max(c(Cplus, H_plus_c_values), na.rm = TRUE) + 2)
  plot(Cplus, ylim = ylim, type = "l", col = "blue",
       main = "CUSUM Control Chart with Cautious Learning",
       xlab = expression(bold("Observations (Phase II)")), ylab = expression(bold(C^"+")),
       cex.axis = 1.2, cex.lab = 0.9, cex.main = 1.2)

  lines(1:N, H_plus_c_values, col = "red", type = "l")
  legend("topright",
         legend = c("CUSUM", "Control Limit"),
         col = c("blue", "red"),
         lwd = c(2, 2),
         lty = c(1, 1),
         cex = 0.9)

  # Summary box
  par(mar = c(1, 2, 1, 2))
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

  rect(0.75, 0.6, 1.3, 1.1, col = "lightgray", border = "black", lwd = 1)
  text(1, 1, "Control Chart Summary", cex = 1, font = 2, adj = 0.5)

  text(0.76, 0.9, sprintf("Initial C.L.: %.2f", H_plus_c_initial), cex = 1, adj = 0)
  text(0.76, 0.8, sprintf("Updated C.L. (Final): %.2f", H_plus_c_final), cex = 1, adj = 0)
  text(0.76, 0.7, sprintf("Phase I Sample Size: %d", n_I), cex = 1, adj = 0)

  text(1.01, 0.9, sprintf("Estimated Alpha: %.2f", alpha0_est), cex = 1, adj = 0)
  text(1.01, 0.8, sprintf("Estimated Beta: %.2f", beta0_est), cex = 1, adj = 0)
  text(1.01, 0.7, sprintf("Value of k_plus: %.4f", k_plus), cex = 1, adj = 0)

  message("Execution completed successfully.")
}
