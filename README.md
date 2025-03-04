---
editor_options: 
  markdown: 
    wrap: sentence
---

# LGCU: Learning Gamma CUSUM

## Overview

**LGCU (Learning Gamma CUSUM)** is an R package for designing and analyzing **CUSUM control charts** applied to **Gamma distributions** with guaranteed performance.
The package provides methods for:

-   **Estimating Average Run Length (ARL)** for both upward and downward CUSUM charts.
-   **Visualizing CUSUM control charts** for process monitoring.
-   **Calibrating control limits** using Monte Carlo simulations.
-   **Applying cautious learning mechanisms** to dynamically update parameters.

The methodologies implemented follow the work of **Madrid-Alvarez, García-Díaz, and Tercero-Gómez (2024)**, ensuring rigorous statistical foundations and practical applications in industrial quality control.

------------------------------------------------------------------------

## Installation

To install the package from GitHub:

``` r
# Install devtools if not already installed
install.packages("devtools")

# Install LGCU from GitHub
devtools::install_github("your_github_username/LGCU")
```

To load the package:

``` r
library(LGCU)
```

------------------------------------------------------------------------

## Functions & Usage

### **1. ARL Estimation for CUSUM Control Charts**

### **Downward CUSUM ARL Calculation (Guaranteed Performance)**

``` r
GICARL_CUSUM_down(alpha, beta, alpha_est, beta_est, beta_ratio, H_minus, H_delta, m)

# Example:
GICARL_CUSUM_down(alpha = 1, beta = 1, alpha_est = 1, beta_est = 1,
                   beta_ratio = 1/2.5, H_minus = -2.792, H_delta = 0, m = 100)
```

### **Upward CUSUM ARL Calculation (Guaranteed Performance)**

``` r
GICARL_CUSUM_up(alpha, beta, alpha_est, beta_est, beta_ratio, H_plus, H_delta, m)

# Example:
GICARL_CUSUM_up(alpha = 1, beta = 1, alpha_est = 1.2, beta_est = 0.8,
                 beta_ratio = 2, H_plus = 6.5081, H_delta = 2.9693, m = 100)
```

### **Downward CUSUM ARL Calculation (With Learning)**

``` r
ARL_Clminus(alpha, beta, alpha0_est, beta0_est, known_alpha, beta_ratio, H_delta, H_minus, n_I, replicates, K_l, delay, tau)

# Example:
ARL_Clminus(
   alpha = 1,
   beta = 1,
   alpha0_est = 1.067,  # alpha = known_alpha
   beta0_est = 0.2760,   # Estimated Beta
   known_alpha = TRUE,
   beta_ratio = 1/2,
   H_delta = 0.6946,
   H_minus = -4.8272,
   n_I = 500,
   replicates = 1000,
   K_l = 0.5,
   delay = 25,
   tau = 1
)
```

### **Upward CUSUM ARL Calculation (With Learning)**

``` r
ARL_Clplus(alpha, beta, alpha0_est, beta0_est, known_alpha, beta_ratio, H_delta, H_plus, n_I, replicates, K_l, delay, tau)

# Example:
ARL_Clplus(
  alpha = 1,
  beta = 1,
  alpha0_est = 1,  # alpha = known_alpha
  beta0_est = 1.1,   # Estimated Beta
  known_alpha = TRUE,
  beta_ratio = 2,
  H_delta = 4.2433,
  H_plus = 8.7434,
  n_I = 200,
  replicates = 100,
  K_l = 2,
  delay = 25,
  tau = 1
)
```

------------------------------------------------------------------------

### **2. CUSUM Control Chart Visualization**

### **Downward CUSUM Chart**

``` r
plot_GICCdown_chart(alpha, beta, beta_ratio, H_delta, H_minus, n_I, n_II, faseI, faseII, known_alpha)

# Example:
plot_GICCdown_chart(alpha = 3, beta = 1, beta_ratio = 1/2, H_delta = 0.9596, H_minus = -4.6901,
                    n_I = 100, n_II = 200, faseI = NULL, faseII = NULL, known_alpha = FALSE)
```

### **Upward CUSUM Chart**

``` r
plot_GICCup_chart(alpha, beta, beta_ratio, H_delta, H_plus, n_I, n_II, faseI, faseII, known_alpha)

# Example:
plot_GICCup_chart(alpha = 1, beta = 1, beta_ratio = 2, H_delta = 0, H_plus = 5.16,
                  n_I = 100, n_II = 200, faseI = NULL, faseII = NULL, known_alpha = TRUE)
```

### **Bidirectional CUSUM Chart**

``` r
plot_GICC_chart2(alpha, beta, beta_ratio_plus, beta_ratio_minus, H_delta_plus, H_plus, H_delta_minus, H_minus, n_I, n_II, faseI, faseII, known_alpha)

# Example:
plot_GICC_chart2(alpha = 1, beta = 1, beta_ratio_plus = 2, beta_ratio_minus = 0.5,
                 H_delta_plus = 2.0, H_plus = 5.0, H_delta_minus = 1.5, H_minus = -4.5,
                 n_I = 100, n_II = 200, faseI = NULL, faseII = NULL, known_alpha = TRUE)
```

------------------------------------------------------------------------

### **3. CUSUM Chart with Learning Mechanism**

### **Downward CUSUM with Learning**

``` r
plot_GICCLdown_Chart(alpha, beta, beta_ratio, H_delta, H_minus, known_alpha, k_l, delay, tau, n_I, n_II, faseI, faseII)

# Example:
plot_GICCLdown_Chart(alpha = 1, beta = 1, beta_ratio = 1/2, H_delta = 4.2433, H_minus= -4.8257,
                     known_alpha = FALSE, k_l = 0.739588, delay = 25, tau = 1,
                     n_I = 200, n_II = 700, faseI = NULL, faseII = NULL)
```

### **Upward CUSUM with Learning**

``` r
plot_GICCLup_Chart(alpha, beta, beta_ratio, H_delta, H_plus, known_alpha, k_l, delay, tau, n_I, n_II, faseI, faseII)

# Example:
plot_GICCLup_Chart(alpha = 1, beta = 1, beta_ratio = 2, H_delta = 2.9819, H_plus = 6.5081,
                    known_alpha = TRUE, k_l = 2, delay = 25, tau = 1,
                    n_I = 200, n_II = 710, faseI = NULL, faseII = NULL)
```

### **Bidirectional CUSUM with Learning**

``` r
plot_GICCL_chart2(alpha, beta, beta_ratio_plus, beta_ratio_minus, H_delta_plus, H_plus, H_delta_minus, H_minus, known_alpha, k_l, delay, tau, n_I, n_II, faseI, faseII)

# Example:
plot_GICCL_chart2(alpha = 1, beta = 1, beta_ratio_plus = 2, beta_ratio_minus = 0.5,
                  H_delta_plus = 3.0, H_plus = 6.5, H_delta_minus = 2.0, H_minus = -5.0,
                  known_alpha = TRUE, k_l = 2, delay = 25, tau = 1,
                  n_I = 200, n_II = 700, faseI = NULL, faseII = NULL)
```

------------------------------------------------------------------------

### **4. Calibration of `H_delta` for Guaranteed Performance**

### **Downward Calibration**

``` r
getDeltaH_down(n_I, alpha, beta, beta_ratio, H_minus, a, b, ARL_esp, m, N_init, N_final, known_alpha)

# Example:
getDeltaH_down(n_I = 100, alpha = 1, beta = 1, beta_ratio = 1/2, H_minus = -4.1497,
               a = 0.1, b = 0.05, ARL_esp = 370, m = 100, N_init = 10, N_final = 1000, known_alpha = TRUE)
```

### **Upward Calibration**

``` r
getDeltaH_up(n_I, alpha, beta, beta_ratio, H_plus, a, b, ARL_esp, m, N_init, N_final, known_alpha)

# Example:
getDeltaH_up(n_I = 100, alpha = 1, beta = 1, beta_ratio = 2, H_plus = 6.8313,
              a = 0.1, b = 0.05, ARL_esp = 370, m = 100, N_init = 10, N_final = 1000, known_alpha = TRUE)
```

------------------------------------------------------------------------

# **5. Calibration of `H_delta` with Learning**

### **Downward Calibration with Learning**

``` r
getDeltaHL_down(n_I, alpha, beta, beta_ratio, H_minus, a, b, ARL_esp, replicates, N_init, N_final, known_alpha, K_l, delay, tau)

# Example:
getDeltaHL_down(n_I = 200, alpha = 1, beta = 1, beta_ratio = 1/1.5,
                H_minus = -6.2913, a = 0.1, b = 0.05, ARL_esp = 370,
                replicates = 10, N_init = 100, N_final = 1000, known_alpha = TRUE,
                K_l = 0.7, delay = 25, tau = 1)
```

### **Upward Calibration with Learning**

``` r
getDeltaHL_up(n_I, alpha, beta, beta_ratio, H_plus, a, b, ARL_esp, replicates, N_init, N_final, known_alpha, K_l, delay, tau)

# Example:
getDeltaHL_up(n_I = 200, alpha = 1, beta = 1, beta_ratio = 2,
              H_plus = 6.8313, a = 0.1, b = 0.05, ARL_esp = 370,
              replicates = 100, N_init = 100, N_final = 500, known_alpha = TRUE,
              K_l = 2, delay = 25, tau = 1)
```

------------------------------------------------------------------------

## References

📖 **Madrid-Alvarez, H. M., García-Díaz, J. C., & Tercero-Gómez, V. G. (2024).** *A CUSUM control chart for gamma distribution with guaranteed performance.* Quality and Reliability Engineering International.

📖 **Madrid-Alvarez, H. M., García-Díaz, J. C., & Tercero-Gómez, V. G. (2024).** *A CUSUM control chart for the Gamma distribution with cautious parameter learning.* Quality Engineering.

------------------------------------------------------------------------

## Contributing

If you would like to contribute, please follow these steps: 1.
Fork the repository.
2.
Create a new branch (`feature-newFunction`).
3.
Commit your changes.
4.
Push to the branch and submit a Pull Request.

------------------------------------------------------------------------

## License

This package is licensed under the **MIT License**, which means you are free to use, modify, and distribute the software, provided that the original copyright and license notice is included in all copies.
This allows both open-source and commercial use.

------------------------------------------------------------------------

## Contact

For any questions or suggestions, please contact: 📩 [**harold.madrid\@unisimon.edu.co**](mailto:harold.madrid@unisimon.edu.co){.email}
