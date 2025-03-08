# LGCU 0.1.2

## 🆕 New Features
- Implementation of CUSUM control charts for the Gamma distribution.
- Monte Carlo simulation methods for process monitoring.
- Dynamic parameter updating with a cautious learning scheme.
- Compatibility with Windows, macOS, and Linux.

## 🔧 Improvements
- Added explicit clarification of the acronym CUSUM in the DESCRIPTION file.
- Included formal references in the package description as requested by CRAN.
- Replaced `cat()` statements with structured `message()` outputs, adding a `verbose` option for clarity.
- Optimization of parameter estimation methods.
- Minor documentation enhancements and formatting improvements.

## 🐞 Bug Fixes
- Fixed issues in the graphical representation of control charts.
- Improved stability and accuracy in Monte Carlo simulations.
- Ensured that `par()` settings are correctly restored using `on.exit()`.
