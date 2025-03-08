# CRAN submission LGCU 0.1.2

## Test environments
- Local Windows 11, R 4.4.3
- win-builder (release and devel)
- macOS Ventura, R Under development
- Ubuntu 22.04 LTS, R 4.4.3

## R CMD check results
- Windows, macOS y Ubuntu: **0 errors ✔ | 0 warnings ✔ | 1 note ✖**
- Pruebas ejecutadas con `devtools::check()` y `rhub::rhub_check()`

## Comments
- This is the first submission to CRAN.
- The **NOTE** detected in CRAN checks refers to the term 
  **"CUSUM"** in the DESCRIPTION file.
- **CUSUM (Cumulative Sum Control Chart)** is a widely used 
  statistical method in Quality Control, and its inclusion in 
  the package description is appropriate.
  Terms flagged as possible misspellings ("CUSUM", "thresholding", "et al.") 
  are correct technical or standard Latin terms used in the statistical literature.
- No other known issues. Ready for review.
