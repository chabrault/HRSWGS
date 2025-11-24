#!/usr/bin/env Rscript

### --- PRE-FLIGHT CHECK FOR PACKAGE DEVELOPMENT --- ###
### Run this before pushing to GitHub or submitting anywhere ###

message("=== 1. Load devtools and dependencies ===")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("roxygen2", quietly = TRUE)) install.packages("roxygen2")
if (!requireNamespace("spelling", quietly = TRUE)) install.packages("spelling")
#if (!requireNamespace("lintr", quietly = TRUE)) install.packages("lintr")
if (!requireNamespace("goodpractice", quietly = TRUE)) install.packages("goodpractice")

library(devtools)

message("=== 2. Re-document package (roxygen) ===")
devtools::document()

message("=== 3. Run tests ===")
devtools::test()

# message("=== 4. Style check (optional but recommended) ===")
# lintr::lint_package()

message("=== 5. Check spelling of documentation ===")
spelling::spell_check_package()

message("=== 6. Check package metadata ===")
devtools::check_man()

message("=== 7. Run CRAN-style R CMD check ===")
devtools::check(args = c("--as-cran"))

message("=== 8. Check if package installs cleanly ===")
devtools::install(upgrade = "never")

message("=== 9. Goodpractice review (optional but useful) ===")
goodpractice::gp()

message("=== 10. Session info ===")
sessionInfo()

message("=== PRE-FLIGHT CHECK COMPLETE ===")
