#' Synthetic dataset
#'
#' A synthetic dataset with the same structure as the GTEx analysis in the paper
#' (three groups playing the role of three studies), provided so that the whole
#' pipeline can be reproduced. The GTEx donor-level data are controlled-access
#' and cannot be redistributed.
#'
#' @format A data frame with 1070 rows and 3 columns:
#' \describe{
#'   \item{study}{integer study/group label (1, 2, 3).}
#'   \item{x}{predictor (genotype dosage, 0/1/2).}
#'   \item{y}{outcome (expression on the analysis scale).}
#' }
#'
#' @details Built from `inst/extdata/synthetic_data.csv` by
#'   `data-raw/synthetic_data.R`. If `data(synthetic_data)` is not yet available,
#'   run that script once (it calls `usethis::use_data()`), or read the CSV
#'   directly with
#'   `read.csv(system.file("extdata", "synthetic_data.csv", package = "RepliBayes"))`.
#'
#' @source Synthetic; see `data-raw/synthetic_data.R`.
"synthetic_data"
