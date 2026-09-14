## data-raw/synthetic_data.R
## Build data/synthetic_data.rda from the shipped CSV.
## Run once, from the package root:  source("data-raw/synthetic_data.R")

## Prefer the CSV in the source tree (works before the package is installed);
## fall back to the installed copy.
csv <- "inst/extdata/synthetic_data.csv"
if (!file.exists(csv))
  csv <- system.file("extdata", "synthetic_data.csv", package = "RepliBayes")
if (!nzchar(csv) || !file.exists(csv))
  stop("synthetic_data.csv not found. Run this from the package root ",
       "(setwd() to the RepliBayes folder).")

synthetic_data <- utils::read.csv(csv)
synthetic_data$study <- as.integer(synthetic_data$study)

## Save as data/synthetic_data.rda (base R; no usethis needed).
if (!dir.exists("data")) dir.create("data")
save(synthetic_data, file = "data/synthetic_data.rda", compress = "bzip2")
message("Wrote data/synthetic_data.rda (", nrow(synthetic_data), " rows).")
