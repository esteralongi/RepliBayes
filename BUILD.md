# Building RepliBayes

The package was scaffolded without an R session, so two generated artefacts are
not yet present: the `man/*.Rd` help files and `data/synthetic_data.rda`. Both
are produced by the standard package-build steps below. Run these once from the
package root.

```r
setwd("RepliBayes")

install.packages(c("devtools", "roxygen2", "usethis"))   # if needed

# 1. Build the dataset object data/synthetic_data.rda from the shipped CSV
source("data-raw/synthetic_data.R")        # calls usethis::use_data()

# 2. Generate man/*.Rd and (re)write NAMESPACE from the roxygen comments
devtools::document()

# 3. Install (compiles the two Stan models on first use)
devtools::install()

# 4. Check
devtools::check()
```

## Notes

- **Stan models.** They live in `inst/stan/` and are compiled at runtime by
  `rstan::stan_model()` the first time a fitting function is called (cached for
  the rest of the session). No `src/` or precompilation is set up; a C++
  toolchain is still required. To precompile at install time instead, migrate to
  `rstantools::rstan_create_package()` and move the models to
  `inst/stan/include` / `src/stanExports_*`.
- **NAMESPACE.** A hand-written `NAMESPACE` is included so the package installs
  before `document()` is run; `devtools::document()` will regenerate it
  identically from the `@export` / `@importFrom` tags.
- **Dependencies.** `rstan`, `posterior`, `tibble`, `dplyr`, `stats`. The
  vignette additionally uses `tidyr`, `ggplot2`, `patchwork`, `BiocStyle`.
