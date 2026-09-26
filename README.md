# RepliBayes

`RepliBayes` provides a concrete, general implementation of the replicability 
framework of Alongi, Altoè & Parmigiani, *The Quantification of Replicability*, 
a set of functions to fit the model and compute the replication probabilities 
(with Monte Carlo standard errors), plus prior-sensitivity and simulation tools.

Methods: Alongi, Altoè and Parmigiani, *The Quantification of Replicability*.

## Installation

`RepliBayes` fits its models with [rstan](https://mc-stan.org/rstan/), so you
need a working C++ toolchain and the `rstan` package before installing.

**1. C++ toolchain** (needed by Stan): Windows →
[Rtools](https://cran.r-project.org/bin/windows/Rtools/); macOS → run
`xcode-select --install` in Terminal; Linux → a system C++ compiler (`g++`).

**2. rstan**:

```r
install.packages("rstan")
```

**3. RepliBayes**:

```r
# install.packages("remotes")
remotes::install_github("esteralongi/RepliBayes")
```

The install above is quick. Add `build_vignettes = TRUE` only if you also want
the tutorial vignette built locally — it compiles the Stan models and takes a
few minutes. The two Stan models are otherwise compiled from source the first
time you fit a model.

## Quick start

```r
library(RepliBayes)

data(synthetic_data)          # three groups acting as three studies

res <- fit_replicability(synthetic_data)   # fits both models, all metrics
res
res$metrics_hierarchical      # study-level + generative metrics, with MCSE
res$metrics_prospective       # a new study from the full model
res$metrics_retrospective     # predicting a left-out study
```

Supply your own priors or threshold:

```r
priors <- default_priors()
priors$mu_beta <- 0
res <- fit_replicability(synthetic_data, priors = priors, eps = 0.86)
```

## The three entry points

```r
# 1. Fit the framework + all metrics
res <- fit_replicability(synthetic_data)

# 2. Prior sensitivity: shift the prior on a chosen heterogeneity
sens <- sensitivity_prior(synthetic_data, target = "tau_beta")
round(sens$p_grid, 3)

# 3. Simulation calibration: choose what to vary in the generative model
sim <- simulate_replicability(res, synthetic_data, vary = "tau_beta", R = 30)
subset(sim, metric == "P_beta")
```

`sensitivity_prior()` and `simulate_replicability()` both accept
`target` / `vary = "tau_beta"` (effect), `"tau_alpha"` (intercept) or `"tau_sig"`
(residual scale).

The metrics are general: they work with any number of studies `S`, any
consensus level(s) `k` (argument `consensus_level`, default 2), and any minimum number of
agreeing studies `m` in the conditional metrics (argument `min_corroborating`, default 1). 
Study effects are passed to the metric functions as `beta`, a list of `S`
posterior-draw matrices (iterations × chains); plain vectors are also accepted.
For the retrospective/prospective analyses, `fit_replicability()` lets you choose the
body of evidence explicitly via `retro_target` / `retro_evidence` (size 1..S-1)
and `pro_evidence` (size 1..S).

## Building blocks

| function | purpose |
|---|---|
| `default_priors()` | default (illustrative) prior hyperparameters |
| `build_stan_data()` | assemble the Stan data list |
| `fit_hierarchical()`, `fit_independence()` | fit a single model |
| `fit_replicability()` | one-call wrapper (fits + all metrics) |
| `sensitivity_prior()` | prior-sensitivity over a heterogeneity's prior |
| `simulate_replicability()` | simulation calibration, varying the generative model |
| `replication_ess()` | metrics with MCSE (hierarchical / independence) |
| `predictive_retro()`, `predictive_pro()` | retrospective / prospective metrics |
| `compute_replication_probs_hier()`, `compute_replication_probs_indep()` | point-estimate metrics (used in the simulation calibration) |
| `replication_threshold_free()` | threshold-free summaries |

## Data

The GTEx donor-level data are controlled-access and are not redistributed. The
package ships a synthetic dataset (`synthetic_data`) with the same structure.

## Reproducing the paper

The analysis scripts live in the companion repository
[The-Quantification-Of-Replicability](https://github.com/esteralongi/The-Quantification-Of-Replicability)
(`analysis/`: `0_motivating_example.R`, `1_fit_replicability.R`,
`2_prior_sensitivity.R`, `3_simulation.R`, `4_figures.R`). The vignette walks
through the whole procedure on the synthetic data.
