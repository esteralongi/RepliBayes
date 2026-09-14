# RepliBayes

Bayesian quantification of replicability. `RepliBayes` fits a Bayesian
hierarchical model to study-level data and reports the replicability of an
effect as posterior probabilities, defined relative to a practical-relevance
threshold `eps` on the scale of the effects (not on p-values). It returns the
empirical, retrospective and prospective replication probabilities, each with a
Monte Carlo standard error.

Methods: Alongi, Altoè and Parmigiani, *The Quantification of Replicability*.

## Installation

```r
# install.packages("devtools")
devtools::install_github("esteralongi/RepliBayes", build_vignettes = TRUE)
```

`RepliBayes` uses [rstan](https://mc-stan.org/rstan/); a working C++ toolchain
is required to compile the two Stan models on first use.

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
priors$mu_a <- 0
res <- fit_replicability(synthetic_data, priors = priors, eps = 0.86)
```

## The three entry points

```r
# 1. Fit the framework + all metrics
res <- fit_replicability(synthetic_data)

# 2. Prior sensitivity: shift the prior on a chosen heterogeneity
sens <- sensitivity_prior(synthetic_data, target = "tau_a")
round(sens$p_grid, 3)

# 3. Simulation calibration: choose what to vary in the generative model
sim <- simulate_replicability(res, synthetic_data, vary = "tau_a", R = 30)
subset(sim, metric == "P_beta")
```

`sensitivity_prior()` and `simulate_replicability()` both accept
`target` / `vary = "tau_a"` (effect), `"tau_alpha"` (intercept) or `"tau_sig"`
(residual scale).

## Building blocks

| function | purpose |
|---|---|
| `default_priors()` | elicited prior hyperparameters |
| `build_stan_data()` | assemble the Stan data list |
| `fit_hierarchical()`, `fit_independence()` | fit a single model |
| `fit_replicability()` | one-call wrapper (fits + all metrics) |
| `sensitivity_prior()` | prior-sensitivity over a heterogeneity's prior |
| `simulate_replicability()` | simulation calibration, varying the generative model |
| `replication_ess()` | metrics with MCSE (hierarchical / independence) |
| `held_mcse()`, `held_pro_mcse()` | retrospective / prospective metrics |
| `compute_replication_probs_hier()`, `compute_replication_probs_indep()` | point-estimate metrics (used in the simulation calibration) |
| `replication_measures_multi()` | threshold-free summaries |

## Data

The GTEx donor-level data are controlled-access and are not redistributed. The
package ships a synthetic dataset (`synthetic_data`) with the same structure.

## Reproducing the paper

The `code/` folder of the source repository contains the scripts used for the
manuscript (`run_models.R`, `simulation.R`, `prior_sensitivity.R`); the vignette
walks through the whole procedure on the synthetic data.
