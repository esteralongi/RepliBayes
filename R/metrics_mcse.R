## metrics_mcse.R -- replication probabilities WITH Monte Carlo standard errors.
##
## MCSE(p_hat) = sqrt( p_hat (1 - p_hat) / ESS ), where
##   - unconditional metrics : ESS = effective sample size of the 0/1 indicator
##                             (posterior::ess_basic, keeps the MCMC autocorrelation);
##   - conditional  metrics  : ESS = number of draws satisfying the conditioning
##                             event.
## The metric DEFINITIONS come from .metric_entries() (metrics.R), so the point
## estimates here are identical to compute_replication_probs_*(). Extract draws
## with permuted = FALSE so the iterations x chains structure (and, for the
## retrospective metrics, the draw-by-draw pairing) is preserved.

## ESS of a 0/1 series [iterations x chains] (NA if constant)
.ess_of <- function(x_ic) {
  x_ic <- as.matrix(x_ic)
  if (length(unique(as.vector(x_ic))) < 2) return(NA_real_)
  da  <- posterior::as_draws_array(array(x_ic, dim = c(nrow(x_ic), ncol(x_ic), 1L)))
  out <- suppressWarnings(posterior::ess_basic(da))
  if (!is.finite(out) || out <= 0) NA_real_ else as.numeric(out)
}

## unconditional metric: p and MCSE = sqrt(p(1-p)/ESS)
.row_uncond <- function(name, u_ic) {
  u_ic <- as.matrix(u_ic); N <- length(u_ic)
  p   <- mean(u_ic); k <- round(p * N)
  ess <- .ess_of(u_ic); if (is.finite(ess)) ess <- min(ess, N)
  mcse <- if (p > 0 && p < 1 && is.finite(ess)) sqrt(p * (1 - p) / ess) else NA_real_
  tibble::tibble(metric = name, p = p, k = k, n = N, ess = ess, mcse = mcse)
}

## conditional metric P(A|B): p = kA/nB, MCSE = sqrt(p(1-p)/nB)
## num_ic must already be A n B; den_ic is the conditioning event B
.row_cond <- function(name, num_ic, den_ic) {
  nB <- sum(den_ic); kA <- sum(num_ic)
  if (nB == 0)
    return(tibble::tibble(metric = name, p = NA_real_, k = 0, n = 0,
                          ess = NA_real_, mcse = NA_real_))
  p    <- kA / nB
  mcse <- if (p > 0 && p < 1) sqrt(p * (1 - p) / nB) else NA_real_
  tibble::tibble(metric = name, p = p, k = kA, n = nB, ess = nB, mcse = mcse)
}

#' Replication probabilities with Monte Carlo standard errors
#'
#' Study-level and (optionally) generative-consistency replication
#' probabilities for a general number of studies `S` and consensus level `k`,
#' each reported with its Monte Carlo standard error. This is the function used
#' to build the supplement tables.
#'
#' @param beta A list of `S` posterior-draw matrices \[iterations x chains\] of
#'   the study-specific effects. Plain vectors are also accepted (treated as a
#'   single chain). Extract with `rstan::extract(fit, permuted = FALSE)` to keep
#'   the chain structure so the ESS reflects the MCMC autocorrelation.
#' @param generative_beta Optional matrix of posterior draws of the generative
#'   effect (`mu_beta`). If `NULL` (default) only the study-level metrics are
#'   returned (independence limit); if supplied, the generative-consistency
#'   metrics are added (hierarchical model).
#' @param eps Practical-relevance threshold on the scale of the effects.
#' @param consensus_level Consensus level(s): at least `consensus_level` of the
#'   `S` studies agree. A vector is allowed (one block of metrics per value).
#'   Default 2.
#' @param min_corroborating For the per-study conditional metrics, the minimum
#'   number of *other* studies required to corroborate the discovery of the
#'   reference study. Default 1.
#'
#' @return A [tibble][tibble::tibble] with one row per metric and columns
#'   `metric`, `p`, `k`, `n`, `ess`, `mcse`. (Here the `k` column is the count of
#'   1s behind each probability, not the consensus level, which is in the metric
#'   name.)
#'
#' @details For a metric with all-zero draws the rule of three gives an upper
#'   bound `< 3/n`; conditional metrics with fewer than 30 conditioning draws
#'   should be flagged as unreliable when reported.
#' @export
replication_ess <- function(beta, generative_beta = NULL, eps, consensus_level = 2, min_corroborating = 1) {
  E <- .metric_entries(beta, generative_beta = generative_beta, eps = eps,
                       consensus_level = consensus_level, min_corroborating = min_corroborating)
  rows <- lapply(E, function(e) {
    if (e$type == "uncond") .row_uncond(e$name, e$ind)
    else                    .row_cond(e$name, e$num, e$den)
  })
  do.call(rbind, rows)
}

#' Retrospective predictive replication probabilities against a left-out study
#'
#' Paired-indicator replication probabilities that compare the observed effect
#' of a left-out study with its predicted effect from a body of evidence, each
#' with a Monte Carlo standard error. This function is agnostic to how many
#' studies formed the body of evidence (that choice is made upstream when
#' `beta_pred` is produced; see [fit_replicability()]).
#'
#' @param beta_obs Matrix \[iterations x chains\] of posterior draws of the
#'   left-out study's effect, fit on that study alone.
#' @param beta_pred Matrix \[iterations x chains\] of predicted effects for the
#'   left-out study, drawn from the model fit on the body of evidence. Paired
#'   draw-by-draw with `beta_obs`.
#' @param eps Practical-relevance threshold on the scale of the effects.
#'
#' @return A [tibble][tibble::tibble] with columns `metric`, `p`, `k`, `n`,
#'   `ess`, `mcse` for the retrospective metrics (`P_overall_l`, `P_null_l`,
#'   `P_non_null_l`, `P_pos_l`, `P_neg_l`).
#'
#' @details `beta_obs` and `beta_pred` come from disjoint data, so the paired
#'   indicator is a valid estimator. Extract both with `permuted = FALSE` (and a
#'   fixed seed for `beta_pred`) so the pairing is reproducible.
#' @export
predictive_retro <- function(beta_obs, beta_pred, eps) {
  beta_obs <- as.matrix(beta_obs); beta_pred <- as.matrix(beta_pred)
  Mo <- abs(beta_obs) > eps; Mr <- abs(beta_pred) > eps
  ind <- list(
    P_overall_l  = (!Mo & !Mr) | (beta_obs >  eps & beta_pred >  eps) | (beta_obs < -eps & beta_pred < -eps),
    P_null_l     = (!Mo & !Mr),
    P_non_null_l = (beta_obs >  eps & beta_pred >  eps) | (beta_obs < -eps & beta_pred < -eps),
    P_pos_l      = (beta_obs >  eps & beta_pred >  eps),
    P_neg_l      = (beta_obs < -eps & beta_pred < -eps)
  )
  do.call(rbind, Map(.row_uncond, names(ind), ind))
}

#' Prospective predictive replication probabilities for a new study
#'
#' Marginal replication probabilities of a single future study-specific effect
#' drawn from a fitted hierarchical model, each with a Monte Carlo standard
#' error. The number of studies in the body of evidence used to fit that model
#' is chosen upstream (see [fit_replicability()]).
#'
#' @param beta_new Matrix \[iterations x chains\] of posterior-predictive draws
#'   of a new study's effect.
#' @param eps Practical-relevance threshold on the scale of the effects.
#'
#' @return A [tibble][tibble::tibble] with columns `metric`, `p`, `k`, `n`,
#'   `ess`, `mcse` for the prospective metrics (`P_non_null_p`, `P_null_p`,
#'   `P_pos_p`, `P_neg_p`).
#' @export
predictive_pro <- function(beta_new, eps) {
  beta_new <- as.matrix(beta_new)
  ind <- list(
    P_non_null_p = abs(beta_new) >  eps,
    P_null_p     = abs(beta_new) <= eps,
    P_pos_p      = beta_new >  eps,
    P_neg_p      = beta_new < -eps
  )
  do.call(rbind, Map(.row_uncond, names(ind), ind))
}
