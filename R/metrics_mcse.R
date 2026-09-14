## metrics_mcse.R -- replication probabilities WITH Monte Carlo standard errors.
##
## MCSE(p_hat) = sqrt( p_hat (1 - p_hat) / ESS ), where
##   - unconditional metrics : ESS = effective sample size of the 0/1 indicator
##                             (posterior::ess_basic, keeps the MCMC autocorrelation);
##   - conditional  metrics  : ESS = number of draws satisfying the conditioning
##                             event.
## These functions use the paper-facing metric names (P_overall_2, P_cond_O1_m1,
## P_gen_pos_2, P_c_pos_2, ...). Extract draws with permuted = FALSE so the
## iterations x chains structure (and, for the retrospective metrics, the
## draw-by-draw pairing) is preserved.

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
#' probabilities, each reported with its Monte Carlo standard error. This is the
#' function used to build the supplement tables.
#'
#' @param a1,a2,a3 Matrices of posterior draws \[iterations x chains\] of the
#'   three study-specific effects. Plain vectors are also accepted (treated as a
#'   single chain). Extract with `rstan::extract(fit, permuted = FALSE)` to keep
#'   the chain structure so the ESS reflects the MCMC autocorrelation.
#' @param mu Optional matrix of posterior draws of the generative effect
#'   (`mu_a`). If `NULL` (default) only the study-level metrics are returned
#'   (independence limit); if supplied, the generative-consistency metrics are
#'   added (hierarchical model).
#' @param eps Practical-relevance threshold on the scale of the effects.
#'
#' @return A [tibble][tibble::tibble] with one row per metric and columns
#'   `metric`, `p`, `k`, `n`, `ess`, `mcse`.
#'
#' @details For a metric with all-zero (`k = 0`) draws the rule of three gives an
#'   upper bound `< 3/n`; conditional metrics with fewer than 30 conditioning
#'   draws should be flagged as unreliable when reported.
#' @export
replication_ess <- function(a1, a2, a3, mu = NULL, eps) {
  a1 <- as.matrix(a1); a2 <- as.matrix(a2); a3 <- as.matrix(a3)
  pos <- (a1 >  eps) + (a2 >  eps) + (a3 >  eps)
  neg <- (a1 < -eps) + (a2 < -eps) + (a3 < -eps)
  nul <- (abs(a1) <= eps) + (abs(a2) <= eps) + (abs(a3) <= eps)
  M1 <- abs(a1) > eps; M2 <- abs(a2) > eps; M3 <- abs(a3) > eps

  rows <- list(
    .row_uncond("P_overall_3",  (nul == 3) | (pos == 3) | (neg == 3)),
    .row_uncond("P_overall_2",  (nul >= 2) | (pos >= 2) | (neg >= 2)),
    .row_uncond("P_non_null_3", (pos == 3) | (neg == 3)),
    .row_uncond("P_non_null_2", (pos >= 2) | (neg >= 2)),
    .row_uncond("P_pos_3",  pos == 3),
    .row_uncond("P_neg_3",  neg == 3),
    .row_uncond("P_null_3", nul == 3),
    .row_uncond("P_pos_2",  pos >= 2),
    .row_uncond("P_neg_2",  neg >= 2),
    .row_uncond("P_null_2", nul >= 2),
    .row_cond("P_cond_O1_m1", (a1 >  eps & (a2 >  eps | a3 >  eps)) |
                              (a1 < -eps & (a2 < -eps | a3 < -eps)), M1),
    .row_cond("P_cond_O2_m1", (a2 >  eps & (a1 >  eps | a3 >  eps)) |
                              (a2 < -eps & (a1 < -eps | a3 < -eps)), M2),
    .row_cond("P_cond_O3_m1", (a3 >  eps & (a1 >  eps | a2 >  eps)) |
                              (a3 < -eps & (a1 < -eps | a2 < -eps)), M3)
  )

  if (!is.null(mu)) {
    mu <- as.matrix(mu)
    mp <- mu >  eps; mn <- mu < -eps; m0 <- abs(mu) <= eps
    rows <- c(rows, list(
      .row_uncond("P_beta",            abs(mu) > eps),
      .row_uncond("P_gen_overall_3",   (pos == 3 & mp) | (neg == 3 & mn) | (nul == 3 & m0)),
      .row_uncond("P_gen_overall_2",   (pos >= 2 & mp) | (neg >= 2 & mn) | (nul >= 2 & m0)),
      .row_uncond("P_gen_non_null_3",  (pos == 3 & mp) | (neg == 3 & mn)),
      .row_uncond("P_gen_non_null_2",  (pos >= 2 & mp) | (neg >= 2 & mn)),
      .row_uncond("P_gen_pos_3",  (pos == 3) & mp),
      .row_uncond("P_gen_neg_3",  (neg == 3) & mn),
      .row_uncond("P_gen_null_3", (nul == 3) & m0),
      .row_uncond("P_gen_pos_2",  (pos >= 2) & mp),
      .row_uncond("P_gen_neg_2",  (neg >= 2) & mn),
      .row_uncond("P_gen_null_2", (nul >= 2) & m0),
      .row_cond("P_c_pos_2",  mp & (pos >= 2), mp),
      .row_cond("P_c_neg_2",  mn & (neg >= 2), mn),
      .row_cond("P_c_null_2", m0 & (nul >= 2), m0)
    ))
  }
  do.call(rbind, rows)
}

#' Retrospective replication probabilities against a left-out study
#'
#' Paired-indicator replication probabilities that compare the observed effect
#' of a left-out study with its predicted effect from the remaining body of
#' evidence, each with a Monte Carlo standard error.
#'
#' @param a_obs_ic Matrix \[iterations x chains\] of posterior draws of the
#'   left-out study's effect, fit on that study alone.
#' @param a_hat_ic Matrix \[iterations x chains\] of predicted effects for the
#'   left-out study, drawn from the model fit on the remaining studies. Paired
#'   draw-by-draw with `a_obs_ic`.
#' @param eps Practical-relevance threshold on the scale of the effects.
#'
#' @return A [tibble][tibble::tibble] with columns `metric`, `p`, `k`, `n`,
#'   `ess`, `mcse` for the retrospective metrics (`P_overall_l`, `P_null_l`,
#'   `P_non_null_l`, `P_pos_l`, `P_neg_l`).
#'
#' @details `a_obs` and `a_hat` come from disjoint data, so the paired indicator
#'   is a valid estimator. Extract both with `permuted = FALSE` (and a fixed
#'   seed for `a_hat`) so the pairing is reproducible.
#' @export
held_mcse <- function(a_obs_ic, a_hat_ic, eps) {
  a_obs <- as.matrix(a_obs_ic); a_hat <- as.matrix(a_hat_ic)
  Mo <- abs(a_obs) > eps; Mr <- abs(a_hat) > eps
  ind <- list(
    P_overall_l  = (!Mo & !Mr) | (a_obs >  eps & a_hat >  eps) | (a_obs < -eps & a_hat < -eps),
    P_null_l     = (!Mo & !Mr),
    P_non_null_l = (a_obs >  eps & a_hat >  eps) | (a_obs < -eps & a_hat < -eps),
    P_pos_l      = (a_obs >  eps & a_hat >  eps),
    P_neg_l      = (a_obs < -eps & a_hat < -eps)
  )
  do.call(rbind, Map(.row_uncond, names(ind), ind))
}

#' Prospective replication probabilities for a new study
#'
#' Marginal replication probabilities of a single future study-specific effect
#' drawn from the full (three-study) hierarchical model, each with a Monte Carlo
#' standard error.
#'
#' @param a_new_ic Matrix \[iterations x chains\] of posterior-predictive draws
#'   of a new study's effect.
#' @param eps Practical-relevance threshold on the scale of the effects.
#'
#' @return A [tibble][tibble::tibble] with columns `metric`, `p`, `k`, `n`,
#'   `ess`, `mcse` for the prospective metrics (`P_non_null_p`, `P_null_p`,
#'   `P_pos_p`, `P_neg_p`).
#' @export
held_pro_mcse <- function(a_new_ic, eps) {
  a_new <- as.matrix(a_new_ic)
  ind <- list(
    P_non_null_p = abs(a_new) >  eps,
    P_null_p     = abs(a_new) <= eps,
    P_pos_p      = a_new >  eps,
    P_neg_p      = a_new < -eps
  )
  do.call(rbind, Map(.row_uncond, names(ind), ind))
}
