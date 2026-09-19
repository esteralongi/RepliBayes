## metrics.R -- indicator definitions and point-estimate replication metrics.
##
## Everything is defined for a general number of studies S and a general
## consensus level (at least `consensus_level` of the S studies agree). A study-specific
## effect is, relative to the threshold eps:
##   positive if a_s > eps, negative if a_s < -eps, null if |a_s| <= eps.
##
## The single source of truth for the metric DEFINITIONS is .metric_entries();
## replication_ess() (metrics_mcse.R) adds Monte Carlo standard errors to the
## same entries, and compute_replication_probs_*() below return the point
## estimates. Metric names carry the consensus level as a suffix, so S = 3,
## consensus_level = 2 reproduces the paper names (P_overall_2, P_cond_O1_m1, P_gen_pos_2,
## P_c_pos_2, ...).

## Conditional probability P(A | B) from 0/1 indicators (A already A n B).
.cond_prob <- function(num, den) {
  d <- sum(den)
  if (d == 0) return(NA_real_)
  sum(num) / d
}

## Ordered list of metric "entries". Each entry is either
##   list(name, type = "uncond", ind)             -- a 0/1 indicator, OR
##   list(name, type = "cond",   num, den)        -- a conditional P(num | den).
##
## a   : list of S draw matrices [iterations x chains] (or vectors).
## generative_beta : optional draws of the generative effect (adds the generative metrics).
## eps : practical-relevance threshold.
## consensus_level   : consensus level(s); one entry block per value.
## min_corroborating   : minimum number of OTHER studies agreeing, for the conditional metrics.
.metric_entries <- function(beta, generative_beta = NULL, eps, consensus_level = 2, min_corroborating = 1) {
  beta <- lapply(beta, as.matrix)
  S <- length(beta)
  consensus_level <- as.integer(consensus_level)
  min_corroborating <- as.integer(min_corroborating)

  posI <- lapply(beta, function(x) x >  eps)
  negI <- lapply(beta, function(x) x < -eps)
  nulI <- lapply(beta, function(x) abs(x) <= eps)
  pos  <- Reduce(`+`, posI)          # count of positive studies (matrix)
  neg  <- Reduce(`+`, negI)
  nul  <- Reduce(`+`, nulI)

  E <- list()
  addu <- function(nm, ind)       E[[length(E) + 1L]] <<- list(name = nm, type = "uncond", ind = ind)
  addc <- function(nm, num, den)  E[[length(E) + 1L]] <<- list(name = nm, type = "cond", num = num, den = den)

  ## study-level agreement, one block per consensus_level value
  for (kk in consensus_level) {
    addu(sprintf("P_overall_%d", kk),  (nul >= kk) | (pos >= kk) | (neg >= kk))
    addu(sprintf("P_non_null_%d", kk), (pos >= kk) | (neg >= kk))
    addu(sprintf("P_pos_%d", kk),      pos >= kk)
    addu(sprintf("P_neg_%d", kk),      neg >= kk)
    addu(sprintf("P_null_%d", kk),     nul >= kk)
  }

  ## per-study conditional: given study i detects, does it agree in direction
  ## with at least min_corroborating other studies?
  for (i in seq_len(S)) {
    others_pos <- pos - posI[[i]]
    others_neg <- neg - negI[[i]]
    Mi <- posI[[i]] | negI[[i]]
    ev <- (posI[[i]] & (others_pos >= min_corroborating)) | (negI[[i]] & (others_neg >= min_corroborating))
    addc(sprintf("P_cond_O%d_m%d", i, min_corroborating), ev, Mi)
  }

  ## generative-consistency metrics (require the generative effect mu)
  if (!is.null(generative_beta)) {
    generative_beta <- as.matrix(generative_beta)
    mp <- generative_beta >  eps; mn <- generative_beta < -eps; m0 <- abs(generative_beta) <= eps
    addu("P_beta", abs(generative_beta) > eps)
    for (kk in consensus_level) {
      addu(sprintf("P_gen_overall_%d", kk),  (pos >= kk & mp) | (neg >= kk & mn) | (nul >= kk & m0))
      addu(sprintf("P_gen_non_null_%d", kk), (pos >= kk & mp) | (neg >= kk & mn))
      addu(sprintf("P_gen_pos_%d", kk),      (pos >= kk) & mp)
      addu(sprintf("P_gen_neg_%d", kk),      (neg >= kk) & mn)
      addu(sprintf("P_gen_null_%d", kk),     (nul >= kk) & m0)
      addc(sprintf("P_c_pos_%d", kk),  mp & (pos >= kk), mp)
      addc(sprintf("P_c_neg_%d", kk),  mn & (neg >= kk), mn)
      addc(sprintf("P_c_null_%d", kk), m0 & (nul >= kk), m0)
    }
  }
  E
}

## point estimate of every entry -> named list
.entries_points <- function(E) {
  out <- lapply(E, function(e)
    if (e$type == "uncond") mean(e$ind) else .cond_prob(e$num, e$den))
  names(out) <- vapply(E, function(e) e$name, character(1))
  out
}

#' Study-level replication probabilities (independence limit)
#'
#' Point estimates of the replication probabilities that use only the
#' study-specific effects (no generative effect). Each probability is the Monte
#' Carlo average of a 0/1 indicator over the posterior draws.
#'
#' @param beta A list of `S` posterior-draw matrices \[iterations x chains\] (or
#'   vectors) of the study-specific effects.
#' @param eps Practical-relevance threshold on the scale of the effects.
#' @param consensus_level Consensus level(s): at least `consensus_level` of the `S` studies agree. A vector
#'   is allowed (one block of metrics per value). Default 2.
#' @param min_corroborating For the per-study conditional metrics, the minimum number of *other*
#'   studies that must agree in direction. Default 1.
#'
#' @return A named list of probabilities, with paper-style names carrying the
#'   consensus level as a suffix (e.g. `P_overall_2`, `P_pos_2`,
#'   `P_cond_O1_m1`).
#'
#' @seealso [replication_ess()] for the same metrics with Monte Carlo standard
#'   errors.
#' @export
compute_replication_probs_indep <- function(beta, eps, consensus_level = 2, min_corroborating = 1) {
  .entries_points(.metric_entries(beta, generative_beta = NULL, eps = eps, consensus_level = consensus_level, min_corroborating = min_corroborating))
}

#' Full replication probabilities (hierarchical model)
#'
#' Point estimates of all replication probabilities: the study-level metrics of
#' [compute_replication_probs_indep()] plus the generative-consistency metrics
#' that compare the study-specific effects with the generative effect
#' `generative_beta`.
#'
#' @inheritParams compute_replication_probs_indep
#' @param generative_beta Posterior draws of the generative effect (`mu_beta`),
#'   paired draw-by-draw with the study effects in `beta`.
#'
#' @return A named list of probabilities (paper-style names). The generative
#'   metrics are `P_beta`, `P_gen_*_k`, `P_c_pos_k`, etc.
#' @export
compute_replication_probs_hier <- function(beta, generative_beta, eps, consensus_level = 2, min_corroborating = 1) {
  .entries_points(.metric_entries(beta, generative_beta = generative_beta, eps = eps, consensus_level = consensus_level, min_corroborating = min_corroborating))
}
