## metrics.R -- point estimates of the replication probabilities.
## Lightweight versions (no ESS, no MCSE); used by the simulation calibration
## and by fit_replicability(). The indicator definitions are identical to those
## in replication_ess() (metrics_mcse.R), kept in sync deliberately.

## A study-specific effect is, relative to the threshold eps:
##   positive if a_s > eps, negative if a_s < -eps, null if |a_s| <= eps.
## Consensus level k = 2 requires at least two of the three studies to agree;
## k = 3 requires all three.

## Conditional probability P(event | conditioning) from 0/1 indicators.
.cond_prob <- function(event, conditioning) {
  d <- sum(conditioning)
  if (d == 0) return(NA_real_)
  sum(event & conditioning) / d
}

#' Study-level replication probabilities (independence limit)
#'
#' Point estimates of the replication probabilities that use only the
#' study-specific effects, i.e. the independence-limit metrics. Each probability
#' is the Monte Carlo average of a 0/1 indicator over the posterior draws.
#'
#' @param a1,a2,a3 Numeric vectors of posterior draws of the three
#'   study-specific effects.
#' @param eps Practical-relevance threshold on the scale of the effects.
#'
#' @return A named list of probabilities (internal metric names, e.g.
#'   `Agr_2`, `P_Rep_2`, `P_pos_2`, `Cond_S1`).
#'
#' @seealso [replication_ess()] for the same metrics with Monte Carlo standard
#'   errors and paper-facing names.
#' @export
compute_replication_probs_indep <- function(a1, a2, a3, eps) {
  pos <- (a1 >  eps) + (a2 >  eps) + (a3 >  eps)   # number of positive studies
  neg <- (a1 < -eps) + (a2 < -eps) + (a3 < -eps)   # number of negative studies
  nul <- (abs(a1) <= eps) + (abs(a2) <= eps) + (abs(a3) <= eps)
  M1  <- abs(a1) > eps; M2 <- abs(a2) > eps; M3 <- abs(a3) > eps  # study detects

  list(
    Agr_3    = mean((nul == 3) | (pos == 3) | (neg == 3)),
    Agr_2    = mean((nul >= 2) | (pos >= 2) | (neg >= 2)),
    P_Rep_3  = mean((pos == 3) | (neg == 3)),
    P_Rep_2  = mean((pos >= 2) | (neg >= 2)),
    P_pos_3  = mean(pos == 3),
    P_neg_3  = mean(neg == 3),
    P_null_3 = mean(nul == 3),
    P_pos_2  = mean(pos >= 2),
    P_neg_2  = mean(neg >= 2),
    P_null_2 = mean(nul >= 2),
    Cond_S1  = .cond_prob((a1 >  eps & (a2 >  eps | a3 >  eps)) |
                            (a1 < -eps & (a2 < -eps | a3 < -eps)), M1),
    Cond_S2  = .cond_prob((a2 >  eps & (a1 >  eps | a3 >  eps)) |
                            (a2 < -eps & (a1 < -eps | a3 < -eps)), M2),
    Cond_S3  = .cond_prob((a3 >  eps & (a1 >  eps | a2 >  eps)) |
                            (a3 < -eps & (a1 < -eps | a2 < -eps)), M3)
  )
}

#' Full replication probabilities (hierarchical model)
#'
#' Point estimates of all replication probabilities: the study-level metrics of
#' [compute_replication_probs_indep()] plus the generative-consistency metrics
#' that compare the study-specific effects with the generative effect `mu`.
#'
#' @param a1,a2,a3 Numeric vectors of posterior draws of the three
#'   study-specific effects.
#' @param mu Numeric vector of posterior draws of the generative effect
#'   (`mu_a`), paired draw-by-draw with `a1`, `a2`, `a3`.
#' @param eps Practical-relevance threshold on the scale of the effects.
#'
#' @return A named list of probabilities (internal metric names). The
#'   generative metrics are `P_beta`, `Agr_Gen_*`, `Pos_Gen_*`, `CondPos_Gen`,
#'   etc.
#' @export
compute_replication_probs_hier <- function(a1, a2, a3, mu, eps) {
  base <- compute_replication_probs_indep(a1, a2, a3, eps)

  pos <- (a1 >  eps) + (a2 >  eps) + (a3 >  eps)
  neg <- (a1 < -eps) + (a2 < -eps) + (a3 < -eps)
  nul <- (abs(a1) <= eps) + (abs(a2) <= eps) + (abs(a3) <= eps)
  mp  <- mu >  eps; mn <- mu < -eps; m0 <- abs(mu) <= eps   # generative effect sign

  gen <- list(
    P_beta             = mean(abs(mu) > eps),
    Agr_Gen_3          = mean((pos == 3 & mp) | (neg == 3 & mn) | (nul == 3 & m0)),
    Agr_Gen_2          = mean((pos >= 2 & mp) | (neg >= 2 & mn) | (nul >= 2 & m0)),
    Agr_Gen_non_null_3 = mean((pos == 3 & mp) | (neg == 3 & mn)),
    Agr_Gen_non_null_2 = mean((pos >= 2 & mp) | (neg >= 2 & mn)),
    Pos_Gen_3  = mean((pos == 3) & mp),
    Neg_Gen_3  = mean((neg == 3) & mn),
    Null_Gen_3 = mean((nul == 3) & m0),
    Pos_Gen_2  = mean((pos >= 2) & mp),
    Neg_Gen_2  = mean((neg >= 2) & mn),
    Null_Gen_2 = mean((nul >= 2) & m0),
    CondPos_Gen  = .cond_prob(pos >= 2, mp),
    CondNeg_Gen  = .cond_prob(neg >= 2, mn),
    CondNull_Gen = .cond_prob(nul >= 2, m0)
  )
  c(base, gen)
}
