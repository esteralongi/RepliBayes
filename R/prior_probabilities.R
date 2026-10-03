## prior_probabilities.R -- replication metrics evaluated under the prior (reference values).
##

#' Replication probabilities under the prior
#'
#' Evaluates the replication metrics on i.i.d. draws from the prior alone, as a
#' reference for the posterior estimates returned by [fit_replicability()]. The
#' metric definitions are identical to the posterior case (the same functions,
#' [replication_ess()], [predictive_retro()] and [predictive_pro()], are used);
#' only the draws differ.
#'
#' For each of `n_draws` i.i.d. replicates the generative mean is drawn from its
#' prior, \eqn{\mu_\beta \sim \mathrm{Student\text{-}t}(\nu, \mathtt{mu\_beta},
#' \mathtt{scale\_beta})}, and:
#' \itemize{
#'   \item \strong{hierarchical}: the between-study heterogeneity is random,
#'     \eqn{\tau_\beta \sim \mathrm{Student\text{-}t}^{+}(\nu,
#'     \mathtt{mu\_tau\_beta}, \mathtt{scale\_tau\_beta})}, and the `S` study
#'     effects are \eqn{\mu_\beta + \tau_\beta\, t_\nu};
#'   \item \strong{independence limit}: `tau_beta` is fixed at the `q_indep`
#'     quantile of its prior, so the `S` study effects are (almost) independent;
#'   \item \strong{retrospective}: a predicted effect is compared with an
#'     observed left-out effect built from \emph{independent} `mu_beta` draws,
#'     mirroring the two separate fits of [fit_replicability()];
#'   \item \strong{prospective}: a single new study effect from the hierarchical
#'     prior predictive.
#' }
#'
#' @param priors Prior hyperparameters (see [default_priors()]). Only the effect
#'   component is used: `nu`, `mu_beta`, `scale_beta`, `mu_tau_beta`,
#'   `scale_tau_beta`.
#' @param eps Practical-relevance threshold on the scale of the effects. Must be
#'   supplied (use the same threshold as the posterior analysis).
#' @param consensus_level Consensus level(s) passed to [replication_ess()].
#'   Default 2.
#' @param min_corroborating Minimum number of other studies required to
#'   corroborate the discovery of the reference study, for the conditional
#'   metrics. Default 1.
#' @param S Number of studies in the super-population; set it to the number of
#'   studies in your data. Default 3.
#' @param n_draws Number of i.i.d. prior replicates. Default 50000.
#' @param q_indep Quantile of the `tau_beta` prior at which the heterogeneity is
#'   fixed in the independence limit. Default 0.999.
#' @param seed Seed for the prior draws.
#'
#' @return An object of class `"RepliBayesPrior"`: a list with the threshold
#'   `eps`, the `priors`, `S`, `n_draws`, and the prior metric tables
#'   `metrics_hierarchical`, `metrics_independence`, `metrics_retrospective` and
#'   `metrics_prospective`. Each is a tibble with the same columns as the
#'   corresponding table of [fit_replicability()] (so they can be joined by
#'   `metric`). For these i.i.d. prior draws the effective sample size `ess`
#'   equals `n_draws` and `mcse` is the (negligible) Monte Carlo error of the
#'   prior value.
#'
#' @seealso [fit_replicability()] for the posterior estimates.
#' @examples
#' \dontrun{
#' pri <- prior_replicability(default_priors(), eps = 0.1, S = 3)
#' pri$metrics_hierarchical
#' }
#' @export
prior_replicability <- function(priors = default_priors(), eps,
                                consensus_level = 2, min_corroborating = 1,
                                S = 3L, n_draws = 50000L, q_indep = 0.999,
                                seed = 42) {
  if (missing(eps))
    stop("'eps' must be supplied (use the same threshold as the posterior analysis).")
  S <- as.integer(S); n_draws <- as.integer(n_draws)
  if (S < 2L) stop("Need at least 2 studies; S = ", S, ".")
  nu <- priors$nu
  set.seed(seed)

  draw_mu   <- function(n) priors$mu_beta + priors$scale_beta * stats::rt(n, df = nu)
  draw_tau  <- function(n) rtrunc_t_pos(n, nu, priors$mu_tau_beta, priors$scale_tau_beta)
  tau_fixed <- qt_trunc_scaled_vec(q_indep, nu, priors$mu_tau_beta, priors$scale_tau_beta)
  col       <- function(x) matrix(x, ncol = 1L)   # one "chain" of i.i.d. draws

  ## hierarchical: shared mu_beta and random tau_beta across the S studies
  mu_h   <- draw_mu(n_draws)
  tau_h  <- draw_tau(n_draws)
  beta_h <- lapply(seq_len(S), function(s) col(mu_h + tau_h * stats::rt(n_draws, df = nu)))
  metrics_hierarchical <- replication_ess(beta_h, generative_beta = col(mu_h), eps = eps,
                                          consensus_level = consensus_level,
                                          min_corroborating = min_corroborating)

  ## independence limit: shared mu_beta, tau_beta fixed at the q_indep quantile
  mu_i   <- draw_mu(n_draws)
  beta_i <- lapply(seq_len(S), function(s) col(mu_i + tau_fixed * stats::rt(n_draws, df = nu)))
  metrics_independence <- replication_ess(beta_i, generative_beta = NULL, eps = eps,
                                          consensus_level = consensus_level,
                                          min_corroborating = min_corroborating)

  ## retrospective: predicted (random tau) vs observed independent left-out (fixed tau),
  mu_pred   <- draw_mu(n_draws); tau_pred <- draw_tau(n_draws)
  beta_pred <- col(mu_pred + tau_pred * stats::rt(n_draws, df = nu))
  mu_obs    <- draw_mu(n_draws)
  beta_obs  <- col(mu_obs + tau_fixed * stats::rt(n_draws, df = nu))
  metrics_retrospective <- predictive_retro(beta_obs, beta_pred, eps = eps)

  ## prospective: a new study from the hierarchical prior predictive (random tau)
  mu_p     <- draw_mu(n_draws); tau_p <- draw_tau(n_draws)
  beta_new <- col(mu_p + tau_p * stats::rt(n_draws, df = nu))
  metrics_prospective <- predictive_pro(beta_new, eps = eps)

  structure(list(
    eps               = eps,
    priors            = priors,
    S                 = S,
    n_draws           = n_draws,
    consensus_level   = consensus_level,
    min_corroborating = min_corroborating,
    metrics_hierarchical  = metrics_hierarchical,
    metrics_independence  = metrics_independence,
    metrics_retrospective = metrics_retrospective,
    metrics_prospective   = metrics_prospective
  ), class = "RepliBayesPrior")
}

#' @export
print.RepliBayesPrior <- function(x, ...) {
  cat("<RepliBayesPrior>\n")
  cat(sprintf("  prior reference values from %d i.i.d. draws;  S = %d,  eps = %.4g\n",
              x$n_draws, x$S, x$eps))
  cat(sprintf("  consensus_level = %s   min_corroborating = %d\n",
              paste(x$consensus_level, collapse = ","), as.integer(x$min_corroborating)))
  cat("  see $metrics_hierarchical, $metrics_independence, $metrics_retrospective, $metrics_prospective\n")
  invisible(x)
}
