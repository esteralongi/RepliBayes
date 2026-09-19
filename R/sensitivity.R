## sensitivity.R -- prior-sensitivity of the replication metrics.

## Fit one setting at the main configuration; if it does not converge, double the
## iterations (and tighten adapt_delta) up to a cap. Convergence criterion:
## no divergences, Rhat < 1.02, min ESS > 100.
.fit_until_converged <- function(model, sd_k, chains, iter_start, iter_max,
                                 seed, adapt_delta, max_treedepth) {
  iter_k <- iter_start
  ad     <- adapt_delta
  repeat {
    fit_k <- rstan::sampling(model, data = sd_k, iter = iter_k, warmup = iter_k / 2,
                             chains = chains, seed = seed, refresh = 0,
                             control = list(adapt_delta = ad, max_treedepth = max_treedepth))
    s       <- rstan::summary(fit_k)$summary
    n_div   <- rstan::get_num_divergent(fit_k)
    max_rh  <- max(s[, "Rhat"],  na.rm = TRUE)
    min_ess <- min(s[, "n_eff"], na.rm = TRUE)
    ok      <- (n_div == 0) & (max_rh < 1.02) & (min_ess > 100)
    if (ok || iter_k >= iter_max)
      return(list(fit = fit_k, ok = ok, iter = iter_k,
                  n_div = n_div, max_rhat = max_rh, min_ess = min_ess))
    iter_k <- iter_k * 2
    ad     <- max(ad, 0.999)
  }
}

#' Prior-sensitivity of the replication metrics
#'
#' Refits the hierarchical model with the data held fixed while the prior on one
#' between-study heterogeneity is shifted across percentiles of its elicited
#' prior (its whole location and scale are rescaled so that its median matches
#' the chosen percentile). Each setting is refit with adaptive convergence and
#' the replication metrics are recomputed with [replication_ess()].
#'
#' @param data A data frame with columns `study`, `x`, `m`.
#' @param priors Prior hyperparameters (see [default_priors()]).
#' @param eps Practical-relevance threshold; if `NULL`, 10\% of the baseline
#'   (predictor = 0) mean outcome.
#' @param target Which heterogeneity's prior to shift: `"tau_beta"` (effect,
#'   default), `"tau_alpha"` (intercept) or `"tau_sig"` (residual scale).
#' @param percentiles Prior percentiles at which to place the median of the
#'   shifted prior (0.5 reproduces the main analysis).
#' @param consensus_level Consensus level(s) passed to [replication_ess()].
#'   Default 2.
#' @param min_corroborating Minimum number of other studies required to
#'   corroborate the discovery of the reference study, for the conditional
#'   metrics. Default 1.
#' @param chains,iter_start,iter_max,seed,adapt_delta,max_treedepth Sampler and
#'   adaptive-convergence settings.
#' @param verbose If `TRUE`, print progress messages.
#'
#' @return A list with:
#'   \describe{
#'     \item{`metrics`}{a tibble stacking [replication_ess()] over the settings,
#'       with extra columns `pctl`, `mult`, `prior_med`, `post_med`,
#'       `iter_used`, `converged`, `max_rhat`, `min_ess`.}
#'     \item{`grid`}{one row per setting with the prior grid and diagnostics.}
#'     \item{`p_grid`}{a metric-by-percentile matrix of the point estimates for
#'       every metric returned.}
#'   }
#'
#' @examples
#' \dontrun{
#' data(synthetic_data)
#' s <- sensitivity_prior(synthetic_data, target = "tau_beta")
#' round(s$p_grid, 3)
#' }
#' @export
sensitivity_prior <- function(data, priors = default_priors(), eps = NULL,
                              target = c("tau_beta", "tau_alpha", "tau_sig"),
                              percentiles = c(0.10, 0.25, 0.50, 0.75, 0.90),
                              consensus_level = 2, min_corroborating = 1,
                              chains = 4, iter_start = 3000, iter_max = 24000,
                              seed = 42, adapt_delta = 0.99, max_treedepth = 15,
                              verbose = TRUE) {
  stopifnot(all(c("study", "x", "m") %in% names(data)))
  target <- match.arg(target)
  if (is.null(eps)) eps <- .default_eps(data)
  nu  <- priors$nu
  S   <- length(unique(data$study))
  say <- function(...) if (verbose) message(...)

  ## prior (mu, scale), the Stan data fields to overwrite, and the matching
  ## posterior tau parameter, for the chosen heterogeneity
  tf <- switch(target,
    tau_beta     = list(mu = priors$mu_tau_beta,     scale = priors$scale_tau_beta,
                     fmu = "prior_mu_tau_beta",     fsc = "prior_scale_tau_beta",     post = "tau_beta"),
    tau_alpha = list(mu = priors$mu_tau_alpha, scale = priors$scale_tau_alpha,
                     fmu = "prior_mu_tau_alpha",  fsc = "prior_scale_tau_alpha",  post = "tau_alpha"),
    tau_sig   = list(mu = priors$mu_tau_sig,   scale = priors$scale_tau_sig,
                     fmu = "prior_mu_tau_sig_m", fsc = "prior_scale_tau_sig_m", post = "tau_sig_m")
  )

  ## grid of prior medians = the requested percentiles of the elicited prior
  target_meds <- qt_trunc_scaled_vec(percentiles, nu, tf$mu, tf$scale)
  med0        <- qt_trunc_scaled_vec(0.5, nu, tf$mu, tf$scale)
  mults       <- target_meds / med0
  grid <- data.frame(pctl = percentiles, mult = mults,
                     mu = mults * tf$mu, scale = mults * tf$scale)
  grid$prior_med <- vapply(seq_len(nrow(grid)),
                           function(ii) qt_trunc_scaled_vec(0.5, nu, grid$mu[ii], grid$scale[ii]),
                           numeric(1))

  sd0   <- build_stan_data(data, priors, "hierarchical")
  model <- .stan_model("model_hierarchical.stan")

  rows <- vector("list", nrow(grid))
  meta <- vector("list", nrow(grid))
  for (gi in seq_len(nrow(grid))) {
    say(sprintf("Setting %d/%d: prior median(%s) = %.4g",
                gi, nrow(grid), target, grid$prior_med[gi]))
    sd_k <- sd0
    sd_k[[tf$fmu]] <- grid$mu[gi]
    sd_k[[tf$fsc]] <- grid$scale[gi]

    res <- .fit_until_converged(model, sd_k, chains = chains, iter_start = iter_start,
                                iter_max = iter_max, seed = seed,
                                adapt_delta = adapt_delta, max_treedepth = max_treedepth)
    arr   <- rstan::extract(res$fit, pars = c("beta", "mu_beta", tf$post), permuted = FALSE)
    a_lst <- lapply(seq_len(S), function(s) arr[, , sprintf("beta[%d]", s)])
    re  <- replication_ess(a_lst, generative_beta = arr[, , "mu_beta"], eps = eps,
                           consensus_level = consensus_level, min_corroborating = min_corroborating)
    re$pctl      <- grid$pctl[gi]
    re$mult      <- grid$mult[gi]
    re$prior_med <- grid$prior_med[gi]
    re$post_med  <- stats::median(as.vector(arr[, , tf$post]))
    re$iter_used <- res$iter
    re$converged <- as.numeric(res$ok)
    re$max_rhat  <- res$max_rhat
    re$min_ess   <- res$min_ess
    rows[[gi]] <- re

    meta[[gi]] <- data.frame(pctl = grid$pctl[gi], prior_med = grid$prior_med[gi],
                             post_med = re$post_med[1], iter_used = res$iter,
                             converged = as.numeric(res$ok),
                             max_rhat = res$max_rhat, min_ess = res$min_ess)
  }
  sens_long <- do.call(rbind, rows)
  meta_df   <- do.call(rbind, meta)
  if (any(meta_df$converged == 0))
    warning("Some settings did not converge at iter_max; raise iter_max or adapt_delta.")

  ## wide metric-by-percentile matrix (every metric returned)
  metric_order <- unique(sens_long$metric)
  sub    <- sens_long[, c("metric", "prior_med", "p")]
  p_grid <- stats::xtabs(p ~ metric + prior_med, data = sub)[metric_order, , drop = FALSE]

  list(metrics = sens_long, grid = meta_df, p_grid = p_grid)
}
