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
                             chains = chains, seed = seed,
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
#' @param target Which heterogeneity's prior to shift: `"tau_a"` (effect,
#'   default), `"tau_alpha"` (intercept) or `"tau_sig"` (residual scale).
#' @param percentiles Prior percentiles at which to place the median of the
#'   shifted prior (0.5 reproduces the main analysis).
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
#'       the k = 2 metrics reported in the paper.}
#'   }
#'
#' @examples
#' \dontrun{
#' data(synthetic_data)
#' s <- sensitivity_prior(synthetic_data, target = "tau_a")
#' round(s$p_grid, 3)
#' }
#' @export
sensitivity_prior <- function(data, priors = default_priors(), eps = NULL,
                              target = c("tau_a", "tau_alpha", "tau_sig"),
                              percentiles = c(0.10, 0.25, 0.50, 0.75, 0.90),
                              chains = 4, iter_start = 3000, iter_max = 24000,
                              seed = 42, adapt_delta = 0.99, max_treedepth = 15,
                              verbose = TRUE) {
  stopifnot(all(c("study", "x", "m") %in% names(data)))
  if (length(unique(data$study)) != 3L)
    stop("RepliBayes's replication metrics are defined for exactly 3 studies.")
  target <- match.arg(target)
  if (is.null(eps)) eps <- .default_eps(data)
  nu  <- priors$nu
  say <- function(...) if (verbose) message(...)

  ## prior (mu, scale), the Stan data fields to overwrite, and the matching
  ## posterior tau parameter, for the chosen heterogeneity
  tf <- switch(target,
    tau_a     = list(mu = priors$mu_tau_a,     scale = priors$scale_tau_a,
                     fmu = "prior_mu_tau_a",     fsc = "prior_scale_tau_a",     post = "tau_a"),
    tau_alpha = list(mu = priors$mu_tau_alpha, scale = priors$scale_tau_alpha,
                     fmu = "prior_mu_tau_int1",  fsc = "prior_scale_tau_int1",  post = "tau_intercept1"),
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
                           function(k) qt_trunc_scaled_vec(0.5, nu, grid$mu[k], grid$scale[k]),
                           numeric(1))

  sd0   <- build_stan_data(data, priors, "hierarchical")
  model <- .stan_model("model_hierarchical.stan")

  rows <- vector("list", nrow(grid))
  meta <- vector("list", nrow(grid))
  for (k in seq_len(nrow(grid))) {
    say(sprintf("Setting %d/%d: prior median(%s) = %.4g",
                k, nrow(grid), target, grid$prior_med[k]))
    sd_k <- sd0
    sd_k[[tf$fmu]] <- grid$mu[k]
    sd_k[[tf$fsc]] <- grid$scale[k]

    res <- .fit_until_converged(model, sd_k, chains = chains, iter_start = iter_start,
                                iter_max = iter_max, seed = seed,
                                adapt_delta = adapt_delta, max_treedepth = max_treedepth)
    arr <- rstan::extract(res$fit, pars = c("a", "mu_a", tf$post), permuted = FALSE)
    re  <- replication_ess(arr[, , "a[1]"], arr[, , "a[2]"], arr[, , "a[3]"],
                           mu = arr[, , "mu_a"], eps = eps)
    re$pctl      <- grid$pctl[k]
    re$mult      <- grid$mult[k]
    re$prior_med <- grid$prior_med[k]
    re$post_med  <- stats::median(as.vector(arr[, , tf$post]))
    re$iter_used <- res$iter
    re$converged <- as.numeric(res$ok)
    re$max_rhat  <- res$max_rhat
    re$min_ess   <- res$min_ess
    rows[[k]] <- re

    meta[[k]] <- data.frame(pctl = grid$pctl[k], prior_med = grid$prior_med[k],
                            post_med = re$post_med[1], iter_used = res$iter,
                            converged = as.numeric(res$ok),
                            max_rhat = res$max_rhat, min_ess = res$min_ess)
  }
  sens_long <- do.call(rbind, rows)
  meta_df   <- do.call(rbind, meta)
  if (any(meta_df$converged == 0))
    warning("Some settings did not converge at iter_max; raise iter_max or adapt_delta.")

  ## wide metric-by-percentile matrix for the k = 2 metrics
  metrics_k2 <- c("P_overall_2", "P_non_null_2", "P_pos_2", "P_null_2", "P_neg_2",
                  "P_cond_O1_m1", "P_cond_O2_m1", "P_cond_O3_m1", "P_beta",
                  "P_gen_overall_2", "P_gen_non_null_2", "P_gen_pos_2", "P_gen_null_2",
                  "P_gen_neg_2", "P_c_pos_2", "P_c_null_2", "P_c_neg_2")
  sub <- sens_long[sens_long$metric %in% metrics_k2, c("metric", "prior_med", "p")]
  p_grid <- stats::xtabs(p ~ metric + prior_med, data = sub)[metrics_k2, , drop = FALSE]

  list(metrics = sens_long, grid = meta_df, p_grid = p_grid)
}
