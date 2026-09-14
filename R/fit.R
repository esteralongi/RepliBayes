## fit.R -- compile/fit the Stan models and the high-level wrapper.

## Compiled-model cache (compile each Stan model once per session).
.model_cache <- new.env(parent = emptyenv())

.stan_model <- function(name) {
  if (!is.null(.model_cache[[name]])) return(.model_cache[[name]])
  f <- system.file("stan", name, package = "RepliBayes")
  if (!nzchar(f)) stop("Stan model '", name, "' not found in the installed package.")
  mod <- rstan::stan_model(file = f)
  assign(name, mod, envir = .model_cache)
  mod
}

#' Fit the hierarchical (partial-pooling) model
#'
#' @param data A data frame with columns `study`, `x`, `m`.
#' @param priors Prior hyperparameters (see [default_priors()]).
#' @param chains,iter,warmup,seed Passed to [rstan::sampling()].
#' @param adapt_delta,max_treedepth NUTS control parameters.
#' @param ... Further arguments passed to [rstan::sampling()].
#'
#' @return A [rstan::stanfit-class] object.
#' @export
fit_hierarchical <- function(data, priors = default_priors(),
                             chains = 4, iter = 3000, warmup = floor(iter / 2),
                             seed = 42, adapt_delta = 0.99, max_treedepth = 15, ...) {
  sdat <- build_stan_data(data, priors, "hierarchical")
  mod  <- .stan_model("model_hierarchical.stan")
  rstan::sampling(mod, data = sdat, chains = chains, iter = iter, warmup = warmup,
                  seed = seed,
                  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
                  ...)
}

#' Fit the independence-limit model
#'
#' Same likelihood as [fit_hierarchical()] but with the between-study
#' heterogeneities fixed at the 99.9th percentile of their priors, so the
#' studies are (almost) independent.
#'
#' @inheritParams fit_hierarchical
#' @return A [rstan::stanfit-class] object.
#' @export
fit_independence <- function(data, priors = default_priors(),
                             chains = 4, iter = 3000, warmup = floor(iter / 2),
                             seed = 42, adapt_delta = 0.99, max_treedepth = 15, ...) {
  sdat <- build_stan_data(data, priors, "independence")
  mod  <- .stan_model("model_independence.stan")
  rstan::sampling(mod, data = sdat, chains = chains, iter = iter, warmup = warmup,
                  seed = seed,
                  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
                  ...)
}

#' Predictive effect of a new study
#'
#' Draws the posterior-predictive effect of a new study from a fitted
#' hierarchical model: for each posterior draw, `mu_a + t_nu * tau_a`.
#'
#' @param fit A hierarchical [rstan::stanfit-class] object.
#' @param nu Degrees of freedom of the study-level Student-t (default 3).
#' @param seed Optional seed for the Student-t draws.
#'
#' @return A numeric vector of posterior-predictive draws of the new effect.
#' @export
predict_new_study <- function(fit, nu = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  post <- rstan::extract(fit, pars = c("mu_a", "tau_a"), permuted = TRUE)
  post$mu_a + stats::rt(length(post$mu_a), df = nu) * post$tau_a
}

## draw the [iterations x chains] predictive effect from a hierarchical fit
.predict_new_study_ic <- function(fit, nu, seed) {
  set.seed(seed)
  mt <- rstan::extract(fit, pars = c("mu_a", "tau_a"), permuted = FALSE)
  di <- dim(mt)                                             # iter x chain x 2
  mt[, , "mu_a"] + matrix(stats::rt(di[1] * di[2], df = nu), di[1], di[2]) * mt[, , "tau_a"]
}

#' Quantify replicability in one call
#'
#' Fits the hierarchical and independence-limit models and returns all
#' replication probabilities (empirical, retrospective and prospective), each
#' with its Monte Carlo standard error.
#'
#' @param data A data frame with columns `study`, `x`, `m`.
#' @param priors Prior hyperparameters (see [default_priors()]).
#' @param eps Practical-relevance threshold on the scale of the effects. If
#'   `NULL` (default) it is set to 10\% of the baseline (predictor = 0) mean
#'   outcome.
#' @param retrospective If `TRUE`, also run the retrospective analysis against a
#'   left-out study (needs at least two studies).
#' @param prospective If `TRUE`, also compute the prospective metrics for a new
#'   study drawn from the full model.
#' @param leave_out The study (value of `data$study`) to hold out in the
#'   retrospective analysis; defaults to the first study.
#' @param chains,iter,warmup,seed,adapt_delta,max_treedepth Sampler settings,
#'   passed to [fit_hierarchical()] / [fit_independence()].
#' @param verbose If `TRUE`, print progress messages.
#' @param ... Further arguments passed to [rstan::sampling()].
#'
#' @return An object of class `"RepliBayes"`: a list with the threshold `eps`,
#'   the `priors`, the fitted models, and the metric tables
#'   `metrics_hierarchical`, `metrics_independence`, `metrics_prospective` and
#'   `metrics_retrospective` (each a tibble from [replication_ess()],
#'   [held_pro_mcse()] or [held_mcse()]).
#'
#' @examples
#' \dontrun{
#' data(synthetic_data)
#' res <- fit_replicability(synthetic_data)
#' res
#' res$metrics_hierarchical
#' }
#' @export
fit_replicability <- function(data, priors = default_priors(), eps = NULL,
                              retrospective = TRUE, prospective = TRUE,
                              leave_out = NULL,
                              chains = 4, iter = 3000, warmup = floor(iter / 2),
                              seed = 42, adapt_delta = 0.99, max_treedepth = 15,
                              verbose = TRUE, ...) {
  stopifnot(all(c("study", "x", "m") %in% names(data)))
  studies <- sort(unique(data$study))
  if (length(studies) != 3L)
    stop("RepliBayes's replication metrics are defined for exactly 3 studies; ",
         "'data$study' has ", length(studies), ".")
  if (is.null(eps)) eps <- .default_eps(data)
  say <- function(...) if (verbose) message(...)

  ## --- full-data fits -------------------------------------------------------
  say("Fitting hierarchical model ...")
  fit_h <- fit_hierarchical(data, priors, chains = chains, iter = iter,
                            warmup = warmup, seed = seed, adapt_delta = adapt_delta,
                            max_treedepth = max_treedepth, ...)
  say("Fitting independence-limit model ...")
  fit_i <- fit_independence(data, priors, chains = chains, iter = iter,
                            warmup = warmup, seed = seed, adapt_delta = adapt_delta,
                            max_treedepth = max_treedepth, ...)

  ## --- empirical metrics (keep chain structure for the ESS) -----------------
  say("Computing replication probabilities ...")
  arr_h <- rstan::extract(fit_h, pars = c("a", "mu_a"), permuted = FALSE)
  report_hier <- replication_ess(arr_h[, , "a[1]"], arr_h[, , "a[2]"], arr_h[, , "a[3]"],
                                 mu = arr_h[, , "mu_a"], eps = eps)
  arr_i <- rstan::extract(fit_i, pars = "a", permuted = FALSE)
  report_indep <- replication_ess(arr_i[, , "a[1]"], arr_i[, , "a[2]"], arr_i[, , "a[3]"],
                                  mu = NULL, eps = eps)

  ## --- prospective: a new study drawn from the full model -------------------
  report_pro <- NULL
  if (prospective) {
    a_pro_ic <- .predict_new_study_ic(fit_h, nu = priors$nu, seed = seed + 1L)
    report_pro <- held_pro_mcse(a_pro_ic, eps)
  }

  ## --- retrospective: predict a left-out study ------------------------------
  report_held <- NULL; fit_held <- NULL; fit_obs <- NULL; lo <- NULL
  if (retrospective && length(studies) >= 2) {
    lo <- if (is.null(leave_out)) studies[1] else leave_out
    say("Retrospective analysis: leaving out study '", lo, "' ...")
    df_held <- data[data$study != lo, , drop = FALSE]
    df_obs  <- data[data$study == lo, , drop = FALSE]
    fit_held <- fit_hierarchical(df_held, priors, chains = chains, iter = iter,
                                 warmup = warmup, seed = seed, adapt_delta = adapt_delta,
                                 max_treedepth = max_treedepth, ...)
    fit_obs  <- fit_independence(df_obs, priors, chains = chains, iter = iter,
                                 warmup = warmup, seed = seed, adapt_delta = adapt_delta,
                                 max_treedepth = max_treedepth, ...)
    a_hat_ic <- .predict_new_study_ic(fit_held, nu = priors$nu, seed = seed)
    a_obs_ic <- rstan::extract(fit_obs, pars = "a", permuted = FALSE)[, , "a[1]"]
    report_held <- held_mcse(a_obs_ic, a_hat_ic, eps)
  }

  structure(list(
    eps       = eps,
    priors    = priors,
    studies   = studies,
    leave_out = lo,
    fit_hierarchical = fit_h,
    fit_independence = fit_i,
    fit_retro_held   = fit_held,
    fit_retro_obs    = fit_obs,
    metrics_hierarchical  = report_hier,
    metrics_independence  = report_indep,
    metrics_prospective   = report_pro,
    metrics_retrospective = report_held
  ), class = "RepliBayes")
}

#' @export
print.RepliBayes <- function(x, ...) {
  cat("<RepliBayes>\n")
  cat(sprintf("  studies: %s\n", paste(x$studies, collapse = ", ")))
  cat(sprintf("  eps    : %.4g\n", x$eps))
  key <- c("P_overall_2", "P_non_null_2", "P_beta", "P_gen_overall_2")
  m <- x$metrics_hierarchical
  m <- m[match(key, m$metric), ]
  m <- m[!is.na(m$metric), ]
  if (nrow(m)) {
    cat("  hierarchical (selected metrics):\n")
    for (i in seq_len(nrow(m)))
      cat(sprintf("    %-16s %.3f +/- %.3f\n", m$metric[i], m$p[i],
                  ifelse(is.na(m$mcse[i]), 0, m$mcse[i])))
  }
  parts <- c(
    if (!is.null(x$metrics_prospective))   "prospective" else NULL,
    if (!is.null(x$metrics_retrospective)) sprintf("retrospective (left-out: %s)", x$leave_out) else NULL
  )
  if (length(parts))
    cat(sprintf("  also computed: %s\n", paste(parts, collapse = ", ")))
  cat("  see $metrics_hierarchical, $metrics_independence, $metrics_prospective, $metrics_retrospective\n")
  invisible(x)
}
