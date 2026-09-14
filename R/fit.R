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
#' hierarchical model: for each posterior draw, `mu_beta + t_nu * tau_beta`.
#'
#' @param fit A hierarchical [rstan::stanfit-class] object.
#' @param nu Degrees of freedom of the study-level Student-t (default 3).
#' @param seed Optional seed for the Student-t draws.
#'
#' @return A numeric vector of posterior-predictive draws of the new effect.
#' @export
predict_new_study <- function(fit, nu = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  post <- rstan::extract(fit, pars = c("mu_beta", "tau_beta"), permuted = TRUE)
  post$mu_beta + stats::rt(length(post$mu_beta), df = nu) * post$tau_beta
}

## draw the [iterations x chains] predictive effect from a hierarchical fit
.predict_new_study_ic <- function(fit, nu, seed) {
  set.seed(seed)
  mt <- rstan::extract(fit, pars = c("mu_beta", "tau_beta"), permuted = FALSE)
  di <- dim(mt)                                             # iter x chain x 2
  mt[, , "mu_beta"] + matrix(stats::rt(di[1] * di[2], df = nu), di[1], di[2]) * mt[, , "tau_beta"]
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
#' @param k Consensus level(s) passed to [replication_ess()]: at least `k` of the
#'   `S` studies agree. A vector is allowed. Default 2.
#' @param m Minimum number of *other* studies agreeing, for the per-study
#'   conditional metrics. Default 1.
#' @param retrospective If `TRUE`, also run the retrospective analysis: predict a
#'   target study from a body of evidence.
#' @param prospective If `TRUE`, also compute the prospective metrics for a new
#'   study drawn from a body of evidence.
#' @param retro_target The study (value of `data$study`) to leave out and predict
#'   in the retrospective analysis; defaults to the first study.
#' @param retro_evidence The studies forming the body of evidence for the
#'   retrospective prediction (must exclude `retro_target`); defaults to all the
#'   other studies. Give a subset of size 1..S-1 to use a smaller evidence base.
#' @param pro_evidence The studies forming the body of evidence for the
#'   prospective new study; defaults to all studies. Give a subset of size 1..S.
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
                              k = 2, m = 1,
                              retrospective = TRUE, prospective = TRUE,
                              retro_target = NULL, retro_evidence = NULL,
                              pro_evidence = NULL,
                              chains = 4, iter = 3000, warmup = floor(iter / 2),
                              seed = 42, adapt_delta = 0.99, max_treedepth = 15,
                              verbose = TRUE, ...) {
  stopifnot(all(c("study", "x", "m") %in% names(data)))
  studies <- sort(unique(data$study))
  S <- length(studies)
  if (S < 2L) stop("Need at least 2 studies; 'data$study' has ", S, ".")
  if (is.null(eps)) eps <- .default_eps(data)
  say <- function(...) if (verbose) message(...)

  ## sampler settings reused for every fit
  fit_h_fun <- function(df) fit_hierarchical(df, priors, chains = chains, iter = iter,
                                             warmup = warmup, seed = seed, adapt_delta = adapt_delta,
                                             max_treedepth = max_treedepth, ...)
  fit_i_fun <- function(df) fit_independence(df, priors, chains = chains, iter = iter,
                                             warmup = warmup, seed = seed, adapt_delta = adapt_delta,
                                             max_treedepth = max_treedepth, ...)
  ## list of per-study effect draw matrices [iter x chain] from a fit of n studies
  a_list_from_fit <- function(fit, n) {
    arr <- rstan::extract(fit, pars = "beta", permuted = FALSE)
    lapply(seq_len(n), function(s) arr[, , sprintf("beta[%d]", s)])
  }

  ## --- full-data fits -------------------------------------------------------
  say("Fitting hierarchical model ...")
  fit_h <- fit_h_fun(data)
  say("Fitting independence-limit model ...")
  fit_i <- fit_i_fun(data)

  ## --- empirical metrics (keep chain structure for the ESS) -----------------
  say("Computing replication probabilities ...")
  mu_h        <- rstan::extract(fit_h, pars = "mu_beta", permuted = FALSE)[, , "mu_beta"]
  report_hier  <- replication_ess(a_list_from_fit(fit_h, S), mu = mu_h, eps = eps, k = k, m = m)
  report_indep <- replication_ess(a_list_from_fit(fit_i, S), mu = NULL,  eps = eps, k = k, m = m)

  ## --- prospective: a new study drawn from a body of evidence ---------------
  report_pro <- NULL; pro_ev <- NULL
  if (prospective) {
    pro_ev <- if (is.null(pro_evidence)) studies else pro_evidence
    stopifnot(all(pro_ev %in% studies), length(pro_ev) >= 1L)
    fit_pro <- if (setequal(pro_ev, studies)) fit_h else
      fit_h_fun(data[data$study %in% pro_ev, , drop = FALSE])
    a_pro_ic <- .predict_new_study_ic(fit_pro, nu = priors$nu, seed = seed + 1L)
    report_pro <- held_pro_mcse(a_pro_ic, eps)
  }

  ## --- retrospective: predict a target study from a body of evidence --------
  report_held <- NULL; fit_ev <- NULL; fit_obs <- NULL; rt <- NULL; re_ev <- NULL
  if (retrospective) {
    rt    <- if (is.null(retro_target))   studies[1]            else retro_target
    re_ev <- if (is.null(retro_evidence)) setdiff(studies, rt)  else retro_evidence
    stopifnot(rt %in% studies, all(re_ev %in% studies),
              !(rt %in% re_ev), length(re_ev) >= 1L)
    say("Retrospective: target '", rt, "', evidence {", paste(re_ev, collapse = ", "), "} ...")
    fit_ev  <- fit_h_fun(data[data$study %in% re_ev, , drop = FALSE])
    fit_obs <- fit_i_fun(data[data$study == rt, , drop = FALSE])
    a_hat_ic <- .predict_new_study_ic(fit_ev, nu = priors$nu, seed = seed)
    a_obs_ic <- rstan::extract(fit_obs, pars = "beta", permuted = FALSE)[, , "beta[1]"]
    report_held <- held_mcse(a_obs_ic, a_hat_ic, eps)
  }

  structure(list(
    eps       = eps,
    priors    = priors,
    studies   = studies,
    S         = S,
    k         = k,
    m         = m,
    retro_target   = rt,
    retro_evidence = re_ev,
    pro_evidence   = pro_ev,
    fit_hierarchical   = fit_h,
    fit_independence   = fit_i,
    fit_retro_evidence = fit_ev,
    fit_retro_obs      = fit_obs,
    metrics_hierarchical  = report_hier,
    metrics_independence  = report_indep,
    metrics_prospective   = report_pro,
    metrics_retrospective = report_held
  ), class = "RepliBayes")
}

#' @export
print.RepliBayes <- function(x, ...) {
  cat("<RepliBayes>\n")
  cat(sprintf("  studies: %s  (S = %d)\n", paste(x$studies, collapse = ", "), x$S))
  cat(sprintf("  eps    : %.4g   k = %s   m = %d\n", x$eps,
              paste(x$k, collapse = ","), as.integer(x$m)))
  k1  <- as.integer(x$k)[1]
  key <- c(sprintf("P_overall_%d", k1), sprintf("P_non_null_%d", k1),
           "P_beta", sprintf("P_gen_overall_%d", k1))
  mm <- x$metrics_hierarchical
  mm <- mm[match(key, mm$metric), ]
  mm <- mm[!is.na(mm$metric), ]
  if (nrow(mm)) {
    cat("  hierarchical (selected metrics):\n")
    for (i in seq_len(nrow(mm)))
      cat(sprintf("    %-18s %.3f +/- %.3f\n", mm$metric[i], mm$p[i],
                  ifelse(is.na(mm$mcse[i]), 0, mm$mcse[i])))
  }
  parts <- c(
    if (!is.null(x$metrics_prospective))
      sprintf("prospective (evidence: %s)", paste(x$pro_evidence, collapse = ", ")) else NULL,
    if (!is.null(x$metrics_retrospective))
      sprintf("retrospective (target: %s; evidence: %s)",
              x$retro_target, paste(x$retro_evidence, collapse = ", ")) else NULL
  )
  if (length(parts))
    cat(sprintf("  also computed: %s\n", paste(parts, collapse = "; ")))
  cat("  see $metrics_hierarchical, $metrics_independence, $metrics_prospective, $metrics_retrospective\n")
  invisible(x)
}
