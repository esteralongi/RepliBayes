## simulation.R -- calibrate the replication metrics by simulation.

## genotype (0/1/2) frequencies of one study
.geno_probs <- function(x) as.numeric(prop.table(table(factor(x, levels = c(0, 1, 2)))))

## simulate one S-study dataset from fixed ground-truth parameters
.simulate_dataset <- function(N_s, geno_probs, a_s, b_s, sig_s) {
  S <- length(N_s)
  df_list <- lapply(seq_len(S), function(s) {
    X <- sample(c(0, 1, 2), N_s[s], replace = TRUE, prob = geno_probs[[s]])
    M <- stats::rnorm(N_s[s], mean = a_s[s] + b_s[s] * X, sd = sig_s[s])
    data.frame(study = s, x = X, m = M)
  })
  do.call(rbind, df_list)
}

#' Calibrate the metrics by simulation, varying the generative model
#'
#' Simulation study in which one component of the generative model is varied
#' while the others are held at their empirical values. For each level, `R`
#' datasets of `S` studies are simulated from fixed ground-truth parameters,
#' each is refit with the hierarchical model, the replication metrics are
#' computed on the posterior, and they are summarised across replicates with a
#' Monte Carlo standard error, \eqn{\mathrm{MCSE}(\bar p) = \mathrm{sd}(\hat
#' p^{(r)})/\sqrt{R_{\mathrm{used}}}}.
#'
#' The empirical ground truth (study sample sizes, genotype frequencies, and the
#' generative locations `mu_beta`, `a`, `mu_alpha`, `alpha`, `mu_sig_m`, `sigma_M`)
#' is taken from a fitted hierarchical model; the chosen heterogeneity is then
#' set to the requested percentiles of its prior.
#'
#' @param fit A fitted hierarchical model: a [rstan::stanfit-class] object or a
#'   `RepliBayes` object (from [fit_replicability()]).
#' @param data A data frame with columns `study`, `x`, `m`, giving the study
#'   sample sizes and genotype frequencies to reuse in the simulation.
#' @param priors Prior hyperparameters used both as the ground-truth prior of
#'   the varied component and as the prior of the fitted model
#'   (see [default_priors()]).
#' @param eps Practical-relevance threshold; if `NULL`, taken from the `RepliBayes`
#'   object or computed as 10\% of the baseline mean outcome.
#' @param vary Which generative component to vary across levels: `"tau_beta"`
#'   (effect heterogeneity, default), `"tau_alpha"` (intercept heterogeneity) or
#'   `"tau_sig"` (residual-scale heterogeneity).
#' @param levels A named numeric vector of prior percentiles defining the
#'   scenarios (default low/medium/high = 0.10/0.50/0.90).
#' @param consensus_level Consensus level(s) for the metrics. Default 2.
#' @param min_corroborating Minimum number of other studies required to
#'   corroborate the discovery of the reference study, for the conditional
#'   metrics. Default 1.
#' @param R Number of simulated datasets per level.
#' @param seed Base seed; replicate `r` uses `seed + r`.
#' @param chains,iter,warmup,adapt_delta,max_treedepth Sampler settings for each
#'   refit.
#' @param verbose If `TRUE`, print progress messages.
#'
#' @return A tibble with one row per (level, metric): columns `level`, `pctl`,
#'   `metric`, `value` (mean across converged replicates), `mcse`, and `n_used`
#'   (number of converged replicates). Metric names are the paper-style names of
#'   [compute_replication_probs_hier()].
#'
#' @details Only replicates that meet the convergence criterion (no divergences,
#'   Rhat < 1.02, min ESS > 100) contribute to the summaries. The predictor `x`
#'   is treated as a genotype dosage in `{0, 1, 2}`.
#'
#' @examples
#' \dontrun{
#' data(synthetic_data)
#' res <- fit_replicability(synthetic_data)
#' sim <- simulate_replicability(res, synthetic_data, vary = "tau_beta", R = 30)
#' subset(sim, metric == "P_beta")
#' }
#' @export
simulate_replicability <- function(fit, data, priors = default_priors(), eps = NULL,
                                   vary = c("tau_beta", "tau_alpha", "tau_sig"),
                                   levels = c(low = 0.10, medium = 0.50, high = 0.90),
                                   consensus_level = 2, min_corroborating = 1,
                                   R = 30, seed = 1000,
                                   chains = 4, iter = 3000, warmup = floor(iter / 2),
                                   adapt_delta = 0.99, max_treedepth = 15,
                                   verbose = TRUE) {
  vary <- match.arg(vary)
  if (inherits(fit, "RepliBayes")) {
    if (is.null(eps)) eps <- fit$eps
    fit <- fit$fit_hierarchical
  }
  stopifnot(all(c("study", "x", "m") %in% names(data)))
  if (is.null(eps)) eps <- .default_eps(data)
  nu  <- priors$nu
  say <- function(...) if (verbose) message(...)

  ## study structure from the data
  groups <- sort(unique(data$study)); S <- length(groups)
  N_s    <- vapply(groups, function(g) sum(data$study == g), integer(1))
  geno   <- lapply(groups, function(g) .geno_probs(data$x[data$study == g]))

  ## empirical ground-truth locations from the fitted model
  post <- rstan::extract(fit)
  mu_b_real     <- mean(post$mu_beta);      b_s_real   <- colMeans(post$beta)
  mu_alpha_real <- mean(post$mu_alpha);  a_s_real   <- colMeans(post$alpha)
  mu_sigma_real <- mean(post$mu_sig_m);  sig_s_real <- colMeans(post$sigma_M)

  ## study positions (fixed) and the varied component's prior
  placement <- seq(0.10, 0.90, length.out = S)
  z_beta    <- stats::qt(placement, df = nu)   # signed effect positions
  u_trunc   <- placement                        # positions for truncated positives
  tp <- switch(vary,
    tau_beta     = list(mu = priors$mu_tau_beta,     scale = priors$scale_tau_beta),
    tau_alpha = list(mu = priors$mu_tau_alpha, scale = priors$scale_tau_alpha),
    tau_sig   = list(mu = priors$mu_tau_sig,   scale = priors$scale_tau_sig)
  )

  model <- .stan_model("model_hierarchical.stan")
  out   <- vector("list", length(levels))

  for (li in seq_along(levels)) {
    lvl_name <- names(levels)[li]; p <- levels[[li]]
    tau_val  <- qt_trunc_scaled_vec(p, nu, tp$mu, tp$scale)
    say(sprintf("Level '%s' (%s at percentile %.2f -> %.4g); %d replicates",
                lvl_name, vary, p, tau_val, R))

    ## ground-truth study parameters for this level
    b_s <- b_s_real; a_s <- a_s_real; sig_s <- sig_s_real
    if (vary == "tau_beta")     b_s   <- mu_b_real + tau_val * z_beta
    if (vary == "tau_alpha") a_s   <- qt_trunc_scaled_vec(u_trunc, nu, mu_alpha_real, tau_val)
    if (vary == "tau_sig")   sig_s <- qt_trunc_scaled_vec(u_trunc, nu, mu_sigma_real, tau_val)

    metrics_list <- vector("list", R); diag_ok <- logical(R); metric_names <- NULL
    for (r in seq_len(R)) {
      set.seed(seed + r)
      df   <- .simulate_dataset(N_s, geno, a_s, b_s, sig_s)
      sdat <- build_stan_data(df, priors, "hierarchical")
      fitr <- rstan::sampling(model, data = sdat, iter = iter, warmup = warmup,
                              chains = chains, seed = seed + r, refresh = 0,
                              control = list(adapt_delta = adapt_delta,
                                             max_treedepth = max_treedepth))
      s   <- rstan::summary(fitr)$summary
      ok  <- rstan::get_num_divergent(fitr) == 0 &&
             max(s[, "Rhat"], na.rm = TRUE) < 1.02 &&
             min(s[, "n_eff"], na.rm = TRUE) > 100
      dr    <- rstan::extract(fitr)
      a_lst <- lapply(seq_len(S), function(s) dr$beta[, s])
      pr    <- compute_replication_probs_hier(a_lst, dr$mu_beta, eps,
                                              consensus_level = consensus_level,
                                              min_corroborating = min_corroborating)
      if (is.null(metric_names)) metric_names <- names(pr)
      metrics_list[[r]] <- as.numeric(pr); diag_ok[r] <- ok
      rm(fitr, dr)
    }

    M      <- do.call(rbind, metrics_list); colnames(M) <- metric_names
    M_ok   <- M[diag_ok, , drop = FALSE]
    n_used <- nrow(M_ok)
    value  <- colMeans(M_ok, na.rm = TRUE)
    mcse   <- apply(M_ok, 2, stats::sd, na.rm = TRUE) / sqrt(n_used)

    out[[li]] <- tibble::tibble(level = lvl_name, pctl = p, metric = metric_names,
                                value = as.numeric(value), mcse = as.numeric(mcse),
                                n_used = n_used)
  }

  res <- do.call(rbind, out)
  res$level <- factor(res$level, levels = names(levels))
  res
}
