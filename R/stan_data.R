## stan_data.R -- assemble the data list passed to the Stan models.

#' Assemble the Stan data list
#'
#' Builds the `data` list for either the hierarchical model or the
#' independence-limit model from a study-level data frame and a set of prior
#' hyperparameters.
#'
#' @param data A data frame with (at least) the columns `study` (group label),
#'   `x` (predictor, e.g. genotype dosage) and `y` (outcome).
#' @param priors A list of prior hyperparameters, as returned by
#'   [default_priors()].
#' @param model Either `"hierarchical"` (partial pooling, the heterogeneities
#'   `tau` are estimated) or `"independence"` (the `tau` are fixed at the
#'   99.9th percentile of their priors, giving the independence limit).
#'
#' @return A named list suitable as the `data` argument of [rstan::sampling()].
#'   For the independence model the entries `prior_tau_beta_fixed`,
#'   `prior_tau_alpha_fixed` and `prior_tau_sig_fixed` are added.
#' @export
build_stan_data <- function(data, priors = default_priors(),
                            model = c("hierarchical", "independence")) {
  model <- match.arg(model)
  stopifnot(all(c("study", "x", "y") %in% names(data)))

  study <- as.integer(as.factor(data$study))
  base <- list(
    N = nrow(data),
    S = length(unique(study)),
    study = study,
    X = as.numeric(data$x),
    Y = as.numeric(data$y),
    nu = priors$nu,
    prior_mu_beta     = priors$mu_beta,     prior_scale_beta     = priors$scale_beta,
    prior_mu_alpha = priors$mu_alpha, prior_scale_alpha = priors$scale_alpha,
    prior_mu_sig = priors$mu_sig,   prior_scale_sig = priors$scale_sig
  )

  if (model == "hierarchical") {
    c(base, list(
      prior_mu_tau_beta     = priors$mu_tau_beta,     prior_scale_tau_beta     = priors$scale_tau_beta,
      prior_mu_tau_alpha  = priors$mu_tau_alpha, prior_scale_tau_alpha  = priors$scale_tau_alpha,
      prior_mu_tau_sig = priors$mu_tau_sig,   prior_scale_tau_sig = priors$scale_tau_sig
    ))
  } else {
    ## independence-limit: tau fixed at the 99.9th percentile of each prior
    c(base, list(
      prior_tau_beta_fixed     = qt_trunc_scaled_vec(0.999, priors$nu, priors$mu_tau_beta,     priors$scale_tau_beta),
      prior_tau_alpha_fixed = qt_trunc_scaled_vec(0.999, priors$nu, priors$mu_tau_alpha, priors$scale_tau_alpha),
      prior_tau_sig_fixed   = qt_trunc_scaled_vec(0.999, priors$nu, priors$mu_tau_sig,   priors$scale_tau_sig)
    ))
  }
}
