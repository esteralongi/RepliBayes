## priors.R -- elicited hyperparameters and truncated-Student-t helpers.

#' Quantile function of a truncated Student-t on (0, Inf)
#'
#' Vectorised quantile function of a Student-t with `df` degrees of freedom,
#' location `mu` and `scale`, truncated to the positive half-line. Used to place
#' the between-study heterogeneities (`tau`) and positive parameters at chosen
#' percentiles of their priors.
#'
#' @param p_vec Numeric vector of probabilities in (0, 1).
#' @param df Degrees of freedom.
#' @param mu Location.
#' @param scale Positive scale.
#'
#' @return Numeric vector of quantiles.
#' @export
qt_trunc_scaled_vec <- function(p_vec, df, mu, scale) {
  p_lower <- stats::pt(-mu / scale, df = df)
  Z <- 1 - p_lower
  mu + scale * stats::qt(p_lower + p_vec * Z, df = df)
}

#' Random draws from a truncated Student-t on (0, Inf)
#'
#' @param n Number of draws.
#' @param nu Degrees of freedom.
#' @param mu Location.
#' @param scale Positive scale.
#'
#' @return Numeric vector of `n` positive draws.
#' @export
rtrunc_t_pos <- function(n, nu, mu, scale) {
  qt_trunc_scaled_vec(stats::runif(n), nu, mu, scale)
}

#' Default (elicited) prior hyperparameters
#'
#' The illustrative empirical-Bayes hyperparameters used in the paper. The three
#' generative means (`mu_a`, `mu_alpha`, `mu_sig`) and their scales, and the
#' three between-study heterogeneities on the log scale (`mu_tau_*`,
#' `scale_tau_*`), define the elicited priors of the hierarchical model. All are
#' Student-t with `nu` degrees of freedom; the intercept, residual scale and the
#' heterogeneities are truncated to the positive half-line.
#'
#' @return A named list of hyperparameters, suitable as the `priors` argument of
#'   [build_stan_data()], [fit_hierarchical()], [fit_independence()] and
#'   [fit_replicability()]. Copy and edit fields to supply your own priors.
#'
#' @examples
#' priors <- default_priors()
#' priors$mu_a <- 0        # centre the effect prior at zero
#' @export
default_priors <- function() {
  list(
    nu           = 3,
    mu_a         =  0.05, scale_a         = 0.60,   # effect beta (untruncated)
    mu_alpha     =  8.5,  scale_alpha     = 3.0,    # intercept alpha (> 0)
    mu_sig       =  4.5,  scale_sig       = 2.0,    # residual sigma (> 0)
    mu_tau_a     = -1.00, scale_tau_a     = 0.30,   # tau_beta  (heterogeneity of the effect)
    mu_tau_alpha = -1.30, scale_tau_alpha = 0.40,   # tau_alpha (heterogeneity of the intercept)
    mu_tau_sig   = -1.00, scale_tau_sig   = 0.28    # tau_sigma (heterogeneity of the residual scale)
  )
}

## Practical-relevance threshold eps = frac * baseline (genotype-0) mean outcome.
.default_eps <- function(data, frac = 0.10) {
  frac * mean(data$m[data$x == 0])
}
