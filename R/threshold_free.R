## threshold_free.R -- summaries of the study-specific effects that do not
## require the practical-relevance threshold eps.

#' Threshold-free replication measures
#'
#' Summaries of agreement among the `S` study-specific effects that do not depend
#' on the practical-relevance threshold `eps`.
#'
#' @param beta A list of `S` numeric vectors of posterior draws of the
#'   study-specific effects.
#' @param dens_n Number of grid points for the density-overlap integral
#'   (default 500).
#'
#' @return A named list with:
#'   \describe{
#'     \item{`mean_max_divergence`, `sd_max_divergence`}{posterior mean and sd
#'       of the range \eqn{\max_s \beta_s - \min_s \beta_s} across studies at
#'       each draw (largest minus smallest study effect).}
#'     \item{`density_overlap`}{the overlapping coefficient: the integral of the
#'       pointwise minimum of the `S` posterior densities (0 = disjoint,
#'       1 = identical distributions).}
#'   }
#'
#' @examples
#' \dontrun{
#' post <- rstan::extract(fit_hier)
#' replication_measures_multi(list(post$beta[, 1], post$beta[, 2], post$beta[, 3]))
#' }
#' @export
replication_measures_multi <- function(beta, dens_n = 500) {
  stopifnot(is.list(beta), length(beta) >= 2L)

  ## draws in a matrix: rows = MCMC draws, columns = the S studies
  draws_mat <- do.call(cbind, beta)

  ## maximum divergence: range (max - min) across the S studies, per draw
  max_diff_per_draw <- apply(draws_mat, 1, max) - apply(draws_mat, 1, min)
  mean_max_diff <- mean(max_diff_per_draw)
  sd_max_diff   <- stats::sd(max_diff_per_draw)

  ## density overlap on a shared grid
  all_draws <- as.vector(draws_mat)
  grid <- seq(min(all_draws), max(all_draws), length.out = dens_n)

  ## approximate a density on the shared grid
  get_dens <- function(x_draws) {
    d   <- stats::density(x_draws)
    app <- stats::approx(d$x, d$y, xout = grid)$y
    app[is.na(app)] <- 0
    app
  }
  dens <- lapply(beta, get_dens)

  ## pointwise minimum across all S curves, integrated
  ## overlap = area shared by all studies' densities (overlapping coefficient)
  min_dens  <- Reduce(pmin, dens)
  overlap   <- sum(min_dens) * (grid[2] - grid[1])

  list(
    mean_max_divergence = mean_max_diff,
    sd_max_divergence   = sd_max_diff,
    density_overlap     = overlap
  )
}
