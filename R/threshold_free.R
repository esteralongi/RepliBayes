## threshold_free.R -- summaries of the study-specific effects that do not
## require the practical-relevance threshold eps.

#' Threshold-free replication measures
#'
#' Summaries of agreement among the three study-specific effects that do not
#' depend on the practical-relevance threshold `eps`.
#'
#' @param d1,d2,d3 Numeric vectors of posterior draws of the three
#'   study-specific effects.
#' @param dens_n Number of grid points for the density-overlap integral
#'   (default 500).
#'
#' @return A named list with:
#'   \describe{
#'     \item{`mean_max_divergence`, `sd_max_divergence`}{posterior mean and sd
#'       of the range \eqn{\max_s \beta_s - \min_s \beta_s} across studies at
#'       each draw (largest minus smallest study effect).}
#'     \item{`density_overlap_3way`}{the overlapping coefficient: the integral of
#'       the pointwise minimum of the three posterior densities (0 = disjoint,
#'       1 = identical distributions).}
#'   }
#'
#' @examples
#' \dontrun{
#' post <- rstan::extract(fit_hier_abc)
#' replication_measures_multi(post$a[, 1], post$a[, 2], post$a[, 3])
#' }
#' @export
replication_measures_multi <- function(d1, d2, d3, dens_n = 500) {

  ## draws in a matrix: rows = MCMC draws, columns = the 3 studies
  draws_mat <- cbind(d1, d2, d3)

  ## maximum divergence: range (max - min) across the 3 studies, per draw
  max_diff_per_draw <- apply(draws_mat, 1, max) - apply(draws_mat, 1, min)
  mean_max_diff <- mean(max_diff_per_draw)
  sd_max_diff   <- stats::sd(max_diff_per_draw)

  ## density overlap on a shared grid
  all_draws <- as.vector(draws_mat)
  grid <- seq(min(all_draws), max(all_draws), length.out = dens_n)

  ## approximate a density on the shared grid
  get_dens <- function(x_draws, shared_grid) {
    dens <- stats::density(x_draws)
    app  <- stats::approx(dens$x, dens$y, xout = shared_grid)$y
    app[is.na(app)] <- 0
    app
  }

  dens1 <- get_dens(d1, grid)
  dens2 <- get_dens(d2, grid)
  dens3 <- get_dens(d3, grid)

  ## pointwise minimum across the 3 curves, integrated
  ## overlap = area shared by all three densities (overlapping coefficient)
  min_dens_3way <- pmin(dens1, dens2, dens3)
  overlap_area  <- sum(min_dens_3way) * (grid[2] - grid[1])

  list(
    mean_max_divergence  = mean_max_diff,
    sd_max_divergence    = sd_max_diff,
    density_overlap_3way = overlap_area
  )
}
