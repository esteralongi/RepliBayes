## Pure-R checks of the metric functions (no Stan needed).

test_that("point-estimate metrics have paper names and valid ranges (S=3, k=2)", {
  set.seed(1)
  a  <- list(rnorm(2000, 1, 0.3), rnorm(2000, 1, 0.3), rnorm(2000, 1, 0.3))
  mu <- rnorm(2000, 1, 0.3)
  eps <- 0.2

  ind  <- compute_replication_probs_indep(a, eps)
  hier <- compute_replication_probs_hier(a, mu, eps)

  expect_true(all(c("P_overall_2", "P_pos_2", "P_cond_O1_m1") %in% names(ind)))
  expect_true(all(c("P_beta", "P_gen_overall_2", "P_gen_pos_2", "P_c_pos_2") %in% names(hier)))

  probs <- unlist(hier); probs <- probs[!is.na(probs)]
  expect_true(all(probs >= 0 & probs <= 1))
  expect_gt(ind$P_pos_2, 0.9)   # clearly positive, homogeneous effect
})

test_that("replication_ess returns the right shape and finite MCSEs", {
  set.seed(2)
  a  <- list(matrix(rnorm(1000 * 4, 1, 0.3), 1000, 4),
             matrix(rnorm(1000 * 4, 1, 0.3), 1000, 4),
             matrix(rnorm(1000 * 4, 1, 0.3), 1000, 4))
  mu <- matrix(rnorm(1000 * 4, 1, 0.3), 1000, 4)

  tab_i <- replication_ess(a, mu = NULL, eps = 0.2)
  tab_h <- replication_ess(a, mu = mu,   eps = 0.2)

  expect_setequal(names(tab_i), c("metric", "p", "k", "n", "ess", "mcse"))
  expect_true(all(c("P_overall_2", "P_pos_2", "P_cond_O1_m1") %in% tab_i$metric))
  expect_true(all(c("P_beta", "P_gen_pos_2", "P_c_pos_2") %in% tab_h$metric))
  expect_gt(nrow(tab_h), nrow(tab_i))

  p <- tab_h$p[!is.na(tab_h$p)]
  expect_true(all(p >= 0 & p <= 1))
})

test_that("point and MCSE versions give identical point estimates", {
  set.seed(3)
  a  <- list(rnorm(5000, 0.5, 0.6), rnorm(5000, 0.5, 0.6), rnorm(5000, 0.5, 0.6))
  mu <- rnorm(5000, 0.5, 0.6); eps <- 0.3

  pt <- compute_replication_probs_hier(a, mu, eps)
  es <- replication_ess(a, mu = mu, eps = eps)
  get <- function(tab, nm) tab$p[tab$metric == nm]

  for (nm in c("P_overall_2", "P_pos_2", "P_beta", "P_gen_pos_2", "P_cond_O1_m1"))
    expect_equal(pt[[nm]], get(es, nm))
})

test_that("generalizes to any S, k, and m", {
  set.seed(4)
  S <- 5
  a <- lapply(seq_len(S), function(i) rnorm(1500, 1, 0.4))
  mu <- rnorm(1500, 1, 0.4)

  ## consensus level k = 3, conditional m = 2, with S = 5 studies
  es <- replication_ess(a, mu = mu, eps = 0.2, k = 3, m = 2)
  expect_true(all(c("P_overall_3", "P_pos_3", "P_gen_pos_3", "P_c_pos_3") %in% es$metric))
  expect_true(all(sprintf("P_cond_O%d_m2", 1:S) %in% es$metric))   # one per study
  expect_false(any(grepl("_2$", es$metric)))                       # k=2 not requested

  ## a vector of consensus levels produces one block each
  es2 <- replication_ess(a, mu = NULL, eps = 0.2, k = c(2, 4))
  expect_true(all(c("P_overall_2", "P_overall_4") %in% es2$metric))
})

test_that("threshold-free overlap is a coefficient in [0, 1] for any S", {
  set.seed(5)
  beta <- list(rnorm(3000, 0, 1), rnorm(3000, 0.2, 1),
               rnorm(3000, -0.1, 1), rnorm(3000, 0.05, 1))   # S = 4
  out <- replication_measures_multi(beta)
  expect_true(out$density_overlap >= 0 && out$density_overlap <= 1)
  expect_gt(out$mean_max_divergence, 0)
})
