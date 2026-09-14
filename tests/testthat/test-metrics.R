## Pure-R checks of the metric functions (no Stan needed).

test_that("point-estimate metrics have the expected names and ranges", {
  set.seed(1)
  a1 <- rnorm(2000, 1, 0.3); a2 <- rnorm(2000, 1, 0.3); a3 <- rnorm(2000, 1, 0.3)
  mu <- rnorm(2000, 1, 0.3)
  eps <- 0.2

  ind <- compute_replication_probs_indep(a1, a2, a3, eps)
  hier <- compute_replication_probs_hier(a1, a2, a3, mu, eps)

  expect_true(all(c("Agr_2", "P_pos_2", "Cond_S1") %in% names(ind)))
  expect_true(all(c("P_beta", "Agr_Gen_2", "Pos_Gen_2", "CondPos_Gen") %in% names(hier)))

  probs <- unlist(hier)
  probs <- probs[!is.na(probs)]
  expect_true(all(probs >= 0 & probs <= 1))

  ## a clearly positive, homogeneous effect -> high positive agreement
  expect_gt(ind$P_pos_2, 0.9)
})

test_that("replication_ess returns the right shape and finite MCSEs", {
  set.seed(2)
  ## fake [iterations x chains] draws
  a1 <- matrix(rnorm(1000 * 4, 1, 0.3), 1000, 4)
  a2 <- matrix(rnorm(1000 * 4, 1, 0.3), 1000, 4)
  a3 <- matrix(rnorm(1000 * 4, 1, 0.3), 1000, 4)
  mu <- matrix(rnorm(1000 * 4, 1, 0.3), 1000, 4)

  tab_i <- replication_ess(a1, a2, a3, mu = NULL, eps = 0.2)
  tab_h <- replication_ess(a1, a2, a3, mu = mu, eps = 0.2)

  expect_setequal(names(tab_i), c("metric", "p", "k", "n", "ess", "mcse"))
  expect_true(all(c("P_overall_2", "P_pos_2", "P_cond_O1_m1") %in% tab_i$metric))
  expect_true(all(c("P_beta", "P_gen_pos_2", "P_c_pos_2") %in% tab_h$metric))
  expect_gt(nrow(tab_h), nrow(tab_i))

  p <- tab_h$p[!is.na(tab_h$p)]
  expect_true(all(p >= 0 & p <= 1))
})

test_that("indicator definitions match between point and MCSE versions", {
  set.seed(3)
  a1 <- rnorm(5000, 0.5, 0.6); a2 <- rnorm(5000, 0.5, 0.6); a3 <- rnorm(5000, 0.5, 0.6)
  mu <- rnorm(5000, 0.5, 0.6); eps <- 0.3

  pt <- compute_replication_probs_hier(a1, a2, a3, mu, eps)
  es <- replication_ess(a1, a2, a3, mu = mu, eps = eps)

  get <- function(tab, nm) tab$p[tab$metric == nm]
  expect_equal(pt$Agr_2,     get(es, "P_overall_2"))
  expect_equal(pt$P_pos_2,   get(es, "P_pos_2"))
  expect_equal(pt$P_beta,    get(es, "P_beta"))
  expect_equal(pt$Pos_Gen_2, get(es, "P_gen_pos_2"))
})

test_that("threshold-free overlap is a coefficient in [0, 1]", {
  set.seed(4)
  d1 <- rnorm(3000, 0, 1); d2 <- rnorm(3000, 0.2, 1); d3 <- rnorm(3000, -0.1, 1)
  out <- replication_measures_multi(d1, d2, d3)
  expect_true(out$density_overlap_3way >= 0 && out$density_overlap_3way <= 1)
  expect_gt(out$mean_max_divergence, 0)
})
