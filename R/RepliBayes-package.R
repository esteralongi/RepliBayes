#' RepliBayes: Bayesian Quantification of Replicability
#'
#' A Bayesian hierarchical framework that quantifies the replicability of an
#' effect across several studies as posterior probabilities, defined relative to
#' a practical-relevance threshold on the scale of the effects (`eps`) rather
#' than on p-values.
#'
#' The package has three top-level entry points:
#' \describe{
#'   \item{[fit_replicability()]}{fits the hierarchical (partial-pooling) and
#'     independence-limit models and returns all replication probabilities
#'     (empirical, retrospective and prospective), each with a Monte Carlo
#'     standard error.}
#'   \item{[sensitivity_prior()]}{prior-sensitivity analysis: refits the model
#'     while shifting the prior on a chosen between-study heterogeneity.}
#'   \item{[simulate_replicability()]}{calibration by simulation, varying one
#'     component of the generative model and summarising the metrics across
#'     replicates.}
#' }
#' Priors are chosen with [default_priors()] (or supplied by the user).
#'
#' The individual building blocks are also exported:
#' [build_stan_data()], [fit_hierarchical()], [fit_independence()],
#' [replication_ess()], [compute_replication_probs_hier()],
#' [compute_replication_probs_indep()], [predictive_retro()], [predictive_pro()],
#' [predict_new_study()] and the threshold-free [replication_measures_multi()].
#'
#' @section GTEx data:
#' The donor-level GTEx data used in the paper are controlled-access and cannot
#' be redistributed. The package ships a synthetic dataset,
#' [synthetic_data], with the same structure so that the whole pipeline can be
#' reproduced.
#'
#' @keywords internal
"_PACKAGE"
