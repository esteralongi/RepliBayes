#' RepliBayes: Bayesian Quantification of Replicability
#'
#' A Bayesian hierarchical framework that quantifies the replicability of an
#' effect across several studies as posterior probabilities, defined relative to
#' a practical-relevance threshold on the scale of the effects (`eps`).
#'
#' There are two ways to use the package, depending on your input.
#'
#' \strong{From raw data.} If your data are in long format, one row per
#' observation, with columns `study` (the study/group label), `x` (the
#' predictor) and `y` (the outcome), everything runs from a single call:
#' \describe{
#'   \item{[fit_replicability()]}{fits the hierarchical (partial-pooling) and
#'     independence-limit models and returns all replication probabilities
#'     (hierarchical, independent, retrospective and prospective), each with a
#'     Monte Carlo standard error.}
#'   \item{[sensitivity_prior()]}{prior-sensitivity analysis: refits the model
#'     while shifting the prior on a chosen between-study heterogeneity.}
#'   \item{[simulate_replicability()]}{calibration by simulation, varying one
#'     component of the generative model and summarising the metrics across
#'     replicates.}
#' }
#' Priors are chosen with [default_priors()] or supplied by the user (with the
#' truncated-Student-t helpers [qt_trunc_scaled_vec()] and [rtrunc_t_pos()]).
#'
#' \strong{From posterior draws.} If you have already fitted your own model and
#' have posterior draws of the study-specific effects, you can compute the
#' metrics directly, without a data frame. These take the draws as a list of `S`
#' \[iterations x chains\] matrices (one per study), and optionally the draws of
#' the generative effect:
#' \describe{
#'   \item{[compute_replication_probs_indep()], [compute_replication_probs_hier()]}{point
#'     estimates of all replication probabilities.}
#'   \item{[replication_ess()]}{the same probabilities, each with a Monte Carlo
#'     standard error.}
#'   \item{[predictive_retro()], [predictive_pro()]}{retrospective and
#'     prospective predictive metrics.}
#'   \item{[replication_threshold_free()]}{threshold-free summaries (divergence
#'     and density overlap), needing no `eps`.}
#' }
#'
#' The lower-level building blocks are also exported:
#' [build_stan_data()], [fit_hierarchical()], [fit_independence()] and
#' [predict_new_study()].
#'
#' @section GTEx data:
#' The donor-level GTEx data used in the paper are controlled-access and cannot
#' be redistributed. The package ships a synthetic dataset,
#' [synthetic_data], with the same structure so that the whole pipeline can be
#' reproduced.
#'
#' @keywords internal
"_PACKAGE"
