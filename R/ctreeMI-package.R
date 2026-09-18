#' ctreeMI: Conditional Inference Trees with Multiple Imputation
#'
#' @description
#' `ctreeMI` implements the **stacked-imputation / Stack / M** workflow for
#' conditional inference trees (ctree) described in Sherlock et al. (2026).
#'
#' ## The problem
#'
#' Missing data are ubiquitous in applied research. Multiple imputation (MI)
#' is the principled solution, but pooling results across imputations is
#' straightforward only for linear-combination models (via Rubin's rules).
#' Conditional inference trees (ctree) cannot be pooled across imputations
#' because structurally different trees -- which split on different variables
#' at different nodes -- define different subgroups, so the targets of
#' inference are no longer aligned across imputations.
#'
#' ## The solution
#'
#' **Stack** the M imputed datasets vertically into one data frame of
#' M x n rows, then fit a **single** ctree. The single tree is coherent and
#' interpretable. The only problem: stacking inflates the nominal sample size
#' by M, so test statistics at each node are similarly inflated and the tree
#' over-splits.
#'
#' **Stack / M correction**: divide each node-level test statistic by M,
#' recompute its p-value from the chi-square reference distribution, reapply
#' the multiplicity adjustment across candidate variables, and compress the
#' tree bottom-up. This is not the same as dividing the significance
#' threshold by M; the two rules coincide only at M = 1, and threshold
#' rescaling under-corrects by an order of magnitude at M = 30. See
#' [ctree_stacked()] for the derivation.
#'
#' ## What the correction does and does not do
#'
#' The correction removes the inflation attributable to the stacked sample
#' size. It does not yield a calibrated node-level test: when the imputation
#' model conditions on the outcome, which is recommended practice, the test
#' rejects more often than its nominal level implies. The procedure recovers
#' known structure and produces reproducible partitions well, but its
#' node-level p-values are not error rates. See the "Node-level calibration"
#' section of [ctree_stacked()].
#'
#' ## Recommended workflow
#'
#' Treat the tree as discovery and confirm it on independent data.
#' [discover_confirm()] splits the sample, imputes each half separately,
#' fits the tree on one half, and tests the resulting partition on the other
#' with a procedure that pools across imputations by Rubin's rules and so has
#' valid error control. [split_holdout()] and [confirm_ctreeMI()] expose the
#' steps individually.
#'
#' ## Main functions
#'
#' [ctree_stacked()] fits the tree. It accepts a `mids` object from
#' [mice::mice()], a list of imputed data frames, or a plain data frame, and
#' returns a fitted tree with full `partykit` compatibility.
#' [confirm_ctreeMI()] tests a fitted tree's partition on held-out data.
#' [node_table()] reports effective sample sizes per terminal node.
#'
#' ## Citation
#'
#' If you use `ctreeMI`, please cite both the package and the methodological
#' paper:
#'
#' ```
#' Sherlock, P., Mansolf, M., Hofheimer, J., Hockett, C. W., O'Connor, T. G.,
#'   Roubinov, D., Graff, J. C., Lai, J.-S., Bush, N. R., Wright, R. J., &
#'   Chiu, Y.-H. M. (2026). Beyond linear risk: A machine learning approach
#'   to understanding perinatal depression in context.
#'   Multivariate Behavioral Research, 1-16.
#'   https://doi.org/10.1080/00273171.2026.2661244
#' ```
#'
#' The underlying ctree algorithm should also be cited:
#'
#' ```
#' Hothorn, T., Hornik, K., & Zeileis, A. (2006). Unbiased recursive
#'   partitioning: A conditional inference framework. Journal of
#'   Computational and Graphical Statistics, 15(3), 651-674.
#'
#' Hothorn, T., & Zeileis, A. (2015). partykit: A modular toolkit for
#'   recursive partitioning in R. Journal of Machine Learning Research,
#'   16, 3905-3909.
#' ```
#'
#' @keywords internal
"_PACKAGE"
