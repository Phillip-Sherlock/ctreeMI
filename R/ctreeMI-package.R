#' ctreeMI: Conditional Inference Trees with Multiple Imputation
#'
#' @description
#' `ctreeMI` implements the **stacked-imputation / Stack ? M** workflow for
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
#' **Stack ? M correction**: Divide the node-level significance threshold by
#' M (`alpha_corrected = alpha / M`) before the pruning decision. This is
#' equivalent to dividing each node's test statistic by M before comparing
#' it to the Bonferroni-corrected critical value. Sherlock et al. (2026)
#' validated this approach via Monte Carlo simulation under MCAR, showing
#' sub-nominal (conservative) type-I error and acceptable power.
#'
#' ## Main function
#'
#' The primary user-facing function is [ctree_stacked()]. It accepts a
#' `mids` object from [mice::mice()], a list of imputed data frames, or a
#' plain data frame, and returns a fitted tree with full `partykit`
#' compatibility.
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
