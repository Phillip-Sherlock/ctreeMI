#' Center Outcomes on Group Means Within Each Imputed Dataset
#'
#' Subtracts each group's mean of the outcome from every observation in that
#' group, separately within each imputed dataset. A two-level factor or
#' logical outcome is coded 0/1 first, so its centered value is the
#' deviation from the group's rate. An ordered factor with three or more
#' levels is replaced by its centered integer score, and an unordered factor
#' with three or more levels by one centered 0/1 indicator per level after
#' the first, named `outcome.level`. Used to grow a tree on within-group
#' variation when the data are nested in a small number of clusters; see
#' [confirm_within()].
#'
#' @param data A `mids` object from [mice::mice()] or a list of imputed data
#'   frames.
#' @param outcome Character vector naming the outcome column(s) to center.
#'   Each must be numeric, logical, or a factor.
#' @param by Name of the grouping column, such as a cohort or site identifier.
#'
#' @details
#' When clusters differ in their average outcome, any predictor that tracks
#' cluster membership looks like a strong splitting variable, and a tree grown
#' on the raw outcome will tend to split on it. Centering removes the
#' between-cluster differences in the outcome, so the splits the tree finds
#' have to reflect associations within clusters. The cost is that subgroups
#' defined by cluster-level variables cannot be discovered; with few clusters
#' they could not be distinguished from cluster differences anyway.
#'
#' Centering changes the outcome's scale and, for a nominal outcome, its
#' shape: a tree grown on a centered binary outcome is a regression tree on
#' deviations from each cluster's rate, and a tree grown on a centered nominal
#' outcome is a multivariate regression tree on deviations from each
#' cluster's level proportions. Node estimates should be taken from an
#' uncentered confirmation. [ctree_stacked()] with `center_by` does the
#' rewriting of the formula that a nominal outcome needs.
#'
#' @return A list of data frames, one per imputation, with the named outcome
#'   column(s) replaced by their centered, numeric versions. The attribute
#'   `"outcome_columns"` is a named list giving, for each outcome, the
#'   column(s) that now carry it.
#'
#' @seealso [confirm_within()], [ctree_stacked()]
#' @examples
#' set.seed(2)
#' imps <- lapply(1:3, function(i) {
#'   g <- factor(rep(c("a", "b", "c"), each = 40))
#'   data.frame(g = g, x = rnorm(120),
#'              y = rnorm(120) + c(a = 0, b = 1, c = 2)[as.character(g)])
#' })
#' cen <- center_outcome(imps, outcome = "y", by = "g")
#' tapply(cen[[1]]$y, cen[[1]]$g, mean)   # all zero
#' @export
center_outcome <- function(data, outcome, by) {
  imps <- .as_imp_list(data)
  nm <- names(imps[[1L]])
  if (!is.character(by) || length(by) != 1L || !by %in% nm)
    stop("`by` must name one column of the data.")
  if (!all(outcome %in% nm))
    stop("Outcome(s) ", paste(setdiff(outcome, nm), collapse = ", "), " not found.")
  cols <- stats::setNames(vector("list", length(outcome)), outcome)
  # levels are taken from the first imputation so every dataset agrees
  for (y in outcome) {
    v <- imps[[1L]][[y]]
    cols[[y]] <- if (is.factor(v) && nlevels(v) > 2L && !is.ordered(v))
      paste(y, levels(v)[-1L], sep = ".") else y
  }
  ctr <- function(v, g) v - stats::ave(v, g, FUN = function(z) mean(z, na.rm = TRUE))
  out <- lapply(imps, function(d) {
    g <- d[[by]]
    for (y in outcome) {
      v <- d[[y]]
      if (is.logical(v)) v <- as.numeric(v)
      if (is.factor(v)) {
        v <- factor(v, levels = levels(imps[[1L]][[y]]), ordered = is.ordered(v))
        if (nlevels(v) == 2L) {
          v <- .as01(v)
        } else if (is.ordered(v)) {
          v <- as.numeric(v)                         # integer score
        } else {
          for (lv in levels(v)[-1L])
            d[[paste(y, lv, sep = ".")]] <- ctr(as.numeric(v == lv), g)
          d[[y]] <- NULL
          next
        }
      }
      if (!is.numeric(v))
        stop("Outcome '", y, "' is not numeric, logical or a factor and cannot be centered.")
      d[[y]] <- ctr(v, g)
    }
    d
  })
  attr(out, "outcome_columns") <- cols
  out
}


#' Confirm a Tree Within Clusters, Three Ways
#'
#' For data nested in a small number of clusters, where the discovery and
#' confirmation halves share clusters, runs the confirmation three ways and
#' places the results side by side: the standard cluster-robust confirmation,
#' the same confirmation with cluster fixed effects, and a separate
#' confirmation within each cluster.
#'
#' @param tree A `ctreeMI` object fitted on the discovery data.
#' @param data The confirmation data as a `mids` object or a list of imputed
#'   data frames, imputed separately from the discovery data.
#' @param cluster Name of the column identifying the cluster, such as a cohort
#'   or site. Required.
#' @param by_cluster Logical. Also confirm within each cluster separately
#'   (default `TRUE`). With many clusters this is the slow part.
#' @param alpha Level at which a split counts as confirmed in the summary.
#'   Default 0.05.
#' @param ... Further arguments passed to [confirm_ctreeMI()], such as
#'   `outcome_type`, `min_node`, `adjust`, `conf.level` and `vcov_type`.
#'
#' @details
#' Splitting whole clusters between the halves (the `cluster` argument of
#' [split_holdout()]) keeps the two halves independent at the cluster level,
#' but leaves half the clusters for confirmation. When that is too few, the
#' alternative is to split participants within every cluster
#' (`strata = cluster` in [split_holdout()]), so the confirmation keeps all
#' clusters. Both halves then contain the same clusters, so a split that
#' merely sorts clusters with different outcome levels can replicate for the
#' wrong reason. The three confirmations are the check:
#'
#' \describe{
#'   \item{robust}{[confirm_ctreeMI()] with `cluster`. The standard test with
#'     cluster-robust variance, on the raw outcome. Node estimates and the
#'     contrasts to report come from here.}
#'   \item{within}{The same test with the cluster identifier added as fixed
#'     effects (`fixed = cluster`), so each split's contrast is estimated from
#'     variation within clusters only. A split that reflected differences
#'     between clusters has a contrast near zero here.}
#'   \item{by_cluster}{[confirm_ctreeMI()] run on each cluster's confirmation
#'     observations alone, without cluster-robust variance. Shows whether a
#'     split points the same way in most clusters or rests on one.}
#' }
#'
#' The `summary` element lines the three up by split. `reading` is
#' `"confirmed within clusters"` when a split passes both the robust and the
#' within-cluster test at `alpha`, `"between-cluster difference"` when it
#' passes the robust test but not the within-cluster one, and
#' `"not confirmed"` otherwise. `clusters_same_sign` counts the clusters whose
#' own contrast shares the sign of the robust contrast, out of
#' `clusters_tested`, the clusters in which the split could be tested at all;
#' small clusters leave some splits untested, which is expected.
#'
#' The claim supported is replication in held-out participants from these
#' clusters. How far it extends to other clusters depends on how consistent
#' the effect is across them, which is what the per-cluster confirmations
#' show.
#'
#' For a binary outcome the robust and within-cluster contrasts are both log
#' odds ratios and are directly comparable. The discovery tree should have
#' been grown on the outcome centered by cluster, with `center_by` in
#' [ctree_stacked()] or [center_outcome()], so that cluster-level variables
#' do not take over the tree; [discover_confirm()] with `within = TRUE` does
#' the whole sequence.
#'
#' @return An object of class `"ctreeMI_confirm_within"`: a list with
#'   elements `robust` and `within` (each a `ctreeMI_confirm` object),
#'   `by_cluster` (a named list of `ctreeMI_confirm` objects, or `NULL` where a
#'   cluster could not be confirmed), `summary` (a data frame with one row per
#'   split and outcome; `contrast_level` is present when a nominal outcome
#'   was confirmed), `cluster`, `n_clusters` and `alpha`.
#'
#' @seealso [confirm_ctreeMI()], [center_outcome()], [discover_confirm()]
#' @examples
#' \dontrun{
#' parts <- split_holdout(d, strata = "cohort", seed = 1)
#' imp_d <- mice::mice(parts$discover, m = 20, printFlag = FALSE, seed = 1)
#' imp_c <- mice::mice(parts$confirm,  m = 20, printFlag = FALSE, seed = 2)
#' tree  <- ctree_stacked(y ~ x1 + x2, data = imp_d, center_by = "cohort")
#' confirm_within(tree, imp_c, cluster = "cohort")
#' }
#' @export
confirm_within <- function(tree, data, cluster, by_cluster = TRUE,
                           alpha = 0.05, ...) {
  if (missing(cluster) || is.null(cluster))
    stop("`cluster` is required: name the column identifying the cluster.")
  dots <- list(...)
  if ("fixed" %in% names(dots))
    stop("`fixed` is set by confirm_within(); pass only `cluster`.")
  imps <- .as_imp_list(data)
  if (!cluster %in% names(imps[[1L]]))
    stop("`cluster` column '", cluster, "' not found in confirmation data.")

  robust <- confirm_ctreeMI(tree, imps, cluster = cluster, ...)
  within <- confirm_ctreeMI(tree, imps, cluster = cluster, fixed = cluster, ...)

  levs <- sort(unique(as.character(imps[[1L]][[cluster]])))
  per <- NULL
  if (isTRUE(by_cluster)) {
    per <- lapply(levs, function(k) {
      dl <- lapply(imps, function(d) d[as.character(d[[cluster]]) == k, , drop = FALSE])
      tryCatch(suppressWarnings(confirm_ctreeMI(tree, dl, ...)),
               error = function(e) NULL)
    })
    names(per) <- levs
  }

  s <- robust$splits
  if (!is.null(s)) {
    w <- within$splits
    key_s <- paste(s$node_id, s$outcome)
    key_w <- paste(w$node_id, w$outcome)
    idx <- match(key_s, key_w)
    s$contrast_within <- w$contrast[idx]
    s$lower_within    <- w$lower[idx]
    s$upper_within    <- w$upper[idx]
    s$p_adj_within    <- w$p_adj[idx]
    if (!is.null(per)) {
      cnt <- t(vapply(seq_len(nrow(s)), function(r) {
        vals <- vapply(per, function(cc) {
          if (is.null(cc) || is.null(cc$splits)) return(NA_real_)
          cs <- cc$splits
          v <- cs$contrast[cs$node_id == s$node_id[r] & cs$outcome == s$outcome[r]]
          if (length(v) == 1L) v else NA_real_
        }, numeric(1L))
        c(tested = sum(!is.na(vals)),
          same   = sum(sign(vals) == sign(s$contrast[r]), na.rm = TRUE))
      }, numeric(2L)))
      s$clusters_tested    <- as.integer(cnt[, "tested"])
      s$clusters_same_sign <- as.integer(cnt[, "same"])
    } else {
      s$clusters_tested <- s$clusters_same_sign <- NA_integer_
    }
    rob_ok <- !is.na(s$p_adj) & s$p_adj < alpha
    wit_ok <- !is.na(s$p_adj_within) & s$p_adj_within < alpha
    s$reading <- ifelse(rob_ok & wit_ok, "confirmed within clusters",
                 ifelse(rob_ok, "between-cluster difference", "not confirmed"))
    s <- s[, c("node_id", "depth", "split_var", "rule", "outcome",
               "contrast", "lower", "upper", "contrast_level", "p_adj",
               "contrast_within", "lower_within", "upper_within", "p_adj_within",
               "n_clusters", "clusters_tested", "clusters_same_sign", "reading")]
    if (all(is.na(s$contrast_level))) s$contrast_level <- NULL
    rownames(s) <- NULL
  }

  structure(list(robust = robust, within = within, by_cluster = per,
                 summary = s, cluster = cluster,
                 n_clusters = length(levs), alpha = alpha),
            class = "ctreeMI_confirm_within")
}


#' @export
print.ctreeMI_confirm_within <- function(x, digits = 3L, ...) {
  cat("Confirmation within clusters, three ways\n")
  cat(sprintf("  %d clusters (%s), all present in the confirmation data\n",
              x$n_clusters, x$cluster))
  cat("  robust: cluster-robust variance on the raw outcome (report these estimates)\n")
  cat("  within: the same test with cluster fixed effects\n")
  if (!is.null(x$by_cluster))
    cat("  by_cluster: each cluster confirmed on its own\n")
  if (is.null(x$summary)) {
    cat("\nNo splits to confirm.\n")
    return(invisible(x))
  }
  s <- x$summary
  out <- data.frame(node = s$node_id, split = s$split_var, outcome = s$outcome,
                    contrast = round(s$contrast, digits),
                    p_adj = signif(s$p_adj, digits),
                    within = round(s$contrast_within, digits),
                    p_adj_within = signif(s$p_adj_within, digits),
                    same_sign = ifelse(is.na(s$clusters_tested), "",
                                       paste0(s$clusters_same_sign, "/", s$clusters_tested)),
                    reading = s$reading, stringsAsFactors = FALSE)
  cat("\n")
  print(out, row.names = FALSE)
  cat(sprintf("\n'confirmed within clusters': adjusted p < %s in both the robust and the within test.\n",
              format(x$alpha)))
  cat("'between-cluster difference': robust test passed, within-cluster test did not.\n")
  cat("same_sign: clusters whose own contrast shares the robust contrast's sign / clusters tested.\n")
  invisible(x)
}
