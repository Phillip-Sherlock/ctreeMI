#' Rescale a Node-Level Test Statistic for Stacking
#'
#' @description
#' Divides a node-level chi-square test statistic by the number of
#' imputations `M` and recomputes the corresponding p-value. This is the
#' Stack / M correction as defined in Sherlock et al. (2026).
#'
#' @details
#' Stacking `M` imputed datasets produces `M * n` rows, and a chi-square
#' test statistic computed on the stacked data scales approximately
#' linearly with `M`. Dividing the statistic by `M` returns it to the
#' scale of a single imputed dataset before the p-value is computed.
#'
#' This is **not** the same as dividing the significance threshold by `M`.
#' The two rules coincide only at `M = 1`:
#'
#' \itemize{
#'   \item Statistic rescaling rejects when `X > M * qchisq(1 - alpha, df)`.
#'   \item Threshold rescaling rejects when `X > qchisq(1 - alpha / M, df)`.
#' }
#'
#' For `df = 1`, `alpha = 0.05` and `M = 30`, the first requires
#' `X > 115.2` and the second only `X > 9.9`. Versions 0.1.0 and 0.2.0 of
#' this package implemented the second rule, which under-corrects severely.
#' See [rescale_alpha()].
#'
#' @param statistic Numeric. The raw test statistic computed on the
#'   stacked data.
#' @param m Integer. Number of imputations.
#' @param df Numeric. Degrees of freedom of the reference chi-square
#'   distribution.
#'
#' @return A list with `statistic_rescaled` and `p_value`.
#'
#' @references
#' Sherlock, P., Mansolf, M., Hofheimer, J., Hockett, C. W., O'Connor, T. G.,
#'   Roubinov, D., Graff, J. C., Lai, J.-S., Bush, N. R., Wright, R. J., &
#'   Chiu, Y.-H. M. (2026). Beyond linear risk: A machine learning approach
#'   to understanding perinatal depression in context.
#'   *Multivariate Behavioral Research*, 1-16.
#'   \doi{10.1080/00273171.2026.2661244}
#'
#' @examples
#' rescale_statistic(115.2, m = 30, df = 1)
#'
#' @export
rescale_statistic <- function(statistic, m, df) {
  if (!is.numeric(statistic) || length(statistic) != 1L)
    stop("`statistic` must be a single numeric value.")
  if (!is.numeric(m) || length(m) != 1L || m < 1)
    stop("`m` must be a single positive integer.")
  if (!is.numeric(df) || length(df) != 1L || df <= 0)
    stop("`df` must be a single positive number.")

  s <- statistic / m
  list(statistic_rescaled = s,
       p_value = stats::pchisq(s, df = df, lower.tail = FALSE))
}


# ---------------------------------------------------------------------------
# Internal: recover degrees of freedom from a (statistic, p-value) pair.
#
# partykit does not store the degrees of freedom of each node-level test, but
# the p-value it stores is a deterministic function of the statistic and df.
# For fixed `statistic`, pchisq(statistic, df, lower.tail = FALSE) is strictly
# increasing in df, so df can be recovered by root-finding and then VERIFIED
# by recomputing the p-value. Verification matters: if a future partykit
# version changes what it stores, this fails loudly instead of silently
# producing a wrongly pruned tree.
# ---------------------------------------------------------------------------

#' @keywords internal
.recover_df <- function(statistic, p, tol = 1e-8) {
  if (!is.finite(statistic) || !is.finite(p) ||
      statistic <= 0 || p <= 0 || p >= 1) return(NA_real_)

  f <- function(d) stats::pchisq(statistic, df = d, lower.tail = FALSE) - p

  lo <- 1e-8
  hi <- 1
  while (f(hi) < 0 && hi < 1e7) hi <- hi * 2
  if (f(hi) < 0 || f(lo) > 0) return(NA_real_)

  out <- try(stats::uniroot(f, c(lo, hi), tol = tol)$root, silent = TRUE)
  if (inherits(out, "try-error")) NA_real_ else out
}

# p-values are occasionally carried on the log scale; detect and undo.
#' @keywords internal
.as_pvalue <- function(v) {
  v <- as.numeric(v)
  fin <- v[is.finite(v)]
  if (length(fin) && all(fin <= 0)) return(exp(v))
  v
}

# Extract the statistic and raw p-value for the split actually taken at a node.
#' @keywords internal
.node_split_test <- function(node, var_names) {
  info <- partykit::info_node(node)
  crit <- info$criterion
  if (is.null(crit) || !is.matrix(crit)) return(NULL)
  if (!all(c("statistic", "p.value") %in% rownames(crit))) return(NULL)

  stat_row <- as.numeric(crit["statistic", ])
  p_row    <- .as_pvalue(crit["p.value", ])
  cols     <- colnames(crit)

  ## Number of candidate variables actually tested at this node -- the
  ## multiplicity the Bonferroni step must account for.
  k <- sum(is.finite(stat_row) & is.finite(p_row))
  if (k < 1L) return(NULL)

  ## Identify the column corresponding to the variable that was split on.
  sp  <- partykit::split_node(node)
  idx <- NA_integer_
  if (!is.null(sp) && !is.null(cols)) {
    vid <- try(as.integer(partykit::varid_split(sp)), silent = TRUE)
    if (!inherits(vid, "try-error") && !is.na(vid) &&
        vid <= length(var_names)) {
      idx <- match(var_names[vid], cols)
    }
  }
  ## Fall back to the most significant candidate, which is what ctree selects.
  if (is.na(idx)) {
    ok <- which(is.finite(p_row))
    if (!length(ok)) return(NULL)
    idx <- ok[which.min(p_row[ok])]
  }

  stat <- stat_row[idx]
  p    <- p_row[idx]
  if (!is.finite(stat) || !is.finite(p)) return(NULL)

  list(statistic = stat, p_value = p, k = k,
       variable = if (!is.null(cols)) cols[idx] else NA_character_)
}


#' Check That Node Statistics Can Be Extracted
#'
#' @description
#' Diagnostic. Fits a small tree and verifies that node-level test
#' statistics and degrees of freedom can be recovered from the fitted
#' object, which the Stack / M correction depends on. Run this once after
#' installing or after upgrading `partykit`.
#'
#' @param verbose Logical. Print a per-node report.
#'
#' @return Invisibly, `TRUE` if extraction succeeded for every inner node.
#'
#' @examples
#' check_stackM_extraction()
#'
#' @export
check_stackM_extraction <- function(verbose = TRUE) {
  set.seed(1)
  n  <- 400
  x1 <- stats::rnorm(n)
  x2 <- factor(sample(c("a", "b", "c"), n, TRUE))
  y  <- 2 * (x1 > 0) + stats::rnorm(n)
  d  <- data.frame(y = y, x1 = x1, x2 = x2)

  fit <- partykit::ctree(
    y ~ x1 + x2, data = d,
    control = partykit::ctree_control(alpha = 0.05,
                                      teststat = "quadratic",
                                      testtype = "Univariate"))

  ids   <- partykit::nodeids(fit)
  inner <- setdiff(ids, partykit::nodeids(fit, terminal = TRUE))
  if (!length(inner)) {
    message("check_stackM_extraction(): the test tree did not split; ",
            "cannot verify extraction.")
    return(invisible(FALSE))
  }

  vn <- names(partykit::data_party(fit))
  ok <- TRUE
  for (id in inner) {
    nd <- .get_node(fit, id)
    tt <- .node_split_test(nd, vn)
    if (is.null(tt)) {
      if (verbose) message(sprintf("  node %-3d EXTRACTION FAILED", id))
      ok <- FALSE; next
    }
    df  <- .recover_df(tt$statistic, tt$p_value)
    chk <- if (is.finite(df))
      stats::pchisq(tt$statistic, df, lower.tail = FALSE) else NA_real_
    good <- is.finite(df) && isTRUE(abs(chk - tt$p_value) < 1e-6)
    if (verbose) {
      message(sprintf(
        "  node %-3d var=%-6s stat=%9.3f  p=%.3e  df=%s  verified=%s",
        id, tt$variable, tt$statistic, tt$p_value,
        if (is.finite(df)) sprintf("%.2f", df) else "NA", good))
    }
    ok <- ok && good
  }

  if (ok) message("check_stackM_extraction(): OK.")
  else    message("check_stackM_extraction(): FAILED -- ",
                  "ctree_stacked() will error rather than return a tree ",
                  "corrected on unreliable statistics.")
  invisible(ok)
}


#' @keywords internal
.get_node <- function(tree, id) {
  partykit::nodeapply(tree, ids = id, FUN = function(n) n)[[1L]]
}

#' @keywords internal
.parent_map <- function(node, parent = NA_integer_, acc = NULL) {
  if (is.null(acc)) acc <- list()
  id <- partykit::id_node(node)
  acc[[as.character(id)]] <- parent
  kids <- partykit::kids_node(node)
  if (!is.null(kids)) for (k in kids) acc <- .parent_map(k, id, acc)
  acc
}


#' Prune a Stacked Tree Using the Stack / M Correction
#'
#' @description
#' Applies the Stack / M correction of Sherlock et al. (2026) to a
#' conditional inference tree fitted on stacked multiply imputed data.
#' Each inner node's test statistic is divided by `M`, the p-value is
#' recomputed, a Bonferroni adjustment for the number of candidate
#' splitting variables tested at that node is reapplied, and nodes whose
#' corrected p-value exceeds `alpha` are collapsed to terminal nodes.
#'
#' @param tree A `party`/`constparty` object fitted on stacked data.
#' @param m Integer. Number of imputations.
#' @param alpha Numeric. Nominal significance threshold.
#' @param verbose Logical. Report how many nodes were pruned.
#'
#' @return A list with the pruned `tree` and a `node_stats` data frame
#'   recording, for every inner node of the unpruned tree, the raw
#'   statistic, recovered degrees of freedom, rescaled statistic, both
#'   p-values, and whether the split was retained.
#'
#' @references
#' Sherlock, P., et al. (2026). Beyond linear risk: A machine learning
#'   approach to understanding perinatal depression in context.
#'   *Multivariate Behavioral Research*, 1-16.
#'   \doi{10.1080/00273171.2026.2661244}
#'
#' @export
prune_stackM <- function(tree, m, alpha = 0.05, verbose = TRUE) {

  ids   <- partykit::nodeids(tree)
  inner <- setdiff(ids, partykit::nodeids(tree, terminal = TRUE))

  empty <- data.frame(
    node_id = integer(0), variable = character(0), n_tested = integer(0),
    statistic = numeric(0), df = numeric(0), statistic_rescaled = numeric(0),
    p_raw = numeric(0), p_rescaled = numeric(0),
    p_rescaled_bonferroni = numeric(0), retained = logical(0),
    stringsAsFactors = FALSE)

  if (!length(inner)) return(list(tree = tree, node_stats = empty))

  vn   <- names(partykit::data_party(tree))
  rows <- vector("list", length(inner))

  for (i in seq_along(inner)) {
    id <- inner[i]
    tt <- .node_split_test(.get_node(tree, id), vn)

    if (is.null(tt)) {
      stop("ctreeMI: could not extract the test statistic at node ", id,
           ". The Stack/M correction cannot be applied. Run ",
           "check_stackM_extraction() to diagnose; this usually indicates ",
           "an incompatible partykit version.")
    }

    df <- .recover_df(tt$statistic, tt$p_value)
    if (!is.finite(df) ||
        abs(stats::pchisq(tt$statistic, df, lower.tail = FALSE) -
            tt$p_value) > 1e-6) {
      stop("ctreeMI: could not recover degrees of freedom at node ", id,
           " (statistic = ", signif(tt$statistic, 6),
           ", p = ", signif(tt$p_value, 6), "). Refusing to prune on ",
           "unverified statistics. Run check_stackM_extraction().")
    }

    rs   <- rescale_statistic(tt$statistic, m = m, df = df)
    ## Reapply the multiplicity adjustment over the variables tested here.
    padj <- 1 - (1 - rs$p_value)^tt$k

    rows[[i]] <- data.frame(
      node_id = id, variable = tt$variable, n_tested = tt$k,
      statistic = tt$statistic, df = df,
      statistic_rescaled = rs$statistic_rescaled,
      p_raw = tt$p_value, p_rescaled = rs$p_value,
      p_rescaled_bonferroni = padj, retained = padj < alpha,
      stringsAsFactors = FALSE)
  }

  node_stats <- do.call(rbind, rows)

  ## Collapse only the SHALLOWEST failing node on each path: once a parent is
  ## pruned its descendants no longer exist, and passing an orphaned id to
  ## nodeprune() would error.
  failing <- node_stats$node_id[!node_stats$retained]
  if (length(failing)) {
    pmap <- .parent_map(partykit::node_party(tree))
    has_failing_ancestor <- vapply(failing, function(id) {
      p <- pmap[[as.character(id)]]
      while (!is.na(p)) {
        if (p %in% failing) return(TRUE)
        p <- pmap[[as.character(p)]]
      }
      FALSE
    }, logical(1))
    to_prune <- failing[!has_failing_ancestor]
    if (length(to_prune)) tree <- partykit::nodeprune(tree, ids = to_prune)
  }

  if (verbose) {
    message(sprintf(
      "[ctreeMI] Stack/M correction: %d of %d splits retained (alpha = %.4f).",
      sum(node_stats$retained), nrow(node_stats), alpha))
  }

  list(tree = tree, node_stats = node_stats)
}
