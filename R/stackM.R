## ---------------------------------------------------------------------------
## The Stack / M correction.
##
## Verified against partykit 1.2-20 internals:
##
##  * With teststat = "quadratic" the node-level statistic for candidate
##    variable j is asymptotically chi-square with
##
##        df_j = dim(g(X_j)) * rank(h(Y))
##
##    dim(g(X_j)) is 1 for a numeric or ordered predictor and L - 1 for an
##    unordered factor with L levels present in that node. rank(h(Y)) is 1 per
##    numeric response component and C - 1 for an unordered factor response
##    with C levels, summed over a multivariate response. So for numeric and
##    ordered predictors df is simply the outcome dimension: 1 univariate,
##    2 bivariate, 3 for three outcomes, and so on.
##
##  * partykit stores p-values on the natural scale, and applies its
##    multiplicity adjustment as p_adj = 1 - (1 - p)^k with k the number of
##    variables it could test at the node, only when testtype = "Bonferroni".
##    Untestable variables are dropped from the stored matrix, so
##    k = ncol(criterion).
##
##  * The stored p-value underflows to exactly 0 once the statistic exceeds
##    roughly 1400 at df = 2, which happens routinely on stacked data. Nothing
##    here may depend on inverting a p-value that has underflowed.
## ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.raw_p <- function(stat, df) stats::pchisq(stat, df = df, lower.tail = FALSE)

## Adjusted p on the natural scale, computed so that it survives underflow:
##
##   p_adj = 1 - (1 - p_raw)^k        the Sidak/Bonferroni step partykit applies
##         = 1 - F^k                  where F = pchisq(stat, df) = 1 - p_raw
##         = -expm1(k * log F)
##
## log F is returned directly by pchisq(lower.tail = TRUE, log.p = TRUE), and it
## is the representation that keeps precision. Once p_raw falls below about
## 1e-17, F rounds to exactly 1, so F^k is 1 and 1 - F^k is exactly 0; whereas
## log F is an ordinary double all the way down (-5.5e-89 at stat = 400, df = 1).
## expm1 then recovers p_adj near zero without cancellation. This is the same
## scale partykit carries its criteria on internally: .extree_node() applies the
## multiplicity step as a multiplication on log(1 - p) and converts back with
## -expm1() before storing.
#' @keywords internal
#' @noRd
.sidak <- function(stat, df, k) {
  if (k <= 1L) return(.raw_p(stat, df))
  -expm1(k * stats::pchisq(stat, df = df, lower.tail = TRUE, log.p = TRUE))
}


# ---------------------------------------------------------------------------
# Degrees of freedom
# ---------------------------------------------------------------------------

#' Recover the degrees of freedom partykit used for a node-level test
#'
#' The stored pair (statistic, p) satisfies p = 1 - (1 - pchisq(statistic, df,
#' lower.tail = FALSE))^k with k known from the fitted control, so df is
#' identified whenever p is interior. df is an integer by construction, so the
#' search is over integers and the match is verified rather than assumed: an
#' inexact match returns NA and the caller falls back on the structural rule.
#'
#' @keywords internal
#' @noRd
.recover_df <- function(statistic, p, k = 1L, df_max = 200L, tol = 1e-6) {
  if (!is.finite(statistic) || !is.finite(p) ||
      statistic <= 0 || p <= 0 || p >= 1) return(NA_real_)
  grid <- seq_len(max(1L, as.integer(df_max)))
  lrec <- log(.sidak(statistic, grid, k))
  lp   <- log(p)
  i    <- which.min(abs(lrec - lp))
  if (!length(i) || !is.finite(lrec[i]) ||
      abs(lrec[i] - lp) > tol * max(1, abs(lp))) return(NA_real_)
  as.numeric(grid[i])
}

#' Rank of the influence function of the response
#' @keywords internal
#' @noRd
.response_rank <- function(tree) {
  y <- try(tree$fitted[["(response)"]], silent = TRUE)
  if (inherits(y, "try-error") || is.null(y)) return(NA_integer_)
  cols <- if (is.data.frame(y)) as.list(y) else list(y)
  sum(vapply(cols, function(v) {
    if (is.factor(v) && !is.ordered(v)) max(nlevels(droplevels(v)) - 1L, 1L) else 1L
  }, integer(1)))
}

#' Structural degrees of freedom, used when the p-value has underflowed
#' @keywords internal
#' @noRd
.structural_df <- function(tree, id, vars, q) {
  na <- stats::setNames(rep(NA_real_, length(vars)), vars)
  if (is.na(q)) return(na)
  nd <- try(partykit::data_party(tree, id), silent = TRUE)
  if (inherits(nd, "try-error")) return(na)
  vapply(vars, function(v) {
    if (!v %in% names(nd)) return(NA_real_)
    x <- nd[[v]]
    dim_x <- if (is.factor(x) && !is.ordered(x)) max(nlevels(droplevels(x)) - 1L, 1L) else 1L
    as.numeric(dim_x * q)
  }, numeric(1), USE.NAMES = TRUE)
}

#' Upper bound for the df search grid, taken from the data
#' @keywords internal
#' @noRd
.df_upper <- function(tree, q) {
  if (is.na(q) || !is.finite(q)) q <- 1L
  mf <- try(tree$data, silent = TRUE)
  lmax <- 1L
  if (!inherits(mf, "try-error") && is.data.frame(mf) && ncol(mf) > 0L) {
    lv <- vapply(mf, function(v) if (is.factor(v) && !is.ordered(v)) nlevels(v) else 2L,
                 integer(1))
    lmax <- max(c(1L, lv - 1L))
  }
  as.integer(max(200L, ceiling(1.25 * q * lmax)))
}


# ---------------------------------------------------------------------------
# Node-level extraction
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.get_node <- function(tree, id) {
  partykit::nodeapply(tree, ids = id, FUN = function(n) n)[[1L]]
}

#' @keywords internal
#' @noRd
.parent_map <- function(node, parent = NA_integer_, acc = NULL) {
  if (is.null(acc)) acc <- list()
  id <- partykit::id_node(node)
  acc[[as.character(id)]] <- parent
  kids <- partykit::kids_node(node)
  if (!is.null(kids)) for (k in kids) acc <- .parent_map(k, id, acc)
  acc
}

#' @keywords internal
#' @noRd
.kid_map <- function(node, acc = NULL) {
  if (is.null(acc)) acc <- list()
  id   <- partykit::id_node(node)
  kids <- partykit::kids_node(node)
  acc[[as.character(id)]] <-
    if (length(kids)) vapply(kids, partykit::id_node, integer(1)) else integer(0)
  if (length(kids)) for (k in kids) acc <- .kid_map(k, acc)
  acc
}

#' Every candidate variable tested at a node, not only the one split on
#'
#' ctree's own rule is to split when *any* candidate variable meets alpha, and
#' that is the rule the correction has to reapply: rescaling can reorder the
#' candidates, because different variables carry different degrees of freedom.
#'
#' @keywords internal
#' @noRd
.node_tests <- function(node, var_names) {
  crit <- partykit::info_node(node)$criterion
  if (is.null(crit) || !is.matrix(crit)) return(NULL)
  if (!all(c("statistic", "p.value") %in% rownames(crit))) return(NULL)

  stat <- as.numeric(crit["statistic", ])
  p    <- as.numeric(crit["p.value", ])
  cols <- colnames(crit)
  keep <- is.finite(stat) & !is.na(p)
  if (!any(keep)) return(NULL)

  ## The variable ctree actually split on, for reporting.
  sp  <- partykit::split_node(node)
  sel <- NA_character_
  if (!is.null(sp) && !is.null(cols)) {
    vid <- try(as.integer(partykit::varid_split(sp)), silent = TRUE)
    if (!inherits(vid, "try-error") && !is.na(vid) && vid <= length(var_names))
      sel <- var_names[vid]
  }
  if (is.na(sel)) sel <- cols[keep][which.min(p[keep])]

  list(variable  = cols[keep],
       statistic = stat[keep],
       p_stored  = p[keep],
       k         = sum(keep),
       selected  = sel)
}


# ---------------------------------------------------------------------------
# Exported
# ---------------------------------------------------------------------------

#' Rescale a Node-Level Test Statistic for Stacking
#'
#' @description
#' Divides a node-level chi-square test statistic by the number of imputations
#' `M` and recomputes the corresponding p-value. This is the Stack / M
#' correction of Sherlock et al. (2026).
#'
#' @details
#' Stacking `M` imputed datasets produces `M * n` rows, and a chi-square
#' statistic computed on the stacked data scales approximately linearly with
#' `M`. Dividing the statistic by `M` returns it to the scale of a single
#' imputed dataset before the p-value is computed.
#'
#' This is **not** the same as dividing the significance threshold by `M`. The
#' two rules coincide only at `M = 1`:
#'
#' \itemize{
#'   \item Statistic rescaling rejects when `X > M * qchisq(1 - alpha, df)`.
#'   \item Threshold rescaling rejects when `X > qchisq(1 - alpha / M, df)`.
#' }
#'
#' For `df = 1`, `alpha = 0.05` and `M = 30`, the first requires `X > 115.2`
#' and the second only `X > 9.9`. Versions 0.1.0 and 0.2.0 implemented the
#' second rule.
#'
#' @param statistic Numeric. The raw test statistic computed on the stacked
#'   data.
#' @param m Integer. Number of imputations.
#' @param df Numeric. Degrees of freedom of the reference chi-square
#'   distribution: the outcome dimension for a numeric or ordered predictor,
#'   `(L - 1)` times the outcome dimension for an unordered factor with `L`
#'   levels. See [prune_stackM()], which derives this for you.
#'
#' @return A list with `statistic_rescaled` and `p_value`.
#'
#' @references
#' Sherlock, P., Mansolf, M., Hofheimer, J., Hockett, C. W., O'Connor, T. G.,
#'   Roubinov, D., Graff, J. C., Lai, J.-S., Bush, N. R., Wright, R. J., &
#'   Chiu, Y.-H. M. (2026). Beyond linear risk: A machine learning approach to
#'   understanding perinatal depression in context.
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
  list(statistic_rescaled = s, p_value = .raw_p(s, df))
}


#' Check That Node Statistics Can Be Extracted
#'
#' @description
#' Diagnostic. Fits small trees under both of the multiplicity settings the
#' correction supports, and verifies that node-level statistics and degrees of
#' freedom can be recovered from the fitted objects. Run it once after
#' installing, or after upgrading `partykit`.
#'
#' @param verbose Logical. Print a per-node report.
#'
#' @return Invisibly, `TRUE` if extraction succeeded everywhere it was tried.
#'
#' @examples
#' check_stackM_extraction(verbose = FALSE)
#'
#' @export
check_stackM_extraction <- function(verbose = TRUE) {

  set.seed(1)
  n  <- 400
  x1 <- stats::rnorm(n)
  x2 <- factor(sample(c("a", "b", "c"), n, TRUE))
  d  <- data.frame(y1 = 2 * (x1 > 0) + stats::rnorm(n),
                   y2 = 1.2 * (x1 > 0) + stats::rnorm(n),
                   x1 = x1, x2 = x2)

  ok <- TRUE
  for (tt in c("Univariate", "Bonferroni")) {
    fit <- partykit::ctree(
      y1 + y2 ~ x1 + x2, data = d,
      control = partykit::ctree_control(alpha = 0.05, teststat = "quadratic",
                                        testtype = tt))
    q     <- .response_rank(fit)
    dfmax <- .df_upper(fit, q)
    bonf  <- isTRUE(fit$info$control$bonferroni)
    ids   <- partykit::nodeids(fit)
    inner <- setdiff(ids, partykit::nodeids(fit, terminal = TRUE))
    if (!length(inner)) {
      if (verbose) message("  testtype = ", tt, ": test tree did not split.")
      next
    }
    vn <- names(partykit::data_party(fit))
    for (id in inner) {
      nt <- .node_tests(.get_node(fit, id), vn)
      if (is.null(nt)) {
        if (verbose) message(sprintf("  [%s] node %-3d EXTRACTION FAILED", tt, id))
        ok <- FALSE; next
      }
      k_stored <- if (bonf) nt$k else 1L
      rec <- vapply(seq_along(nt$variable), function(j)
        .recover_df(nt$statistic[j], nt$p_stored[j], k_stored, dfmax), numeric(1))
      str <- .structural_df(fit, id, nt$variable, q)
      cmp <- is.finite(rec) & is.finite(str)
      good <- !any(cmp) || all(rec[cmp] == str[cmp])
      if (verbose)
        message(sprintf(
          "  [%-10s] node %-3d k=%d  recovered df: %s  structural: %s  agree=%s",
          tt, id, nt$k,
          paste(ifelse(is.finite(rec), format(rec), "NA"), collapse = ","),
          paste(ifelse(is.finite(str), format(str), "NA"), collapse = ","),
          good))
      ok <- ok && good
    }
  }

  if (ok) message("check_stackM_extraction(): OK.")
  else    message("check_stackM_extraction(): FAILED -- the recovered and ",
                  "structural degrees of freedom disagree. Please report this ",
                  "at https://github.com/Phillip-Sherlock/ctreeMI/issues .")
  invisible(ok)
}


#' Prune a Stacked Tree Using the Stack / M Correction
#'
#' @description
#' Applies the Stack / M correction of Sherlock et al. (2026) to a conditional
#' inference tree fitted on stacked multiply imputed data. Every node-level
#' test statistic is divided by `M`, the p-values are recomputed from the
#' chi-square reference distribution `ctree` uses, the multiplicity adjustment
#' over candidate splitting variables is reapplied, and the tree is compressed
#' bottom-up.
#'
#' @details
#' # What is tested
#'
#' `ctree` splits a node when *any* candidate variable meets `alpha`, so that
#' is the rule reapplied here: all candidates tested at the node are rescaled
#' and the smallest corrected p-value decides. Testing only the variable that
#' was split on would not be equivalent, because candidates carry different
#' degrees of freedom and rescaling can reorder them.
#'
#' # How the tree is compressed
#'
#' Pruning is bottom-up, as in the analysis for the paper. An internal node is
#' collapsed only once all of its own internal descendants have been collapsed,
#' so a node whose descendant survives the correction keeps its split.
#'
#' # Degrees of freedom
#'
#' With `teststat = "quadratic"` the node-level statistic is chi-square with
#' `df` equal to the rank of the covariance matrix of the linear statistic:
#' the outcome dimension for a numeric or ordered predictor (1 for a
#' univariate outcome, 2 for a bivariate outcome, and so on), and `(L - 1)`
#' times the outcome dimension for an unordered factor with `L` levels present
#' in the node. These are derived for every node and every candidate variable,
#' by inverting the (statistic, p-value) pairs `partykit` stored and falling
#' back on the structural rule where the stored p-value has underflowed to
#' zero, which happens routinely on stacked data.
#'
#' @param tree A `party`/`constparty` object fitted by [partykit::ctree()] on
#'   stacked data, with `teststat = "quadratic"` and `testtype` either
#'   `"Univariate"` or `"Bonferroni"`. Both are handled: whether the stored
#'   p-values already carry `partykit`'s multiplicity adjustment is read from
#'   the fitted control rather than assumed.
#' @param m Integer. Number of imputations.
#' @param alpha Numeric. Nominal significance threshold.
#' @param verbose Logical. Report how many splits were retained.
#'
#' @return A list with the pruned `tree` and a `node_stats` data frame with one
#'   row per internal node of the unpruned tree, recording the variable split
#'   on, its raw and rescaled statistics, the degrees of freedom and how they
#'   were obtained, the p-values before and after correction, the smallest
#'   corrected p-value over all candidates at that node, and whether the split
#'   was retained. The full per-candidate table is attached as
#'   `attr(node_stats, "candidates")`.
#'
#' @references
#' Sherlock, P., et al. (2026). Beyond linear risk: A machine learning
#'   approach to understanding perinatal depression in context.
#'   *Multivariate Behavioral Research*, 1-16.
#'   \doi{10.1080/00273171.2026.2661244}
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' d <- data.frame(x = stats::rnorm(n))
#' d$y1 <- stats::rnorm(n) + 1.2 * (d$x > 0)
#' d$y2 <- stats::rnorm(n) + 0.9 * (d$x > 0)
#' grown <- partykit::ctree(y1 + y2 ~ x, data = d)
#' out <- prune_stackM(grown, m = 5, verbose = FALSE)
#' out$node_stats
#'
#' @export
prune_stackM <- function(tree, m, alpha = 0.05, verbose = TRUE) {
  .prune_stackM(tree, m = m, alpha = alpha, df = NULL, verbose = verbose)
}

## Internal worker. The `df` argument is used by the test suite to hold the
## degrees of freedom constant while checking the pruning logic in isolation.
#' @keywords internal
#' @noRd
.prune_stackM <- function(tree, m, alpha = 0.05, df = NULL, verbose = TRUE) {

  if (!inherits(tree, "party"))
    stop("`tree` must be a party object from partykit::ctree().")
  if (!is.numeric(m) || length(m) != 1L || !is.finite(m) || m < 1)
    stop("`m` must be a single number >= 1.")
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be strictly between 0 and 1.")
  if (!is.null(df) && !identical(df, "outcome") && !is.numeric(df))
    stop("`df` must be NULL, \"outcome\", or a single positive number.")

  ctrl <- tree$info$control
  if (!is.null(ctrl)) {
    if (!identical(as.character(ctrl$teststat)[1L], "quadratic"))
      stop("The Stack/M correction rescales a chi-square statistic and so ",
           "requires teststat = \"quadratic\"; this tree used \"",
           as.character(ctrl$teststat)[1L], "\".")
    if (isTRUE(ctrl$splittest))
      stop("The Stack/M correction is not defined for splittest = TRUE, ",
           "which uses a Monte Carlo reference distribution.")
    tt <- as.character(ctrl$testtype)[1L]
    if (!tt %in% c("Bonferroni", "Univariate"))
      stop("The Stack/M correction requires testtype = \"Bonferroni\" or ",
           "\"Univariate\"; this tree used \"", tt, "\".")
  }
  bonf <- isTRUE(ctrl$bonferroni)

  ids   <- partykit::nodeids(tree)
  inner <- setdiff(ids, partykit::nodeids(tree, terminal = TRUE))

  empty_cand <- data.frame(
    node_id = integer(0), variable = character(0), statistic = numeric(0),
    df = numeric(0), df_source = character(0), statistic_rescaled = numeric(0),
    p_stored = numeric(0), p_corrected = numeric(0), stringsAsFactors = FALSE)
  empty <- data.frame(
    node_id = integer(0), variable = character(0), n_tested = integer(0),
    statistic = numeric(0), df = numeric(0), df_source = character(0),
    statistic_rescaled = numeric(0), p_raw = numeric(0),
    p_rescaled = numeric(0), p_rescaled_bonferroni = numeric(0),
    p_min = numeric(0), variable_min = character(0), retained = logical(0),
    stringsAsFactors = FALSE)
  attr(empty, "candidates") <- empty_cand

  if (!length(inner)) return(list(tree = tree, node_stats = empty))

  q     <- .response_rank(tree)
  dfmax <- .df_upper(tree, q)
  vn    <- names(partykit::data_party(tree))
  rows  <- vector("list", length(inner))
  cands <- vector("list", length(inner))

  for (i in seq_along(inner)) {
    id <- inner[i]
    nt <- .node_tests(.get_node(tree, id), vn)
    if (is.null(nt))
      stop("ctreeMI: could not extract the test statistics at node ", id,
           ". Run check_stackM_extraction() to diagnose; this usually means ",
           "an incompatible partykit version, or a tree fitted with ",
           "ctree_control(saveinfo = FALSE).")

    nv       <- length(nt$variable)
    k_stored <- if (bonf) nt$k else 1L

    ## ---- degrees of freedom, per candidate variable ----------------------
    if (is.numeric(df)) {
      dfv <- rep(as.numeric(df[1L]), nv); src <- rep("fixed", nv)
    } else if (identical(df, "outcome")) {
      if (is.na(q))
        stop("ctreeMI: could not determine the outcome dimension; supply `df`.")
      dfv <- rep(as.numeric(q), nv); src <- rep("outcome", nv)
    } else {
      dfv <- vapply(seq_len(nv), function(j)
        .recover_df(nt$statistic[j], nt$p_stored[j], k_stored, dfmax), numeric(1))
      src <- ifelse(is.finite(dfv), "recovered", NA_character_)
      if (anyNA(dfv)) {
        str <- .structural_df(tree, id, nt$variable, q)
        take <- is.na(dfv) & is.finite(str)
        dfv[take] <- str[take]; src[take] <- "structural"
      }
      if (anyNA(dfv) && !is.na(q)) {
        take <- is.na(dfv); dfv[take] <- as.numeric(q); src[take] <- "outcome"
      }
      if (anyNA(dfv))
        stop("ctreeMI: could not determine the degrees of freedom at node ", id,
             ". Supply `df` explicitly, or run check_stackM_extraction().")
    }

    ## ---- rescale and reapply the multiplicity adjustment -----------------
    stat_m <- nt$statistic / m
    p_corr <- .sidak(stat_m, dfv, nt$k)

    cands[[i]] <- data.frame(
      node_id = id, variable = nt$variable, statistic = nt$statistic,
      df = dfv, df_source = src, statistic_rescaled = stat_m,
      p_stored = nt$p_stored, p_corrected = p_corr,
      stringsAsFactors = FALSE)

    jmin <- which.min(p_corr)
    jsel <- match(nt$selected, nt$variable)
    if (is.na(jsel)) jsel <- jmin

    rows[[i]] <- data.frame(
      node_id = id, variable = nt$variable[jsel], n_tested = nt$k,
      statistic = nt$statistic[jsel], df = dfv[jsel], df_source = src[jsel],
      statistic_rescaled = stat_m[jsel],
      p_raw = .raw_p(nt$statistic[jsel], dfv[jsel]),
      p_rescaled = .raw_p(stat_m[jsel], dfv[jsel]),
      p_rescaled_bonferroni = p_corr[jsel],
      p_min = p_corr[jmin], variable_min = nt$variable[jmin],
      retained = p_corr[jmin] < alpha,
      stringsAsFactors = FALSE)
  }

  node_stats <- do.call(rbind, rows)
  cand_tab   <- do.call(rbind, cands)

  ## ---- bottom-up compression -------------------------------------------
  ## keep(node) = passes(node) OR any(keep(child)). A node whose descendant
  ## survives is never reached by the bottom-up sweep, so it keeps its split.
  kids   <- .kid_map(partykit::node_party(tree))
  passes <- stats::setNames(node_stats$retained, as.character(node_stats$node_id))
  keep   <- new.env(parent = emptyenv())
  resolve <- function(id) {
    key <- as.character(id)
    kd  <- kids[[key]]
    if (!length(kd)) return(FALSE)                       # terminal node
    kept_below <- any(vapply(kd, resolve, logical(1)))
    val <- isTRUE(unname(passes[key])) || kept_below
    assign(key, val, envir = keep)
    val
  }
  resolve(partykit::id_node(partykit::node_party(tree)))

  keep_id <- vapply(inner, function(id)
    isTRUE(get0(as.character(id), envir = keep, ifnotfound = FALSE)), logical(1))
  node_stats$retained <- keep_id
  attr(node_stats, "candidates") <- cand_tab

  drop <- inner[!keep_id]
  if (length(drop)) {
    pmap <- .parent_map(partykit::node_party(tree))
    topmost <- vapply(drop, function(id) {
      p <- pmap[[as.character(id)]]
      while (!is.na(p)) { if (p %in% drop) return(FALSE); p <- pmap[[as.character(p)]] }
      TRUE
    }, logical(1))
    tree <- partykit::nodeprune(tree, ids = drop[topmost])
  }

  if (verbose)
    message(sprintf(
      "[ctreeMI] Stack/M correction: %d of %d splits retained (alpha = %.4f).",
      sum(keep_id), length(keep_id), alpha))

  list(tree = tree, node_stats = node_stats)
}
