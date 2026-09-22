#' Split Data Into Discovery and Confirmation Sets
#'
#' Partitions a data frame into two disjoint subsets of original
#' observations, before imputation. Each subset should then be imputed
#' separately, so that neither half's outcomes inform the other's imputed
#' predictor values.
#'
#' @param data A data frame, typically containing missing values.
#' @param prop Proportion of observations assigned to the discovery set.
#'   Default 0.5.
#' @param strata Optional name of a column on which to stratify the split,
#'   so that its distribution is preserved in both halves. Ordinarily the
#'   outcome, if categorical. Rows with a missing value in this column are
#'   allocated at random. When `cluster` is given, `strata` must be constant
#'   within each cluster, and the clusters are stratified.
#' @param cluster Optional name of a column identifying a higher-level unit
#'   such as a school, district or site. When given, whole clusters are
#'   assigned to one half or the other, so that no cluster contributes
#'   observations to both. `prop` then refers to the proportion of clusters.
#' @param seed Optional integer seed.
#'
#' @details
#' The split is performed on original observations, not on stacked rows,
#' and before any imputation. This matters. Imputing the full data once and
#' then splitting leaves each half's imputed values informed by the other
#' half's outcomes, which compromises the independence the confirmation
#' step relies on. Impute the two halves separately.
#'
#' ## Clustered data
#'
#' When observations are nested in higher-level units, splitting rows at
#' random places members of the same unit on both sides, and the
#' confirmation sample is then not independent of the discovery sample in
#' the way the test assumes. Passing the unit identifier as `cluster` keeps
#' every unit intact on one side. The confirmation test in
#' [confirm_ctreeMI()] should then use a cluster-robust variance by passing
#' the same `cluster` argument there.
#'
#' With few clusters the two halves may be unequal in rows even when equal
#' in clusters. The returned `n_clusters` reports the count on each side.
#'
#' @return A list of class `"ctreeMI_split"` with elements `discover` and
#'   `confirm`, each a data frame; `index`, the row indices of `data`
#'   assigned to the discovery set; `cluster`, the cluster column name or
#'   `NULL`; and, when clustered, `n_clusters` for each half.
#'
#' @seealso [confirm_ctreeMI()], [discover_confirm()]
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(200), y = rnorm(200))
#' d$x[sample(200, 30)] <- NA
#' parts <- split_holdout(d, prop = 0.5, seed = 1)
#' nrow(parts$discover); nrow(parts$confirm)
#' @export
split_holdout <- function(data, prop = 0.5, strata = NULL, cluster = NULL,
                          seed = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  if (!is.numeric(prop) || prop <= 0 || prop >= 1)
    stop("`prop` must be strictly between 0 and 1.")
  n <- nrow(data)
  if (n < 4L) stop("Too few observations to split.")
  if (!is.null(seed)) set.seed(seed)

  if (!is.null(cluster)) {
    # ---- sample whole clusters, so a level-2 unit never straddles the split
    if (!cluster %in% names(data))
      stop("`cluster` column '", cluster, "' not found in `data`.")
    cl <- data[[cluster]]
    if (anyNA(cl)) stop("`cluster` contains missing values; every row must belong to a cluster.")
    cl    <- as.character(cl)
    units <- unique(cl)
    if (length(units) < 4L) stop("Fewer than four clusters; cannot split by cluster.")
    if (is.null(strata)) {
      k   <- max(1L, min(length(units) - 1L, round(prop * length(units))))
      sel <- sample(units, size = k)
    } else {
      if (!strata %in% names(data))
        stop("`strata` column '", strata, "' not found in `data`.")
      # stratification must be at the cluster level
      cs <- unique(data.frame(cl = cl, st = as.character(data[[strata]]),
                              stringsAsFactors = FALSE))
      if (nrow(cs) != length(units))
        stop("`strata` varies within `cluster`; when splitting by cluster the ",
             "stratifying variable must be constant within each cluster.")
      cs$st[is.na(cs$st)] <- ".NA."
      sel <- unlist(lapply(split(cs$cl, cs$st), function(u) {
        k <- round(prop * length(u))
        if (length(u) == 1L) return(if (stats::runif(1) < prop) u else character(0))
        sample(u, size = k)
      }), use.names = FALSE)
    }
    idx <- which(cl %in% sel)
  } else if (is.null(strata)) {
    idx <- sort(sample.int(n, size = round(prop * n)))
  } else {
    if (!strata %in% names(data))
      stop("`strata` column '", strata, "' not found in `data`.")
    g   <- data[[strata]]
    g   <- as.character(g)
    g[is.na(g)] <- ".NA."
    idx <- unlist(lapply(split(seq_len(n), g), function(rows) {
      k <- round(prop * length(rows))
      if (length(rows) == 1L) return(if (stats::runif(1) < prop) rows else integer(0))
      sort(sample(rows, size = k))
    }), use.names = FALSE)
    idx <- sort(idx)
  }
  if (length(idx) == 0L || length(idx) == n)
    stop("Split produced an empty half; adjust `prop` or `strata`.")

  out <- list(discover = data[idx, , drop = FALSE],
              confirm  = data[-idx, , drop = FALSE],
              index    = idx,
              cluster  = cluster)
  if (!is.null(cluster)) {
    out$n_clusters <- c(discover = length(unique(data[[cluster]][idx])),
                        confirm  = length(unique(data[[cluster]][-idx])))
  }
  structure(out, class = "ctreeMI_split")
}


#' Confirm a Discovered Tree on Independent Data
#'
#' Applies the terminal-node rules of a tree fitted by [ctree_stacked()] to
#' an independently imputed confirmation sample, and tests whether the
#' outcome differs across those subgroups using a test with valid error
#' control under multiple imputation.
#'
#' @param tree An object of class `ctreeMI`, fitted on the discovery data.
#' @param data The confirmation data as a `mids` object from
#'   [mice::mice()] or a list of imputed data frames. Must have been
#'   imputed separately from the discovery data; see [split_holdout()].
#' @param outcome_type Either `"auto"` (default), in which case numeric
#'   outcomes are modelled by linear regression and two-level factors by
#'   logistic regression, or a character vector of the same length as the
#'   number of outcomes with entries `"continuous"` or `"binary"`.
#' @param min_node Minimum number of confirmation observations, per
#'   imputation, that a terminal node must receive. Nodes below this are
#'   reported but excluded from the test. Default 5.
#' @param adjust Method for adjusting the per-split p-values across the
#'   internal nodes of the tree, passed to [stats::p.adjust()]. Default
#'   `"holm"`. Use `"none"` to report unadjusted values.
#' @param conf.level Confidence level for the pooled contrast reported with
#'   each split. Default 0.95.
#' @param cluster Optional name of a column in the confirmation data
#'   identifying a higher-level unit such as a school or district. When
#'   given, every model's covariance matrix is replaced by a cluster-robust
#'   estimate from [sandwich::vcovCL()], so that within-cluster correlation
#'   is accounted for in every test and interval. Requires the `sandwich`
#'   package.
#' @param vcov_type The small-sample adjustment passed to
#'   [sandwich::vcovCL()] when `cluster` is given. Default `"HC1"`.
#'
#' @details
#' ## Why a separate confirmation step
#'
#' A tree uses the data to choose which variables to split on and where.
#' Testing the resulting subgroups on the same data is invalid regardless of
#' how the missing values were handled, because the partition was selected
#' to maximise exactly the separation being tested. Independently of that,
#' the node-level test of [ctree_stacked()] is not calibrated under
#' outcome-conditioned imputation; see the "Node-level calibration" section
#' of that help page.
#'
#' Both problems are avoided by treating the discovery tree as a hypothesis
#' about subgroup structure and testing it on data the tree has not seen.
#' This function implements that step. The subgroups are fixed by the
#' discovery tree, so the multiplicity of the search is not carried into
#' the test, and the confirmation test pools across imputations by Rubin's
#' rules rather than stacking, so the between-imputation variance that
#' miscalibrates the stacked test is properly accounted for.
#'
#' ## What is tested
#'
#' Two families of test are run, both pooled across imputations by the
#' multivariate Wald procedure of Li, Raghunathan and Rubin (1991). The
#' pooling is implemented within this package rather than through
#' [mice::D1()] so that a cluster-robust covariance can be supplied; with
#' `cluster = NULL` the two agree.
#'
#' The **omnibus test** fits, for each outcome, a model with terminal-node
#' membership as the only predictor, and tests whether the outcome differs
#' across nodes at all. This establishes that the partition carries
#' information, not that every split is real: a tree with one genuine split
#' and one spurious one will still yield a small \emph{p} value.
#'
#' The **per-split tests** address that directly. For each internal node of
#' the discovery tree, the confirmation observations falling within that
#' node are divided by the node's own split rule, and the outcome is tested
#' for a difference between the resulting children. Because the parent node
#' and its rule are fixed by the discovery tree, each is a pre-specified
#' contrast, and the test is valid. A split whose children do not differ in
#' independent data is one the discovery tree should not have made. The
#' per-split \emph{p} values are adjusted across internal nodes by the
#' method in `adjust`.
#'
#' Each split is also reported with the pooled difference between its
#' children and a confidence interval, on the scale of the outcome for a
#' continuous response and as a log odds ratio for a binary one. These are
#' the quantities to report: the p-value says whether the split survived,
#' the contrast says how much the subgroups differ and how precisely that is
#' known. [prune_unconfirmed()] collapses the splits that did not survive.
#'
#' A split below a non-confirmed split has no clear interpretation, since
#' the partition it refines was not itself supported; the `depth` column
#' allows this to be read off. Pooled node-level estimates are returned so
#' that the pattern of differences can be examined directly.
#'
#' ## Clustered data
#'
#' When observations are nested in units, the members of a unit are not
#' independent and a model-based variance understates uncertainty. Passing
#' the unit identifier as `cluster` substitutes a cluster-robust covariance
#' for every fitted model before pooling, so that the omnibus test, the
#' per-split tests, the contrasts and the node estimates all reflect the
#' effective sample size. The discovery split should have kept units
#' intact; see [split_holdout()].
#'
#' The reference distributions carry the usual small-sample adjustments: the
#' Reiter (2007) adjustment for the multivariate test and the Barnard-Rubin
#' (1999) adjustment for the scalar contrasts, both of which cap the
#' degrees of freedom by the complete-data degrees of freedom. Unclustered,
#' that is the residual degrees of freedom, and the tests then agree with
#' [mice::D1()] and [mice::pool()]. Clustered, the complete-data degrees of
#' freedom are taken as the number of clusters minus one, which is the
#' conventional choice for cluster-robust inference and is what makes the
#' tests more conservative as clusters become few.
#'
#' Cluster-robust variance is nonetheless unreliable with very few clusters,
#' and a warning is issued below twenty. The number of clusters behind each
#' node and each split is reported so that this can be judged.
#'
#' ## Node assignment under imputation
#'
#' Each imputed confirmation dataset is passed through the discovery tree's
#' split rules. Because imputed predictor values vary across imputations,
#' an observation whose imputed value straddles a split point may be
#' assigned to different nodes in different imputations. This is not an
#' error; it reflects imputation uncertainty, and pooling by Rubin's rules
#' accounts for it.
#'
#' @return An object of class `"ctreeMI_confirm"`: a list with elements
#'   \describe{
#'     \item{`test`}{A data frame with one row per outcome: `outcome`,
#'       `F`, `df1`, `df2`, `p`, and `riv`, the relative increase in
#'       variance due to nonresponse.}
#'     \item{`splits`}{A data frame with one row per internal node and
#'       outcome: `node_id`, `depth`, `split_var`, `rule` (the split as
#'       written), `outcome`, `contrast` with `lower` and `upper` (the
#'       pooled difference between the children and its confidence
#'       interval), `F`, `df1`, `df2`, `p`, `p_adj`, and `n_confirm` (mean
#'       observations in the parent node per imputation). `NA` where a child
#'       node fell below `min_node` in some imputation.}
#'     \item{`nodes`}{A data frame of pooled per-node estimates: `node_id`,
#'       `outcome`, `estimate` (mean or proportion), `se`, `n_confirm`
#'       (mean observations per imputation), `n_clusters`, and `tested`.}
#'     \item{`m`}{Number of imputations used.}
#'     \item{`n_confirm`}{Number of confirmation observations.}
#'     \item{`cluster`, `n_clusters`, `vcov_type`}{The cluster column,
#'       the number of clusters in the confirmation data, and the variance
#'       estimator used (`"model"` when unclustered).}
#'     \item{`excluded`}{Terminal-node ids excluded for having fewer than
#'       `min_node` observations, if any.}
#'   }
#'
#' @references
#' Barnard, J., and Rubin, D. B. (1999). Small-sample degrees of freedom
#' with multiple imputation. \emph{Biometrika}, 86, 948--955.
#'
#' Li, K. H., Raghunathan, T. E., and Rubin, D. B. (1991). Large-sample
#' significance levels from multiply imputed data using moment-based
#' statistics and an F reference distribution. \emph{Journal of the
#' American Statistical Association}, 86, 1065--1073.
#'
#' Reiter, J. P. (2007). Small-sample degrees of freedom for multi-component
#' significance tests with multiple imputation for missing data.
#' \emph{Biometrika}, 94, 502--508.
#'
#' @seealso [split_holdout()], [discover_confirm()], [ctree_stacked()]
#' @examples
#' \dontrun{
#' library(mice)
#' set.seed(7)
#' n <- 600
#' d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' d$y <- rnorm(n) + 1.2 * (d$x1 > 0)
#' d$x1[sample(n, 90)] <- NA
#'
#' parts <- split_holdout(d, seed = 7)
#' imp_d <- mice(parts$discover, m = 20, printFlag = FALSE, seed = 1)
#' imp_c <- mice(parts$confirm,  m = 20, printFlag = FALSE, seed = 2)
#'
#' tree <- ctree_stacked(y ~ x1 + x2, data = imp_d, verbose = FALSE)
#' confirm_ctreeMI(tree, imp_c)
#' }
#' @export
confirm_ctreeMI <- function(tree, data, outcome_type = "auto", min_node = 5L,
                            adjust = "holm", conf.level = 0.95,
                            cluster = NULL, vcov_type = "HC1") {
  if (!inherits(tree, "ctreeMI"))
    stop("`tree` must be a 'ctreeMI' object from ctree_stacked().")
  info <- attr(tree, "ctreeMI_info")
  if (is.null(info$formula)) stop("Tree carries no formula; cannot identify the outcome.")

  imps <- .as_imp_list(data)
  M    <- length(imps)
  if (M < 2L) stop("Confirmation data must contain at least two imputations.")

  outcomes <- all.vars(info$formula[[2L]])
  if (!all(outcomes %in% names(imps[[1L]])))
    stop("Outcome(s) ", paste(setdiff(outcomes, names(imps[[1L]])), collapse = ", "),
         " not found in confirmation data.")

  # ---- clustering ---------------------------------------------------------
  if (!is.null(cluster)) {
    if (!requireNamespace("sandwich", quietly = TRUE))
      stop("Package 'sandwich' is required for cluster-robust variance; ",
           "install it with install.packages(\"sandwich\").")
    if (!cluster %in% names(imps[[1L]]))
      stop("`cluster` column '", cluster, "' not found in confirmation data.")
    n_cl <- length(unique(imps[[1L]][[cluster]]))
    if (n_cl < 20L)
      warning("Only ", n_cl, " clusters in the confirmation data. Cluster-robust ",
              "variance is unreliable with few clusters, and the reference ",
              "distributions, whose degrees of freedom are bounded by the number ",
              "of clusters, will have little power. Treat p-values as approximate ",
              "and rely on the contrasts and their intervals.", call. = FALSE)
  } else {
    n_cl <- NA_integer_
  }

  term <- partykit::nodeids(tree, terminal = TRUE)
  if (length(term) < 2L)
    stop("The discovery tree has a single terminal node; there is no partition to confirm.")

  types <- .resolve_types(imps[[1L]], outcomes, outcome_type)

  # ---- node assignment per imputation -------------------------------------
  assigned <- lapply(imps, function(d)
    factor(stats::predict(tree, newdata = d, type = "node"), levels = term))

  counts <- vapply(assigned, function(a) as.numeric(table(a)), numeric(length(term)))
  dimnames(counts) <- list(as.character(term), NULL)
  size <- rowMeans(counts)
  keep <- rownames(counts)[apply(counts, 1L, min) >= min_node]
  drop <- setdiff(as.character(term), keep)
  if (length(keep) < 2L)
    stop("Fewer than two terminal nodes receive at least `min_node` confirmation observations.")

  # clusters per terminal node (mean over imputations)
  node_cl <- if (is.null(cluster)) rep(NA_real_, length(term)) else
    vapply(as.character(term), function(nd) mean(vapply(seq_len(M), function(i)
      length(unique(imps[[i]][[cluster]][assigned[[i]] == nd])), numeric(1L))),
      numeric(1L))

  # ---- omnibus test + pooled node estimates, per outcome --------------------
  test_rows <- vector("list", length(outcomes))
  node_rows <- vector("list", length(outcomes))

  for (j in seq_along(outcomes)) {
    y    <- outcomes[j]
    type <- types[j]

    # omnibus: y ~ .node with intercept; test non-intercept coefficients = 0
    om <- lapply(seq_len(M), function(i) {
      d <- imps[[i]]; d$.node <- assigned[[i]]
      dk <- d[d$.node %in% keep, , drop = FALSE]
      dk$.node <- factor(as.character(dk$.node), levels = keep)
      .fit_extract(stats::reformulate(".node", y), dk, type,
                   cl = if (is.null(cluster)) NULL else dk[[cluster]],
                   vcov_type = vcov_type, drop_intercept = TRUE)
    })
    d1 <- .pool_d1(om)
    test_rows[[j]] <- data.frame(outcome = y,
                                 F   = if (is.null(d1)) NA_real_ else d1$F,
                                 df1 = if (is.null(d1)) NA_real_ else d1$df1,
                                 df2 = if (is.null(d1)) NA_real_ else d1$df2,
                                 p   = if (is.null(d1)) NA_real_ else d1$p,
                                 riv = if (is.null(d1)) NA_real_ else d1$riv)
    if (is.null(d1))
      warning("Omnibus test failed for outcome '", y, "'.", call. = FALSE)

    # node estimates: y ~ 0 + .node so coefficients are node means
    # (linear probability for a binary outcome, giving proportions)
    nm <- lapply(seq_len(M), function(i) {
      d <- imps[[i]]; d$.node <- droplevels(assigned[[i]])   # empty nodes: absent, not NA
      if (type == "binary") d[[y]] <- .as01(d[[y]])
      .fit_extract(stats::reformulate("0 + .node", y), d, "continuous",
                   cl = if (is.null(cluster)) NULL else d[[cluster]],
                   vcov_type = vcov_type, drop_intercept = FALSE)
    })
    est <- se <- rep(NA_real_, length(term)); names(est) <- names(se) <- as.character(term)
    ok_nm <- !vapply(nm, is.null, logical(1L))
    if (any(ok_nm)) {
      nm <- nm[ok_nm]
      for (nd in as.character(term)) {
        key <- paste0(".node", nd)
        q <- vapply(nm, function(f) if (key %in% names(f$coef)) f$coef[[key]] else NA_real_, numeric(1L))
        u <- vapply(nm, function(f) if (key %in% rownames(f$vcov)) f$vcov[key, key] else NA_real_, numeric(1L))
        good <- !is.na(q) & !is.na(u)
        if (sum(good) >= 2L) {
          dfc <- min(vapply(nm[good], function(f) f$dfcom, numeric(1L)))
          ps  <- .pool_scalar(q[good], u[good], dfcom = dfc)
          est[nd] <- ps$estimate; se[nd] <- ps$se
        }
      }
    }
    node_rows[[j]] <- data.frame(node_id    = as.integer(term),
                                 outcome    = y,
                                 estimate   = unname(est),
                                 se         = unname(se),
                                 n_confirm  = as.numeric(size[as.character(term)]),
                                 n_clusters = unname(node_cl),
                                 tested     = as.character(term) %in% keep,
                                 row.names  = NULL)
  }

  # ---- per-split tests ------------------------------------------------------
  internal <- .internal_nodes(tree)
  under <- lapply(stats::setNames(partykit::nodeids(tree), partykit::nodeids(tree)),
                  function(id) partykit::nodeids(tree, from = id, terminal = TRUE))
  split_rows <- list()
  tq <- stats::qnorm(1 - (1 - conf.level) / 2)

  for (nd in internal) {
    kid_ids <- nd$kids
    side <- lapply(assigned, function(a) {
      t <- as.integer(as.character(a))
      out <- rep(NA_character_, length(t))
      for (k in kid_ids) out[t %in% under[[as.character(k)]]] <- as.character(k)
      factor(out, levels = as.character(kid_ids))
    })
    kc <- vapply(side, function(sd) as.numeric(table(sd)), numeric(length(kid_ids)))
    kc <- matrix(kc, nrow = length(kid_ids))
    ok <- all(apply(kc, 1L, min) >= min_node)
    n_parent <- mean(colSums(kc))
    cl_parent <- if (is.null(cluster)) NA_real_ else
      mean(vapply(seq_len(M), function(i)
        length(unique(imps[[i]][[cluster]][!is.na(side[[i]])])), numeric(1L)))

    for (j in seq_along(outcomes)) {
      y <- outcomes[j]; type <- types[j]
      row <- data.frame(node_id = nd$node_id, depth = nd$depth,
                        split_var = nd$split_var, rule = nd$rule, outcome = y,
                        contrast = NA_real_, lower = NA_real_, upper = NA_real_,
                        F = NA_real_, df1 = NA_real_, df2 = NA_real_,
                        p = NA_real_, p_adj = NA_real_,
                        n_confirm = n_parent, n_clusters = cl_parent,
                        stringsAsFactors = FALSE)
      if (ok) {
        fx <- lapply(seq_len(M), function(i) {
          d <- imps[[i]]; d$.side <- side[[i]]
          dk <- d[!is.na(d$.side), , drop = FALSE]
          .fit_extract(stats::reformulate(".side", y), dk, type,
                       cl = if (is.null(cluster)) NULL else dk[[cluster]],
                       vcov_type = vcov_type, drop_intercept = TRUE)
        })
        d1 <- .pool_d1(fx)
        if (!is.null(d1)) {
          row$F <- d1$F; row$df1 <- d1$df1; row$df2 <- d1$df2; row$p <- d1$p
        }
        # pooled contrast: the coefficient with the largest |estimate|
        fxo <- fx[!vapply(fx, is.null, logical(1L))]
        if (length(fxo) >= 2L) {
          cn <- names(fxo[[1L]]$coef)
          qm <- t(vapply(fxo, function(f) f$coef[cn], numeric(length(cn))))
          # with one coefficient vapply simplifies to a vector and t() gives
          # 1 x M; reshape so rows are always imputations and columns coefficients
          if (nrow(qm) != length(fxo)) qm <- matrix(qm, nrow = length(fxo))
          k  <- which.max(abs(colMeans(qm)))
          q  <- qm[, k]
          u  <- vapply(fxo, function(f) f$vcov[cn[k], cn[k]], numeric(1L))
          ps <- .pool_scalar(q, u, dfcom = min(vapply(fxo, function(f) f$dfcom, numeric(1L))))
          crit <- if (is.finite(ps$df)) stats::qt(1 - (1 - conf.level) / 2, ps$df) else tq
          row$contrast <- ps$estimate
          row$lower <- ps$estimate - crit * ps$se
          row$upper <- ps$estimate + crit * ps$se
        }
      }
      split_rows[[length(split_rows) + 1L]] <- row
    }
  }
  splits <- if (length(split_rows)) do.call(rbind, split_rows) else NULL
  if (!is.null(splits)) {
    for (y in outcomes) {
      w <- splits$outcome == y & !is.na(splits$p)
      splits$p_adj[w] <- stats::p.adjust(splits$p[w], method = adjust)
    }
    splits <- splits[order(splits$outcome, splits$depth, splits$node_id), ]
    rownames(splits) <- NULL
  }

  structure(list(test       = do.call(rbind, test_rows),
                 splits     = splits,
                 nodes      = do.call(rbind, node_rows),
                 m          = M,
                 n_confirm  = nrow(imps[[1L]]),
                 cluster    = cluster,
                 n_clusters = n_cl,
                 vcov_type  = if (is.null(cluster)) "model" else vcov_type,
                 excluded   = as.integer(drop),
                 outcomes   = outcomes,
                 types      = types),
            class = "ctreeMI_confirm", adjust = adjust, conf.level = conf.level)
}


#' Discover a Tree and Confirm It on Held-Out Data
#'
#' A single call that splits the data, imputes each half separately, fits
#' the Stack/M tree on the discovery half, and tests the resulting
#' partition on the confirmation half. Equivalent to calling
#' [split_holdout()], [mice::mice()] twice, [ctree_stacked()] and
#' [confirm_ctreeMI()] in sequence.
#'
#' @param formula Model formula, as for [ctree_stacked()].
#' @param data A data frame containing missing values.
#' @param prop Proportion assigned to discovery. Default 0.5.
#' @param m Number of imputations for each half. Default 30.
#' @param strata Optional stratification column for the split.
#' @param cluster Optional name of a column identifying a higher-level unit.
#'   When given, the split keeps every unit intact on one side, the unit
#'   identifier is excluded from the imputation model, and the confirmation
#'   test uses a cluster-robust variance. See [split_holdout()] and
#'   [confirm_ctreeMI()].
#' @param seed Optional integer seed. The confirmation imputation uses
#'   `seed + 1L` so that the two halves are imputed with different streams.
#' @param alpha Nominal level for the discovery tree. Default 0.05.
#' @param mice_args A list of additional arguments passed to
#'   [mice::mice()] for both halves. If it contains a `predictorMatrix`,
#'   that matrix is used as given and `cluster` is not removed from it.
#' @param vcov_type Passed to [confirm_ctreeMI()]; ignored unless `cluster`
#'   is given.
#' @param ... Further arguments passed to [ctree_stacked()].
#'
#' @details
#' Imputing the halves separately is deliberate. A single imputation of the
#' full data lets each half's outcomes inform the other's imputed
#' predictors, and the confirmation test then no longer sees independent
#' data. The cost is that each imputation model is fitted to half the
#' observations. For most designs this is a reasonable trade; for small
#' samples it may not be, and the user should judge.
#'
#' @return An object of class `"ctreeMI_dc"` containing `tree` (the
#'   `ctreeMI` fit), `confirmation` (a `ctreeMI_confirm` object), `split`
#'   (the `ctreeMI_split` object), and the two `mids` objects as
#'   `imp_discover` and `imp_confirm`.
#'
#' @seealso [split_holdout()], [confirm_ctreeMI()], [ctree_stacked()]
#' @examples
#' \dontrun{
#' set.seed(11)
#' n <- 800
#' d <- data.frame(x1 = rnorm(n), x2 = factor(sample(letters[1:3], n, TRUE)))
#' d$y <- rnorm(n) + 1.5 * (d$x1 > 0.3) + (d$x2 == "a")
#' d$x1[sample(n, 120)] <- NA
#' d$x2[sample(n, 80)]  <- NA
#'
#' res <- discover_confirm(y ~ x1 + x2, data = d, m = 20, seed = 11)
#' res$tree
#' res$confirmation
#' }
#' @export
discover_confirm <- function(formula, data, prop = 0.5, m = 30L, strata = NULL,
                             cluster = NULL, seed = NULL, alpha = 0.05,
                             mice_args = list(), vcov_type = "HC1", ...) {
  if (!is.null(cluster) && cluster %in% all.vars(formula[[3L]]))
    warning("`cluster` variable '", cluster, "' also appears as a predictor in ",
            "`formula`. It will be used to define the split and the variance, ",
            "and the tree will be free to split on it; that is rarely intended.",
            call. = FALSE)

  parts <- split_holdout(data, prop = prop, strata = strata, cluster = cluster,
                         seed = seed)

  args_d <- c(list(data = parts$discover, m = m, printFlag = FALSE), mice_args)
  args_c <- c(list(data = parts$confirm,  m = m, printFlag = FALSE), mice_args)
  if (!is.null(seed)) {                # mice() rejects seed = NULL outright
    args_d$seed <- seed
    args_c$seed <- seed + 1L
  }
  # A cluster identifier is an identifier, not a predictor: keep it out of
  # the imputation model unless the user supplied their own predictorMatrix.
  if (!is.null(cluster) && is.null(mice_args$predictorMatrix)) {
    drop_cluster <- function(a) {
      pm <- mice::make.predictorMatrix(a$data)
      pm[, cluster] <- 0L
      pm[cluster, ] <- 0L
      a$predictorMatrix <- pm
      a
    }
    args_d <- drop_cluster(args_d)
    args_c <- drop_cluster(args_c)
  }
  imp_d <- do.call(mice::mice, args_d)
  imp_c <- do.call(mice::mice, args_c)

  tree <- ctree_stacked(formula, data = imp_d, alpha = alpha, verbose = FALSE, ...)

  conf <- if (length(partykit::nodeids(tree, terminal = TRUE)) >= 2L) {
    confirm_ctreeMI(tree, imp_c, cluster = cluster, vcov_type = vcov_type)
  } else {
    message("Discovery tree has a single node; no partition to confirm.")
    NULL
  }

  structure(list(tree = tree, confirmation = conf, split = parts,
                 imp_discover = imp_d, imp_confirm = imp_c),
            class = "ctreeMI_dc")
}


#' @export
print.ctreeMI_confirm <- function(x, digits = 3L, ...) {
  cat("Confirmation of a Stack/M partition on independent data\n")
  cat(sprintf("  %d confirmation observations, %d imputations, pooled by Rubin's rules\n",
              x$n_confirm, x$m))
  if (!is.null(x$cluster))
    cat(sprintf("  %d clusters (%s); cluster-robust variance, %s\n",
                x$n_clusters, x$cluster, x$vcov_type))
  if (length(x$excluded))
    cat("  Node(s) excluded for insufficient confirmation observations: ",
        paste(x$excluded, collapse = ", "), "\n", sep = "")
  cat("\nOmnibus test that the outcome differs across terminal nodes (Li-Raghunathan-Rubin D1):\n")
  t <- x$test
  t$F   <- round(t$F, digits); t$df2 <- round(t$df2, 1)
  t$riv <- round(t$riv, digits)
  t$p   <- format.pval(t$p, digits = digits)
  print(t, row.names = FALSE)
  if (!is.null(x$splits)) {
    lev <- attr(x, "conf.level"); if (is.null(lev)) lev <- 0.95
    cat("\nPer-split tests (children of each internal node compared within it):\n")
    sp <- x$splits
    ci <- ifelse(is.na(sp$contrast), "",
                 sprintf("%.*f [%.*f, %.*f]", digits, sp$contrast,
                         digits, sp$lower, digits, sp$upper))
    tab <- data.frame(node = sp$node_id, depth = sp$depth, split = sp$rule,
                      outcome = sp$outcome, difference = ci,
                      p_adj = ifelse(is.na(sp$p_adj), NA,
                                     format.pval(sp$p_adj, digits = digits)),
                      n = round(sp$n_confirm, 0),
                      stringsAsFactors = FALSE)
    if (!is.null(x$cluster)) tab$clusters <- round(sp$n_clusters, 0)
    print(tab, row.names = FALSE)
    cat("  difference: pooled contrast between children with ",
        format(100 * lev), "% interval; p_adj: ", attr(x, "adjust"),
        "-adjusted across splits.\n",
        "  A split beneath a non-confirmed split is not interpretable.\n", sep = "")
  }
  cat("\nPooled node estimates:\n")
  nd <- x$nodes
  nd$estimate  <- round(nd$estimate, digits)
  nd$se        <- round(nd$se, digits)
  nd$n_confirm <- round(nd$n_confirm, 1)
  if (is.null(x$cluster)) nd$n_clusters <- NULL else nd$n_clusters <- round(nd$n_clusters, 0)
  print(nd, row.names = FALSE)
  invisible(x)
}

#' @export
print.ctreeMI_dc <- function(x, ...) {
  cat("Discover-then-confirm Stack/M analysis\n")
  cat(sprintf("  discovery n = %d, confirmation n = %d\n\n",
              nrow(x$split$discover), nrow(x$split$confirm)))
  cat("Discovery tree:\n")
  print(x$tree, ...)
  cat("\n")
  if (is.null(x$confirmation)) {
    cat("No partition to confirm (single-node tree).\n")
  } else {
    print(x$confirmation, ...)
  }
  invisible(x)
}


# ---- internal helpers -------------------------------------------------------

# Fit one model to one imputed dataset and return its coefficient vector and
# covariance matrix, model-based or cluster-robust. Returns NULL if the fit
# fails or is rank-deficient, so callers can drop it.
.fit_extract <- function(formula, d, type, cl = NULL, vcov_type = "HC1",
                         drop_intercept = TRUE) {
  fit <- tryCatch(
    if (type == "binary") stats::glm(formula, data = d, family = stats::binomial())
    else stats::lm(formula, data = d),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  b <- stats::coef(fit)
  if (anyNA(b)) return(NULL)
  if (!is.null(cl)) {
    # The cluster vector must align with the rows the model actually used.
    # lm/glm drop rows with missing values; mirror that here, then refuse
    # to proceed if the lengths still disagree rather than misalign silently.
    cl <- as.character(cl)
    if (!is.null(fit$na.action)) cl <- cl[-fit$na.action]
    if (length(cl) != stats::nobs(fit)) return(NULL)
  }
  V <- tryCatch(
    if (is.null(cl)) stats::vcov(fit)
    else sandwich::vcovCL(fit, cluster = cl, type = vcov_type),
    error = function(e) NULL)
  if (is.null(V) || anyNA(V)) return(NULL)
  if (drop_intercept && "(Intercept)" %in% names(b)) {
    keep <- names(b) != "(Intercept)"
    b <- b[keep]; V <- V[keep, keep, drop = FALSE]
  }
  # Complete-data degrees of freedom for the small-sample adjustment. With
  # independent observations this is the residual df. With clustering the
  # effective number of independent units is the number of clusters, and
  # G - 1 is the conventional choice for cluster-robust inference.
  dfcom <- if (is.null(cl)) stats::df.residual(fit) else length(unique(cl)) - 1L
  list(coef = b, vcov = V, dfcom = dfcom)
}

# Multivariate Wald test pooled across imputations (Li, Raghunathan and
# Rubin, 1991), testing that the coefficient vector is zero. Takes a list of
# .fit_extract() results and tolerates NULL entries.
.pool_d1 <- function(fits) {
  fits <- fits[!vapply(fits, is.null, logical(1L))]
  m <- length(fits)
  if (m < 2L) return(NULL)
  cn <- names(fits[[1L]]$coef)
  if (!all(vapply(fits, function(f) identical(names(f$coef), cn), logical(1L))))
    return(NULL)
  k     <- length(cn)
  coefs <- t(vapply(fits, function(f) f$coef, numeric(k)))
  if (is.null(dim(coefs)) || nrow(coefs) != m) coefs <- matrix(coefs, nrow = m)
  qbar  <- colMeans(coefs)
  ubar  <- Reduce(`+`, lapply(fits, function(f) f$vcov)) / m
  b     <- stats::cov(coefs)
  ubar_inv <- tryCatch(solve(ubar), error = function(e) NULL)
  if (is.null(ubar_inv)) return(NULL)
  r  <- (1 + 1 / m) * sum(diag(b %*% ubar_inv)) / k
  d1 <- as.numeric(t(qbar) %*% ubar_inv %*% qbar) / (k * (1 + r))
  t  <- k * (m - 1)
  # large-sample reference df (Li, Raghunathan and Rubin, 1991)
  nu_lrr <- if (r <= 0) Inf
            else if (t > 4) 4 + (t - 4) * (1 + (1 - 2 / t) / r)^2
            else t * (1 + 1 / k) * (1 + 1 / r)^2 / 2
  # small-sample adjustment (Reiter, 2007, eqs 1-2), as in mitml / mice::D1
  dfcom <- min(vapply(fits, function(f) f$dfcom, numeric(1L)))
  nu <- nu_lrr
  if (is.finite(dfcom) && dfcom > 0 && r > 0 && t > 4) {
    a     <- r * t / (t - 2)
    vstar <- ((dfcom + 1) / (dfcom + 3)) * dfcom
    c0 <- 1 / (t - 4)
    c1 <- vstar - 2 * (1 + a)
    c2 <- vstar - 4 * (1 + a)
    z  <- 1 / c2 +
          c0 * (a^2 * c1 / ((1 + a)^2 * c2)) +
          c0 * (8 * a^2 * c1 / ((1 + a) * c2^2) + 4 * a^2 / ((1 + a) * c2)) +
          c0 * (4 * a^2 / (c2 * c1) + 16 * a^2 * c1 / c2^3) +
          c0 * (8 * a^2 / c2^2)
    v <- 4 + 1 / z
    # the adjustment can only reduce df; if it misbehaves (small dfcom with
    # large r drives c2 negative) fall back to the cap it is meant to impose
    nu <- if (is.finite(v) && v >= 4 && v <= nu_lrr) v else min(nu_lrr, dfcom)
  } else if (is.finite(dfcom) && dfcom > 0) {
    nu <- min(nu_lrr, dfcom)
  }
  list(F = d1, df1 = k, df2 = nu,
       p = stats::pf(d1, k, nu, lower.tail = FALSE), riv = r, dfcom = dfcom)
}

# Rubin's rules for a scalar: pooled estimate, total standard error, and the
# degrees of freedom of the t reference, with the Barnard-Rubin (1999)
# small-sample adjustment when complete-data df are supplied.
.pool_scalar <- function(q, u, dfcom = NULL) {
  m    <- length(q)
  qbar <- mean(q)
  ubar <- mean(u)
  b    <- if (m > 1L) stats::var(q) else 0
  tot  <- ubar + (1 + 1 / m) * b
  r    <- if (ubar > 0) (1 + 1 / m) * b / ubar else 0
  nu   <- if (r <= 0) Inf else (m - 1) * (1 + 1 / r)^2
  if (!is.null(dfcom) && is.finite(dfcom) && dfcom > 0 && tot > 0) {
    lambda <- (1 + 1 / m) * b / tot
    nu_obs <- ((dfcom + 1) / (dfcom + 3)) * dfcom * (1 - lambda)
    nu     <- 1 / (1 / nu + 1 / nu_obs)
  }
  list(estimate = qbar, se = sqrt(tot), df = nu, riv = r)
}

# 0/1 coding for a binary outcome given as factor, logical, or numeric.
.as01 <- function(v) {
  if (is.factor(v))  return(as.numeric(v == levels(v)[2L]))
  if (is.logical(v)) return(as.numeric(v))
  as.numeric(v)
}

# One entry per internal node: id, depth, split variable, child ids.
.internal_nodes <- function(tree) {
  dat      <- partykit::data_party(tree)
  varnames <- names(dat)
  out <- list()
  walk <- function(node, depth) {
    kids <- partykit::kids_node(node)
    if (length(kids) == 0L) return(invisible(NULL))
    sp    <- partykit::split_node(node)
    vname <- varnames[sp$varid]
    labs  <- tryCatch(split_labels(sp, vname, dat[[sp$varid]], NULL),
                      error = function(e) rep(NA_character_, length(kids)))
    out[[length(out) + 1L]] <<- list(
      node_id   = partykit::id_node(node),
      depth     = depth,
      split_var = vname,
      # A binary split is fully described by one side of it, so print only
      # the first condition; multiway splits need every branch listed.
      rule      = if (length(labs) == 2L && !anyNA(labs)) labs[1L]
                  else paste(labs, collapse = " | "),
      kid_rules = labs,
      kids      = vapply(kids, partykit::id_node, integer(1L)))
    for (k in kids) walk(k, depth + 1L)
    invisible(NULL)
  }
  walk(partykit::node_party(tree), 1L)
  out
}

.as_imp_list <- function(data) {
  if (inherits(data, "mids")) {
    lapply(seq_len(data$m), function(i) mice::complete(data, i))
  } else if (is.list(data) && all(vapply(data, is.data.frame, logical(1L)))) {
    data
  } else {
    stop("`data` must be a 'mids' object or a list of imputed data frames.")
  }
}

.resolve_types <- function(d, outcomes, outcome_type) {
  if (identical(outcome_type, "auto")) {
    vapply(outcomes, function(y) {
      v <- d[[y]]
      if (is.numeric(v)) return("continuous")
      if (is.factor(v) && nlevels(v) == 2L) return("binary")
      if (is.logical(v)) return("binary")
      stop("Outcome '", y, "' is neither numeric nor a two-level factor; ",
           "specify `outcome_type` explicitly.")
    }, character(1L))
  } else {
    if (length(outcome_type) != length(outcomes))
      stop("`outcome_type` must have one entry per outcome.")
    if (!all(outcome_type %in% c("continuous", "binary")))
      stop("`outcome_type` entries must be 'continuous' or 'binary'.")
    stats::setNames(outcome_type, outcomes)
  }
}


#' Collapse Splits That Did Not Confirm
#'
#' Takes a discovery tree and its confirmation, and returns the tree with
#' every split that failed to confirm collapsed into a terminal node.
#'
#' @param tree The `ctreeMI` object that was confirmed.
#' @param confirmation A `ctreeMI_confirm` object from [confirm_ctreeMI()].
#' @param alpha Level at which a split counts as confirmed. Default 0.05.
#' @param which Which p-value to use: `"adjusted"` (default) or `"raw"`.
#' @param outcome For a multivariate outcome, the name of the outcome whose
#'   tests govern pruning, or `"any"` (default) to retain a split confirmed
#'   for any outcome, or `"all"` to require every outcome.
#'
#' @details
#' Pruning is bottom-up and a split is collapsed only once all of its own
#' internal descendants have been collapsed, matching [prune_stackM()]. A
#' split whose p-value could not be computed, because a child node fell
#' below `min_node` in some imputation, is treated as unconfirmed; the tree
#' offers no evidence for it on the confirmation data.
#'
#' The result is the partition the confirmation data supports. Its terminal
#' nodes can be reported without the caveat that attaches to the discovery
#' tree, since the splits defining them were tested on data the tree had not
#' seen.
#'
#' @return A `ctreeMI` object with unconfirmed splits collapsed. Its
#'   `ctreeMI_info` attribute gains `confirmed`, a data frame recording which
#'   splits were retained and why.
#'
#' @seealso [confirm_ctreeMI()], [discover_confirm()]
#' @examples
#' \dontrun{
#' res <- discover_confirm(y ~ ., data = d, m = 30, seed = 1)
#' pruned <- prune_unconfirmed(res$tree, res$confirmation)
#' plot(pruned)
#' }
#' @export
prune_unconfirmed <- function(tree, confirmation, alpha = 0.05,
                              which = c("adjusted", "raw"),
                              outcome = "any") {
  if (!inherits(tree, "ctreeMI"))
    stop("`tree` must be a 'ctreeMI' object.")
  if (!inherits(confirmation, "ctreeMI_confirm"))
    stop("`confirmation` must be a 'ctreeMI_confirm' object from confirm_ctreeMI().")
  which <- match.arg(which)
  sp <- confirmation$splits
  if (is.null(sp) || !nrow(sp)) return(tree)

  pcol <- if (which == "adjusted") "p_adj" else "p"

  # a split is confirmed if its p-value meets alpha; NA counts as unconfirmed
  ok_by_node <- split(sp, sp$node_id)
  passes <- vapply(ok_by_node, function(g) {
    if (!identical(outcome, "any") && !identical(outcome, "all"))
      g <- g[g$outcome == outcome, , drop = FALSE]
    pv <- g[[pcol]]
    if (!length(pv)) return(FALSE)
    hit <- !is.na(pv) & pv < alpha
    if (identical(outcome, "all")) all(hit) else any(hit)
  }, logical(1L))
  names(passes) <- names(ok_by_node)

  kids_of <- function(node) partykit::kids_node(node)

  # bottom-up: keep a split if it confirmed, or if any descendant split is kept
  keep <- new.env(parent = emptyenv())
  resolve <- function(node) {
    id  <- as.character(partykit::id_node(node))
    kd  <- kids_of(node)
    if (!length(kd)) return(FALSE)
    below <- any(vapply(kd, resolve, logical(1L)))
    val   <- isTRUE(unname(passes[id])) || below
    assign(id, val, envir = keep)
    val
  }
  resolve(partykit::node_party(tree))

  # Collapsing changes which nodes are terminal, so the fitted node
  # assignments must be remapped: every observation in a node that has been
  # absorbed now belongs to the node that absorbed it. partykit validates
  # this, and a party object whose fitted ids name nodes that no longer
  # exist is rejected.
  remap <- integer(0)
  collapse <- function(node) {
    idn <- as.integer(partykit::id_node(node))
    id  <- as.character(idn)
    kd  <- kids_of(node)
    if (!length(kd)) {                       # already terminal: maps to itself
      remap[id] <<- idn
      return(node)
    }
    if (!isTRUE(get0(id, envir = keep, ifnotfound = FALSE))) {
      for (t in partykit::nodeids(tree, from = idn, terminal = TRUE))
        remap[as.character(t)] <<- idn
      return(partykit::partynode(idn))
    }
    node$kids <- lapply(kd, collapse)
    node
  }

  new_root <- collapse(partykit::node_party(tree))

  fit <- tree$fitted
  old_ids <- as.character(fit[["(fitted)"]])
  new_ids <- unname(remap[old_ids])
  if (anyNA(new_ids))
    stop("Internal error remapping fitted nodes during pruning.")
  fit[["(fitted)"]] <- new_ids

  out <- partykit::party(new_root,
                         data = partykit::data_party(tree),
                         fitted = fit,
                         terms = tree$terms)
  class(out) <- class(tree)
  info <- attr(tree, "ctreeMI_info")
  info$confirmed <- data.frame(
    node_id   = as.integer(names(passes)),
    confirmed = unname(passes),
    retained  = vapply(names(passes),
                       function(i) isTRUE(get0(i, envir = keep, ifnotfound = FALSE)),
                       logical(1L)),
    row.names = NULL)
  info$pruned_by <- paste0("confirmation, ", which, " p < ", alpha)
  attr(out, "ctreeMI_info") <- info
  out
}


#' A Methods Paragraph For a Discover-Then-Confirm Analysis
#'
#' Generates a paragraph describing a discover-then-confirm analysis in the
#' form usually required by a methods section: how the sample was divided,
#' how each half was imputed, the tree that was discovered, and which of its
#' splits survived testing on the confirmation data.
#'
#' @param object A `ctreeMI_dc` object from [discover_confirm()], or a
#'   `ctreeMI_confirm` object from [confirm_ctreeMI()].
#' @param tree The discovery tree, required when `object` is a
#'   `ctreeMI_confirm` object and ignored otherwise.
#' @param digits Number of digits used in the reported quantities.
#'
#' @return An object of class `"ctreeMI_report"`, as returned by
#'   [report_ctreeMI()]: a list whose `text` element is the paragraph.
#'
#' @seealso [discover_confirm()], [report_ctreeMI()]
#' @examples
#' \dontrun{
#' res <- discover_confirm(y ~ ., data = d, m = 30, seed = 1)
#' report_confirm(res)
#' }
#' @export
report_confirm <- function(object, tree = NULL, digits = 3) {
  if (inherits(object, "ctreeMI_dc")) {
    tree <- object$tree
    cf   <- object$confirmation
    n_d  <- nrow(object$split$discover)
    n_c  <- nrow(object$split$confirm)
  } else if (inherits(object, "ctreeMI_confirm")) {
    if (is.null(tree)) stop("`tree` is required when `object` is a ctreeMI_confirm.")
    cf  <- object
    n_d <- attr(tree, "ctreeMI_info")$n_original
    n_c <- object$n_confirm
  } else {
    stop("`object` must be a 'ctreeMI_dc' or 'ctreeMI_confirm' object.")
  }
  if (is.null(cf)) stop("No confirmation to report.")

  info <- attr(tree, "ctreeMI_info")
  sp   <- cf$splits
  n_sp <- if (is.null(sp)) 0L else length(unique(sp$node_id))
  conf_ids <- if (is.null(sp)) integer(0) else
    unique(sp$node_id[!is.na(sp$p_adj) & sp$p_adj < 0.05])
  n_conf <- length(conf_ids)
  n_term <- length(partykit::nodeids(tree, terminal = TRUE))

  clustered <- !is.null(cf$cluster)
  split_desc <- if (clustered) {
    ncl <- if (inherits(object, "ctreeMI_dc") && !is.null(object$split$n_clusters))
      object$split$n_clusters else c(NA, cf$n_clusters)
    paste0("divided at the level of ", cf$cluster, " into a discovery set of ",
           n_d, " observations", if (!is.na(ncl[1L])) paste0(" in ", ncl[1L], " ", cf$cluster, "s") else "",
           " and a confirmation set of ", n_c, " observations",
           if (!is.na(ncl[2L])) paste0(" in ", ncl[2L], " ", cf$cluster, "s") else "",
           ", so that no ", cf$cluster, " contributed to both")
  } else {
    paste0("divided at random into a discovery set of ", n_d,
           " and a confirmation set of ", n_c)
  }
  txt <- paste0(
    "The sample of ", n_d + n_c, " observations was ", split_desc,
    ". Each was multiply imputed separately, with M = ", info$m,
    " imputations, so that neither set's outcomes informed the other's ",
    "imputed predictor values. A conditional inference tree was fitted to the ",
    "stacked discovery imputations with the Stack/M correction at alpha = ",
    format(info$alpha, digits = digits), ", yielding ", n_term,
    " terminal node", if (n_term == 1L) "" else "s", " defined by ", n_sp,
    " split", if (n_sp == 1L) "" else "s", ". Because the node-level test of ",
    "that procedure is not calibrated to its nominal level under ",
    "outcome-conditioned imputation, and because testing a partition on the ",
    "data used to select it is invalid in any case, each split was then tested ",
    "on the confirmation set. The confirmation observations falling within each ",
    "internal node were divided by that node's own rule and the resulting ",
    "children compared, pooling across imputations by the method of Li, ",
    "Raghunathan and Rubin (1991)",
    if (clustered) paste0(" with variance estimated robust to clustering by ",
                          cf$cluster) else "",
    " and adjusting across splits by Holm's method. ",
    n_conf, " of ", n_sp, " split", if (n_sp == 1L) "" else "s",
    " met the 0.05 level on the confirmation data",
    if (n_conf < n_sp) paste0(
      "; the remainder", if (n_sp - n_conf == 1L) " was" else " were",
      " not supported and the corresponding node",
      if (n_sp - n_conf == 1L) " was" else "s were", " collapsed") else "", ".")

  structure(list(text = txt, n_discover = n_d, n_confirm = n_c,
                 m = info$m, alpha = info$alpha,
                 n_splits = n_sp, n_confirmed = n_conf,
                 n_terminal = n_term),
            class = "ctreeMI_report")
}
