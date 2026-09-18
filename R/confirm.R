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
#'   allocated at random.
#' @param seed Optional integer seed.
#'
#' @details
#' The split is performed on original observations, not on stacked rows,
#' and before any imputation. This matters. Imputing the full data once and
#' then splitting leaves each half's imputed values informed by the other
#' half's outcomes, which compromises the independence the confirmation
#' step relies on. Impute the two halves separately.
#'
#' @return A list of class `"ctreeMI_split"` with elements `discover` and
#'   `confirm`, each a data frame, and `index`, the row indices of `data`
#'   assigned to the discovery set.
#'
#' @seealso [confirm_ctreeMI()], [discover_confirm()]
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(200), y = rnorm(200))
#' d$x[sample(200, 30)] <- NA
#' parts <- split_holdout(d, prop = 0.5, seed = 1)
#' nrow(parts$discover); nrow(parts$confirm)
#' @export
split_holdout <- function(data, prop = 0.5, strata = NULL, seed = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  if (!is.numeric(prop) || prop <= 0 || prop >= 1)
    stop("`prop` must be strictly between 0 and 1.")
  n <- nrow(data)
  if (n < 4L) stop("Too few observations to split.")
  if (!is.null(seed)) set.seed(seed)

  if (is.null(strata)) {
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

  structure(list(discover = data[idx, , drop = FALSE],
                 confirm  = data[-idx, , drop = FALSE],
                 index    = idx),
            class = "ctreeMI_split")
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
#' multivariate Wald procedure of Li, Raghunathan and Rubin (1991),
#' implemented in [mice::D1()].
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
#'       (mean observations per imputation), and `tested`.}
#'     \item{`m`}{Number of imputations used.}
#'     \item{`n_confirm`}{Number of confirmation observations.}
#'     \item{`excluded`}{Terminal-node ids excluded for having fewer than
#'       `min_node` observations, if any.}
#'   }
#'
#' @references
#' Li, K. H., Raghunathan, T. E., and Rubin, D. B. (1991). Large-sample
#' significance levels from multiply imputed data using moment-based
#' statistics and an F reference distribution. \emph{Journal of the
#' American Statistical Association}, 86, 1065--1073.
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
                            adjust = "holm", conf.level = 0.95) {
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

  term <- partykit::nodeids(tree, terminal = TRUE)
  if (length(term) < 2L)
    stop("The discovery tree has a single terminal node; there is no partition to confirm.")

  # ---- outcome types ----------------------------------------------------
  types <- .resolve_types(imps[[1L]], outcomes, outcome_type)

  # ---- node assignment per imputation -----------------------------------
  assigned <- lapply(imps, function(d)
    factor(stats::predict(tree, newdata = d, type = "node"), levels = term))

  # Nodes too thin to test. A node must reach `min_node` in EVERY imputation,
  # not merely on average: an empty node in one imputation would give that
  # imputation's model a different coefficient vector and break the pooled test.
  counts <- vapply(assigned, function(a) as.numeric(table(a)), numeric(length(term)))
  dimnames(counts) <- list(as.character(term), NULL)
  size <- rowMeans(counts)
  keep <- rownames(counts)[apply(counts, 1L, min) >= min_node]
  drop <- setdiff(as.character(term), keep)
  if (length(keep) < 2L)
    stop("Fewer than two terminal nodes receive at least `min_node` confirmation observations.")

  # ---- per-outcome test and pooled node estimates -----------------------
  test_rows <- vector("list", length(outcomes))
  node_rows <- vector("list", length(outcomes))

  for (j in seq_along(outcomes)) {
    y    <- outcomes[j]
    type <- types[j]

    fits  <- vector("list", M)
    means <- matrix(NA_real_, M, length(term), dimnames = list(NULL, as.character(term)))
    ses   <- means

    for (i in seq_len(M)) {
      d <- imps[[i]]
      d$.node <- assigned[[i]]
      dk <- d[d$.node %in% keep, , drop = FALSE]
      dk$.node <- factor(as.character(dk$.node), levels = keep)

      fits[[i]] <- if (type == "binary") {
        stats::glm(stats::reformulate(".node", y), data = dk, family = stats::binomial())
      } else {
        stats::lm(stats::reformulate(".node", y), data = dk)
      }

      for (nd in as.character(term)) {
        v <- d[[y]][d$.node == nd]
        if (type == "binary") {
          v <- if (is.factor(v)) as.numeric(v == levels(d[[y]])[2L]) else as.numeric(v)
        }
        k <- length(v)
        means[i, nd] <- if (k) mean(v) else NA_real_
        ses[i, nd]   <- if (k > 1L) stats::sd(v) / sqrt(k) else NA_real_
      }
    }

    d1 <- tryCatch(mice::D1(mice::as.mira(fits)), error = function(e) NULL)
    if (is.null(d1)) {
      test_rows[[j]] <- data.frame(outcome = y, F = NA_real_, df1 = NA_real_,
                                   df2 = NA_real_, p = NA_real_, riv = NA_real_)
      warning("D1 test failed for outcome '", y, "'; check for separation or empty nodes.")
    } else {
      r  <- d1$result
      cn <- colnames(r)
      pick <- function(nm, pos) {
        if (!is.null(cn) && nm %in% cn) unname(r[1L, nm]) else unname(r[1L, pos])
      }
      test_rows[[j]] <- data.frame(outcome = y,
                                   F   = pick("F.value", 1L),
                                   df1 = pick("df1",     2L),
                                   df2 = pick("df2",     3L),
                                   p   = pick("P(>F)",   4L),
                                   riv = pick("RIV",     5L))
    }

    # Rubin's rules per node: mean of means; W-bar + (1 + 1/M) B
    pooled <- vapply(as.character(term), function(nd) {
      mu <- means[, nd]; s <- ses[, nd]
      if (all(is.na(mu))) return(c(NA_real_, NA_real_))
      qbar <- mean(mu, na.rm = TRUE)
      W    <- mean(s^2, na.rm = TRUE)
      B    <- if (sum(!is.na(mu)) > 1L) stats::var(mu, na.rm = TRUE) else 0
      c(qbar, sqrt(W + (1 + 1 / M) * B))
    }, numeric(2L))

    node_rows[[j]] <- data.frame(node_id   = as.integer(term),
                                 outcome   = y,
                                 estimate  = pooled[1L, ],
                                 se        = pooled[2L, ],
                                 n_confirm = as.numeric(size[as.character(term)]),
                                 tested    = as.character(term) %in% keep,
                                 row.names = NULL)
  }

  # ---- per-split tests ---------------------------------------------------
  # For every internal node: take confirmation observations inside it, split
  # them by that node's own rule, test children against one another.
  internal <- .internal_nodes(tree)
  # terminal ids beneath each node, so membership can be read from `assigned`
  under <- lapply(stats::setNames(partykit::nodeids(tree), partykit::nodeids(tree)),
                  function(id) partykit::nodeids(tree, from = id, terminal = TRUE))
  split_rows <- list()

  for (nd in internal) {
    kid_ids <- nd$kids
    # side label per observation per imputation; NA if outside this node
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

    for (j in seq_along(outcomes)) {
      y <- outcomes[j]; type <- types[j]
      row <- data.frame(node_id = nd$node_id, depth = nd$depth,
                        split_var = nd$split_var, rule = nd$rule, outcome = y,
                        contrast = NA_real_, lower = NA_real_, upper = NA_real_,
                        F = NA_real_, df1 = NA_real_, df2 = NA_real_,
                        p = NA_real_, p_adj = NA_real_, n_confirm = n_parent,
                        stringsAsFactors = FALSE)
      if (ok) {
        fits <- vector("list", M)
        for (i in seq_len(M)) {
          d <- imps[[i]]; d$.side <- side[[i]]
          dk <- d[!is.na(d$.side), , drop = FALSE]
          fits[[i]] <- if (type == "binary") {
            stats::glm(stats::reformulate(".side", y), data = dk, family = stats::binomial())
          } else {
            stats::lm(stats::reformulate(".side", y), data = dk)
          }
        }
        d1 <- tryCatch(mice::D1(mice::as.mira(fits)), error = function(e) NULL)
        if (!is.null(d1)) {
          r <- d1$result; cn <- colnames(r)
          pick <- function(nm, pos) if (!is.null(cn) && nm %in% cn) unname(r[1L, nm]) else unname(r[1L, pos])
          row$F <- pick("F.value", 1L); row$df1 <- pick("df1", 2L)
          row$df2 <- pick("df2", 3L);   row$p   <- pick("P(>F)", 4L)
        }
        # Pooled contrast between the children, with a confidence interval.
        # For a two-way split this is the difference the split asserts; for a
        # multiway split it is the largest coefficient against the reference.
        ct <- tryCatch({
          pl <- summary(mice::pool(mice::as.mira(fits)), conf.int = TRUE,
                        conf.level = conf.level)
          pl <- pl[pl$term != "(Intercept)", , drop = FALSE]
          if (!nrow(pl)) NULL else {
            k  <- which.max(abs(pl$estimate))
            # mice names the interval columns by percentile, e.g. "2.5 %";
            # they are the final two columns regardless of the level used.
            nm <- names(pl)
            c(pl$estimate[k], pl[[nm[length(nm) - 1L]]][k], pl[[nm[length(nm)]]][k])
          }
        }, error = function(e) NULL)
        if (!is.null(ct)) {
          row$contrast <- ct[1L]; row$lower <- ct[2L]; row$upper <- ct[3L]
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

  structure(list(test      = do.call(rbind, test_rows),
                 splits    = splits,
                 nodes     = do.call(rbind, node_rows),
                 m         = M,
                 n_confirm = nrow(imps[[1L]]),
                 excluded  = as.integer(drop),
                 outcomes  = outcomes,
                 types     = types),
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
#' @param seed Optional integer seed. The confirmation imputation uses
#'   `seed + 1L` so that the two halves are imputed with different streams.
#' @param alpha Nominal level for the discovery tree. Default 0.05.
#' @param mice_args A list of additional arguments passed to
#'   [mice::mice()] for both halves.
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
                             seed = NULL, alpha = 0.05, mice_args = list(), ...) {
  parts <- split_holdout(data, prop = prop, strata = strata, seed = seed)

  args_d <- c(list(data = parts$discover, m = m, printFlag = FALSE), mice_args)
  args_c <- c(list(data = parts$confirm,  m = m, printFlag = FALSE), mice_args)
  if (!is.null(seed)) {                # mice() rejects seed = NULL outright
    args_d$seed <- seed
    args_c$seed <- seed + 1L
  }
  imp_d <- do.call(mice::mice, args_d)
  imp_c <- do.call(mice::mice, args_c)

  tree <- ctree_stacked(formula, data = imp_d, alpha = alpha, verbose = FALSE, ...)

  conf <- if (length(partykit::nodeids(tree, terminal = TRUE)) >= 2L) {
    confirm_ctreeMI(tree, imp_c)
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

  txt <- paste0(
    "The sample of ", n_d + n_c, " observations was divided at random into a ",
    "discovery set of ", n_d, " and a confirmation set of ", n_c,
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
    "Raghunathan and Rubin (1991) and adjusting across splits by Holm's method. ",
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
