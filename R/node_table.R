## ---------------------------------------------------------------------------
## ctreeMI 0.2.0: node_table() and report_ctreeMI()
## ---------------------------------------------------------------------------

#' @importFrom stats terms formula
NULL

## ---------------------------------------------------------------------------
## node_table(): terminal-node summary with effective sample size
## ---------------------------------------------------------------------------

#' Terminal-Node Summary With Effective Sample Sizes
#'
#' @description
#' Summarises the nodes of a fitted tree: the split path leading to each node,
#' its size in the stacked data, its effective size in original observations
#' (stacked size divided by `M`), and node-level outcome summaries.
#'
#' @details
#' Terminal-node sample sizes printed by `partykit` refer to the stacked data,
#' in which every observation appears `M` times. `effective_n` divides by `M`
#' to give the number of original observations behind each node. Node-level
#' means and proportions need no such adjustment: they are already averages
#' over the imputations.
#'
#' @param object A `ctreeMI` object from [ctree_stacked()], or any `party`
#'   object (in which case `M` is taken to be 1).
#' @param terminal_only Logical. Report terminal nodes only (default), or every
#'   node including the root.
#' @param max_levels Integer or `NULL`. Truncate factor level sets longer than
#'   this in the printed split path.
#' @param digits Number of digits for the outcome summaries.
#'
#' @return A data frame of class `"ctreeMI_nodes"` with one row per node:
#'   `node_id`, `depth`, `n_stacked`, `effective_n`, any outcome summaries,
#'   `path`, and a `conditions` list column holding the split conditions
#'   individually.
#'
#' @seealso [ctree_stacked()], [report_ctreeMI()]
#'
#' @examples
#' set.seed(1)
#' imps <- lapply(1:5, function(i) {
#'   x <- stats::rnorm(200)
#'   data.frame(x = x, y = stats::rnorm(200) + 1.5 * (x > 0))
#' })
#' fit <- ctree_stacked(y ~ x, data = imps, verbose = FALSE)
#' node_table(fit)
#'
#' @export
node_table <- function(object, terminal_only = TRUE, max_levels = NULL,
                       digits = 3) {

  if (!inherits(object, "party")) {
    stop("`object` must be a party/ctreeMI object from ctree_stacked().")
  }

  info <- attr(object, "ctreeMI_info")
  m    <- if (is.null(info)) 1L else info$m

  ## ---- collect paths by walking the tree ---------------------------------
  dat      <- partykit::data_party(object)
  varnames <- names(dat)

  paths <- new.env(parent = emptyenv())

  walk <- function(node, acc) {
    id <- partykit::id_node(node)
    assign(as.character(id), acc, envir = paths)

    kids <- partykit::kids_node(node)
    if (length(kids) == 0L) return(invisible(NULL))

    sp    <- partykit::split_node(node)
    vname <- varnames[sp$varid]
    conds <- split_labels(sp, vname, dat[[sp$varid]], max_levels)

    for (k in seq_along(kids)) {
      lab <- if (k <= length(conds)) conds[k] else NA_character_
      walk(kids[[k]], c(acc, lab))
    }
    invisible(NULL)
  }

  walk(partykit::node_party(object), character(0))

  ## ---- which nodes to report ---------------------------------------------
  ids <- if (isTRUE(terminal_only)) {
    partykit::nodeids(object, terminal = TRUE)
  } else {
    partykit::nodeids(object)
  }

  cond_list <- lapply(ids, function(i) get(as.character(i), envir = paths))

  ## ---- node sizes ---------------------------------------------------------
  fit_ids   <- partykit::predict.party(object, type = "node")
  n_stacked <- as.integer(sapply(ids, function(i) {
    sum(node_members(object, i, fit_ids))
  }))

  out <- data.frame(
    node_id     = as.integer(ids),
    depth       = as.integer(lengths(cond_list)),
    n_stacked   = n_stacked,
    effective_n = round(n_stacked / m, 1),
    path        = vapply(cond_list, function(z) {
      if (length(z) == 0L) "root" else paste(z, collapse = " & ")
    }, character(1)),
    stringsAsFactors = FALSE
  )
  out$conditions <- cond_list

  ## ---- outcome summaries --------------------------------------------------
  resp <- response_frame(object)
  if (!is.null(resp)) {
    summ <- outcome_summary(resp, fit_ids, ids, digits)
    if (!is.null(summ)) out <- cbind(out, summ)
  }

  ord <- c("node_id", "depth", "n_stacked", "effective_n",
           setdiff(names(out), c("node_id", "depth", "n_stacked",
                                 "effective_n", "path", "conditions")),
           "path", "conditions")
  out <- out[, ord, drop = FALSE]
  rownames(out) <- NULL

  attr(out, "m")            <- m
  attr(out, "ctreeMI_info") <- info
  class(out) <- c("ctreeMI_nodes", "data.frame")
  out
}


## ---- internal: which stacked rows fall in node i --------------------------
node_members <- function(object, id, fit_ids) {
  kids <- partykit::nodeids(object, from = id, terminal = TRUE)
  fit_ids %in% kids
}


## ---- internal: readable labels for each kid of a split --------------------
split_labels <- function(sp, vname, xvar, max_levels) {

  ## numeric / ordered split
  if (!is.null(sp$breaks)) {
    br    <- sp$breaks
    right <- if (is.null(sp$right)) TRUE else isTRUE(sp$right)
    brs   <- format(br, digits = 5, trim = TRUE)

    if (length(br) == 1L) {
      if (right) {
        return(c(paste0(vname, " <= ", brs), paste0(vname, " > ", brs)))
      } else {
        return(c(paste0(vname, " < ", brs), paste0(vname, " >= ", brs)))
      }
    }
    labs <- c(paste0(vname, " <= ", brs[1]))
    for (j in seq_len(length(br) - 1L)) {
      labs <- c(labs, paste0(brs[j], " < ", vname, " <= ", brs[j + 1L]))
    }
    labs <- c(labs, paste0(vname, " > ", brs[length(br)]))
    return(labs)
  }

  ## categorical split via index
  if (!is.null(sp$index)) {
    lv  <- if (is.factor(xvar)) levels(xvar) else sort(unique(xvar))
    idx <- sp$index
    nk  <- max(idx, na.rm = TRUE)
    vapply(seq_len(nk), function(k) {
      sel <- lv[which(idx == k)]
      paste0(vname, " in ", brace_set(sel, max_levels))
    }, character(1))
  } else {
    rep(NA_character_, 2L)
  }
}


## ---- internal: format a level set, optionally truncated --------------------
brace_set <- function(levs, max_levels) {
  levs <- as.character(levs)
  if (!is.null(max_levels) && length(levs) > max_levels) {
    shown <- levs[seq_len(max_levels)]
    extra <- length(levs) - max_levels
    return(paste0("{", paste(shown, collapse = ", "),
                  ", +", extra, " more}"))
  }
  paste0("{", paste(levs, collapse = ", "), "}")
}


## ---- internal: pull the response as a data frame ---------------------------
response_frame <- function(object) {
  dat <- partykit::data_party(object)

  ## names of the response variable(s), from the model formula
  rn <- NULL
  info <- attr(object, "ctreeMI_info")
  if (!is.null(info$formula)) {
    lhs <- info$formula[[2L]]
    rn  <- all.vars(lhs)
  }
  if (is.null(rn) || length(rn) == 0L) {
    tt <- try(stats::terms(object), silent = TRUE)
    if (!inherits(tt, "try-error")) {
      rn <- all.vars(stats::formula(tt))[1L]
    }
  }

  fitcol <- dat[["(response)"]]
  if (!is.null(fitcol)) {
    if (is.data.frame(fitcol)) return(fitcol)
    nm <- if (length(rn) == 1L) rn else "response"
    out <- data.frame(fitcol, stringsAsFactors = FALSE)
    names(out) <- nm
    return(out)
  }

  if (is.null(rn) || !all(rn %in% names(dat))) return(NULL)
  dat[, rn, drop = FALSE]
}


## ---- internal: per-node outcome summaries ----------------------------------
outcome_summary <- function(resp, fit_ids, ids, digits) {

  members <- lapply(ids, function(i) fit_ids == i)

  cols <- list()
  for (nm in names(resp)) {
    v <- resp[[nm]]

    if (is.numeric(v)) {
      cols[[paste0(nm, "_mean")]] <- round(vapply(members, function(sel) {
        if (!any(sel)) NA_real_ else mean(v[sel], na.rm = TRUE)
      }, numeric(1)), digits)

    } else if (is.factor(v) && nlevels(v) == 2L) {
      pos <- levels(v)[2]
      cols[[paste0(nm, "_prop_", pos)]] <- round(vapply(members, function(sel) {
        if (!any(sel)) NA_real_ else mean(v[sel] == pos, na.rm = TRUE)
      }, numeric(1)), digits)

    } else if (is.factor(v)) {
      cols[[paste0(nm, "_mode")]] <- vapply(members, function(sel) {
        if (!any(sel)) return(NA_character_)
        tb <- table(v[sel])
        names(tb)[which.max(tb)]
      }, character(1))
    }
  }

  if (length(cols) == 0L) return(NULL)
  as.data.frame(cols, stringsAsFactors = FALSE)
}


## ---- print method for the node table --------------------------------------

#' Print a Node Table
#'
#' @param x A `"ctreeMI_nodes"` object from [node_table()].
#' @param ... Passed to [print.data.frame()].
#'
#' @return `x`, invisibly.
#' @export
print.ctreeMI_nodes <- function(x, ...) {
  m <- attr(x, "m")
  cat(sprintf(
    "-- ctreeMI terminal nodes (M = %s; effective_n = n_stacked / M) --\n",
    if (is.null(m)) "?" else m))
  y <- as.data.frame(x)
  y$conditions <- NULL
  print(y, row.names = FALSE, ...)
  invisible(x)
}


## ---------------------------------------------------------------------------
## report_ctreeMI(): study-level methods paragraph
## ---------------------------------------------------------------------------

#' A Methods Paragraph For a ctreeMI Analysis
#'
#' @description
#' Generates a paragraph describing the analysis in the form usually required
#' by a methods section: the number of imputations, the size of the stacked
#' dataset, the model, the Stack / M correction and the threshold applied, how
#' many candidate splits survived it, and the shape of the resulting tree.
#'
#' @param object A `ctreeMI` object from [ctree_stacked()].
#' @param digits Number of digits used in the node summaries.
#'
#' @return An object of class `"ctreeMI_report"`: a list whose `text` element
#'   is the paragraph, alongside the quantities it reports.
#'
#' @seealso [ctree_stacked()], [node_table()]
#'
#' @examples
#' set.seed(1)
#' imps <- lapply(1:5, function(i) {
#'   x <- stats::rnorm(200)
#'   data.frame(x = x, y = stats::rnorm(200) + 1.5 * (x > 0))
#' })
#' fit <- ctree_stacked(y ~ x, data = imps, verbose = FALSE)
#' report_ctreeMI(fit)
#'
#' @export
report_ctreeMI <- function(object, digits = 3) {

  if (!inherits(object, "ctreeMI")) {
    stop("`object` must be a ctreeMI object from ctree_stacked().")
  }

  info <- attr(object, "ctreeMI_info")
  tab  <- node_table(object, terminal_only = TRUE)

  n_term <- nrow(tab)
  depth  <- if (n_term) max(tab$depth) else 0L
  fvar   <- paste(deparse(info$formula), collapse = " ")

  txt <- paste0(
    "Missing data were handled using multiple imputation with M = ",
    info$m, " imputed datasets. Following the stacked-imputation workflow ",
    "for recursive partitioning (Rodgers et al., 2021), the imputed ",
    "datasets were concatenated vertically, yielding a stacked dataset of ",
    format(info$n_stacked, big.mark = ","), " rows (",
    format(info$n_original, big.mark = ","), " observations x ", info$m,
    " imputations). A single conditional inference tree (Hothorn, Hornik, ",
    "& Zeileis, 2006) was then fit to the stacked data using the model ",
    fvar, ". Because stacking inflates the nominal sample size by a factor ",
    "of M and correspondingly inflates node-level test statistics, each ",
    "node-level test statistic was divided by M and its p value recomputed ",
    "before node selection (the Stack/M correction; Sherlock et al., 2026), ",
    "with splits retained at alpha = ", format(info$alpha), ". ",
    "Of ", info$n_splits_before, " candidate split",
    if (identical(info$n_splits_before, 1L)) "" else "s",
    " identified on the stacked data, ", info$n_splits_after,
    " met the corrected criterion. The resulting tree ",
    "contained ", n_term, " terminal node",
    if (n_term == 1L) "" else "s",
    " with a maximum depth of ", depth, " split",
    if (depth == 1L) "" else "s",
    ". Terminal-node sample sizes are reported on the original scale ",
    "(stacked node size divided by M). Analyses were conducted in R using ",
    "the ctreeMI package."
  )

  out <- list(
    text          = txt,
    m             = info$m,
    n_original    = info$n_original,
    n_stacked     = info$n_stacked,
    alpha           = info$alpha,
    correction      = info$correction,
    n_splits_before = info$n_splits_before,
    n_splits_after  = info$n_splits_after,
    node_stats      = info$node_stats,
    n_terminal    = n_term,
    max_depth     = depth,
    formula       = info$formula
  )
  class(out) <- "ctreeMI_report"
  out
}


#' Print a Methods Paragraph
#'
#' @param x A `"ctreeMI_report"` object from [report_ctreeMI()].
#' @param width Wrapping width in characters.
#' @param ... Currently unused.
#'
#' @return `x`, invisibly.
#' @export
print.ctreeMI_report <- function(x, width = 76, ...) {
  cat(paste(strwrap(x$text, width = width), collapse = "\n"), "\n")
  invisible(x)
}
