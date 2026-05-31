#' Print Method for ctreeMI Objects
#'
#' @description
#' Prints a summary of the stacked-imputation settings used to fit the
#' tree, followed by the standard [partykit::ctree()] output.
#'
#' @param x An object of class `"ctreeMI"`, as returned by
#'   [ctree_stacked()].
#' @param ... Further arguments passed to the partykit print method.
#'
#' @return `x`, invisibly.
#'
#' @export
print.ctreeMI <- function(x, ...) {
  info <- attr(x, "ctreeMI_info")

  cat("-- ctreeMI: Stacked-Imputation Conditional Inference Tree --\n")
  if (!is.null(info)) {
    cat(sprintf("  Imputations (M)      : %d\n",   info$m))
    cat(sprintf("  Rows per imputation  : %d\n",   info$n_original))
    cat(sprintf("  Stacked rows         : %d\n",   info$n_stacked))
    cat(sprintf("  Nominal alpha        : %.4f\n", info$alpha_nominal))
    cat(sprintf("  Corrected alpha (a/M): %.6f\n", info$alpha_applied))
    cat(sprintf("  Formula              : %s\n",
                deparse(info$formula, width.cutoff = 60L)))
  }
  cat("------------------------------------------------------------\n\n")

  # Delegate to partykit's print method (strip ctreeMI class first)
  x_party <- x
  class(x_party) <- class(x_party)[class(x_party) != "ctreeMI"]
  print(x_party, ...)

  invisible(x)
}


#' Summary Method for ctreeMI Objects
#'
#' @description
#' Returns a named list with details about the stacked-imputation fit
#' and the resulting tree structure (number of terminal nodes, tree depth).
#'
#' @param object An object of class `"ctreeMI"`, as returned by
#'   [ctree_stacked()].
#' @param ... Currently unused.
#'
#' @return A list (invisibly) with components:
#'   \describe{
#'     \item{`ctreeMI_info`}{Stacking metadata from [ctree_stacked()].}
#'     \item{`n_terminal_nodes`}{Number of terminal nodes in the fitted tree.}
#'     \item{`depth`}{Maximum depth of the fitted tree.}
#'   }
#'
#' @export
summary.ctreeMI <- function(object, ...) {

  info <- attr(object, "ctreeMI_info")

  node_ids     <- partykit::nodeids(object)
  terminal_ids <- partykit::nodeids(object, terminal = TRUE)
  n_terminal   <- length(terminal_ids)
  depth        <- tree_depth(object)

  cat("-- ctreeMI Summary --------------------------------------------------\n")
  if (!is.null(info)) {
    cat(sprintf("  Imputations (M)      : %d\n",   info$m))
    cat(sprintf("  Rows per imputation  : %d\n",   info$n_original))
    cat(sprintf("  Stacked rows         : %d\n",   info$n_stacked))
    cat(sprintf("  Nominal alpha        : %.4f\n", info$alpha_nominal))
    cat(sprintf("  Corrected alpha (a/M): %.6f\n", info$alpha_applied))
    cat(sprintf("  Formula              : %s\n",
                deparse(info$formula, width.cutoff = 60L)))
  }
  cat("-- Tree structure ---------------------------------------------------\n")
  cat(sprintf("  Total nodes          : %d\n", length(node_ids)))
  cat(sprintf("  Terminal nodes       : %d\n", n_terminal))
  cat(sprintf("  Maximum depth        : %d\n", depth))
  cat("---------------------------------------------------------------------\n")

  out <- list(
    ctreeMI_info     = info,
    n_terminal_nodes = n_terminal,
    depth            = depth
  )
  invisible(out)
}


# ## Internal: compute tree depth############################################-

#' @keywords internal
tree_depth <- function(tree) {
  # Recursively find max depth of a party/constparty object
  node <- partykit::node_party(tree)
  depth_node(node)
}

#' @keywords internal
depth_node <- function(node) {
  kids <- partykit::kids_node(node)
  if (length(kids) == 0L) return(0L)
  1L + max(sapply(kids, depth_node))
}
