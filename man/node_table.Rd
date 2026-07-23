\name{node_table}
\alias{node_table}
\alias{print.ctreeMI_nodes}
\title{Terminal-Node Table with Effective Sample Sizes}
\description{
  Summarizes the nodes of a fitted \code{\link{ctree_stacked}} tree as a
  data frame, giving for each node the sequence of splits that defines it,
  its size in the stacked data, and its \emph{effective} sample size on the
  original scale (stacked size divided by the number of imputations
  \eqn{M}).

  Because a stacked dataset repeats every observation \eqn{M} times, node
  sizes printed by \code{print()} and \code{plot()} are \eqn{M} times
  larger than the number of observations they represent. Reporting those
  numbers directly overstates precision. This function performs the
  division so that node sizes can be reported and interpreted correctly.
}
\usage{
node_table(object, terminal_only = TRUE, max_levels = NULL, digits = 3)
}
\arguments{
  \item{object}{An object of class \code{"ctreeMI"} as returned by
    \code{\link{ctree_stacked}}.}
  \item{terminal_only}{Logical. If \code{TRUE} (default) only terminal
    nodes are returned. If \code{FALSE}, inner nodes are included as well,
    which is useful for describing the tree structure.}
  \item{max_levels}{Integer or \code{NULL}. For splits on unordered
    factors with many levels, the number of levels to print before
    truncating (e.g. \code{{A, B, +3 more}}). \code{NULL} (default) prints
    the full level set.}
  \item{digits}{Integer. Number of decimal places for outcome summaries.}
}
\value{
  A data frame of class \code{"ctreeMI_nodes"} with one row per node and
  the following columns:
  \describe{
    \item{\code{node_id}}{Node identifier, matching \code{print()} and
      \code{plot()} output.}
    \item{\code{depth}}{Number of splits from the root to this node.}
    \item{\code{n_stacked}}{Number of rows in the stacked data.}
    \item{\code{effective_n}}{\code{n_stacked / M}: the node size on the
      scale of the original data. This is the number to report.}
    \item{outcome summaries}{Node means for numeric responses, or the
      proportion in the second level for binary factor responses. One
      column per response variable, so multivariate outcomes are
      supported.}
    \item{\code{path}}{The split conditions leading to the node, joined
      by \code{" & "}, suitable for pasting into a manuscript table.}
    \item{\code{conditions}}{A list-column holding the same conditions as
      a character vector, for programmatic use.}
  }
}
\details{
  Note that \code{effective_n} is exact at the root and is an
  observation-scale rescaling elsewhere. Because imputed values differ
  across imputations, an individual's \eqn{M} rows are not guaranteed to
  follow the same path down the tree, so \code{effective_n} should be read
  as the node's size expressed on the original scale rather than as a
  count of distinct individuals.
}
\references{
  Sherlock, P., Mansolf, M., Hofheimer, J., Hockett, C. W.,
  O'Connor, T. G., Roubinov, D., Graff, J. C., Lai, J.-S.,
  Bush, N. R., Wright, R. J., & Chiu, Y.-H. M. (2026).
  Beyond linear risk: A machine learning approach to understanding
  perinatal depression in context.
  \emph{Multivariate Behavioral Research}, 1--16.
  \doi{10.1080/00273171.2026.2661244}
}
\seealso{
  \code{\link{ctree_stacked}}, \code{\link{report_ctreeMI}}
}
\examples{
set.seed(1)
make_df <- function(i) {
  set.seed(i)
  n  <- 200
  x1 <- rnorm(n)
  x2 <- factor(sample(c("A", "B", "C"), n, TRUE))
  y  <- 2 * (x1 > 0) + 1.5 * (x2 == "A") + rnorm(n)
  data.frame(y = y, x1 = x1, x2 = x2)
}
imp_list <- lapply(1:10, make_df)
fit <- ctree_stacked(y ~ x1 + x2, data = imp_list, verbose = FALSE)

# Terminal nodes with effective sample sizes
node_table(fit)

# Truncate long factor level sets
node_table(fit, max_levels = 2)

# Include inner nodes
node_table(fit, terminal_only = FALSE)
}
