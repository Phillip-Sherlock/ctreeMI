\name{summary.ctreeMI}
\alias{summary.ctreeMI}
\title{Summary Method for ctreeMI Objects}
\description{
  Prints and returns a summary of the stacked-imputation fit and the
  resulting tree structure (number of nodes, maximum depth).
}
\usage{
\method{summary}{ctreeMI}(object, \dots)
}
\arguments{
  \item{object}{An object of class \code{"ctreeMI"} as returned by
    \code{\link{ctree_stacked}}.}
  \item{\dots}{Currently unused.}
}
\value{
  A list (returned invisibly) with components:
  \describe{
    \item{\code{ctreeMI_info}}{Stacking metadata from
      \code{\link{ctree_stacked}}.}
    \item{\code{n_terminal_nodes}}{Number of terminal nodes.}
    \item{\code{depth}}{Maximum depth of the fitted tree.}
  }
}
\seealso{
  \code{\link{ctree_stacked}}, \code{\link{print.ctreeMI}}
}
\examples{
set.seed(1)
imp_list <- lapply(1:5, function(i) {
  set.seed(i)
  data.frame(y = rnorm(80), x = rnorm(80))
})
fit <- ctree_stacked(y ~ x, data = imp_list, verbose = FALSE)
summary(fit)
}
