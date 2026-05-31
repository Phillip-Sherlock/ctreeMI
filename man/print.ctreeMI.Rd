\name{print.ctreeMI}
\alias{print.ctreeMI}
\title{Print Method for ctreeMI Objects}
\description{
  Prints a header summarizing the stacked-imputation settings, followed
  by the standard \pkg{partykit} tree output.
}
\usage{
\method{print}{ctreeMI}(x, \dots)
}
\arguments{
  \item{x}{An object of class \code{"ctreeMI"} as returned by
    \code{\link{ctree_stacked}}.}
  \item{\dots}{Further arguments passed to the \pkg{partykit} print
    method.}
}
\value{
  \code{x}, invisibly.
}
\seealso{
  \code{\link{ctree_stacked}}, \code{\link{summary.ctreeMI}}
}
\examples{
set.seed(1)
imp_list <- lapply(1:5, function(i) {
  set.seed(i)
  data.frame(y = rnorm(80), x = rnorm(80))
})
fit <- ctree_stacked(y ~ x, data = imp_list, verbose = FALSE)
print(fit)
}
