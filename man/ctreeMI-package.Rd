\name{ctreeMI-package}
\alias{ctreeMI-package}
\alias{ctreeMI}
\docType{package}
\title{
  \packageTitle{ctreeMI}
}
\description{
  \packageDescription{ctreeMI}
}
\details{
  The main function is \code{\link{ctree_stacked}}, which accepts a
  \code{mids} object from \pkg{mice}, a list of imputed data frames, or a
  plain data frame. It returns a \code{ctreeMI} object that inherits from
  \pkg{partykit}'s \code{constparty} class, so all standard methods
  (\code{plot}, \code{predict}, \code{\link[partykit]{nodeids}}, etc.)
  work without modification.

  Two utility functions are also exported:
  \code{\link{stack_imputations}} (stack a list of data frames) and
  \code{\link{rescale_alpha}} (compute \code{alpha / M}).
}
\author{
  \packageAuthor{ctreeMI}

  Maintainer: \packageMaintainer{ctreeMI}
}
\references{
  Sherlock, P., Mansolf, M., Hofheimer, J., Hockett, C. W.,
  O'Connor, T. G., Roubinov, D., Graff, J. C., Lai, J.-S.,
  Bush, N. R., Wright, R. J., & Chiu, Y.-H. M. (2026).
  Beyond linear risk: A machine learning approach to understanding
  perinatal depression in context.
  \emph{Multivariate Behavioral Research}, 1--16.
  \doi{10.1080/00273171.2026.2661244}

  Hothorn, T., Hornik, K., & Zeileis, A. (2006).
  Unbiased recursive partitioning: A conditional inference framework.
  \emph{Journal of Computational and Graphical Statistics}, 15(3), 651--674.
  \doi{10.1198/106186006X133933}

  Hothorn, T., & Zeileis, A. (2015).
  \pkg{partykit}: A modular toolkit for recursive partitioning in R.
  \emph{Journal of Machine Learning Research}, 16, 3905--3909.
}
\keyword{package}
\seealso{
  \code{\link{ctree_stacked}},
  \code{\link{stack_imputations}},
  \code{\link{rescale_alpha}},
  \code{\link[partykit]{ctree}},
  \code{\link[mice]{mice}}
}
