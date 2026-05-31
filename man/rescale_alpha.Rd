\name{rescale_alpha}
\alias{rescale_alpha}
\title{Rescale Significance Threshold for the Stack / M Correction}
\description{
  Computes the adjusted significance threshold \code{alpha / M} used in
  the Stack / M correction of Sherlock et al. (2026). Dividing the
  nominal alpha by the number of imputations \eqn{M} counteracts the
  inflated test statistics that arise from stacking \eqn{M} copies of
  the data.
}
\usage{
rescale_alpha(alpha = 0.05, m)
}
\arguments{
  \item{alpha}{Numeric. Nominal significance level (default \code{0.05}).
    Must be strictly between 0 and 1.}
  \item{m}{A single positive integer. Number of imputations.}
}
\value{
  A single numeric value: \code{alpha / m}.
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
  \code{\link{ctree_stacked}}, \code{\link{stack_imputations}}
}
\examples{
rescale_alpha(0.05, 30)   # 0.001666...
rescale_alpha(0.05, 10)   # 0.005
rescale_alpha(0.01, 5)    # 0.002
}
