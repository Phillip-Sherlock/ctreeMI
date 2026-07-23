\name{report_ctreeMI}
\alias{report_ctreeMI}
\alias{print.ctreeMI_report}
\title{Methods Paragraph for a Stacked-Imputation Tree}
\description{
  Generates a short methods paragraph describing how a
  \code{\link{ctree_stacked}} model was fit, populated with the actual
  values used: the number of imputations, the original and stacked sample
  sizes, the nominal and applied significance thresholds, and the size and
  depth of the resulting tree. The paragraph is intended to be edited and
  pasted into a manuscript, and it standardizes how the Stack/M correction
  is described across studies.
}
\usage{
report_ctreeMI(object, digits = 3)
}
\arguments{
  \item{object}{An object of class \code{"ctreeMI"} as returned by
    \code{\link{ctree_stacked}}.}
  \item{digits}{Integer. Significant digits used when reporting the
    applied significance threshold.}
}
\value{
  An object of class \code{"ctreeMI_report"}: a list whose
  \code{text} element is the methods paragraph as a single character
  string, alongside the underlying quantities (\code{m},
  \code{n_original}, \code{n_stacked}, \code{alpha_nominal},
  \code{alpha_applied}, \code{n_terminal}, \code{max_depth},
  \code{formula}) so they can be used programmatically. The print method
  wraps and displays the paragraph.
}
\details{
  The generated text describes the imputation and stacking workflow and
  the Stack/M threshold correction, and states that node sizes are
  reported on the original scale. It does not describe the imputation
  model itself, which is specific to the analysis and must be added by
  the author: the number of iterations, the variables included, and the
  imputation methods used should be reported alongside this paragraph.
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
  \code{\link{ctree_stacked}}, \code{\link{node_table}}
}
\examples{
set.seed(1)
make_df <- function(i) {
  set.seed(i)
  n  <- 200
  x1 <- rnorm(n)
  y  <- 2 * (x1 > 0) + rnorm(n)
  data.frame(y = y, x1 = x1)
}
imp_list <- lapply(1:10, make_df)
fit <- ctree_stacked(y ~ x1, data = imp_list, verbose = FALSE)

# Print the methods paragraph
report_ctreeMI(fit)

# Access the raw text or the underlying quantities
rep <- report_ctreeMI(fit)
rep$n_terminal
substr(rep$text, 1, 60)
}
