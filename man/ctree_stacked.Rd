\name{ctree_stacked}
\alias{ctree_stacked}
\title{Conditional Inference Tree on Stacked Multiply Imputed Data}
\description{
  Fits a conditional inference tree (\code{ctree}) on stacked multiply
  imputed datasets using the Stack / M rescaling procedure described in
  Sherlock et al. (2026). Multiply imputed datasets are concatenated
  vertically ("stacked"), and the significance threshold used for
  node-level pruning is divided by the number of imputations \eqn{M} to
  counteract the artificially inflated sample size. This yields a single,
  coherent, interpretable tree that incorporates imputation variability
  without requiring the pooling of structurally different trees.
}
\usage{
ctree_stacked(formula, data, m = NULL, alpha = 0.05, verbose = TRUE, ...)
}
\arguments{
  \item{formula}{A model formula, passed to \code{\link[partykit]{ctree}}.}
  \item{data}{A \code{mids} object from \code{\link[mice]{mice}}, a list
    of imputed data frames, or a single complete data frame. If a single
    complete data frame is supplied, the function falls back to a standard
    \code{\link[partykit]{ctree}} call with a warning.}
  \item{m}{Integer. Number of imputations to use. Defaults to all
    available datasets. Ignored when \code{data} is a plain data frame.}
  \item{alpha}{Numeric. Nominal significance threshold for node-level
    splitting (default \code{0.05}). The corrected threshold actually
    applied is \code{alpha / m}. Must be strictly between 0 and 1.}
  \item{verbose}{Logical. If \code{TRUE} (default), prints a message
    summarising the stacking and correction applied.}
  \item{...}{Additional arguments passed to
    \code{\link[partykit]{ctree_control}}.}
}
\value{
  An object of class \code{c("ctreeMI", "constparty", "party")}.
  All \pkg{partykit} methods (\code{plot}, \code{predict}, etc.) work on
  this object. An additional \code{ctreeMI_info} attribute carries:
  \describe{
    \item{\code{m}}{Number of imputations used.}
    \item{\code{n_original}}{Rows in one imputed dataset.}
    \item{\code{n_stacked}}{Total rows in the stacked dataset (\eqn{M \times n}).}
    \item{\code{alpha_nominal}}{Nominal alpha supplied by the user.}
    \item{\code{alpha_applied}}{Corrected alpha applied (\code{alpha / m}).}
    \item{\code{formula}}{The model formula.}
    \item{\code{call}}{The matched call.}
  }
  Node-level sample sizes reported by \code{print()} and \code{plot()}
  reflect the stacked dataset. Divide by \eqn{M} to obtain effective
  per-node counts in the original data.
}
\details{
  \strong{Methodological background}

  Conditional inference trees (Hothorn, Hornik & Zeileis, 2006) use
  permutation-based significance tests to select splits, providing
  built-in protection against spurious partitioning. When data are
  multiply imputed, pooling trees fitted to separate imputations is
  infeasible: structurally different trees define different subgroups,
  making the targets of inference incomparable across imputations.

  Rodgers et al. (2021) proposed stacking the \eqn{M} imputed datasets
  and fitting a single tree to the combined data. This produces one
  coherent, interpretable tree. The complication is that stacking inflates
  the nominal sample size by \eqn{M}, causing test statistics at each
  node to be similarly inflated.

  Sherlock et al. (2026) proposed and validated the \strong{Stack / M
  correction}: use a significance threshold of \code{alpha / M}. Monte
  Carlo simulations under MCAR confirmed sub-nominal (conservative)
  type-I error and acceptable power, making this approach well-suited for
  exploratory analyses where interpretability is prioritised.
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

  Rodgers, J., Khoo, S.-T., & Ludtke, O. (2021).
  Handling missing data in structural equation models using multiple
  imputation and stacking.
  \emph{Structural Equation Modeling}, 28(6), 915--930.
  \doi{10.1080/10705511.2021.1916925}
}
\seealso{
  \code{\link[partykit]{ctree}},
  \code{\link[partykit]{ctree_control}},
  \code{\link[mice]{mice}},
  \code{\link{stack_imputations}},
  \code{\link{rescale_alpha}},
  \code{\link{print.ctreeMI}},
  \code{\link{summary.ctreeMI}}
}
\examples{
\dontrun{
  library(mice)

  # Introduce missingness into the airquality dataset
  set.seed(42)
  aq <- airquality
  aq$Ozone[sample(nrow(aq), 20)]   <- NA
  aq$Solar.R[sample(nrow(aq), 15)] <- NA

  # Impute (M = 20)
  imp <- mice(aq, m = 20, printFlag = FALSE)

  # Fit ctree with Stack/M correction
  fit <- ctree_stacked(Ozone ~ Solar.R + Wind + Temp + Month,
                       data  = imp,
                       alpha = 0.05)
  print(fit)
  plot(fit)
}

# Example using a list of data frames (no mice required)
set.seed(1)
make_df <- function(i) {
  set.seed(i)
  n  <- 100
  x1 <- rnorm(n)
  y  <- x1 + rnorm(n)
  data.frame(y = y, x1 = x1)
}
imp_list <- lapply(1:10, make_df)
fit <- ctree_stacked(y ~ x1, data = imp_list, alpha = 0.05, verbose = FALSE)
print(fit)
}
