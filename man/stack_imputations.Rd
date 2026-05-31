\name{stack_imputations}
\alias{stack_imputations}
\title{Stack Multiply Imputed Datasets}
\description{
  Concatenates a list of imputed data frames into a single stacked data
  frame. An imputation-index column (\code{.imp}) is added to identify
  which imputed dataset each row originated from.
}
\usage{
stack_imputations(data_list, imp_col = ".imp")
}
\arguments{
  \item{data_list}{A list of data frames, all with the same dimensions
    and column names, representing \eqn{M} imputed versions of the same
    dataset.}
  \item{imp_col}{Character string. Name of the imputation-index column
    added to the stacked data (default \code{".imp"}). Set to \code{NULL}
    to suppress the column.}
}
\value{
  A single data frame with \eqn{M \times n} rows, where \eqn{n} is the
  number of rows in each imputed dataset. If \code{imp_col} is not
  \code{NULL}, an integer column recording the imputation index is
  appended.
}
\references{
  Sherlock, P., Mansolf, M., Hofheimer, J., Hockett, C. W.,
  O'Connor, T. G., Roubinov, D., Graff, J. C., Lai, J.-S.,
  Bush, N. R., Wright, R. J., & Chiu, Y.-H. M. (2026).
  Beyond linear risk: A machine learning approach to understanding
  perinatal depression in context.
  \emph{Multivariate Behavioral Research}, 1--16.
  \doi{10.1080/00273171.2026.2661244}

  Rodgers, J., Khoo, S.-T., & Ludtke, O. (2021).
  Handling missing data in structural equation models using multiple
  imputation and stacking.
  \emph{Structural Equation Modeling}, 28(6), 915--930.
  \doi{10.1080/10705511.2021.1916925}
}
\seealso{
  \code{\link{ctree_stacked}}, \code{\link{rescale_alpha}}
}
\examples{
df1 <- data.frame(x = 1:5, y = c(2, 4, 6, 8, 10))
df2 <- data.frame(x = 1:5, y = c(2, 3, 6, 9, 10))
df3 <- data.frame(x = 1:5, y = c(1, 4, 5, 8, 11))

stacked <- stack_imputations(list(df1, df2, df3))
nrow(stacked)           # 15
table(stacked$.imp)     # 5 rows per imputation
}
