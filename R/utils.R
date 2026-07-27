#' Stack Multiply Imputed Datasets
#'
#' @description
#' Concatenates a list of imputed data frames into a single "stacked" data
#' frame. An imputation index column (`.imp`) is added to identify which
#' imputed dataset each row originated from. This is the first step of the
#' stacked-imputation workflow described in Rodgers et al. (2021) and
#' applied to conditional inference trees in Sherlock et al. (2026).
#'
#' @param data_list A list of data frames, each the same dimensions,
#'   representing M imputed versions of the same dataset.
#' @param imp_col Character string. Name of the imputation-index column
#'   added to the stacked data (default `".imp"`). Set to `NULL` to
#'   suppress this column.
#'
#' @return A single data frame with `M * n` rows, where `n` is the number
#'   of rows in each imputed dataset.
#'
#' @references
#' Sherlock, P., et al. (2026). Beyond linear risk: A machine learning
#'   approach to understanding perinatal depression in context.
#'   *Multivariate Behavioral Research*, 1-16.
#'   \doi{10.1080/00273171.2026.2661244}
#'
#' Rodgers, D. M., Jacobucci, R., & Grimm, K. J. (2021). A multiple
#'   imputation approach for handling missing data in classification and
#'   regression trees. *Journal of Behavioral Data Science*, 1(1), 127-153.
#'   \doi{10.35566/jbds/v1n1/p6}
#'
#' @examples
#' df1 <- data.frame(x = 1:5, y = c(2, 4, 6, 8, 10))
#' df2 <- data.frame(x = 1:5, y = c(2, 3, 6, 9, 10))
#' stacked <- stack_imputations(list(df1, df2))
#' nrow(stacked) # 10
#' table(stacked$.imp)
#'
#' @export
stack_imputations <- function(data_list, imp_col = ".imp") {

  if (!is.list(data_list) || length(data_list) < 1) {
    stop("`data_list` must be a non-empty list of data frames.")
  }
  if (!all(sapply(data_list, is.data.frame))) {
    stop("All elements of `data_list` must be data frames.")
  }

  # Check consistent dimensions and column names
  ref_dim  <- dim(data_list[[1]])
  ref_cols <- colnames(data_list[[1]])

  for (i in seq_along(data_list)) {
    if (!identical(dim(data_list[[i]]), ref_dim)) {
      stop(sprintf(
        "Imputed dataset %d has different dimensions (%d x %d) from ",
        "dataset 1 (%d x %d).",
        i, nrow(data_list[[i]]), ncol(data_list[[i]]),
        ref_dim[1], ref_dim[2]
      ))
    }
    if (!identical(colnames(data_list[[i]]), ref_cols)) {
      stop(sprintf(
        "Imputed dataset %d has different column names from dataset 1.", i
      ))
    }
  }

  # Add imputation index if requested
  if (!is.null(imp_col)) {
    data_list <- lapply(seq_along(data_list), function(i) {
      d <- data_list[[i]]
      d[[imp_col]] <- i
      d
    })
  }

  do.call(rbind, data_list)
}
