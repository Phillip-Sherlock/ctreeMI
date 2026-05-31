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
#' Rodgers, J., Khoo, S.-T., & L?dtke, O. (2021). Handling missing data in
#'   structural equation models using multiple imputation and stacking.
#'   *Structural Equation Modeling*, 28(6), 915-930.
#'   \doi{10.1080/10705511.2021.1916925}
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


#' Rescale Significance Threshold for Stack / M Correction
#'
#' @description
#' Computes the adjusted significance threshold `alpha / M` used in the
#' Stack / M correction of Sherlock et al. (2026). Dividing the nominal
#' alpha by the number of imputations M counteracts the inflated test
#' statistics that arise from stacking M copies of the data, thereby
#' preserving approximately nominal type-I error control (conservative)
#' in the node-level pruning step of [partykit::ctree()].
#'
#' @param alpha Numeric. Nominal significance level (default 0.05).
#' @param m Integer. Number of imputations.
#'
#' @return A single numeric value: `alpha / m`.
#'
#' @references
#' Sherlock, P., et al. (2026). Beyond linear risk: A machine learning
#'   approach to understanding perinatal depression in context.
#'   *Multivariate Behavioral Research*, 1-16.
#'   \doi{10.1080/00273171.2026.2661244}
#'
#' @examples
#' rescale_alpha(0.05, 30)   # 0.001667
#' rescale_alpha(0.05, 10)   # 0.005
#'
#' @export
rescale_alpha <- function(alpha = 0.05, m) {
  if (missing(m) || !is.numeric(m) || length(m) != 1 || m < 1) {
    stop("`m` must be a single positive integer.")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be strictly between 0 and 1.")
  }
  alpha / m
}
