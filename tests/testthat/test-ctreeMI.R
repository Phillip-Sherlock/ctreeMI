test_that("rescale_alpha works correctly", {
  expect_equal(rescale_alpha(0.05, 30), 0.05 / 30)
  expect_equal(rescale_alpha(0.05, 10), 0.05 / 10)
  expect_equal(rescale_alpha(0.01, 5),  0.01 / 5)
  expect_error(rescale_alpha(0.05, 0))
  expect_error(rescale_alpha(1.0,  10))
  expect_error(rescale_alpha(0.0,  10))
  expect_error(rescale_alpha(0.05))
})

test_that("stack_imputations works on a list of data frames", {
  df1 <- data.frame(x = 1:5, y = c(2, 4, 6, 8, 10))
  df2 <- data.frame(x = 1:5, y = c(2, 3, 6, 9, 10))
  df3 <- data.frame(x = 1:5, y = c(1, 4, 5, 8, 11))

  stacked <- stack_imputations(list(df1, df2, df3))

  expect_equal(nrow(stacked), 15)
  expect_equal(ncol(stacked), 3)           # x, y, .imp
  expect_true(".imp" %in% names(stacked))
  expect_equal(sort(unique(stacked$.imp)), 1:3)
})

test_that("stack_imputations adds no index column when imp_col = NULL", {
  df1 <- data.frame(x = 1:3, y = 4:6)
  df2 <- data.frame(x = 1:3, y = 5:7)
  stacked <- stack_imputations(list(df1, df2), imp_col = NULL)
  expect_equal(ncol(stacked), 2)
  expect_false(".imp" %in% names(stacked))
})

test_that("stack_imputations errors on mismatched dims", {
  df1 <- data.frame(x = 1:5, y = 1:5)
  df2 <- data.frame(x = 1:4, y = 1:4)
  expect_error(stack_imputations(list(df1, df2)))
})

test_that("stack_imputations errors on mismatched columns", {
  df1 <- data.frame(x = 1:5, y = 1:5)
  df2 <- data.frame(a = 1:5, b = 1:5)
  expect_error(stack_imputations(list(df1, df2)))
})

test_that("ctree_stacked runs on a list of data frames and returns ctreeMI", {
  skip_if_not_installed("partykit")

  set.seed(1)
  make_df <- function() {
    n <- 100
    x1 <- rnorm(n)
    x2 <- sample(c("A", "B"), n, replace = TRUE)
    y  <- x1 + rnorm(n)
    data.frame(y = y, x1 = x1, x2 = factor(x2))
  }
  imp_list <- lapply(1:5, function(i) { set.seed(i); make_df() })

  fit <- ctree_stacked(y ~ x1 + x2, data = imp_list, alpha = 0.05,
                       verbose = FALSE)

  expect_s3_class(fit, "ctreeMI")
  expect_s3_class(fit, "constparty")

  info <- attr(fit, "ctreeMI_info")
  expect_equal(info$m, 5)
  expect_equal(info$n_original, 100)
  expect_equal(info$n_stacked, 500)
  expect_equal(info$alpha_nominal, 0.05)
  expect_equal(info$alpha_applied, 0.05 / 5)
})

test_that("ctree_stacked warns and falls back for single data frame", {
  skip_if_not_installed("partykit")

  set.seed(99)
  df <- data.frame(y = rnorm(50), x = rnorm(50))
  expect_warning(
    fit <- ctree_stacked(y ~ x, data = df, verbose = FALSE),
    "Only one dataset"
  )
  expect_false(inherits(fit, "ctreeMI"))
})

test_that("print.ctreeMI runs without error", {
  skip_if_not_installed("partykit")

  set.seed(2)
  imp_list <- lapply(1:3, function(i) {
    data.frame(y = rnorm(60), x = rnorm(60))
  })
  fit <- ctree_stacked(y ~ x, data = imp_list, verbose = FALSE)
  expect_output(print(fit), "ctreeMI")
})

test_that("summary.ctreeMI runs without error", {
  skip_if_not_installed("partykit")

  set.seed(3)
  imp_list <- lapply(1:3, function(i) {
    data.frame(y = rnorm(60), x = rnorm(60))
  })
  fit <- ctree_stacked(y ~ x, data = imp_list, verbose = FALSE)
  # Call the S3 method directly to capture both output and return value
  expect_output(
    out <- summary.ctreeMI(fit),
    "ctreeMI Summary"
  )
  expect_type(out, "list")
  expect_true("n_terminal_nodes" %in% names(out))
  expect_true("depth" %in% names(out))
})
