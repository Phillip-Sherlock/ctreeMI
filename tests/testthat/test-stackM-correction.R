## Regression tests for the Stack/M correctness fix in 0.3.0.
##
## The bug these guard against: 0.1.0-0.2.0 divided alpha by M instead of
## dividing the node statistic by M, which under-corrects severely.

test_that("statistic rescaling is NOT equivalent to threshold rescaling", {
  alpha <- 0.05
  for (df in c(1, 2, 5)) {
    ## M = 1 is the only case where the two rules agree
    expect_equal(1 * qchisq(1 - alpha, df),
                 qchisq(1 - alpha / 1, df), tolerance = 1e-10)

    for (m in c(5, 10, 30, 50)) {
      stat_rule <- m * qchisq(1 - alpha, df)   # correct
      alpha_rule <- qchisq(1 - alpha / m, df)  # the old, wrong rule
      expect_gt(stat_rule, alpha_rule)
    }
  }
})

test_that("rescale_statistic divides the statistic and recomputes p", {
  r <- rescale_statistic(115.2, m = 30, df = 1)
  expect_equal(r$statistic_rescaled, 115.2 / 30)
  expect_equal(r$p_value, pchisq(115.2 / 30, df = 1, lower.tail = FALSE))

  ## a statistic exactly at the corrected boundary yields p == alpha
  alpha <- 0.05
  boundary <- 30 * qchisq(1 - alpha, df = 1)
  expect_equal(rescale_statistic(boundary, m = 30, df = 1)$p_value,
               alpha, tolerance = 1e-8)

  expect_error(rescale_statistic(1, m = 0, df = 1))
  expect_error(rescale_statistic(1, m = 5, df = 0))
})

test_that("degrees of freedom can be recovered and verified", {
  for (df in c(1, 2, 3, 7)) {
    for (stat in c(5, 20, 100)) {
      p  <- pchisq(stat, df = df, lower.tail = FALSE)
      df_hat <- ctreeMI:::.recover_df(stat, p)
      expect_equal(df_hat, df, tolerance = 1e-4)
    }
  }
})

test_that("the correction removes splits relative to no correction", {
  skip_on_cran()
  set.seed(42)
  m <- 10
  n <- 300
  ## Pure null: no relationship between predictors and outcome.
  imp <- lapply(seq_len(m), function(i) {
    set.seed(100 + i)
    data.frame(y  = rnorm(n),
               x1 = rnorm(n),
               x2 = rnorm(n))
  })

  fit <- ctree_stacked(y ~ x1 + x2, data = imp, alpha = 0.05,
                       verbose = FALSE)
  info <- attr(fit, "ctreeMI_info")

  expect_identical(info$correction, "statistic/M")
  expect_equal(info$m, m)
  expect_true(info$n_splits_after <= info$n_splits_before)

  ## Under a pure null the corrected tree should be small.
  expect_lt(info$n_splits_after, 3)

  ## node_stats must line up with what was retained
  if (nrow(info$node_stats)) {
    expect_true(all(info$node_stats$statistic_rescaled ==
                      info$node_stats$statistic / m))
    expect_true(all(info$node_stats$p_rescaled >= info$node_stats$p_raw))
  }
})

test_that("ctreeMI_info no longer advertises an applied alpha", {
  skip_on_cran()
  set.seed(1)
  imp <- lapply(1:5, function(i) {
    set.seed(i)
    data.frame(y = rnorm(120), x = rnorm(120))
  })
  fit  <- ctree_stacked(y ~ x, data = imp, verbose = FALSE)
  info <- attr(fit, "ctreeMI_info")

  expect_null(info$alpha_applied)
  expect_null(info$alpha_nominal)
  expect_equal(info$alpha, 0.05)
})
