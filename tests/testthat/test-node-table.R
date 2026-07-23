make_fit <- function(m = 10, n = 200, seed_base = 0) {
  mk <- function(i) {
    set.seed(seed_base + i)
    x1 <- rnorm(n)
    x2 <- factor(sample(c("A", "B", "C"), n, TRUE))
    y  <- 2 * (x1 > 0) + 1.5 * (x2 == "A") + rnorm(n)
    data.frame(y = y, x1 = x1, x2 = x2)
  }
  ctree_stacked(y ~ x1 + x2, data = lapply(seq_len(m), mk), verbose = FALSE)
}

test_that("node_table returns expected structure", {
  skip_if_not_installed("partykit")
  fit <- make_fit()
  nt  <- node_table(fit)

  expect_s3_class(nt, "ctreeMI_nodes")
  expect_s3_class(nt, "data.frame")
  for (col in c("node_id", "depth", "n_stacked", "effective_n",
                "path", "conditions")) {
    expect_true(col %in% names(nt))
  }
  expect_type(nt$conditions, "list")
})

test_that("effective_n equals n_stacked / M and totals recover original n", {
  skip_if_not_installed("partykit")
  m   <- 10
  n   <- 200
  fit <- make_fit(m = m, n = n)
  nt  <- node_table(fit)

  expect_equal(nt$effective_n, round(nt$n_stacked / m, 1))
  expect_equal(sum(nt$n_stacked), m * n)
  expect_equal(sum(nt$effective_n), n, tolerance = 0.5)
})

test_that("node ids match partykit terminal node ids", {
  skip_if_not_installed("partykit")
  fit <- make_fit()
  expect_equal(node_table(fit)$node_id,
               as.integer(partykit::nodeids(fit, terminal = TRUE)))
})

test_that("terminal_only = FALSE includes inner nodes and root", {
  skip_if_not_installed("partykit")
  fit  <- make_fit()
  all_ <- node_table(fit, terminal_only = FALSE)

  expect_equal(all_$node_id, as.integer(partykit::nodeids(fit)))
  expect_true(1L %in% all_$node_id)
  expect_equal(all_$depth[all_$node_id == 1L], 0L)
  expect_equal(all_$path[all_$node_id == 1L], "root")
  # root holds every stacked row
  expect_equal(max(all_$n_stacked), attr(fit, "ctreeMI_info")$n_stacked)
})

test_that("depth equals number of conditions in the path", {
  skip_if_not_installed("partykit")
  fit <- make_fit()
  nt  <- node_table(fit)
  expect_equal(nt$depth, as.integer(lengths(nt$conditions)))
})

test_that("max_levels truncates long factor level sets", {
  skip_if_not_installed("partykit")
  mk <- function(i) {
    set.seed(i)
    n <- 300
    g <- factor(sample(paste0("G", 1:6), n, TRUE))
    x <- rnorm(n)
    y <- 2 * (g %in% c("G1", "G2")) + x + rnorm(n)
    data.frame(y = y, x = x, g = g)
  }
  fit <- ctree_stacked(y ~ x + g, data = lapply(1:10, mk), verbose = FALSE)

  full  <- node_table(fit)
  trunc <- node_table(fit, max_levels = 1)
  expect_true(any(grepl("more\\}", trunc$path)))
  expect_false(any(grepl("more\\}", full$path)))
})

test_that("binary factor outcomes are summarized as proportions", {
  skip_if_not_installed("partykit")
  mk <- function(i) {
    set.seed(i)
    n <- 300
    x <- rnorm(n)
    y <- factor(rbinom(n, 1, plogis(1.5 * x)),
                levels = c(0, 1), labels = c("No", "Yes"))
    data.frame(y = y, x = x)
  }
  fit <- ctree_stacked(y ~ x, data = lapply(1:10, mk), verbose = FALSE)
  nt  <- node_table(fit)

  prop_col <- grep("^y_prop_", names(nt), value = TRUE)
  expect_length(prop_col, 1L)
  expect_true(all(nt[[prop_col]] >= 0 & nt[[prop_col]] <= 1))
})

test_that("node_table handles a tree with no splits", {
  skip_if_not_installed("partykit")
  flat <- lapply(1:5, function(i) {
    set.seed(i)
    data.frame(y = rnorm(60), x = rnorm(60))
  })
  fit <- ctree_stacked(y ~ x, data = flat, verbose = FALSE)
  nt  <- node_table(fit)

  expect_equal(nrow(nt), 1L)
  expect_equal(nt$depth, 0L)
  expect_equal(nt$path, "root")
  expect_equal(nt$effective_n, 60)
})

test_that("node_table rejects non-party input", {
  expect_error(node_table(data.frame(x = 1)))
})

test_that("report_ctreeMI returns text and matching quantities", {
  skip_if_not_installed("partykit")
  m   <- 10
  n   <- 200
  fit <- make_fit(m = m, n = n)
  rep <- report_ctreeMI(fit)
  nt  <- node_table(fit)

  expect_s3_class(rep, "ctreeMI_report")
  expect_type(rep$text, "character")
  expect_length(rep$text, 1L)

  expect_equal(rep$m, m)
  expect_equal(rep$n_original, n)
  expect_equal(rep$n_stacked, m * n)
  expect_equal(rep$alpha_applied, 0.05 / m)
  expect_equal(rep$n_terminal, nrow(nt))
  expect_equal(rep$max_depth, max(nt$depth))

  # key quantities appear in the prose
  expect_true(grepl(paste0("M = ", m), rep$text, fixed = TRUE))
  expect_true(grepl("Stack/M", rep$text, fixed = TRUE))
})

test_that("report_ctreeMI rejects non-ctreeMI input", {
  skip_if_not_installed("partykit")
  expect_error(report_ctreeMI(data.frame(x = 1)))
})

test_that("print methods run without error", {
  skip_if_not_installed("partykit")
  fit <- make_fit()
  expect_output(print(node_table(fit)), "effective_n")
  expect_output(print(report_ctreeMI(fit)), "Stack/M")
})
