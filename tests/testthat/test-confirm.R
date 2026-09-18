make_confirm_data <- function(n = 600, seed = 3, signal = 1.2) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                  x3 = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$y <- rnorm(n) + signal * (d$x1 > 0)
  d$x1[sample(n, round(0.15 * n))] <- NA
  d$x2[sample(n, round(0.10 * n))] <- NA
  d
}

impute_list <- function(d, m = 8, seed = 1) {
  imp <- mice::mice(d, m = m, printFlag = FALSE, seed = seed)
  lapply(seq_len(m), function(i) mice::complete(imp, i))
}

test_that("split_holdout partitions rows disjointly and respects prop", {
  d <- make_confirm_data(n = 200)
  p <- split_holdout(d, prop = 0.5, seed = 1)
  expect_s3_class(p, "ctreeMI_split")
  expect_equal(nrow(p$discover) + nrow(p$confirm), nrow(d))
  expect_equal(nrow(p$discover), 100)
  expect_length(intersect(rownames(p$discover), rownames(p$confirm)), 0)
  expect_equal(p$index, sort(p$index))
})

test_that("split_holdout stratifies on a factor", {
  d <- make_confirm_data(n = 300)
  d$g <- factor(rep(c("u", "v", "w"), each = 100))
  p <- split_holdout(d, prop = 0.5, strata = "g", seed = 2)
  expect_equal(as.vector(table(p$discover$g)), c(50, 50, 50))
})

test_that("split_holdout rejects bad inputs", {
  d <- make_confirm_data(n = 50)
  expect_error(split_holdout(d, prop = 1.2))
  expect_error(split_holdout(d, strata = "nope"))
  expect_error(split_holdout(as.matrix(d[, 1:2])))
})

test_that("confirm_ctreeMI returns expected structure", {
  skip_if_not_installed("mice")
  skip_if_not_installed("partykit")
  d  <- make_confirm_data(n = 600)
  p  <- split_holdout(d, seed = 5)
  id <- impute_list(p$discover, seed = 10)
  ic <- impute_list(p$confirm,  seed = 20)
  tree <- ctree_stacked(y ~ x1 + x2 + x3, data = id, verbose = FALSE)
  skip_if(length(partykit::nodeids(tree, terminal = TRUE)) < 2,
          "discovery tree did not split; nothing to confirm")

  cf <- confirm_ctreeMI(tree, ic)
  expect_s3_class(cf, "ctreeMI_confirm")
  expect_true(all(c("test", "nodes", "m", "n_confirm") %in% names(cf)))
  expect_equal(cf$m, 8)
  expect_equal(cf$n_confirm, nrow(p$confirm))
  expect_true(all(c("outcome", "F", "df1", "df2", "p", "riv") %in% names(cf$test)))
  expect_equal(nrow(cf$test), 1)
  expect_true(all(c("node_id", "estimate", "se", "n_confirm", "tested") %in%
                  names(cf$nodes)))
  expect_equal(nrow(cf$nodes), length(partykit::nodeids(tree, terminal = TRUE)))
})

test_that("confirm_ctreeMI returns per-split tests with the right structure", {
  skip_if_not_installed("mice")
  skip_if_not_installed("partykit")
  d  <- make_confirm_data(n = 700, signal = 1.5)
  p  <- split_holdout(d, seed = 4)
  id <- impute_list(p$discover, seed = 15)
  ic <- impute_list(p$confirm,  seed = 25)
  tree <- ctree_stacked(y ~ x1 + x2 + x3, data = id, verbose = FALSE)
  skip_if(length(partykit::nodeids(tree, terminal = TRUE)) < 2)

  cf <- confirm_ctreeMI(tree, ic)
  n_internal <- length(setdiff(partykit::nodeids(tree),
                               partykit::nodeids(tree, terminal = TRUE)))
  expect_false(is.null(cf$splits))
  expect_equal(nrow(cf$splits), n_internal)
  expect_true(all(c("node_id", "depth", "split_var", "F", "p", "p_adj",
                    "n_confirm") %in% names(cf$splits)))
  expect_equal(cf$splits$depth[cf$splits$node_id == 1L], 1L)
  # root split is on the signal variable and should confirm
  root <- cf$splits[cf$splits$node_id == 1L, ]
  expect_equal(root$split_var, "x1")
  expect_lt(root$p_adj, 0.05)
  # Holm-adjusted p is never smaller than raw
  ok <- !is.na(cf$splits$p)
  expect_true(all(cf$splits$p_adj[ok] >= cf$splits$p[ok] - 1e-12))
})

test_that("a spurious downstream split does not confirm", {
  skip_if_not_installed("mice")
  skip_if_not_installed("partykit")
  # Signal only through x1. Any split on x2 or x3 the discovery tree makes
  # is spurious, and its per-split test on independent data should not
  # reject at the adjusted level.
  d  <- make_confirm_data(n = 1000, signal = 1.5, seed = 31)
  p  <- split_holdout(d, seed = 31)
  id <- impute_list(p$discover, seed = 35)
  ic <- impute_list(p$confirm,  seed = 45)
  tree <- ctree_stacked(y ~ x1 + x2 + x3, data = id, verbose = FALSE)
  cf <- confirm_ctreeMI(tree, ic)
  skip_if(is.null(cf$splits) || !any(cf$splits$split_var %in% c("x2", "x3")),
          "discovery tree made no spurious split to test")
  spur <- cf$splits[cf$splits$split_var %in% c("x2", "x3") & !is.na(cf$splits$p_adj), ]
  expect_true(all(spur$p_adj > 0.05))
})

test_that("confirm_ctreeMI detects a real partition", {
  skip_if_not_installed("mice")
  skip_if_not_installed("partykit")
  d  <- make_confirm_data(n = 800, signal = 1.5)
  p  <- split_holdout(d, seed = 6)
  id <- impute_list(p$discover, seed = 30)
  ic <- impute_list(p$confirm,  seed = 40)
  tree <- ctree_stacked(y ~ x1 + x2 + x3, data = id, verbose = FALSE)
  skip_if(length(partykit::nodeids(tree, terminal = TRUE)) < 2)

  cf <- confirm_ctreeMI(tree, ic)
  expect_lt(cf$test$p, 0.01)
  # node estimates should differ in the direction of the signal
  expect_gt(diff(range(cf$nodes$estimate, na.rm = TRUE)), 0.5)
})

test_that("confirm_ctreeMI does not reject on a null partition", {
  skip_if_not_installed("mice")
  skip_if_not_installed("partykit")
  # Build a tree on signal data, then confirm on data with NO signal.
  # The confirmation test should not reject at a high rate.
  d  <- make_confirm_data(n = 800, signal = 1.5, seed = 8)
  p  <- split_holdout(d, seed = 8)
  id <- impute_list(p$discover, seed = 50)
  tree <- ctree_stacked(y ~ x1 + x2 + x3, data = id, verbose = FALSE)
  skip_if(length(partykit::nodeids(tree, terminal = TRUE)) < 2)

  null_conf <- make_confirm_data(n = 400, signal = 0, seed = 9)
  ic <- impute_list(null_conf, seed = 60)
  cf <- confirm_ctreeMI(tree, ic)
  expect_true(is.finite(cf$test$p))
  expect_gt(cf$test$p, 0.001)   # a null partition should not be overwhelming
  # and the node estimates should not separate: both near the overall null mean
  expect_lt(diff(range(cf$nodes$estimate, na.rm = TRUE)), 0.5)
})

test_that("confirm_ctreeMI handles a binary outcome", {
  skip_if_not_installed("mice")
  skip_if_not_installed("partykit")
  set.seed(12)
  n <- 700
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- factor(rbinom(n, 1, plogis(-0.5 + 1.8 * (d$x1 > 0))))
  d$x1[sample(n, 100)] <- NA
  p  <- split_holdout(d, seed = 12)
  id <- impute_list(p$discover, seed = 70)
  ic <- impute_list(p$confirm,  seed = 80)
  tree <- ctree_stacked(y ~ x1 + x2, data = id, verbose = FALSE)
  skip_if(length(partykit::nodeids(tree, terminal = TRUE)) < 2)

  cf <- confirm_ctreeMI(tree, ic)
  expect_equal(unname(cf$types), "binary")
  expect_true(all(cf$nodes$estimate >= 0 & cf$nodes$estimate <= 1, na.rm = TRUE))
})

test_that("confirm_ctreeMI refuses a single-node tree", {
  skip_if_not_installed("mice")
  skip_if_not_installed("partykit")
  d  <- make_confirm_data(n = 300, signal = 0, seed = 14)
  id <- impute_list(d, m = 5, seed = 90)
  tree <- ctree_stacked(y ~ x1 + x2 + x3, data = id, verbose = FALSE)
  skip_if(length(partykit::nodeids(tree, terminal = TRUE)) >= 2,
          "null tree happened to split; cannot test single-node refusal")
  expect_error(confirm_ctreeMI(tree, id), "single terminal node")
})

test_that("discover_confirm runs end to end", {
  skip_if_not_installed("mice")
  skip_if_not_installed("partykit")
  d <- make_confirm_data(n = 600, signal = 1.5)
  res <- discover_confirm(y ~ x1 + x2 + x3, data = d, m = 6, seed = 21)
  expect_s3_class(res, "ctreeMI_dc")
  expect_s3_class(res$tree, "ctreeMI")
  expect_s3_class(res$split, "ctreeMI_split")
  expect_s3_class(res$imp_discover, "mids")
  expect_s3_class(res$imp_confirm, "mids")
  expect_equal(res$imp_discover$m, 6)
  if (!is.null(res$confirmation)) {
    expect_s3_class(res$confirmation, "ctreeMI_confirm")
    expect_output(print(res), "Discover-then-confirm")
  }
})

test_that("confirmation imputations are independent of discovery imputations", {
  skip_if_not_installed("mice")
  d <- make_confirm_data(n = 300)
  res <- discover_confirm(y ~ x1 + x2 + x3, data = d, m = 4, seed = 33)
  # different seeds -> different mids objects; rows do not overlap
  expect_false(identical(res$imp_discover$seed, res$imp_confirm$seed))
  expect_length(intersect(rownames(res$split$discover),
                          rownames(res$split$confirm)), 0)
})

test_that("per-split output carries the rule and a pooled contrast", {
  skip_if_not_installed("mice")
  skip_if_not_installed("partykit")
  d  <- make_confirm_data(n = 700, signal = 1.5)
  p  <- split_holdout(d, seed = 55)
  id <- impute_list(p$discover, seed = 55)
  ic <- impute_list(p$confirm,  seed = 65)
  tree <- ctree_stacked(y ~ x1 + x2 + x3, data = id, verbose = FALSE)
  skip_if(length(partykit::nodeids(tree, terminal = TRUE)) < 2)

  cf <- confirm_ctreeMI(tree, ic)
  expect_true(all(c("rule", "contrast", "lower", "upper") %in% names(cf$splits)))
  root <- cf$splits[cf$splits$node_id == 1L, ]
  expect_match(root$rule, "x1")           # rule names the variable and threshold
  expect_true(is.finite(root$contrast))
  expect_lt(root$lower, root$upper)        # interval is ordered
  expect_true(root$lower > 0 || root$upper < 0)  # excludes zero for a real split
  expect_output(print(cf), "difference")
})

test_that("prune_unconfirmed keeps confirmed splits and drops the rest", {
  skip_if_not_installed("mice")
  skip_if_not_installed("partykit")
  d  <- make_confirm_data(n = 800, signal = 1.5, seed = 71)
  p  <- split_holdout(d, seed = 71)
  id <- impute_list(p$discover, seed = 71)
  ic <- impute_list(p$confirm,  seed = 81)
  tree <- ctree_stacked(y ~ x1 + x2 + x3, data = id, verbose = FALSE)
  skip_if(length(partykit::nodeids(tree, terminal = TRUE)) < 2)

  cf  <- confirm_ctreeMI(tree, ic)
  out <- prune_unconfirmed(tree, cf)
  expect_s3_class(out, "ctreeMI")
  # a real signal should survive, so the pruned tree keeps at least one split
  expect_gte(length(partykit::nodeids(out, terminal = TRUE)), 2)
  # never larger than the tree it came from
  expect_lte(length(partykit::nodeids(out, terminal = TRUE)),
             length(partykit::nodeids(tree, terminal = TRUE)))
  expect_false(is.null(attr(out, "ctreeMI_info")$confirmed))
})

test_that("prune_unconfirmed collapses everything when nothing confirms", {
  skip_if_not_installed("mice")
  skip_if_not_installed("partykit")
  d  <- make_confirm_data(n = 800, signal = 1.5, seed = 73)
  p  <- split_holdout(d, seed = 73)
  id <- impute_list(p$discover, seed = 73)
  tree <- ctree_stacked(y ~ x1 + x2 + x3, data = id, verbose = FALSE)
  skip_if(length(partykit::nodeids(tree, terminal = TRUE)) < 2)

  # confirm against data with no structure at all
  null_c <- make_confirm_data(n = 400, signal = 0, seed = 74)
  cf  <- confirm_ctreeMI(tree, impute_list(null_c, seed = 84))
  out <- prune_unconfirmed(tree, cf)
  expect_equal(length(partykit::nodeids(out, terminal = TRUE)), 1L)
})

test_that("report_confirm produces a methods paragraph", {
  skip_if_not_installed("mice")
  d   <- make_confirm_data(n = 600, signal = 1.5)
  res <- discover_confirm(y ~ x1 + x2 + x3, data = d, m = 6, seed = 77)
  skip_if(is.null(res$confirmation))
  rp <- report_confirm(res)
  expect_s3_class(rp, "ctreeMI_report")
  expect_type(rp$text, "character")
  expect_match(rp$text, "discovery set")
  expect_match(rp$text, "confirmation set")
  expect_match(rp$text, "Li")
  expect_equal(rp$n_discover, nrow(res$split$discover))
  expect_output(print(rp), "discovery set")
})
