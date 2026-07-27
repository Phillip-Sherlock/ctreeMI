## Regression tests for the 1.0.0 fixes.
##
## What these guard against, all confirmed against partykit 1.2-20:
##
##  1. The stored p-value underflows to exactly 0 on stacked data, which made
##     0.3.0 abort with "could not recover degrees of freedom".
##  2. Degrees of freedom must follow the outcome dimension for numeric and
##     ordered predictors, and (L - 1) x outcome dim for unordered factors.
##  3. prune_stackM() must handle a tree fitted with the partykit default
##     testtype = "Bonferroni", where the stored p-values already carry the
##     multiplicity adjustment.
##  4. The correction must reproduce the bottom-up procedure used for the
##     paper when the degrees of freedom are held fixed.

library(partykit)

# ---------------------------------------------------------------------------
# 1. underflow
# ---------------------------------------------------------------------------

test_that("large stacked samples do not abort the correction", {
  skip_on_cran()
  set.seed(1)
  n <- 1200; m <- 30
  imps <- lapply(seq_len(m), function(i) {
    x <- rnorm(n); g <- factor(sample(letters[1:4], n, TRUE))
    data.frame(x = x, g = g,
               y1 = rbinom(n, 1, plogis(-1 + 1.1 * (x > 0.5) + 0.7 * (g == "a"))),
               y2 = rbinom(n, 1, plogis(-1 + 0.9 * (x > 0.5))))
  })
  fit <- expect_no_error(
    ctree_stacked(y1 + y2 ~ x + g, data = imps, alpha = 0.05, verbose = FALSE))
  ns <- attr(fit, "ctreeMI_info")$node_stats
  expect_gt(nrow(ns), 0)
  ## partykit's own p-value underflows here; the corrected one must not
  expect_true(any(ns$statistic > 1000))
  expect_true(all(is.finite(ns$p_rescaled_bonferroni)))
})

test_that("the corrected p-value stays positive where partykit's underflows", {
  set.seed(2); n <- 6000
  d <- data.frame(x = rnorm(n))
  d$y1 <- rnorm(n) + 2.5 * (d$x > 0)
  d$y2 <- rnorm(n) + 2.0 * (d$x > 0)
  g <- ctree(y1 + y2 ~ x, data = d,
             control = ctree_control(maxdepth = 1, testtype = "Univariate"))
  cr <- nodeapply(g, 1, function(nd) info_node(nd)$criterion)[[1]]
  expect_gt(cr["statistic", 1], 1000)
  expect_equal(cr["p.value", 1], 0)              # partykit underflows
  out <- prune_stackM(g, m = 20, verbose = FALSE)
  expect_gt(out$node_stats$p_rescaled_bonferroni[1], 0)
  expect_true(is.finite(out$node_stats$p_rescaled_bonferroni[1]))
  expect_identical(out$node_stats$df_source[1], "structural")
})

# ---------------------------------------------------------------------------
# 2. degrees of freedom
# ---------------------------------------------------------------------------

test_that("df follows the outcome dimension for numeric predictors", {
  set.seed(3); n <- 1200
  d <- data.frame(x = rnorm(n))
  d$y1 <- rnorm(n) + 0.45 * d$x
  d$y2 <- rnorm(n) + 0.40 * d$x
  d$y3 <- rnorm(n) + 0.35 * d$x
  dfs <- vapply(list(y1 ~ x, y1 + y2 ~ x, y1 + y2 + y3 ~ x), function(f) {
    g <- ctree(f, data = d, control = ctree_control(maxdepth = 1,
                                                    testtype = "Univariate"))
    prune_stackM(g, m = 2, verbose = FALSE)$node_stats$df[1]
  }, numeric(1))
  expect_equal(dfs, c(1, 2, 3))
})

test_that("df for an unordered factor is (levels - 1) x outcome dimension", {
  set.seed(4); L <- 6; n <- 3000
  d <- data.frame(site = factor(sample(seq_len(L), n, TRUE)), x = rnorm(n))
  b <- rnorm(L, 0, 0.5)
  d$y1 <- rnorm(n) + b[as.integer(d$site)]
  d$y2 <- rnorm(n) + 0.5 * b[as.integer(d$site)]
  g <- ctree(y1 + y2 ~ site + x, data = d,
             control = ctree_control(maxdepth = 1, testtype = "Univariate"))
  cand <- attr(prune_stackM(g, m = 5, verbose = FALSE)$node_stats, "candidates")
  root <- cand[cand$node_id == 1, ]
  expect_equal(root$df[root$variable == "site"], (L - 1) * 2)
  expect_equal(root$df[root$variable == "x"], 2)
})

test_that("the df search grid covers what the data can produce", {
  set.seed(5); n <- 500
  d <- data.frame(g = factor(sample(letters, n, TRUE)), x = rnorm(n))
  d$y1 <- rnorm(n); d$y2 <- rnorm(n); d$y3 <- rnorm(n)
  g <- ctree(y1 + y2 + y3 ~ g + x, data = d, control = ctree_control(maxdepth = 1))
  expect_gte(ctreeMI:::.df_upper(g, 3), (26 - 1) * 3)
})

# ---------------------------------------------------------------------------
# 3. both multiplicity settings
# ---------------------------------------------------------------------------

test_that("df is recovered correctly under Bonferroni as well as Univariate", {
  set.seed(6); n <- 3000
  d <- data.frame(x = rnorm(n), g = factor(sample(letters[1:5], n, TRUE)))
  d$y1 <- rnorm(n) + 0.35 * (d$x > 0)
  d$y2 <- rnorm(n) + 0.30 * (d$x > 0)
  for (tt in c("Univariate", "Bonferroni")) {
    g <- ctree(y1 + y2 ~ x + g, data = d,
               control = ctree_control(maxdepth = 1, testtype = tt))
    cand <- attr(prune_stackM(g, m = 3, verbose = FALSE)$node_stats, "candidates")
    root <- cand[cand$node_id == 1, ]
    expect_equal(root$df[root$variable == "x"], 2, info = tt)
    expect_equal(root$df[root$variable == "g"], 4 * 2, info = tt)
    expect_true(all(root$df_source == "recovered"), info = tt)
  }
})

test_that("check_stackM_extraction passes under both settings", {
  expect_true(check_stackM_extraction(verbose = FALSE))
})

# ---------------------------------------------------------------------------
# 4. agreement with the procedure used for the paper
# ---------------------------------------------------------------------------

test_that("fixed df reproduces the published bottom-up pruning", {
  skip_on_cran()
  set.seed(7); n <- 4000
  d <- data.frame(x = rnorm(n), z = rnorm(n), w = runif(n))
  d$y1 <- rnorm(n) + 0.55 * (d$x > 0) + 0.35 * d$z
  d$y2 <- rnorm(n) + 0.45 * (d$x > 0)
  grown <- ctree(y1 + y2 ~ x + z + w, data = d,
                 control = ctree_control(testtype = "Univariate"))

  ## The analysis code for the paper: divide the statistic by m, recompute the
  ## adjusted p with df = 2, prune terminal pairs, repeat.
  reference <- function(x, alpha, m) {
    repeat {
      term  <- nodeids(x, terminal = TRUE)
      inner <- setdiff(nodeids(x), term)
      sct <- nodeapply(x, nodeids(x), function(n) info_node(n)$criterion)
      names(sct) <- nodeids(x)
      kids <- nodeapply(x, inner, function(n) vapply(kids_node(n), id_node, 1L))
      tp   <- Filter(function(p) all(p %in% term), kids)
      n2p <- c()
      for (p in tp) {
        cr <- sct[[min(p) - 1]]
        st <- cr["statistic", ] / m
        pv <- 1 - pchisq(st, 2)^ncol(cr)
        if (!any(pv < alpha)) n2p <- c(n2p, min(p) - 1)
      }
      if (length(n2p) == 0) break else x <- nodeprune(x, n2p)
    }
    x
  }
  for (m in c(2, 5, 10, 30)) {
    a <- reference(grown, alpha = 0.05, m = m)
    b <- ctreeMI:::.prune_stackM(grown, m = m, alpha = 0.05, df = 2,
                                 verbose = FALSE)$tree
    expect_identical(capture.output(print(a)), capture.output(print(b)),
                     info = paste("m =", m))
  }
})

test_that("all candidates are tested, not only the variable split on", {
  set.seed(8); n <- 2000
  d <- data.frame(x = rnorm(n), g = factor(sample(letters[1:4], n, TRUE)))
  d$y1 <- rnorm(n) + 0.30 * d$x + 0.30 * (d$g == "a")
  d$y2 <- rnorm(n) + 0.25 * d$x
  g <- ctree(y1 + y2 ~ x + g, data = d,
             control = ctree_control(maxdepth = 1, testtype = "Univariate"))
  ns   <- prune_stackM(g, m = 4, verbose = FALSE)$node_stats
  cand <- attr(ns, "candidates")
  expect_true(all(c("p_min", "variable_min") %in% names(ns)))
  expect_equal(ns$p_min[1], min(cand$p_corrected[cand$node_id == 1]))
  expect_true(ns$p_min[1] <= ns$p_rescaled_bonferroni[1])
})

test_that("pruning is monotone in m", {
  skip_on_cran()
  set.seed(9); n <- 3000
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  d$y1 <- rnorm(n) + 0.5 * (d$x > 0) + 0.3 * d$z
  d$y2 <- rnorm(n) + 0.4 * (d$x > 0)
  grown <- ctree(y1 + y2 ~ x + z, data = d,
                 control = ctree_control(testtype = "Univariate"))
  w <- vapply(c(1, 2, 5, 10, 40), function(m)
    width(prune_stackM(grown, m = m, verbose = FALSE)$tree), numeric(1))
  expect_equal(w[1], width(grown))
  expect_true(all(diff(w) <= 0))
})

test_that("a node keeps its split when a descendant survives the correction", {
  ## bottom-up compression must not delete a failing node that sits above a
  ## surviving one; that is what distinguishes it from deleting every failure
  skip_on_cran()
  set.seed(10); n <- 4000
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  d$y <- rnorm(n) + 0.5 * (d$x > 0) + 0.25 * (d$z > 1)
  grown <- ctree(y ~ x + z, data = d, control = ctree_control(testtype = "Univariate"))
  out <- prune_stackM(grown, m = 3, verbose = FALSE)
  ns  <- out$node_stats
  kept <- ns$node_id[ns$retained]
  ## every retained node's parent must also be retained
  pm <- ctreeMI:::.parent_map(node_party(grown))
  for (id in kept) {
    p <- pm[[as.character(id)]]
    if (!is.na(p)) expect_true(p %in% kept, info = paste("node", id))
  }
})

# ---------------------------------------------------------------------------
# 5. stacking hygiene
# ---------------------------------------------------------------------------

test_that("the imputation index is never a candidate splitting variable", {
  set.seed(11); n <- 200
  imps <- lapply(1:6, function(i) data.frame(y = rnorm(n), x = rnorm(n)))
  fit <- ctree_stacked(y ~ ., data = imps, verbose = FALSE)
  expect_false(".imp" %in% names(partykit::data_party(fit)))
  cand <- attr(attr(fit, "ctreeMI_info")$node_stats, "candidates")
  if (!is.null(cand) && nrow(cand)) expect_false(any(cand$variable == ".imp"))
})

test_that("minsplit and minbucket refer to original observations", {
  set.seed(12); n <- 400
  imps <- lapply(1:5, function(i) {
    x <- rnorm(n); data.frame(x = x, y = rnorm(n) + 0.8 * (x > 0))
  })
  fit <- ctree_stacked(y ~ x, data = imps, verbose = FALSE, minbucket = 20)
  expect_equal(attr(fit, "ctreeMI_info")$minbucket, 100)
  fit2 <- ctree_stacked(y ~ x, data = imps, verbose = FALSE, minbucket = 20,
                        scale_minsize = FALSE)
  expect_equal(attr(fit2, "ctreeMI_info")$minbucket, 20)
})

test_that("prune_stackM rejects reference distributions it cannot rescale", {
  set.seed(13); n <- 400
  d <- data.frame(x = rnorm(n)); d$y <- rnorm(n) + (d$x > 0)
  expect_error(
    prune_stackM(ctree(y ~ x, data = d,
                       control = ctree_control(teststat = "maximum")), m = 2),
    "quadratic")
  expect_error(
    prune_stackM(ctree(y ~ x, data = d,
                       control = ctree_control(testtype = "MonteCarlo",
                                               nresample = 99)), m = 2),
    "Bonferroni")
  expect_error(prune_stackM(lm(y ~ x, data = d), m = 2), "party")
  fit <- ctree(y ~ x, data = d)
  expect_error(prune_stackM(fit, m = 0), "m` must be")
  expect_error(prune_stackM(fit, m = 2, alpha = 1.5), "alpha")
})
