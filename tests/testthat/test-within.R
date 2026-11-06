make_cluster_data <- function(G = 6, n_per = 60, seed = 5, binary = FALSE) {
  set.seed(seed)
  n <- G * n_per
  cohort <- factor(rep(sprintf("c%02d", seq_len(G)), each = n_per))
  u  <- rnorm(G, 0, 0.6)[as.integer(cohort)]
  x1 <- rnorm(n); x2 <- rnorm(n)
  x3 <- u + rnorm(n, 0, 0.3)               # tracks the cluster, no within-cluster effect
  eta <- u + 1.5 * (x1 > 0)
  y <- if (binary) factor(rbinom(n, 1, plogis(eta - 0.5)), labels = c("no", "yes"))
       else eta + rnorm(n)
  d <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, cohort = cohort)
  d$x1[sample(n, round(0.15 * n))] <- NA
  d$x2[sample(n, round(0.10 * n))] <- NA
  d
}

impute_halves <- function(parts, m = 3) {
  list(d = mice::mice(parts$discover, m = m, printFlag = FALSE, seed = 11),
       c = mice::mice(parts$confirm,  m = m, printFlag = FALSE, seed = 12))
}

test_that("center_outcome removes group means within each imputed dataset", {
  imps <- lapply(1:3, function(i) {
    g <- factor(rep(c("a", "b", "c"), each = 20))
    data.frame(g = g, x = rnorm(60), y = rnorm(60) + c(a = 0, b = 1, c = 2)[as.character(g)])
  })
  cen <- center_outcome(imps, outcome = "y", by = "g")
  expect_length(cen, 3)
  for (d in cen) expect_equal(as.vector(tapply(d$y, d$g, mean)), c(0, 0, 0))
  expect_equal(cen[[1]]$x, imps[[1]]$x)   # other columns untouched

  # two-level factor outcome is coded 0/1 before centering
  imps_b <- lapply(imps, function(d) { d$y <- factor(d$y > 1, labels = c("no", "yes")); d })
  cen_b <- center_outcome(imps_b, outcome = "y", by = "g")
  expect_true(is.numeric(cen_b[[1]]$y))
  expect_equal(as.vector(tapply(cen_b[[1]]$y, cen_b[[1]]$g, mean)), c(0, 0, 0))

  # logical outcomes are coded 0/1, as in confirm_ctreeMI()
  imps_l <- lapply(imps, function(d) { d$y <- d$y > 1; d })
  cen_l <- center_outcome(imps_l, outcome = "y", by = "g")
  expect_equal(as.vector(tapply(cen_l[[1]]$y, cen_l[[1]]$g, mean)), c(0, 0, 0))

  # several outcomes at once
  imps_2 <- lapply(imps, function(d) { d$y2 <- d$y * 2 + 1; d })
  cen_2 <- center_outcome(imps_2, outcome = c("y", "y2"), by = "g")
  expect_equal(as.vector(tapply(cen_2[[1]]$y2, cen_2[[1]]$g, mean)), c(0, 0, 0))

  expect_error(center_outcome(imps, outcome = "nope", by = "g"))
  expect_error(center_outcome(imps, outcome = "y", by = "nope"))
  # a three-level unordered factor is expanded into centered indicators
  imps_3 <- lapply(imps, function(d) { d$y <- factor(sample(c("p", "q", "r"), 60, TRUE)); d })
  cen_3 <- center_outcome(imps_3, outcome = "y", by = "g")
  expect_identical(attr(cen_3, "outcome_columns")$y, c("y.q", "y.r"))
  expect_false("y" %in% names(cen_3[[1]]))
  expect_equal(as.vector(tapply(cen_3[[1]]$y.q, cen_3[[1]]$g, mean)), c(0, 0, 0))
})

test_that("ctree_stacked centers the outcome when center_by is given", {
  d <- make_cluster_data()
  parts <- split_holdout(d, strata = "cohort", seed = 1)
  imp <- mice::mice(parts$discover, m = 3, printFlag = FALSE, seed = 11)
  tree <- ctree_stacked(y ~ x1 + x2 + x3, data = imp, center_by = "cohort",
                        verbose = FALSE)
  expect_s3_class(tree, "ctreeMI")
  expect_identical(attr(tree, "ctreeMI_info")$center_by, "cohort")
  expect_warning(
    ctree_stacked(y ~ x1 + cohort, data = imp, center_by = "cohort", verbose = FALSE),
    "also appears as a predictor")
  # a `.` formula must not be able to split on the grouping column
  tree_dot <- ctree_stacked(y ~ ., data = imp, center_by = "cohort", verbose = FALSE)
  expect_false("cohort" %in% names(partykit::data_party(tree_dot)))
})

test_that("confirm_ctreeMI accepts fixed effects and tests only the contrast", {
  d <- make_cluster_data()
  parts <- split_holdout(d, strata = "cohort", seed = 2)
  imp <- impute_halves(parts)
  tree <- ctree_stacked(y ~ x1 + x2 + x3, data = imp$d, center_by = "cohort",
                        verbose = FALSE)
  skip_if(length(partykit::nodeids(tree, terminal = TRUE)) < 2L,
          "discovery tree did not split in this replicate")
  cf <- suppressWarnings(confirm_ctreeMI(tree, imp$c, cluster = "cohort", fixed = "cohort"))
  expect_s3_class(cf, "ctreeMI_confirm")
  expect_identical(cf$fixed, "cohort")
  tested <- !is.na(cf$splits$p)
  expect_true(any(tested))
  expect_true(all(cf$splits$df1[tested] == 1))          # one contrast per split, not the dummies
  expect_true(all(is.finite(cf$splits$contrast) | is.na(cf$splits$p)))
  expect_true(all(cf$splits$p[!is.na(cf$splits$p)] >= 0 & cf$splits$p[!is.na(cf$splits$p)] <= 1))
  expect_error(confirm_ctreeMI(tree, imp$c, fixed = "nope"), "not found")
  expect_error(confirm_ctreeMI(tree, imp$c, fixed = c("a", "b")), "one column")
})

test_that("confirm_within lines up the three confirmations", {
  d <- make_cluster_data()
  parts <- split_holdout(d, strata = "cohort", seed = 3)
  imp <- impute_halves(parts)
  tree <- ctree_stacked(y ~ x1 + x2 + x3, data = imp$d, center_by = "cohort",
                        verbose = FALSE)
  skip_if(length(partykit::nodeids(tree, terminal = TRUE)) < 2L,
          "discovery tree did not split in this replicate")
  cw <- suppressWarnings(confirm_within(tree, imp$c, cluster = "cohort"))
  expect_s3_class(cw, "ctreeMI_confirm_within")
  expect_s3_class(cw$robust, "ctreeMI_confirm")
  expect_s3_class(cw$within, "ctreeMI_confirm")
  expect_null(cw$robust$fixed)
  expect_identical(cw$within$fixed, "cohort")
  expect_equal(cw$n_clusters, 6)
  expect_length(cw$by_cluster, 6)
  s <- cw$summary
  expect_equal(nrow(s), nrow(cw$robust$splits))
  expect_true(all(s$reading %in% c("confirmed within clusters",
                                   "between-cluster difference", "not confirmed")))
  expect_true(all(s$clusters_tested <= 6))
  expect_true(all(s$clusters_same_sign <= s$clusters_tested))
  expect_output(print(cw), "three ways")

  cw2 <- suppressWarnings(confirm_within(tree, imp$c, cluster = "cohort", by_cluster = FALSE))
  expect_null(cw2$by_cluster)
  expect_true(all(is.na(cw2$summary$clusters_tested)))
  expect_error(confirm_within(tree, imp$c), "required")
})

test_that("confirm_within works for a binary outcome", {
  d <- make_cluster_data(binary = TRUE, seed = 8)
  parts <- split_holdout(d, strata = "cohort", seed = 4)
  imp <- impute_halves(parts)
  tree <- ctree_stacked(y ~ x1 + x2 + x3, data = imp$d, center_by = "cohort",
                        verbose = FALSE)
  skip_if(length(partykit::nodeids(tree, terminal = TRUE)) < 2L,
          "discovery tree did not split in this replicate")
  cw <- suppressWarnings(confirm_within(tree, imp$c, cluster = "cohort"))
  expect_s3_class(cw, "ctreeMI_confirm_within")
  expect_identical(unname(cw$robust$types), "binary")
  expect_identical(unname(cw$within$types), "binary")
})

test_that("discover_confirm runs the within-cluster design end to end", {
  d <- make_cluster_data()
  res <- suppressWarnings(
    discover_confirm(y ~ x1 + x2 + x3, data = d, m = 3, cluster = "cohort",
                     within = TRUE, seed = 21))
  expect_s3_class(res, "ctreeMI_dc")
  # every cluster on both sides of the split
  expect_setequal(unique(res$split$discover$cohort), levels(d$cohort))
  expect_setequal(unique(res$split$confirm$cohort),  levels(d$cohort))
  # cluster kept in the imputation model
  expect_true(any(res$imp_discover$predictorMatrix[, "cohort"] == 1))
  expect_identical(attr(res$tree, "ctreeMI_info")$center_by, "cohort")
  if (!is.null(res$confirmation)) {
    expect_s3_class(res$within, "ctreeMI_confirm_within")
    expect_identical(res$confirmation, res$within$robust)
    expect_output(print(res), "three ways")
  }
  expect_error(discover_confirm(y ~ x1, data = d, m = 3, within = TRUE, seed = 1),
               "requires `cluster`")
  expect_error(discover_confirm(y ~ x1, data = d, m = 3, cluster = "cohort",
                                strata = "x2", within = TRUE, seed = 1),
               "must be NULL or equal")
})

test_that("report_confirm describes the within-cluster design", {
  d <- make_cluster_data()
  res <- suppressWarnings(
    discover_confirm(y ~ x1 + x2 + x3, data = d, m = 3, cluster = "cohort",
                     within = TRUE, seed = 21))
  skip_if(is.null(res$confirmation), "discovery tree did not split in this replicate")
  rp <- report_confirm(res)
  expect_s3_class(rp, "ctreeMI_report")
  expect_match(rp$text, "within each cohort")
  expect_match(rp$text, "centered on its cohort mean")
  expect_match(rp$text, "fixed effects")
  expect_no_match(rp$text, "no cohort contributed to both")
  rp2 <- report_confirm(res$within, tree = res$tree)
  expect_identical(rp2$text, rp$text)
})

test_that("confirm_within handles several outcomes of mixed type", {
  d <- make_cluster_data(seed = 9)
  d$yb <- factor(rbinom(nrow(d), 1, plogis(d$y - 1)), labels = c("no", "yes"))
  parts <- split_holdout(d, strata = "cohort", seed = 6)
  imp <- impute_halves(parts)
  tree <- ctree_stacked(y + yb ~ x1 + x2 + x3, data = imp$d, center_by = "cohort",
                        verbose = FALSE)
  skip_if(length(partykit::nodeids(tree, terminal = TRUE)) < 2L,
          "discovery tree did not split in this replicate")
  cw <- suppressWarnings(confirm_within(tree, imp$c, cluster = "cohort"))
  s <- cw$summary
  expect_setequal(unique(s$outcome), c("y", "yb"))
  expect_equal(nrow(s), 2 * length(unique(s$node_id)))
  expect_identical(unname(cw$robust$types), c("continuous", "binary"))
  expect_identical(unname(cw$within$types), c("continuous", "binary"))
  expect_true(all(s$clusters_same_sign <= s$clusters_tested, na.rm = TRUE))
})
