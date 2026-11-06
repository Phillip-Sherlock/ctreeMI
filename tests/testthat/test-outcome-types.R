make_typed_data <- function(G = 6, n_per = 80, seed = 13) {
  set.seed(seed)
  n <- G * n_per
  cohort <- factor(rep(sprintf("c%02d", seq_len(G)), each = n_per))
  u  <- rnorm(G, 0, 0.5)[as.integer(cohort)]
  x1 <- rnorm(n); x2 <- rnorm(n)
  eta <- u + 1.6 * (x1 > 0)
  lat <- eta + rnorm(n)
  y_ord <- cut(lat, c(-Inf, -0.3, 0.9, Inf), labels = c("low", "mid", "high"), ordered_result = TRUE)
  # nominal: level "b" more likely when x1 > 0, "c" more likely in high-u cohorts
  pb <- exp(1.5 * (x1 > 0)); pc <- exp(u); pa <- 1
  y_nom <- factor(apply(cbind(pa, pb, pc), 1L, function(w) sample(c("a", "b", "c"), 1L, prob = w)))
  d <- data.frame(y_ord = y_ord, y_nom = y_nom, x1 = x1, x2 = x2, cohort = cohort)
  d$x1[sample(n, round(0.15 * n))] <- NA
  d
}

impute_both <- function(d, m = 3) {
  parts <- split_holdout(d, strata = "cohort", seed = 2)
  list(parts = parts,
       d = mice::mice(parts$discover, m = m, printFlag = FALSE, seed = 31),
       c = mice::mice(parts$confirm,  m = m, printFlag = FALSE, seed = 32))
}

test_that("outcome types are resolved for ordered and unordered factors", {
  d <- make_typed_data()
  tp <- ctreeMI:::.resolve_types(d, c("y_ord", "y_nom", "x1"), "auto")
  expect_identical(unname(tp), c("ordinal", "nominal", "continuous"))
  expect_error(ctreeMI:::.resolve_types(d, "y_ord", "count"), "entries must be")
})

test_that("cluster-robust covariance for multinom agrees with sandwich on a two-level problem", {
  skip_if_not_installed("nnet"); skip_if_not_installed("sandwich")
  set.seed(4)
  n <- 400; g <- factor(rep(1:8, each = 50)); x <- rnorm(n); s <- factor(x > 0)
  y <- factor(rbinom(n, 1, plogis(-0.3 + 1.2 * (x > 0) + rnorm(8)[as.integer(g)])),
              labels = c("no", "yes"))
  d <- data.frame(y = y, s = s, g = g)
  fm <- nnet::multinom(y ~ s, data = d, trace = FALSE, Hess = TRUE)
  fg <- stats::glm(y ~ s, data = d, family = binomial())
  for (type in c("HC0", "HC1")) {
    Vm <- ctreeMI:::.vcov_multinom(fm, d, cl = d$g, vcov_type = type)
    Vg <- sandwich::vcovCL(fg, cluster = d$g, type = type)
    expect_equal(unname(Vm), unname(Vg), tolerance = 1e-4)
  }
  # model-based covariance passes through unchanged
  expect_equal(ctreeMI:::.vcov_multinom(fm, d), stats::vcov(fm))
  expect_error(ctreeMI:::.vcov_multinom(fm, d, cl = d$g, vcov_type = "HC3"), "HC0")
})

test_that("confirm_ctreeMI handles an ordinal outcome", {
  skip_if_not_installed("MASS")
  d  <- make_typed_data()
  im <- impute_both(d)
  tree <- ctree_stacked(y_ord ~ x1 + x2, data = im$d, verbose = FALSE)
  skip_if(length(partykit::nodeids(tree, terminal = TRUE)) < 2L, "no split")
  cf <- confirm_ctreeMI(tree, im$c)
  expect_identical(unname(cf$types), "ordinal")
  expect_true(all(cf$splits$df1[!is.na(cf$splits$p)] == 1))
  expect_true(all(cf$nodes$level %in% c("low", "mid", "high")))
  expect_equal(nrow(cf$nodes), 3 * length(partykit::nodeids(tree, terminal = TRUE)))
  # proportions within a node sum to one
  sums <- tapply(cf$nodes$estimate, cf$nodes$node_id, sum)
  expect_equal(as.vector(sums[!is.na(sums)]), rep(1, sum(!is.na(sums))), tolerance = 1e-6)
  expect_output(print(cf), "cumulative odds")
  cfc <- suppressWarnings(confirm_ctreeMI(tree, im$c, cluster = "cohort"))
  expect_true(all(is.finite(cfc$splits$contrast[!is.na(cfc$splits$p)])))
  cff <- suppressWarnings(confirm_ctreeMI(tree, im$c, cluster = "cohort", fixed = "cohort"))
  expect_true(all(cff$splits$df1[!is.na(cff$splits$p)] == 1))
  pr <- prune_unconfirmed(tree, cf)
  expect_s3_class(pr, "ctreeMI")
})

test_that("confirm_ctreeMI handles a nominal outcome", {
  skip_if_not_installed("nnet")
  d  <- make_typed_data()
  im <- impute_both(d)
  tree <- ctree_stacked(y_nom ~ x1 + x2, data = im$d, verbose = FALSE)
  skip_if(length(partykit::nodeids(tree, terminal = TRUE)) < 2L, "no split")
  cf <- confirm_ctreeMI(tree, im$c)
  expect_identical(unname(cf$types), "nominal")
  tested <- !is.na(cf$splits$p)
  expect_true(all(cf$splits$df1[tested] == 2))               # two levels against the reference
  expect_true(all(cf$splits$contrast_level[tested] %in% c("b", "c")))
  expect_equal(nrow(cf$nodes), 3 * length(partykit::nodeids(tree, terminal = TRUE)))
  expect_output(print(cf), "nominal")
  cfc <- suppressWarnings(confirm_ctreeMI(tree, im$c, cluster = "cohort"))
  expect_true(all(is.finite(cfc$splits$contrast[!is.na(cfc$splits$p)])))
  expect_error(suppressWarnings(confirm_ctreeMI(tree, im$c, cluster = "cohort", vcov_type = "HC3")),
               NA)   # a failed fit is dropped, not fatal
  cff <- suppressWarnings(confirm_ctreeMI(tree, im$c, cluster = "cohort", fixed = "cohort"))
  expect_true(all(cff$splits$df1[!is.na(cff$splits$p)] == 2))
})

test_that("center_outcome expands a nominal outcome and scores an ordinal one", {
  d  <- make_typed_data()
  imps <- lapply(1:2, function(i) d[complete.cases(d), ])
  cen <- center_outcome(imps, outcome = c("y_ord", "y_nom"), by = "cohort")
  cols <- attr(cen, "outcome_columns")
  expect_identical(cols$y_ord, "y_ord")
  expect_identical(cols$y_nom, c("y_nom.b", "y_nom.c"))
  expect_true(is.numeric(cen[[1]]$y_ord))
  expect_false("y_nom" %in% names(cen[[1]]))
  for (v in c("y_ord", "y_nom.b", "y_nom.c"))
    expect_equal(as.vector(tapply(cen[[1]][[v]], cen[[1]]$cohort, mean)), rep(0, 6))
})

test_that("the within-cluster design runs for ordinal and nominal outcomes", {
  skip_if_not_installed("MASS"); skip_if_not_installed("nnet")
  d <- make_typed_data()
  for (f in list(y_ord ~ x1 + x2, y_nom ~ x1 + x2)) {
    res <- suppressWarnings(
      discover_confirm(f, data = d, m = 3, cluster = "cohort", within = TRUE, seed = 7))
    info <- attr(res$tree, "ctreeMI_info")
    expect_identical(info$formula, f)
    if (!is.null(res$confirmation)) {
      expect_s3_class(res$within, "ctreeMI_confirm_within")
      expect_identical(unique(res$within$summary$outcome), all.vars(f[[2L]]))
      expect_output(print(res$within), "three ways")
    }
  }
  # the nominal tree was grown on indicator columns but keeps the user's formula
  res_n <- suppressWarnings(
    discover_confirm(y_nom ~ x1 + x2, data = d, m = 3, cluster = "cohort", within = TRUE, seed = 8))
  expect_identical(all.vars(attr(res_n$tree, "ctreeMI_info")$discovery_formula[[2L]]),
                   c("y_nom.b", "y_nom.c"))
})
