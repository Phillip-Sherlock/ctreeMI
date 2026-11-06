# Simulation study for the within-cluster design (ctreeMI 1.3.0)
#
# Usage:  Rscript within_cluster_simulation.R [reps] [outfile] [truths]
#         truths is a comma-separated subset of null,real,weak (default all)
#
# Factors
#   G        number of cohorts: 5, 8, 12, 20
#   sizes    equal (100 per cohort) or unequal (sizes vary about 5:1,
#            same total)
#   truth    null (cohort intercepts only), real (within-cohort effect of
#            x1 of 1.2 SD, plus cohort intercepts) or weak (0.5 SD)
#
# Every data set has cohort intercepts with SD 0.6, a cohort-proxy
# predictor x4 = cohort intercept + noise, which has no within-cohort
# relation to the outcome but a real between-cohort one, and a cohort-level
# noise predictor x5 = a cohort constant + noise, unrelated to the outcome.
# x2 and x3 are individual-level noise. Predictors are 15 to 30% MAR, the
# outcome 10%.
#
# Designs compared on the same data
#   A  within-cluster design as shipped: strata = cohort, tree grown on the
#      cohort-centered outcome, confirmed robust / with fixed effects /
#      per cohort  (discover_confirm(cluster = , within = TRUE))
#   B  the same split, but the tree grown on the raw outcome and confirmed
#      with cluster-robust variance only (what 1.2.0 allowed)
#   C  whole cohorts split between halves, G >= 8 only
#      (discover_confirm(cluster = ))
#
# Outcomes recorded per replicate and design
#   splits       number of splits in the discovery tree
#   proxy        the tree split on x4 (cohort proxy) anywhere
#   noise        the tree split on x2, x3 or x5 anywhere
#   root_x1      the root split is on x1 (real structure only)
#   x1_robust    an x1 split was confirmed by the robust test
#   x1_within    ... and by the fixed-effects test (design A only)
#   x1_same      share of cohorts whose own contrast had the sign of the
#                robust contrast, for the root x1 split (design A only)
#   x4_robust    an x4 split was confirmed by the robust test: a real
#                between-cohort effect, which design C may legitimately
#                confirm and design A should not count as within-cluster
#   x4_within    an x4 split was 'confirmed within clusters' (A): an error
#   spur_robust  a split on x2, x3 or x5 was confirmed by the robust test
#   spur_within  a split on x2, x3 or x5 was 'confirmed within clusters' (A)
#   For design A, the *_within outcomes are 0 when the tree had no splits;
#   for B and C they are NA (not applicable).
#
# Results are appended to the CSV one replicate at a time, so a partial
# run is usable.

suppressPackageStartupMessages({ library(ctreeMI); library(mice) })
args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args) >= 1) as.integer(args[1]) else 25L
OUT  <- if (length(args) >= 2) args[2] else "within_cluster_simulation.csv"
TRUTHS <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else c("null", "real", "weak")
M <- 5L

make_data <- function(seed, G, sizes, truth) {
  set.seed(seed)
  n_per <- if (sizes == "equal") rep(100L, G) else {
    w <- exp(stats::rnorm(G, 0, 0.55)); w <- w / sum(w)
    pmax(30L, round(w * 100 * G))
  }
  n <- sum(n_per)
  cohort <- factor(rep(sprintf("c%02d", seq_len(G)), times = n_per))
  u  <- stats::rnorm(G, 0, 0.6)[as.integer(cohort)]
  x1 <- stats::rnorm(n); x2 <- stats::rnorm(n); x3 <- stats::rnorm(n)
  x4 <- u + stats::rnorm(n, 0, 0.3)
  x5 <- stats::rnorm(G)[as.integer(cohort)] + stats::rnorm(n, 0, 0.3)
  eff <- switch(truth, real = 1.2, weak = 0.5, 0)
  y  <- u + eff * (x1 > 0) + stats::rnorm(n)
  d  <- data.frame(y, x1, x2, x3, x4, x5, cohort)
  ys <- as.numeric(scale(y))
  for (v in c("x1", "x2", "x3", "x4", "x5"))
    d[[v]][stats::runif(n) < stats::plogis(-1.0 + 0.5 * ys)] <- NA
  d$y[stats::runif(n) < 0.10] <- NA
  d
}

split_vars <- function(tree) {
  nodes <- ctreeMI:::.internal_nodes(tree)
  vapply(nodes, `[[`, "", "split_var")
}

summarize_confirm <- function(tree, cf, within = NULL, has_within = FALSE) {
  sv <- split_vars(tree)
  noise <- c("x2", "x3", "x5")
  na_w <- if (has_within) 0 else NA
  out <- c(splits = length(sv), proxy = as.numeric("x4" %in% sv),
           noise = as.numeric(any(sv %in% noise)),
           root_x1 = as.numeric(length(sv) > 0 && sv[1] == "x1"),
           x1_robust = 0, x1_within = na_w, x1_same = NA,
           x4_robust = 0, x4_within = na_w,
           spur_robust = 0, spur_within = na_w)
  if (is.null(cf) || is.null(cf$splits)) return(out)
  s <- cf$splits
  ok <- !is.na(s$p_adj) & s$p_adj < 0.05
  out["x1_robust"]   <- as.numeric(any(ok & s$split_var == "x1"))
  out["x4_robust"]   <- as.numeric(any(ok & s$split_var == "x4"))
  out["spur_robust"] <- as.numeric(any(ok & s$split_var %in% noise))
  if (!is.null(within)) {
    ws <- within$summary
    cw <- ws$reading == "confirmed within clusters"
    out["x1_within"]   <- as.numeric(any(cw & ws$split_var == "x1"))
    out["x4_within"]   <- as.numeric(any(cw & ws$split_var == "x4"))
    out["spur_within"] <- as.numeric(any(cw & ws$split_var %in% noise))
    r <- which(ws$split_var == "x1" & ws$depth == 1)
    if (length(r) == 1 && !is.na(ws$clusters_tested[r]) && ws$clusters_tested[r] > 0)
      out["x1_same"] <- ws$clusters_same_sign[r] / ws$clusters_tested[r]
  }
  out
}

one_rep <- function(r, G, sizes, truth) {
  d <- make_data(seed = 10000 * G + 1000 * (sizes == "unequal") +
                   100 * c(null = 0, real = 1, weak = 2)[[truth]] + r,
                 G = G, sizes = sizes, truth = truth)
  f <- y ~ x1 + x2 + x3 + x4 + x5
  rows <- list()

  # A and B share the split and the imputations
  parts <- split_holdout(d, strata = "cohort", seed = r)
  imp_d <- mice(parts$discover, m = M, printFlag = FALSE, seed = r + 1)
  imp_c <- mice(parts$confirm,  m = M, printFlag = FALSE, seed = r + 2)

  # A confirmation that cannot be run (a terminal node too small to test in
  # the confirmation half) counts as nothing confirmed.
  try_confirm <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) NULL)

  tree_a <- ctree_stacked(f, data = imp_d, center_by = "cohort", verbose = FALSE)
  cw_a <- if (length(partykit::nodeids(tree_a, terminal = TRUE)) > 1)
    try_confirm(confirm_within(tree_a, imp_c, cluster = "cohort")) else NULL
  rows$A <- summarize_confirm(tree_a, if (is.null(cw_a)) NULL else cw_a$robust, cw_a,
                              has_within = TRUE)

  tree_b <- ctree_stacked(f, data = imp_d, verbose = FALSE)
  cf_b <- if (length(partykit::nodeids(tree_b, terminal = TRUE)) > 1)
    try_confirm(confirm_ctreeMI(tree_b, imp_c, cluster = "cohort")) else NULL
  rows$B <- summarize_confirm(tree_b, cf_b)

  if (G >= 8) {
    # whole cohorts split between halves, as discover_confirm(cluster = ) does,
    # with the cohort identifier kept out of the imputation model
    parts_c <- split_holdout(d, cluster = "cohort", seed = r)
    pm <- mice::make.predictorMatrix(parts_c$discover)
    pm[, "cohort"] <- 0L; pm["cohort", ] <- 0L
    imp_cd <- mice(parts_c$discover, m = M, printFlag = FALSE, seed = r,      predictorMatrix = pm)
    imp_cc <- mice(parts_c$confirm,  m = M, printFlag = FALSE, seed = r + 1L, predictorMatrix = pm)
    tree_c <- ctree_stacked(f, data = imp_cd, verbose = FALSE)
    cf_c <- if (length(partykit::nodeids(tree_c, terminal = TRUE)) > 1)
      try_confirm(confirm_ctreeMI(tree_c, imp_cc, cluster = "cohort")) else NULL
    rows$C <- summarize_confirm(tree_c, cf_c)
  }
  do.call(rbind, lapply(names(rows), function(k)
    data.frame(G = G, sizes = sizes, truth = truth, rep = r, design = k,
               as.list(rows[[k]]), stringsAsFactors = FALSE)))
}

grid <- expand.grid(G = c(5L, 8L, 12L, 20L), sizes = c("equal", "unequal"),
                    truth = TRUTHS, stringsAsFactors = FALSE)
done <- if (file.exists(OUT)) utils::read.csv(OUT, stringsAsFactors = FALSE) else NULL
for (r in seq_len(REPS)) {
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    if (!is.null(done) && any(done$G == g$G & done$sizes == g$sizes &
                              done$truth == g$truth & done$rep == r)) next
    res <- tryCatch(one_rep(r, g$G, g$sizes, g$truth),
                    error = function(e) {
                      message(sprintf("rep %d G=%d %s %s: %s", r, g$G, g$sizes, g$truth,
                                      conditionMessage(e))); NULL })
    if (!is.null(res))
      utils::write.table(res, OUT, sep = ",", row.names = FALSE,
                         col.names = !file.exists(OUT), append = file.exists(OUT))
  }
  message(sprintf("[%s] replicate %d of %d done", format(Sys.time(), "%H:%M:%S"), r, REPS))
}

# ---- summary shipped with the package ---------------------------------------
# One row per design, cohort-size condition, cohort count and truth, with the
# mean of each outcome over replicates (NA where not applicable). This is
# what the "Few clusters" vignette reads; the full results stay in the
# repository under simulations/results/.
res <- utils::read.csv(OUT, stringsAsFactors = FALSE)
res <- res[!duplicated(res[c("G", "sizes", "truth", "rep", "design")]), ]
metrics <- c("splits", "proxy", "noise", "root_x1", "x1_robust", "x1_within",
             "x1_same", "x4_robust", "x4_within", "spur_robust", "spur_within")
rate <- function(x) if (all(is.na(x))) NA_real_ else round(mean(x, na.rm = TRUE), 3)
summ <- aggregate(res[metrics], by = res[c("truth", "design", "sizes", "G")], FUN = rate)
summ$n_reps <- aggregate(res$rep, by = res[c("truth", "design", "sizes", "G")],
                         FUN = function(r) length(unique(r)))$x
summ <- summ[order(summ$truth, summ$design, summ$sizes, summ$G),
             c("truth", "design", "sizes", "G", "n_reps", metrics)]
utils::write.csv(summ, sub("\\.csv$", "_summary.csv", OUT), row.names = FALSE)
