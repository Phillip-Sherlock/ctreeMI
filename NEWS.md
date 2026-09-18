# ctreeMI 1.1.0

## New: discover-then-confirm workflow

Three functions implement the workflow recommended in the documentation of
`ctree_stacked()`: treat the tree as discovery, and test its partition on
independent data with a procedure that has valid error control.

* `split_holdout()` partitions a data frame into discovery and confirmation
  sets on original observations, before imputation, so that each half can
  be imputed separately.
* `confirm_ctreeMI()` applies the terminal-node rules of a fitted tree to an
  independently imputed confirmation sample and tests whether the outcome
  differs across the resulting subgroups. The test pools across imputations
  by the Li-Raghunathan-Rubin procedure (`mice::D1()`) rather than
  stacking, so the between-imputation variance that miscalibrates the
  stacked node-level test is properly accounted for. Two families of test
  are returned: an omnibus test that the outcome differs across terminal
  nodes at all, and a per-split test for every internal node, in which the
  confirmation observations inside that node are divided by its own rule
  and the children compared. The per-split tests identify which individual
  splits survive independent data, and are Holm-adjusted across internal
  nodes by default. Pooled per-node estimates are returned alongside.
* `discover_confirm()` runs the full sequence from a raw data frame.
* `prune_unconfirmed()` returns the tree with every split that failed to
  confirm collapsed, so the partition that survived can be reported and
  plotted directly.
* `report_confirm()` generates a methods paragraph describing the design,
  in the idiom of `report_ctreeMI()`.

Each split is reported with the pooled difference between its children and
a confidence interval, alongside the test. The package argues that
node-level p-values should not be read as error rates; it would be
inconsistent to report only p-values here.

The confirmation test addresses two problems at once: the selection effect
of testing subgroups on the data used to find them, and the miscalibration
of the stacked node-level test under outcome-conditioned imputation.

## Documentation now renders as written

The roxygen comments throughout the package use markdown syntax for
emphasis, code, cross-references and section headings, but markdown
processing had never been enabled, so these rendered literally on CRAN.
It is now enabled. A percent sign in `?ctree_stacked` that Rd read as a
comment marker, silently truncating a sentence in the 1.0.1 manual, is now
escaped.

## Corrected package-level documentation

`?ctreeMI` described the correction as dividing the significance threshold
by M. That is threshold rescaling, which versions 0.1.0 and 0.2.0
implemented and which under-corrects by an order of magnitude at M = 30;
the package has applied the correction to the statistic since 0.3.0, and
`?ctree_stacked` says so. The package-level page now describes the
mechanism correctly, no longer characterizes the node-level test as
conservative, and points to the confirmation workflow.

# ctreeMI 1.0.1

Documentation only. No function has changed and results are identical to
1.0.0.

## Corrected description of node-level calibration

`?ctree_stacked` reported sub-nominal, and therefore conservative, type-I
error under MCAR. That holds under the marginal imputation model used in the
simulations of Sherlock et al. (2026), which conditioned on neither the
outcome nor the remaining predictors. It does not hold under
outcome-conditioned imputation, which is recommended practice and is what
`mice()` does by default when the outcome is present in the data frame. In
that setting the node-level test rejects more often than `alpha` implies,
increasingly so as the missingness rate rises.

The documentation now states this, notes that omitting the outcome from the
imputation model reverses the direction rather than restoring calibration,
and points to `node_table()` and to refitting on independent imputations as
the checks to use in place of node-level p-values.

The package DESCRIPTION previously described the result as "a conservative
but interpretable single tree". The word "conservative" has been removed for
the same reason.

## Added a scope section

`?ctree_stacked` now states that the procedure assumes missingness at random,
and that recovery degrades under MNAR for every imputation-based approach.

Supporting simulations are archived at
<https://doi.org/10.5281/zenodo.21939940>.

# ctreeMI 1.0.0

0.3.0 introduced the right correction. This release makes it survive contact
with real stacked data, and derives the degrees of freedom instead of relying
on a p-value that has usually underflowed by the time it is needed.
**Results will differ from 0.3.0.**

## `ctree_stacked()` no longer aborts on large stacked datasets

`partykit` stores node-level p-values on the natural scale, and they underflow
to exactly 0 once the statistic passes roughly 1400 at `df = 2`. On stacked
data that is the normal case, not an edge case: 3,000 rows with M = 30 already
produces a root statistic above 5,000. 0.3.0 recovered the degrees of freedom
by inverting the stored p-value and, failing that, stopped with

    could not recover degrees of freedom at node 1 (statistic = 5187.76,
    p = 0). Refusing to prune on unverified statistics.

so the package could not fit the kind of dataset it was written for. 1.0.0
recovers `df` where the stored p-value is interior and falls back on the
structural rule (below) where it is not, and computes the corrected p-value as
`-expm1(k * pchisq(stat, df, lower.tail = TRUE, log.p = TRUE))`, which stays
accurate to about 1e-300 rather than collapsing to 0.

## Degrees of freedom are derived per node and per candidate variable

With `teststat = "quadratic"` the node-level statistic is chi-square with `df`
equal to the rank of the covariance matrix of the linear statistic:

* numeric or ordered predictor: `df` = the outcome dimension, so 1 for a
  univariate outcome, 2 for a bivariate outcome, 3 for three outcomes;
* unordered factor with `L` levels **present in that node**: `(L - 1)` times
  the outcome dimension, which changes from node to node as levels drop out of
  a branch.

0.3.0 recovered a single continuous `df` for the variable that was split on.
1.0.0 recovers an integer `df` for every candidate variable at every node and
cross-checks it against the structural rule. `df = 2` (or any number) holds it
fixed, and `df = "outcome"` uses the outcome dimension throughout.

## All candidate variables are tested, not only the one split on

`ctree` splits a node when *any* candidate variable meets `alpha`, and that is
the rule the correction must reapply. Rescaling can reorder the candidates,
because they carry different degrees of freedom, so the variable with the
smallest p-value before the correction is not always the smallest after.
`node_stats` gains `p_min` and `variable_min` recording the decisive candidate,
and the full per-candidate table is attached as
`attr(node_stats, "candidates")`.

## Pruning is bottom-up

An internal node is now collapsed only once all of its own internal
descendants have been collapsed. A node whose descendant survives the
correction keeps its split even if it would fail on its own. 0.3.0 deleted
every failing node outright, taking any surviving descendants with it, which
is more aggressive than the procedure used for the paper. With `df` held
fixed, `prune_stackM()` now reproduces that procedure exactly; this is
asserted in the test suite for M = 2, 5, 10 and 30.

## `prune_stackM()` handles a default `ctree()` fit

`ctree_stacked()` grows with `testtype = "Univariate"` so that the stored
p-values are raw, but `prune_stackM()` is documented as accepting any tree
fitted on stacked data, and `partykit`'s default is `testtype = "Bonferroni"`,
where the stored p-values already carry the multiplicity adjustment. Inverting
those as if they were raw returned a plausible but wrong `df` (2.32 instead of
2; 12.05 instead of 8 in testing) and the internal verification could not
detect it, because it re-solved the same equation it had just inverted.
1.0.0 reads the multiplicity setting from the fitted control and inverts the
matching relationship, and `check_stackM_extraction()` now exercises both
settings and checks the recovered `df` against the structural rule rather than
against itself.

## Other changes

* `minsplit` and `minbucket` are multiplied by M (`scale_minsize = TRUE`), so
  they refer to original rather than stacked observations. The `partykit`
  defaults previously allowed terminal nodes holding fewer than one original
  observation.
* The imputation index is no longer added to the model frame before fitting,
  so a `y ~ .` formula can no longer split on imputation number. The
  documented usage in 0.3.0 invited exactly that.
* `NAMESPACE` was out of date: `prune_stackM()`, `rescale_statistic()` and
  `check_stackM_extraction()` had `@export` tags but were not exported, so
  they were unreachable despite being advertised in `DESCRIPTION`.
* `ctree_stacked()` gains `df` and `scale_minsize` arguments. `formula`,
  `data`, `m`, `alpha` and `verbose` are unchanged.
* Trees fitted with `teststat = "maximum"`, `testtype = "MonteCarlo"` or
  `splittest = TRUE` are rejected with an explanatory error rather than
  silently rescaled: those settings do not use a chi-square reference
  distribution.
* Fixed a test that called `summary.ctreeMI()` by name; the method is
  registered for S3 dispatch but not exported, so the call failed and
  `R CMD check` did not pass on 0.3.0.
* Corrected the Rodgers et al. (2021) reference to the paper actually cited in
  Sherlock et al. (2026): Rodgers, Jacobucci & Grimm, *Journal of Behavioral
  Data Science*, 1(1), 127-153.
* `README.md` still described the `alpha / M` rule from 0.1.0-0.2.0 and
  pointed at a repository URL that does not exist; both corrected.
* Removed `R/zzz.R`. The `.onLoad()` / `registerS3method()` workaround is
  unnecessary now that `NAMESPACE` carries the `S3method()` entries, and it
  registered only two of the four methods.

# ctreeMI 0.3.0

## Correctness fix -- please refit any tree produced by 0.1.0 or 0.2.0

Versions 0.1.0 and 0.2.0 did not implement the Stack/M correction of
Sherlock et al. (2026). They divided the significance *threshold* by M and
passed `alpha / M` to `ctree_control()`. The published and validated method
divides the node-level chi-square *statistic* by M and recomputes the
p-value.

These are not equivalent. Writing `q(p, df)` for the chi-square quantile
function, the statistic rule rejects when `X > M * q(1 - alpha, df)` while
the threshold rule rejects when `X > q(1 - alpha / M, df)`. They coincide
only at `M = 1`. At `df = 1` and `alpha = 0.05`:

| M  | statistic rule | threshold rule |
|----|----------------|----------------|
| 5  | X > 19.2       | X > 6.63       |
| 10 | X > 38.4       | X > 7.88       |
| 30 | X > 115.2      | X > 9.88       |
| 50 | X > 192.1      | X > 10.83      |

Threshold rescaling therefore under-corrects severely, and trees fitted
with 0.1.0-0.2.0 contain more splits than the published method supports.
**Refit any such tree with 0.3.0 and re-check any reported subgroups.**

## Changes

* `ctree_stacked()` now grows the tree at the nominal `alpha` and then
  applies the correction post hoc: node statistics are divided by M,
  p-values recomputed, the multiplicity adjustment reapplied, and failing
  nodes collapsed. Split selection in partykit does not depend on `alpha`,
  and the correction is strictly stricter than the nominal threshold, so
  this yields the same tree as growing under the corrected criterion.
* New `rescale_statistic()`: divides a statistic by M and recomputes its
  p-value.
* New `prune_stackM()`: applies the correction to any `party` object
  fitted on stacked data; returns the pruned tree and a per-node table of
  raw and rescaled statistics.
* New `check_stackM_extraction()`: diagnostic confirming node statistics
  and degrees of freedom can be recovered from the installed `partykit`.
* `rescale_alpha()` is deprecated and now warns. It never implemented the
  published correction.
* `ctreeMI_info` replaces `alpha_nominal`/`alpha_applied` with `alpha`,
  `correction`, `node_stats`, `n_splits_before` and `n_splits_after`.
* `report_ctreeMI()` now describes the statistic-based correction and
  reports how many candidate splits survived it.

# ctreeMI 0.2.0

* New `node_table()`: returns a data frame of tree nodes with the split
  path defining each node, its size in the stacked data, and its
  effective sample size on the original scale (stacked size / M).
* New `report_ctreeMI()`: generates a methods paragraph describing the
  fitted model, populated with the actual M, sample sizes, significance
  thresholds, and tree size.

# ctreeMI 0.1.0

* Initial CRAN release.
* `ctree_stacked()`: fit a conditional inference tree on stacked multiply
  imputed data with the Stack/M significance-threshold correction
  (Sherlock et al., 2026).
* `stack_imputations()`: stack a list of imputed data frames vertically.
* `rescale_alpha()`: compute the Stack/M corrected significance threshold.
* `print.ctreeMI()` and `summary.ctreeMI()` S3 methods.
* Accepts `mids` objects from `mice`, lists of data frames, or plain data
  frames (falls back to standard `ctree` with a warning).
