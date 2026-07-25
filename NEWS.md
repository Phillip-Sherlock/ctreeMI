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
