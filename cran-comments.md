## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* Ubuntu 24.04 (local), R 4.3.3
* GitHub Actions: macOS release, Windows release, Ubuntu release / devel / oldrel-1
* win-builder (devel and release)

## Reason for this submission

This is a correctness release and it changes the numerical results of
`ctree_stacked()` relative to 0.3.0. Two problems are fixed:

1. `ctree_stacked()` aborted on realistic input. The correction recovered the
   degrees of freedom of each node-level test by inverting the p-value stored
   by partykit, and that p-value underflows to exactly 0 on stacked data
   (3,000 rows x M = 30 already yields a root statistic above 5,000). The
   package therefore stopped with an error on the datasets it is intended for.
   Degrees of freedom are now derived structurally where the stored p-value has
   underflowed, and the corrected p-value is computed on the log scale.

2. The correction tested only the variable each node was split on, and pruned
   every failing node outright. It now tests all candidate variables, as ctree
   itself does, and compresses the tree bottom-up, which is the procedure
   described in the reference paper.

NEWS.md documents these changes and advises refitting trees produced by
earlier versions. The exported API is unchanged apart from two new arguments
to `ctree_stacked()`; `prune_stackM()`, `rescale_statistic()` and
`check_stackM_extraction()` are newly reachable because they were documented
as exported in 0.3.0 but omitted from NAMESPACE.

## Downstream dependencies

None.

## Methodology

This package implements the stacked-imputation Stack/M correction for
conditional inference trees described in:

  Sherlock et al. (2026). Beyond linear risk: A machine learning approach to
  understanding perinatal depression in context. Multivariate Behavioral
  Research. doi:10.1080/00273171.2026.2661244

It wraps partykit::ctree() and is a methodological extension of, not a
replacement for, partykit. Both partykit and the methodological paper are
cited in DESCRIPTION and in the function documentation.
