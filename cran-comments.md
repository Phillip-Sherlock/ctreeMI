## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* Ubuntu 24.04 (local), R 4.3.3
* GitHub Actions: macOS release, Windows release, Ubuntu release / devel / oldrel-1
* win-builder (devel and release)

## Reason for this submission

Documentation only. No code has changed and results are identical to 1.0.0.

The DESCRIPTION and ?ctree_stacked described the type-I error behaviour of the
procedure as sub-nominal. That is accurate for the imputation model used in
the paper the method comes from, which conditioned on neither the outcome nor
the remaining predictors, but not for the outcome-conditioned imputation
recommended in practice and produced by mice() defaults. The documentation now
states the behaviour under both, and the change is recorded in NEWS.md.
