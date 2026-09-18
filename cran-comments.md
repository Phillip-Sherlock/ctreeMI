## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* macOS Sonoma 14.8.4 (local), R 4.4.2, aarch64-apple-darwin20
* win-builder (R-devel and R-release)

## Reason for this submission

New functionality and a documentation correction.

Three new exported functions (split_holdout, confirm_ctreeMI,
discover_confirm) implement a discover-then-confirm workflow in which a
tree fitted by ctree_stacked() is tested on independently imputed held-out
data using mice::D1(). No existing function has changed.

The package-level help page described the correction as dividing the
significance threshold by M, which is not what the package does and is
contradicted by ?ctree_stacked. It now describes the mechanism correctly.

## Downstream dependencies

None.
