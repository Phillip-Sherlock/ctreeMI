## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* macOS Sonoma 14.8.4 (local), R 4.4.2, aarch64-apple-darwin20
* win-builder (R-devel and R-release)

## Reason for this submission

New functionality. The discover-then-confirm workflow introduced in 1.1.0
now supports clustered data: split_holdout() and discover_confirm() can
assign whole clusters to one half of the split, and confirm_ctreeMI() can
use a cluster-robust variance from sandwich::vcovCL() in every pooled test.
The Rubin pooling is now implemented internally so that an arbitrary
covariance matrix can be supplied; a test confirms agreement with mice::D1()
in the unclustered case.

sandwich is added to Suggests. No existing function's default behaviour has
changed.

## Downstream dependencies

None.
