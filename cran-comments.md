## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* macOS Sonoma 14.8.4 (local), R 4.4.2, aarch64-apple-darwin20
* win-builder (R-devel and R-release)

## Reason for this submission

New functionality. The discover-then-confirm workflow now has a design for
data nested in few clusters, where splitting whole clusters would leave too
few to confirm on: participants are split within every cluster, the tree is
grown on the outcome centered by cluster (center_outcome(), and center_by in
ctree_stacked()), and the partition is confirmed three ways
(confirm_within()), including with cluster fixed effects (the new `fixed`
argument of confirm_ctreeMI()). discover_confirm() gains `within` to run
the sequence.

confirm_ctreeMI() also gains ordinal and nominal outcomes, through
MASS::polr() and nnet::multinom(); MASS and nnet are added to Suggests and
are only required when such an outcome is confirmed.

All new arguments default to their previous behavior. No existing
function's default behavior has changed.
