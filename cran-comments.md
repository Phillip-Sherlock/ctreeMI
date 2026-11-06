## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Number of updates in past 6 months: 7

Version 1.2.0 was released on 12 October 2026. This submission is not a
bug-fix resubmission: it adds a design for data nested in few clusters,
ordinal and nominal outcomes in the confirmation step, and two vignettes.

## Test environments

* macOS, local, R 4.4.2 (aarch64-apple-darwin20)
* win-builder, R-release 4.6.2 (x86_64-w64-mingw32)
* win-builder, R-devel (2026-11-05 r90840)

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
