## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new submission.

## Downstream dependencies

None — this is the first release, so there are no reverse dependencies
to check.

## Notes on test environments

Tested on:
* macOS (local): R 4.x, no errors/warnings/notes
* Ubuntu 22.04 (GitHub Actions): R release, R devel, R oldrel-1
* Windows (GitHub Actions): R release

## Notes on methodology

This package implements the stacked-imputation / Stack/M correction for
conditional inference trees, as described in the peer-reviewed paper:

  Sherlock et al. (2026). Beyond linear risk: A machine learning
  approach to understanding perinatal depression in context.
  Multivariate Behavioral Research. doi:10.1080/00273171.2026.2661244

The package wraps partykit::ctree() and is designed as a methodological
extension, not a replacement, of partykit. Both partykit and the
methodological paper are cited in DESCRIPTION, all function
documentation, and the vignette.
