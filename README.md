# ctreeMI

<!-- badges: start -->
[![R-CMD-check](https://github.com/Phillip-Sherlock/ctreeMI/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Phillip-Sherlock/ctreeMI/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/ctreeMI)](https://CRAN.R-project.org/package=ctreeMI)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

Conditional inference trees for multiply imputed data, using the
stacked-imputation **Stack / M** workflow of
[Sherlock et al. (2026)](https://doi.org/10.1080/00273171.2026.2661244).

## The problem

Multiple imputation is the principled way to handle missing data, but its
results cannot be pooled across imputations for trees. Trees are hierarchical
and path dependent: once the root split differs between two imputations, every
node below it describes a different subgroup, and there is nothing comparable
to average.

## The method

Stack the M imputed datasets vertically (M × n rows) and grow **one** tree on
the combined data, so the algorithm sees the variability that imputation
introduced. Stacking multiplies the nominal sample size by M, which inflates
every node-level chi-square statistic, so before deciding which splits to keep:

1. divide each node-level statistic by M;
2. recompute the node-level p-values from the chi-square reference
   distribution `ctree` uses, and reapply the multiplicity adjustment across
   candidate splitting variables;
3. compress the tree bottom-up, collapsing an internal node once all of its
   own internal descendants have been collapsed and no candidate variable
   still meets `alpha`.

This is **not** the same as dividing `alpha` by M. Writing `q(p, df)` for the
chi-square quantile function, statistic rescaling rejects when
`X > M * q(1 - alpha, df)` and threshold rescaling when
`X > q(1 - alpha / M, df)`. At `df = 1`, `alpha = 0.05`, `M = 30` the first
requires `X > 115.2` and the second only `X > 9.9`.

The Monte Carlo study in the supplement to the paper shows the correction is
conservative: node-level type-I error of about 0.034, 0.018 and 0.007 under
10%, 30% and 50% MCAR missingness with M = 30, at some cost in power for weak
effects under heavy missingness.

## Installation

```r
install.packages("ctreeMI")

# development version
# remotes::install_github("Phillip-Sherlock/ctreeMI")
```

## Usage

```r
library(ctreeMI)
library(mice)

imp <- mice(my_data, m = 30, printFlag = FALSE)

# univariate outcome
fit <- ctree_stacked(depression ~ ., data = imp, alpha = 0.05)

# bivariate outcome, as in the paper
fit <- ctree_stacked(PND + PPD ~ cohort + education + marital + bmi,
                     data = imp, alpha = 0.05)

plot(fit)
predict(fit, newdata = new_data)

node_table(fit)      # terminal nodes with effective (unstacked) sample sizes
report_ctreeMI(fit)  # a methods paragraph you can paste into a manuscript
```

Already grew your own tree? Apply the correction post hoc, which is how the
analysis in the paper was done:

```r
grown  <- partykit::ctree(PND + PPD ~ ., data = stacked)
out    <- prune_stackM(grown, m = 30, alpha = 0.05)
out$tree
out$node_stats                          # what the correction did, node by node
attr(out$node_stats, "candidates")      # every candidate variable at every node
```

## Degrees of freedom

With `teststat = "quadratic"` the node-level statistic is chi-square with
degrees of freedom equal to the rank of the covariance matrix of the linear
statistic:

| Predictor | `df` |
| --- | --- |
| numeric or ordered | outcome dimension (1 univariate, 2 bivariate, ...) |
| unordered factor, `L` levels present in the node | `(L - 1) ×` outcome dimension |

These are derived per node and per candidate variable, so univariate,
bivariate and higher-dimensional outcomes all work, and factor predictors are
handled correctly even though their `df` changes as levels drop out of a
branch.

## Key functions

| Function | Description |
| --- | --- |
| `ctree_stacked()` | Stack, grow and correct in one call |
| `prune_stackM()` | Apply the correction post hoc to a tree you grew yourself |
| `rescale_statistic()` | The arithmetic of the correction for a single statistic |
| `stack_imputations()` | Stack a list of imputed data frames |
| `node_table()` | Terminal nodes with effective sample sizes and split paths |
| `report_ctreeMI()` | A methods paragraph describing the analysis |
| `check_stackM_extraction()` | Diagnostic; run after upgrading `partykit` |

## Reading the output

* Terminal-node sample sizes printed by `partykit` refer to the **stacked**
  data. Divide by M for the effective number of original observations, or use
  `node_table()`, which does it for you.
* Node-level means and proportions need no adjustment: they are already
  averages over the imputations.
* `minsplit` and `minbucket` are multiplied by M internally, so you specify
  them in units of original observations.

## Citation

```r
citation("ctreeMI")
```

Sherlock, P., Mansolf, M., Hofheimer, J., Hockett, C. W., O'Connor, T. G.,
Roubinov, D., Graff, J. C., Lai, J.-S., Bush, N. R., Wright, R. J., & Chiu,
Y.-H. M. (2026). Beyond linear risk: A machine learning approach to
understanding perinatal depression in context. *Multivariate Behavioral
Research*, 1–16. <https://doi.org/10.1080/00273171.2026.2661244>

The underlying algorithm is `partykit::ctree()` (Hothorn, Hornik & Zeileis,
2006; Hothorn & Zeileis, 2015); please cite it as well.
