# ctreeMI <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/phillipsherlock/ctreeMI/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/phillipsherlock/ctreeMI/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/ctreeMI)](https://CRAN.R-project.org/package=ctreeMI)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## Overview

`ctreeMI` implements the **stacked-imputation / Stack ÷ M** workflow for
conditional inference trees (ctree), as described in **Sherlock et al.
(2026)**. It provides a principled, easy-to-use way to apply `ctree` when
your data contain missing values requiring multiple imputation.

## The method in one paragraph

When data are multiply imputed, trees cannot be pooled across imputations
because structurally different trees define different subgroups —
incompatible targets of inference. The solution is to **stack** the M
imputed datasets vertically (M × n rows) and fit **one tree** to the
combined data. The catch: stacking inflates test statistics by M.
The fix: apply a significance threshold of `alpha / M`
(**Stack ÷ M correction**), which is equivalent to dividing each
node-level statistic by M before pruning. Validated via Monte Carlo
simulation in Sherlock et al. (2026), this yields a conservative,
interpretable single tree.

## Installation

```r
# From CRAN (once released):
install.packages("ctreeMI")

# Development version from GitHub:
# install.packages("remotes")
remotes::install_github("phillipsherlock/ctreeMI")
```

## Quick start

```r
library(ctreeMI)
library(mice)

# 1. Impute your data
imp <- mice(my_data, m = 30, printFlag = FALSE)

# 2. Fit ctree with Stack/M correction — one function call
fit <- ctree_stacked(outcome ~ pred1 + pred2 + pred3,
                     data  = imp,
                     alpha = 0.05)

# 3. Everything you'd do with a regular ctree works
print(fit)
plot(fit)
predict(fit, newdata = new_data)
```

## Key functions

| Function            | Description                                               |
|---------------------|-----------------------------------------------------------|
| `ctree_stacked()`   | Main function: fit ctree on stacked MI data               |
| `stack_imputations()` | Stack a list of imputed data frames                     |
| `rescale_alpha()`   | Compute the Stack ÷ M corrected significance threshold    |

## Interpreting output

- The returned object is a standard `partykit` `constparty`/`party` object
  with an added `"ctreeMI"` class — all `partykit` methods work.
- Node sample sizes in `print()` and `plot()` reflect the **stacked**
  dataset. **Divide by M** to get the effective per-node n in your
  original data.
- Node-level estimates (means, proportions) are averaged over imputations
  and are the correct summaries to report.

## Citation

If you use `ctreeMI`, please cite both:

**Methods paper (Stack ÷ M approach):**

> Sherlock, P., Mansolf, M., Hofheimer, J., Hockett, C. W., O'Connor, T. G.,
> Roubinov, D., Graff, J. C., Lai, J.-S., Bush, N. R., Wright, R. J., &
> Chiu, Y.-H. M. (2026). Beyond linear risk: A machine learning approach to
> understanding perinatal depression in context.
> *Multivariate Behavioral Research*, 1–16.
> https://doi.org/10.1080/00273171.2026.2661244

**ctree algorithm:**

> Hothorn, T., Hornik, K., & Zeileis, A. (2006). Unbiased recursive
> partitioning: A conditional inference framework. *Journal of Computational
> and Graphical Statistics*, 15(3), 651–674.

> Hothorn, T., & Zeileis, A. (2015). partykit: A modular toolkit for
> recursive partitioning in R. *Journal of Machine Learning Research*,
> 16, 3905–3909.

In R:

```r
citation("ctreeMI")
```

## License

GPL (>= 3)
