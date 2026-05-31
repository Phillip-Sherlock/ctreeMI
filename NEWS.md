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
