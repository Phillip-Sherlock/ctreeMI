# Force registration of S3 methods at load time so they take
# precedence over methods from imported packages (e.g. partykit).
.onLoad <- function(libname, pkgname) {
  registerS3method("print",   "ctreeMI", print.ctreeMI,   envir = asNamespace(pkgname))
  registerS3method("summary", "ctreeMI", summary.ctreeMI, envir = asNamespace(pkgname))
}
