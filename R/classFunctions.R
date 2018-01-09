print.CoxCACE <- function(x, ...) {
  print("CACE = ", x$CACE, "\n")
  invisible(x)
}


coef.CoxCACE <- function(x, ...) {
  x$CACE
}

