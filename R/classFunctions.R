print.CoxCACE <- function(x, ...) {
  print(paste("CACE = ", signif(x$CACE,4)))
  invisible(x)
}


coef.CoxCACE <- function(x, ...) {
  x$CACE
}

summary.CoxCACE <- function(object, ...) {
  cat("Complier Average Causal Effect Cox Proportional Hazards\n\n")

  cMat <- matrix(NA, 1, 4)
  cMat[1,1] <- 1object$CACE
  cMat[1,2] <- NA
  cMat[1,3] <- NA
  cMat[1,4] <- NA
  colnames(cMat) <- c("Estimate", "Std. Error", "z value", "P(>|z|)")
  rownames(cMat) <- c("CACE")

  cat("Coefficients:\n")
  print.default(cMat)
  cat("\n\n")

  invisible(object)
}
