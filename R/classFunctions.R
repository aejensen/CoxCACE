is.CoxCACE <- function(x) {
  inherits(x, "CoxCACE")
}

coef.CoxCACE <- function(object, ...) {
  object$CACE
}

summary.CoxCACE <- function(object, ...) {
  cat("Complier Average Causal Effect for Cox Proportional Hazards\n\n")

  cMat <- matrix(NA, 1, 4)
  cMat[1,1] <- object$CACE
  cMat[1,2] <- NA
  cMat[1,3] <- NA
  cMat[1,4] <- NA
  colnames(cMat) <- c("Estimate", "Std. Error", "z value", "P(>|z|)")
  rownames(cMat) <- c("CACE")

  cat("Coefficients:\n")
  print.default(cMat)
  cat("\n\n")

  cat("Value of estimating function at estimate:", object$U, "\n")
  cat("Number of iterations:", object$iter, "\n")

  invisible(object)
}

print.CoxCACE <- function(object, ...) {
  summary(object)
}

cumhaz <- function(m) {
  out <- cbind(m$Lambda[, "time"],
               m$Lambda[, "LambdaN"],
               m$Lambda[, "LambdaC"])
  colnames(out) <- c("time", "LambdaN", "LambdaC")
  out
}
