print.CoxCACE <- function(x, ...) {
  summary(x)
}


coef.CoxCACE <- function(x, ...) {
  x$CACE
}

summary.CoxCACE <- function(object, ...) {
  cat("Complier Average Causal Effect Cox Proportional Hazards\n\n")

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

plot.CoxCACE <- function(object, type, ...) {
  if(type == "cumhaz") {
    yMax <- max(object$Lambda[, c("Lambda1", "Lambda2")])
    plot(object$Lambda[, "time"], object$Lambda[, "Lambda1"],
         type="s", ylim=c(0, yMax), xlab="t", ylab=expression(Lambda(t)), bty="n")
    lines(object$Lambda[, "time"], object$Lambda[, "Lambda2"], type="s")
  } else if(type == "survival") {
    stop("not implemented yet")
  } else {
    stop("Wrong type as argument. Use either cumhaz or survival")
  }
}
