plot.CoxCACE <- function(object, type, ...) {
  if(type == "cumhaz") {
    time <- object$Lambda[, "time"]
    LN <- object$Lambda[, "LambdaN"]
    LC <- object$Lambda[, "LambdaC"]
    
    yMax <- max(c(LN, LC))
    
    plot(time, LN, type="s", ylim=c(0, yMax), xlab="Time", 
         ylab="Cumulative baseline hazards", bty="n", lty=2, lwd=2)
    lines(time, LC, type="s", lty=1, lwd=2)
    legend("topleft", c(expression(Lambda^N), expression(Lambda^C)), lty=1:2, bty="n", lwd=2)
  } else if(type == "survival") {
    time <- object$Lambda[, "time"]
    sN <- exp(-object$Lambda[, "LambdaN"])
    sC <- exp(-object$Lambda[, "LambdaC"])
    
    plot(time, sN, ylim=c(0,1), lty=2, lwd=2, type="s", bty="n")
    lines(time, sC, lty=1, lwd=2, type="s")
    legend("topright", c(expression(S^N), expression(S^C)), lty=1:2, bty="n", lwd=2)
  } else if(type == "CACE") {
    sC0 <- exp(-object$Lambda[, "LambdaC"])
    plot(0,0)
  } else {
    stop("Wrong type as argument. Use either cumhaz, survival or CACE.")
  }
}


  