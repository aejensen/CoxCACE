CCCP <- function(data, lower = -10, upper = 10) {
  time = data$time
  status = data$status
  R <- data$R
  C <- data$C

  out <- tryCatch({
      rs <- uniroot(function(b) {
              U(time, status, R, C, b)
            }, interval=c(lower, upper), extendInt = "yes", check.conv=TRUE)
      list(CACE = rs$root,
           U    = rs$f.root,
           iter = rs$iter,
           convergence = TRUE)
    },
    error = function(c) {
      list(CACE = NA,
           U    = NA,
           iter = NA,
           convergence = FALSE)
    }
  )

  if(out$convergence == TRUE) {
    #If we found the optimal beta, get Lambda and dLambda at that value
    lamOpt <- Lambda(time, status, R, C, out$CACE)

    LambdaMat <- cbind(lamOpt$time, t(lamOpt$Lambda))
    colnames(LambdaMat) <- c("time", "LambdaN", "LambdaC")

    dLambdaMat <- cbind(lamOpt$time, t(lamOpt$dLambda))
    colnames(dLambdaMat) <- c("time", "dLambdaN", "dLambdaC")

    out$Lambda <- LambdaMat
    out$dLambda <- dLambdaMat
  } else {
    message("Convergence failed")
    opt$Lambda <- NULL
    out$dLambda <- NULL
  }

  class(out) <- append(class(out), "CoxCACE")
  out
}

