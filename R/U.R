U <- function(time, status, R, C, psi, eps=0.001, verbose=TRUE) {
  #Estimating function for psi
  est <- Lambda(time, status, R, C, psi, eps, verbose)
  sum(est$dU)
}

UCR <- function(time, status, R, C, psi, maxConditionNumber=300, verbose=TRUE) {
  #Estimating function for psi
  est <- LambdaCR(time, status, R, C, psi, maxConditionNumber, verbose)
  apply(est$dU, 1, sum)
}

