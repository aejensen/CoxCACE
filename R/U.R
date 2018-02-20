U <- function(time, status, R, C, psi) {
  #Estimating function for psi
  est <- Lambda(time, status, R, C, psi)
  sum(est$dU)
}

UCR <- function(time, status, R, C, psi) {
  #Estimating function for psi
  est <- LambdaCR(time, status, R, C, psi)
  apply(est$dU, 1, sum)
}

