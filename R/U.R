U = function(time, status, R, C, psi) {
  #Estimating function for psi
  est <- Lambda(time, status, R, C, psi)
  sum(est$dU)
}
