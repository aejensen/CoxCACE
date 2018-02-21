estimateCR <- function(time, status, R, C, psiStart, maxConditionNumber=300, verbose=TRUE) {
  rootSearch <- rootSolve::multiroot(function(psi) {
    UCR(time, status, R, C, psi, maxConditionNumber, verbose)
  }, start=psiStart, maxiter=250)

  roots <- rootSearch$root
  values <- rootSearch$f.root
  iter <- rootSearch$iter
  
  lambda <- LambdaCR(time, status, R, C, roots, maxConditionNumber, verbose)
  
  
  list(psi = roots, value = values, norm = sqrt(sum(values^2)),
       iter = iter, lambda = lambda)
}

