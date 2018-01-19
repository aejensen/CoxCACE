UU <- function(time, status, R, C, beta) {
  apply(rbind(beta), 1, function(x) Ubeta(time, status, R, C, x))
}

dU1 <- function(time, status, R, C, beta, h = 1e-100) {
  pp <- rbind(beta) %x% rbind(1)
  J <- diag(1)
  Im(unlist(UU(time, status, R, C, pp + h * 1i * J))) / h
}

findRoot = function(d, b, alpha = 3/4, tol = 0.001, iterMax = 20) {
  delta <- 1
  i <- 1
  beta = b
  res1 = Ubeta(d, b)

  while ((delta > tol) & (i < iterMax)) {
    dU = dU1(d, beta)
    betaNew = beta - alpha * solve(dU) %*% res1
    res1 = Ubeta(d, betaNew)
    beta = betaNew
    delta = sqrt(res1^2)
    i <- i + 1
  }

  list(beta = beta, U = res1,
       converge = as.numeric(i < iterMax),
       iter = i - 1)
}

