simulateComplianceData <- function(n, psi, lambdaC = 1/80, lambdaN = 1/40) {
  p.tr <- 0.5  # Probability of assigned to treatment
  pc.tr <- 0.6 # Probability of compliance

  lambda00 <- lambdaN #non-complier baseline hazard
  lambda01 <- lambdaC #complier baseline hazard

  #Simulate treatment assignements
  R <- rbinom(n, 1, p.tr)

  #Simulate compliance
  E1 <- rbinom(n, 1, pc.tr)

  n00 <- sum((R == 0) & (E1 == 0))
  n01 <- sum((R == 0) & (E1 == 1))
  n10 <- sum((R == 1) & (E1 == 0))
  n11 <- sum((R == 1) & (E1 == 1))

  #Simulate event times
  T00.tmp0 <- rexp(n00) / lambda00
  T01.tmp0 <- rexp(n01) / lambda01
  T10.tmp0 <- rexp(n10) / lambda00
  T11.tmp0 <- rexp(n11) / (lambda01 * exp(psi))

  #Simulate censoring times
  cen1 <- 60
  cen2 <- 100
  entryt2 <- 75

  entry.time.00 <- runif(n00, 0, entryt2)
  entry.time.01 <- runif(n01, 0, entryt2)
  entry.time.10 <- runif(n10, 0, entryt2)
  entry.time.11 <- runif(n11, 0, entryt2)
  cen2.00 <- apply(cbind(100 - entry.time.00, cen1), 1, min)
  cen2.01 <- apply(cbind(100 - entry.time.01, cen1), 1, min)
  cen2.10 <- apply(cbind(100 - entry.time.10, cen1), 1, min)
  cen2.11 <- apply(cbind(100 - entry.time.11, cen1), 1, min)

  #Censor the event times
  T00 <- apply(cbind(T00.tmp0, cen2.00), 1, min)
  T01 <- apply(cbind(T01.tmp0, cen2.01), 1, min)
  T10 <- apply(cbind(T10.tmp0, cen2.10), 1, min)
  T11 <- apply(cbind(T11.tmp0, cen2.11), 1, min)

  status00 <- 0 * T00
  status00[T00 == T00.tmp0] <- 1

  status01 <- 0 * T01
  status01[T01 == T01.tmp0] <- 1

  status10 <- 0 * T10
  status10[T10 == T10.tmp0] <- 1

  status11 <- 0 * T11
  status11[T11 == T11.tmp0] <- 1

  d10 <- cbind(T10, status10, E1[(R == 1) & (E1 == 0)], R[(R == 1) & (E1 == 0)])
  d11 <- cbind(T11, status11, E1[(R == 1) & (E1 == 1)], R[(R == 1) & (E1 == 1)])
  d0 <- cbind(c(T00, T01), c(status00, status01), rep(0, n00 + n01), rep(0, n00 + n01))

  d <- data.frame(rbind(d10, d11, d0))
  names(d) <- c("time", "status", "C", "R")

  d
}
