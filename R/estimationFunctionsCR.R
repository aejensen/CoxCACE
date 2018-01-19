LambdaEstCR = function(time, status, R, C, psi, eps=0.01) {
  J <- 2 #number of competing risks
  n <- length(R)

  stime <- sort(time[status != 0]) #sorted event times
  k <- length(stime)

  #Allocate output objects
  dU <- matrix(0, 2, k)           #increments for U(psi)
  Lambda <- matrix(0, 2 * J, k)   #Cumulative hazards, Lambda[1,] = Lambda_1^N,
                                  #                  , Lambda[2,] = Lambda_1^C,
                                  #                  , Lambda[3,] = Lambda_2^N,
                                  #                  , Lambda[4‚] = Lambda_2^C
  dLambda <- matrix(0, 2 * J, k)  #Cumulative hazards increments

  #Estimate compliance probability
  n11 <- sum(R == 1 & C == 1)
  n10 <- sum(R == 1 & C == 0)
  pc <- n11 / (n11 + n10)

  ##########################################
  #### Estimation at the first jump time ###
  ##########################################
  dN1 <- as.numeric(time == stime[1] & status == 1)   #jumps for cause 1
  dN2 <- as.numeric(time == stime[1] & status == 2)   #jumps for cause 2
  dN <- matrix(c(dN1, dN2), ncol = 1)   #

  risk <- as.numeric(time >= stime[1])  #at-risk indicator

  #Calculate Hn and Hc functions
  Hc <- pc
  Hn <- 1 - pc

  #Setup matrices
  X1 <- cbind(risk * (R * (1 - C) + (1 - R) * Hn),
              risk * (R * C * exp(psi[1]) + (1 - R) * Hc))

  X2 <- cbind(risk * (R * (1 - C) + (1 - R) * Hn),
              risk * (R * C * exp(psi[2]) + (1 - R) * Hc))
  X <- as.matrix(Matrix::bdiag(X1, X2))

  W1 <- diag(exp(-psi[1] * R * C), n)
  W2 <- diag(exp(-psi[2] * R * C), n)
  W <- as.matrix(Matrix::bdiag(W1, W2))

  crossMat <- t(X) %*% W %*% X

  #Calculate Lambda increment and update Lambda
  dLambda[, 1] <- matrix(0, 2, 1)
  if (det(crossMat) > eps) {
    dLambda[, 1] <- solve(crossMat) %*% (t(X) %*% W %*% dN)
  } else {
    message("Singularity error")
  }
  Lambda[, 1] <- dLambda[, 1] + 0

  #Calculate increment for U(psi)
  V <- matrix(NA, 2*n, 2)
  V[,1] <- c(R * C, R * C)
  V[,2] <- c(R * C, R * C)
  dU[,1] <- t(V) %*% (dN - X %*% dLambda[, 1])

  ###################################
  #### Loop through all the jumps ###
  ###################################
  for (j in 2:k) {
    dN1 <- as.numeric(time == stime[j] & status == 1)   #jumps for cause 1
    dN2 <- as.numeric(time == stime[j] & status == 2)   #jumps for cause 2
    dN <- matrix(c(dN1, dN2), ncol = 1)   #

    risk <- as.numeric(time >= stime[j]) #at-risk indicator

    #Calculate Hn and Hc functions
    Hnumerator <- pc * exp(-Lambda[2, j - 1] - Lambda[4, j - 1]) +
                 (1 - pc) * exp(-Lambda[1, j - 1] - Lambda[3, j - 1])
    Hn <- ((1 - pc) * exp(-Lambda[1, j - 1] - Lambda[3, j - 1])) / Hnumerator
    Hc <- (pc * exp(-Lambda[2, j - 1] - Lambda[4, j - 1])) / Hnumerator

    #Setup matrices
    X1 <- cbind(risk * (R * (1 - C) + (1 - R) * Hn),
                risk * (R * C * exp(psi[1]) + (1 - R) * Hc))

    X2 <- cbind(risk * (R * (1 - C) + (1 - R) * Hn),
                risk * (R * C * exp(psi[2]) + (1 - R) * Hc))

    X <- as.matrix(Matrix::bdiag(X1, X2))

    W1 <- diag(exp(-psi[1] * R * C), n)
    W2 <- diag(exp(-psi[2] * R * C), n)
    W <- as.matrix(Matrix::bdiag(W1, W2))

    crossMat <- t(X) %*% W %*% X

    #Calculate Lambda increment and update Lambda
    if (det(crossMat) > eps) {
      dLambda[, j] <- solve(crossMat) %*% (t(X) %*% W %*% matrix(dN, ncol = 1))
    } else {
      message("Singularity error")
    }
    Lambda[, j] <- dLambda[, j] + Lambda[, j - 1]

    #Calculate increment for U(beta)
    V.t <- matrix(NA, 2, 2*n)
    V.t[1,] <- c(R * C, R * C)
    V.t[2,] <- c(R * C, R * C)
    #V <- matrix(NA, 2*n, 2)
    #V[,1] <- c(risk * R * C, risk * R * C)
    #V[,2] <- c(risk * R * C, risk * R * C)
    dU[,j] <- V.t %*% (dN - X %*% dLambda[, j])
  }

  list(time = stime, Lam = Lambda, dLam = dLambda, dU = dU)
}

UbetaCR = function(time, status, R, C, psi) {
  #Estimating function for psi
  est <- LambdaEstCR(time, status, R, C, psi)
  apply(est$dU, 1, sum)
}

