Lambda = function(time, status, R, C, psi, eps=0.01, verbose=TRUE) {
  n <- length(R)
  stime <- sort(time[status != 0])
  k <- length(stime)

  #Allocate output objects
  dU <- rep(0, k)             #increments for U(psi)
  Lambda <- matrix(0, 2, k)   #Cum hazards, Lambda[1,] = Lambda^N(t), Lambda[2,] = Lambda^C(t)
  dLambda <- matrix(0, 2, k)  #Cum hazards increments

  #Estimate compliance probability
  n11 <- sum(R == 1 & C == 1)
  n10 <- sum(R == 1 & C == 0)
  pc <- n11 / (n11 + n10)

  ##########################################
  #### Estimation at the first jump time ###
  ##########################################
  dN <- as.numeric(time == stime[1])   #jumps
  risk <- as.numeric(time >= stime[1]) #at-risk indicator

  #Calculate Hn and Hc functions
  Hc <- pc
  Hn <- 1 - pc

  #Setup matrices
  X <- cbind(risk * (R * (1 - C) + (1 - R) * Hn),
             risk * (R * C * exp(psi) + (1 - R) * Hc))
  W <- diag(exp(-psi * R * C), n)

  crossMat <- t(X) %*% W %*% X

  #Calculate increments and update Lambda
  dLambda[, 1] <- matrix(0, 2, 1)
  if (det(crossMat) > eps) {
    dLambda[, 1] <- solve(crossMat) %*% (t(X) %*% W %*% matrix(dN, ncol = 1))
  } else {
    if(vebose) {
      message("Ill-conditioned matrix at time = ", stime[1])
    }
  }
  Lambda[, 1] <- dLambda[, 1] + 0

  #Calculate increment for U(psi)
  v <- risk * R * C
  dU[1] <- sum(v * dN) - matrix(v, nrow = 1) %*% X %*% dLambda[, 1]

  ###################################
  #### Loop through all the jumps ###
  ###################################
  for (j in 2:k) {
    dN <- as.numeric(time == stime[j])   #jumps
    risk <- as.numeric(time >= stime[j]) #at-risk indicator

    #Calculate Hn and Hc functions
    Hnumerator <- pc * exp(-Lambda[1, j - 1]) + (1 - pc) * exp(-Lambda[2, j - 1])
    Hn <- ((1 - pc) * exp(-Lambda[1, j - 1])) / Hnumerator
    Hc <- (pc * exp(-Lambda[2, j - 1])) / Hnumerator

    #Setup matrices
    X <- cbind(risk * (R * (1 - C) + (1 - R) * Hn),
               risk * (R * C * exp(psi) + (1 - R) * Hc))
    W <- diag(exp(-psi * R * C), n)
    crossMat <- t(X) %*% W %*% X

    #Calculate increments and update Lambda
    dLambda[, j] <- matrix(0, 2, 1)
    if (det(crossMat) > eps) {
      dLambda[, j] <- solve(crossMat) %*% (t(X) %*% W %*% matrix(dN, ncol = 1))
    } else {
      if(verbose) {
        message("Ill-conditioned matrix at time = ", stime[j])
      }
    }
    Lambda[, j] <- dLambda[, j] + Lambda[, j - 1]

    #Calculate increment for U(beta)
    v <- risk * R * C
    dU[j] <- sum(v * dN) - matrix(v, nrow = 1) %*% X %*% dLambda[, j]
  }

  list(time = stime, Lam = Lambda, dLam = dLambda, dU = dU)
}
