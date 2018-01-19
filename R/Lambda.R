Lambda = function(time, status, R, C, psi, eps=0.001, verbose=TRUE) {
  n <- length(R)
  stime <- sort(time[status != 0])
  k <- length(stime)

  #Allocate output objects
  dU <- rep(0, k)             #increments for U(psi)
  Lambda <- matrix(0, 2, k)   #Cumulative hazards:
                              #Lambda[1,] = Lambda^N(t)
                              #Lambda[2,] = Lambda^C(t)
  dLambda <- matrix(0, 2, k)  #Cumulative hazards increments

  #Estimate compliance probability
  n11 <- sum(R == 1 & C == 1)
  n10 <- sum(R == 1 & C == 0)
  pc <- n11 / (n11 + n10)

  ##########################################
  #### Estimation at the first jump time ###
  ##########################################
  dN <- as.numeric(time == stime[1])  #jumps
  Y <- as.numeric(time >= stime[1])   #at-risk indicator

  #Calculate Hn and Hc functions
  Hc <- pc
  Hn <- 1 - pc

  #Setup design matrix
  X <- cbind(Y * (R * (1 - C) + (1 - R) * Hn),
             Y * (R * C * exp(psi) + (1 - R) * Hc))

  #Setup weight matrix
  W <- diag(exp(-psi * R * C), n)

  #Calculate increments and update Lambda
  crossMat <- t(X) %*% W %*% X
  if (det(crossMat) > eps) {
    dLambda[, 1] <- solve(crossMat) %*% (t(X) %*% W %*% matrix(dN, ncol = 1))
  } else {
    if(verbose) {
      message("Ill-conditioned design matrix at time = ", stime[1])
    }
  }
  Lambda[, 1] <- dLambda[, 1] + 0

  #Calculate increment for U(psi)
  v <- matrix(Y * R * C, ncol = 1)
  dU[1] <- t(v) %*% (matrix(dN, ncol=1) - X %*% dLambda[, 1])

  ###################################
  #### Loop through all the jumps ###
  ###################################
  for (j in 2:k) {
    dN <- as.numeric(time == stime[j]) #jumps
    Y <- as.numeric(time >= stime[j])  #at-risk indicator

    #Calculate Hn and Hc functions
    Hdenom <- pc * exp(-Lambda[2, j - 1]) + (1 - pc) * exp(-Lambda[1, j - 1])
    Hn <- ((1 - pc) * exp(-Lambda[1, j - 1])) / Hdenom
    Hc <- (pc * exp(-Lambda[2, j - 1])) / Hdenom

    #Setup design matrix
    X <- cbind(Y * (R * (1 - C) + (1 - R) * Hn),
               Y * (R * C * exp(psi) + (1 - R) * Hc))

    #Setup weight matrix
    W <- diag(exp(-psi * R * C), n)

    #Calculate increments and update Lambda
    crossMat <- t(X) %*% W %*% X
    if (det(crossMat) > eps) {
      dLambda[, j] <- solve(crossMat) %*% (t(X) %*% W %*% matrix(dN, ncol = 1))
    } else {
      if(verbose) {
        message("Ill-conditioned design matrix at time = ", stime[j])
      }
    }
    Lambda[, j] <- dLambda[, j] + Lambda[, j - 1]

    #Calculate increment for U(psi)
    v <- matrix(Y * R * C, ncol = 1)
    dU[j] <- t(v) %*% (matrix(dN, ncol=1) - X %*% dLambda[, j])
  }

  list(time = stime, Lambda = Lambda, dLambda = dLambda, dU = dU)
}

