LambdaCR <- function(time, status, R, C, psi, eps=0.001, verbose=TRUE) {
  n <- length(R)
  stime <- sort(time[status != 0])
  k <- length(stime)

  psi1 <- psi[1]
  psi2 <- psi[2]

  #Allocate output objects
  dU <- matrix(0, 2, k)       #increments for U(psi)
  Lambda <- matrix(0, 4, k)   #Cumulative hazards:
                              #Lambda[1,] = Lambda^N_1(t)
                              #Lambda[2,] = Lambda^C_1(t)
                              #Lambda[3,] = Lambda^N_2(t)
                              #Lambda[4,] = Lambda^C_2(t)
  dLambda <- matrix(0, 4, k)  #Cumulative hazards increments

  #Estimate compliance probability
  n11 <- sum(R == 1 & C == 1)
  n10 <- sum(R == 1 & C == 0)
  pc <- n11 / (n11 + n10)

  ##########################################
  #### Estimation at the first jump time ###
  ##########################################
  dN1 <- as.numeric(time == stime[1] & status == 1)  #jumps for cause 1
  dN2 <- as.numeric(time == stime[1] & status == 2)  #jumps for cause 2
  dN <- matrix(c(dN1, dN2), ncol = 1)                #2n x 1 vector of jumps

  Y <- as.numeric(time >= stime[1])                  #at-risk indicator

  #Calculate Hn and Hc functions
  Hc <- pc
  Hn <- 1 - pc

  #Setup design matrix
  X1 <- cbind(Y * (R * (1 - C) + (1 - R) * Hn),
              Y * (R * C * exp(psi1) + (1 - R) * Hc))

  X2 <- cbind(Y * (R * (1 - C) + (1 - R) * Hn),
              Y * (R * C * exp(psi2) + (1 - R) * Hc))

  X <- as.matrix(Matrix::bdiag(X1, X2))

  #Setup weight matrix
  W1 <- diag(exp(-psi1 * R * C), n)
  W2 <- diag(exp(-psi2 * R * C), n)
  W <- as.matrix(Matrix::bdiag(W1, W2))

  #Calculate increments and update Lambda
  crossMat <- t(X) %*% W %*% X
  #determinant <- det(crossMat)
  #if (determinant > eps) {
    dLambda[, 1] <- solve(crossMat) %*% (t(X) %*% W %*% dN)
  #} else {
  #  if(verbose) {
  #    message("Ill-conditioned design matrix at time = ", stime[1], ", determinant = ", determinant)
  #  }
  #}
  Lambda[, 1] <- dLambda[, 1] + 0

  #Calculate increment for U(psi)
  V <- matrix(NA, 2*n , 2)
  V[,1] <- c(Y * R * C, rep(0, n))
  V[,2] <- c(rep(0, n), Y * R * C)

  dU[, 1] <- t(V) %*% (dN - X %*% dLambda[, 1])

  ###################################
  #### Loop through all the jumps ###
  ###################################
  for (j in 2:k) {
    dN1 <- as.numeric(time == stime[j] & status == 1) #jumps for cause 1
    dN2 <- as.numeric(time == stime[j] & status == 2) #jumps for cause 2
    dN <- matrix(c(dN1, dN2), ncol = 1)               #2n x 1 vector of jumps

    Y <- as.numeric(time >= stime[j])                  #at-risk indicator

    #Calculate Hn and Hc functions
    Hdenom <- pc * exp(-Lambda[2, j - 1] - Lambda[4, j - 1]) + (1 - pc) * exp(-Lambda[1, j - 1] - Lambda[3, j - 1])
    Hn <- ((1 - pc) * exp(-Lambda[1, j - 1] - Lambda[3, j - 1])) / Hdenom
    Hc <- (pc * exp(-Lambda[2, j - 1] - Lambda[4, j - 1])) / Hdenom

    #Setup design matrix
    X1 <- cbind(Y * (R * (1 - C) + (1 - R) * Hn),
                Y * (R * C * exp(psi1) + (1 - R) * Hc))

    X2 <- cbind(Y * (R * (1 - C) + (1 - R) * Hn),
                Y * (R * C * exp(psi2) + (1 - R) * Hc))

    X <- as.matrix(Matrix::bdiag(X1, X2))

    #Setup weight matrix
    W1 <- diag(exp(-psi1 * R * C), n)
    W2 <- diag(exp(-psi2 * R * C), n)
    W <- as.matrix(Matrix::bdiag(W1, W2))

    #Calculate increments and update Lambda
    crossMat <- t(X) %*% W %*% X
    increment <- tryCatch(solve(crossMat) %*% (t(X) %*% W %*% dN), error = function(e) e)
    if(any(class(out) == "error")) {
      dLambda[, j] <- rep(0, 4)
    } else {
      dLambda[, j] <- increment
    }
    #determinant <- det(crossMat)
    #if (determinant > eps) {
    #  dLambda[, j] <- solve(crossMat) %*% (t(X) %*% W %*% dN)
    #} else {
    #  if(verbose) {
    #    message("Ill-conditioned design matrix at time = ", stime[1], ", determinant = ", determinant)
      #}
    #}
    Lambda[, j] <- dLambda[, j] + Lambda[, j - 1]

    #Calculate increment for U(psi)
    V <- matrix(NA, 2*n , 2)
    V[,1] <- c(Y * R * C, rep(0, n))
    V[,2] <- c(rep(0, n), Y * R * C)

    dU[, j] <- t(V) %*% (dN - X %*% dLambda[, j])
  }

  list(time = stime, Lambda = Lambda, dLambda = dLambda, dU = dU)
}

