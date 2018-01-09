coxCACE <- function(data, lower = -10, upper = 10) {
  time = data$time
  status = data$status
  R <- data$R
  C <- data$C

  estimate <- tryCatch(uniroot(function(b) Ubeta(time, status, R, C, b), interval=c(lower, upper))$root,
    error = function(c) NA,
    warning = function(c) "warning",
    message = function(c) "message"
  )

  out <- list(CACE = estimate,
              convergence = 1)
  class(out) <- append(class(out), "CoxCACE")
  out
}

Ubeta = function(time, status, R, C, beta) {
  #Estimating function for beta
  est <- LambdaEstimator(time, status, R, C, beta)
  sum(est$dU)
}

LambdaEstimator = function(time, status, R, C, beta, eps=0.01) {
  stime = sort(time[status == 1])
  k = length(stime)

  #Allocate output objects
  dU <- rep(0, k)            #increments for U(beta)
  Lambda = matrix(0, 2, k)   #Cumulative hazards
  dLambda = matrix(0, 2, k)  #Cumulative hazards increments

  #Estimate compliance probability
  n11 <- sum(R == 1 & C == 1)
  n10 <- sum(R == 1 & C == 0)
  pc = n11 / (n11 + n10)

  ##########################################
  #### Estimation at the first jump time ###
  ##########################################
  dN = as.numeric(time == stime[1])   #jumps
  risk = as.numeric(time >= stime[1]) #at-risk indicator

  #Calculate Hn and Hc functions
  Hc = pc
  Hn = 1 - pc

  #Setup matrices
  X = cbind(risk * (R * C + (1 - R) * Hc),
            risk * (R * (1 - C) + (1 - R) * Hn))
  Y =  cbind(exp(R * C * beta) * risk * (R * C + (1 - R) * Hc),
             exp(R * C * beta) * risk * (R * (1 - C) + (1 - R) * Hn))
  crossMat = t(X) %*% Y

  #Calculate Lambda increment and update Lambda
  dLambda[, 1] = matrix(0, 2, 1)
  if (determinant2(crossMat) > eps) {
    dLambda[, 1] = solve(crossMat) %*% (t(X) %*% matrix(dN, ncol = 1))
  }
  Lambda[, 1] = dLambda[, 1] + 0

  #Calculate increment for U(beta)
  v = risk * R * C
  dU[1] = sum(v * dN) - matrix(v, nrow = 1) %*% Y %*% dLambda[, 1]

  ###################################
  #### Loop through all the jumps ###
  ###################################
  for (j in 2:k) {
    dN = as.numeric(time == stime[j])   #jumps
    risk = as.numeric(time >= stime[j]) #at-risk indicator

    #Calculate Hn and Hc functions
    Hnumerator = exp(-Lambda[2, j - 1]) * (1 - pc) + exp(-Lambda[1, j - 1]) * pc
    Hc = (exp(-Lambda[1, j - 1]) * pc) / Hnumerator
    Hn = (exp(-Lambda[2, j - 1]) * (1 - pc)) / Hnumerator

    #Setup matrices
    X = cbind(risk * (R * C + (1 - R) * Hc),
              risk * (R * (1 - C) + (1 - R) * Hn))
    Y =  cbind(exp(R * C * beta) * risk * (R * C + (1 - R) * Hc),
               exp(R * C * beta) * risk * (R * (1 - C) + (1 - R) * Hn))
    crossMat = t(X) %*% Y

    #Calculate Lambda increment and update Lambda
    dLambda[, j] = matrix(0, 2, 1)
    if (determinant2(crossMat) > eps) {
      dLambda[, j] = solve(crossMat) %*% (t(X) %*% matrix(dN, ncol = 1))
    }
    Lambda[, j] = dLambda[, j] + Lambda[, j - 1]

    #Calculate increment for U(beta)
    v = risk * R * C
    vY = matrix(c(exp(beta) * sum(risk * R * C), 0), 1, 2)
    dU[j] = (sum(v * dN) - vY %*% dLambda[, j])
  }

  list(Lam = Lambda, dLam = dLambda, dU = dU)
}

determinant2 <- function(m) {
  abs(m[1, 1] * m[2, 2] - m[1, 2]^2)
}

