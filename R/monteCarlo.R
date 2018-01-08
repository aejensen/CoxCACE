rm(list=ls())

library(devtools)
install_github('aejensen/CoxCACE')

library(parallel)

set.seed(123456789)

monteCarloSim <- function(n, psi, B = 1000) {
  unlist(mclapply(1:B, function(sim) {
    dat <- simulateData(n = n, psi = psi)
    coxCACE(dat, lower=-10, upper=10)
  }, mc.cores=8))
}

sim1.150 <- monteCarloSim(n = 150, psi = 0, B = 1000)
sim2.150 <- monteCarloSim(n = 150, psi = 0.405, B = 1000)
sim3.150 <- monteCarloSim(n = 150, psi = -0.405, B = 1000)

sim1.300 <- monteCarloSim(n = 300, psi = 0, B = 1000)
sim2.300 <- monteCarloSim(n = 300, psi = 0.405, B = 1000)
sim3.300 <- monteCarloSim(n = 300, psi = -0.405, B = 1000)

sim1.600 <- monteCarloSim(n = 600, psi = 0, B = 1000)
sim2.600 <- monteCarloSim(n = 600, psi = 0.405, B = 1000)
sim3.600 <- monteCarloSim(n = 600, psi = -0.405, B = 1000)

sim1.1200 <- monteCarloSim(n = 1200, psi = 0, B = 1000)
sim2.1200 <- monteCarloSim(n = 1200, psi = 0.405, B = 1000)
sim3.1200 <- monteCarloSim(n = 1200, psi = -0.405, B = 1000)


boxplot(cbind(sim1.150, sim1.300, sim1.600, sim1.1200))


boxplot(sim1.150)
boxplot(sim2.150)
boxplot(sim3.150)
boxplot(sim4.150)
boxplot(sim5.150)

