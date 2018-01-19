# Complier Average Causal Effect for Cox Proportional Hazards

To install the development version of CoxCACE run the following commands from within R

```{r}
library(devtools)
install_github('aejensen/CoxCACE')
```

# Example
```{r}
library(CoxCACE)

set.seed(12345)
d <- simulateComplianceData(400, 2)
m <- CCCP(d)
summary(m)

plot(m, type="cumhaz")
abline(0, 1/80) #true complier cumulative baseline hazard
abline(0, 1/40, lty=2) #true non-complier cumulative baseline hazard

plot(m, type="survival")
curve(exp(-x / 80), 0, 60, add=TRUE)
curve(exp(-x / 40), 0, 60, add=TRUE)
```
