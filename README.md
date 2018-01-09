To install the development version of CoxCACE run the following commands from within R

```{r}
library(devtools)
install_github('aejensen/CoxCACE')
```

# Example
```{r}
library(CoxCACE)

dat <- simulateData(n = 150, psi = 0)
m <- coxCACE(dat)
summary(m)
```

**C**omplier Average **C**ausal Effect for **C**ox **P**roportional Hazards (CCCP)
