# Complier Average Causal Effect for Cox Proportional Hazards

To install the development version of CoxCACE run the following commands from within R

```{r}
library(devtools)
install_github('aejensen/CoxCACE')
```

# Example
```{r}
library(CoxCACE)

set.seed(123456789)
dat <- simulateData(n = 150, psi = 0)
m <- CCCP(dat)
summary(m)
plot(m, type="cumhaz")
```
