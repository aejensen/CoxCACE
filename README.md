To install the development version of TruncComp run the following commands from within R

```{r}
library(devtools)
install_github('aejensen/CoxCACE')
```

# Example
```{r}
library(CoxCACE)

dat <- simulateData(n = 150, psi = 0)
coxCACE(dat, lower=-3, upper=3)
```
