# goodness-of-fit test of the uniform distribution via steck's determinant
## Kolmogorov–Smirnov statistic
```
alpha0 <- 1
x <- runif(10)
stat <- KSStat(x = x, alpha0 = alpha0)
stat
p <- KSPvalue(stat)
p

## Compared with the existing R function
## Same results
ks.test(x,punif)
```

The parameter `alpha0` is a tuning parameter controlling the detection range, that is, the range of x that the differences between F(x) and F_hat(x) will be considered.

## The higher criticism statistic
```
alpha0 <- 1
x <- runif(10)
stat <- HCStat(x = x, alpha0 = alpha0)
stat
p <- HCPvalue(stat)
p
```
## Berk-Jone statistic
```
alpha0 <- 1
x <- runif(10)
stat <- BJStat(x = x, alpha0 = alpha0)
stat
p <- BJPvalue(stat)
p
```
The package also provide functions to compute the one-sided Berk-Jone statistics, namely BJPlusStat and BJMinusStat. 


# Advanced usage
You can also precisely specify the detection range via the argument `index`. It determines which ordered statistic will be used in the test. For example
```
index <- c(2,3,4)
x <- runif(10)
stat <- HCStat(x = x, index = index)
stat
p <- HCPvalue(stat)
p
```
Will compute the higher criticism statistic using only 2th, 3th and 4th ordered samples.

