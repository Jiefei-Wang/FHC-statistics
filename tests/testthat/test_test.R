context("Test statistic")
set.seed(123)
N=1000
alpha_true <- 0.5

std_err <- sqrt(alpha_true*(1-alpha_true)/N)
test_that("HC",{
  record=c()
  critical=HCCritical(alpha_true,5,index = 2:3)
  for(i in 1:N){
    p=runif(5)
    stat=HCStat(p,index=2:3)
    record=c(record,stat>critical)
  }
  alpha=mean(record)
  expect_true(abs(alpha-alpha_true)<std_err*qnorm(0.9999))
})

test_that("BJ+",{
  record=c()
  critical=BJPlusCritical(alpha_true,5,index = 2:3)
  for(i in 1:N){
    p=runif(5)
    stat=BJPlusStat(p,index=2:3)
    record=c(record,stat<critical)
  }
  alpha=mean(record)
  expect_true(abs(alpha-alpha_true)<std_err*qnorm(0.9999))
})


test_that("BJ-",{
  record=c()
  critical=BJMinusCritical(alpha_true,5,index = 2:3)
  for(i in 1:N){
    p=runif(5)
    stat=BJMinusStat(p,index=2:3)
    record=c(record,stat<critical)
  }
  alpha=mean(record)
  expect_true(abs(alpha-alpha_true)<std_err*qnorm(0.9999))
})

test_that("BJ",{
  record=c()
  critical=BJCritical(alpha_true,5,index = 2:3)
  for(i in 1:N){
    p=runif(5)
    stat=BJStat(p,index=2:3)
    record=c(record,stat<critical)
  }
  alpha=mean(record)
  expect_true(abs(alpha-alpha_true)<std_err*qnorm(0.9999))
})

test_that("KS",{
  record=c()
  critical=KSCritical(alpha_true,5,index = 2:3)
  for(i in 1:N){
    p=runif(5)
    stat=KSStat(p,index=2:3)
    record=c(record,stat>critical)
  }
  alpha=mean(record)
  expect_true(abs(alpha-alpha_true)<std_err*qnorm(0.9999))
})
