context("Test statistic")

N=4000
test_that("HC",{
  record=c()
  critical=HCCritical(0.05,5,index = 2:3)
  for(i in 1:N){
    p=runif(5)
    stat=HCStat(p,index=2:3)
    record=c(record,stat>critical)
  }
  alpha=mean(record)
  expect_true(alpha<0.06&&alpha>0.04)
})

test_that("BJ+",{
  record=c()
  critical=BJPlusCritical(0.05,5,index = 2:3)
  for(i in 1:N){
    p=runif(5)
    stat=BJPlusStat(p,index=2:3)
    record=c(record,stat<critical)
  }
  alpha=mean(record)
  expect_true(alpha<0.06&&alpha>0.04)
})


test_that("BJ-",{
  record=c()
  critical=BJMinusCritical(0.05,5,index = 2:3)
  for(i in 1:N){
    p=runif(5)
    stat=BJMinusStat(p,index=2:3)
    record=c(record,stat<critical)
  }
  alpha=mean(record)
  expect_true(alpha<0.06&&alpha>0.04)
})

test_that("BJ",{
  record=c()
  critical=BJCritical(0.05,5,index = 2:3)
  for(i in 1:N){
    p=runif(5)
    stat=BJStat(p,index=2:3)
    record=c(record,stat<critical)
  }
  alpha=mean(record)
  expect_true(alpha<0.06&&alpha>0.04)
})

test_that("KS",{
  record=c()
  critical=KSCritical(0.05,5,index = 2:3)
  for(i in 1:N){
    p=runif(5)
    stat=KSStat(p,index=2:3)
    record=c(record,stat>critical)
  }
  alpha=mean(record)
  expect_true(alpha<0.06&&alpha>0.04)
})
