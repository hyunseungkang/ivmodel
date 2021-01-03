context("Check accuracy of CLR using Card data")


test_that("CLR test without exogeneous covariates and intercept, single-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,"nearc4"]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,intercept=FALSE)
  teststat_out = CLR(foo)
  expect_equal(as.numeric(teststat_out$ci),c(0.462792094452, 0.470402558674))
  expect_equal(as.numeric(teststat_out$p.value),0)
  expect_equal(as.numeric(teststat_out$test.stat),6679.92118721)
})

test_that("CLR test without exogeneous covariates and intercept, mult-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,c("nearc4","nearc2")]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,intercept=FALSE)
  teststat_out = CLR(foo)
  expect_equal(as.numeric(teststat_out$ci),c(0.463889965605, 0.471283479854))
  expect_equal(as.numeric(teststat_out$p.value),0)
  expect_equal(as.numeric(teststat_out$test.stat),8327.04383792)
})

test_that("CLR test with intercept, but without exogeneous covariates, single-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,"nearc4"]
  
  foo = ivmodel(Y=Y,D=D,Z=Z)
  teststat_out = CLR(foo)
  expect_equal(as.numeric(teststat_out$ci),c(0.143037537021, 0.250862570435))
  expect_equal(as.numeric(teststat_out$p.value),0)
  expect_equal(as.numeric(teststat_out$test.stat),82.7445324192)
})


test_that("CLR test with intercept, but without exogeneous covariates, mult-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,c("nearc4","nearc2")]
  
  foo = ivmodel(Y=Y,D=D,Z=Z)
  teststat_out = CLR(foo)
  expect_equal(as.numeric(teststat_out$ci),c(0.158839802213, 0.27417906789))
  expect_equal(as.numeric(teststat_out$p.value),0)
  expect_equal(as.numeric(teststat_out$test.stat),99.3797747842)
})


test_that("CLR test with exogeneous covariates, but no intercept, single-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,"nearc4"]
  Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
          "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
          "reg668", "smsa66")
  X=card.data[,Xname]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,X=X,intercept=FALSE)
  teststat_out = CLR(foo)
  expect_equal(as.numeric(teststat_out$ci),c(0.280764970946, 0.346153931564))
  expect_equal(as.numeric(teststat_out$p.value),0)
  expect_equal(as.numeric(teststat_out$test.stat),113.307424185)
})

test_that("CLR test with exogeneous covariates, but no intercept, multi-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,c("nearc4","nearc2")]
  Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
          "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
          "reg668", "smsa66")
  X=card.data[,Xname]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,X=X,intercept=FALSE)
  teststat_out = CLR(foo)
  expect_equal(as.numeric(teststat_out$ci),c(0.309132826837, 0.356869350444))
  expect_equal(as.numeric(teststat_out$p.value),0)
  expect_equal(as.numeric(teststat_out$test.stat),267.050612418)
})


test_that("CLR test with exogeneous covariates and intercept, single-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,"nearc4"]
  Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
          "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
          "reg668", "smsa66")
  X=card.data[,Xname]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,X=X)
  teststat_out = CLR(foo)
  expect_equal(as.numeric(teststat_out$ci),c(0.0248043722948, 0.284824550722))
  expect_equal(as.numeric(teststat_out$p.value),0.0200276297596)
  expect_equal(as.numeric(teststat_out$test.stat),5.41527923822)
})

test_that("CLR test with exogeneous covariates, but no intercept, multi-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,c("nearc4","nearc2")]
  Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
          "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
          "reg668", "smsa66")
  X=card.data[,Xname]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,X=X)
  teststat_out = CLR(foo)
  expect_equal(as.numeric(teststat_out$ci),c(0.0621199910211, 0.336180869927))
  expect_equal(as.numeric(teststat_out$p.value),0.00346295807184)
  expect_equal(as.numeric(teststat_out$test.stat),9.26245429367)
})