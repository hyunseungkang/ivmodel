context("Check accuracy of AR using Card data")


test_that("AR test without exogeneous covariates and intercept, single-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,"nearc4"]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,intercept=FALSE)
  teststat_out = AR.test(foo)
  expect_equal(as.numeric(teststat_out$ci),c(0.462792094733,0.470402558386))
  expect_equal(teststat_out$p.value,0)
  expect_equal(teststat_out$df,c(1,3009))
  expect_equal(teststat_out$Fstat,6679.92118721)
})

test_that("AR test without exogeneous covariates and intercept, mult-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,c("nearc4","nearc2")]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,intercept=FALSE)
  teststat_out = AR.test(foo)
  expect_equal(as.numeric(teststat_out$ci),c(0.465107588912,0.470042589959))
  expect_equal(teststat_out$p.value,0)
  expect_equal(teststat_out$df,c(2,3008))
  expect_equal(teststat_out$Fstat,4165.66455412)
})

test_that("AR test with intercept, but without exogeneous covariates, single-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,"nearc4"]
  
  foo = ivmodel(Y=Y,D=D,Z=Z)
  teststat_out = AR.test(foo)
  expect_equal(as.numeric(teststat_out$ci),c(0.143037504074,0.250862634531))
  expect_equal(teststat_out$p.value,0)
  expect_equal(teststat_out$df,c(1,3008))
  expect_equal(teststat_out$Fstat,82.7445324192)
})


test_that("AR test with intercept, but without exogeneous covariates, mult-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,c("nearc4","nearc2")]
  
  foo = ivmodel(Y=Y,D=D,Z=Z)
  teststat_out = AR.test(foo)
  expect_equal(as.numeric(teststat_out$ci),c(0.166084633236,0.260253933988))
  expect_equal(teststat_out$p.value,0)
  expect_equal(teststat_out$df,c(2,3007))
  expect_equal(teststat_out$Fstat,51.3582485426)
})


test_that("AR test with exogeneous covariates, but no intercept, single-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,"nearc4"]
  Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
          "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
          "reg668", "smsa66")
  X=card.data[,Xname]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,X=X,intercept=FALSE)
  teststat_out = AR.test(foo)
  expect_equal(as.numeric(teststat_out$ci),c(0.280764970615,0.346153931966))
  expect_equal(teststat_out$p.value,0)
  expect_equal(teststat_out$df,c(1,2995))
  expect_equal(teststat_out$Fstat,113.307424185)
})

test_that("AR test with exogeneous covariates, but no intercept, multi-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,c("nearc4","nearc2")]
  Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
          "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
          "reg668", "smsa66")
  X=card.data[,Xname]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,X=X,intercept=FALSE)
  teststat_out = AR.test(foo)
  expect_equal(as.numeric(teststat_out$ci),c(0.310404897446,0.355337343207))
  expect_equal(teststat_out$p.value,0)
  expect_equal(teststat_out$df,c(2,2994))
  expect_equal(teststat_out$Fstat,134.812593279)
})


test_that("AR test with exogeneous covariates and intercept, single-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,"nearc4"]
  Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
          "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
          "reg668", "smsa66")
  X=card.data[,Xname]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,X=X)
  teststat_out = AR.test(foo)
  expect_equal(as.numeric(teststat_out$ci),c(0.0248048359651, 0.284823593339))
  expect_equal(teststat_out$p.value,0.0200276297596)
  expect_equal(teststat_out$df,c(1,2994))
  expect_equal(teststat_out$Fstat,5.41527923822)
})

test_that("AR test with exogeneous covariates, but no intercept, multi-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,c("nearc4","nearc2")]
  Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
          "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
          "reg668", "smsa66")
  X=card.data[,Xname]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,X=X)
  teststat_out = AR.test(foo)
  expect_equal(as.numeric(teststat_out$ci),c(0.0536002610089, 0.361980791255))
  expect_equal(teststat_out$p.value,.00532805613556)
  expect_equal(teststat_out$df,c(2,2993))
  expect_equal(teststat_out$Fstat,5.24393512598)
})




