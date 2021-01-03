context("Check accuracy of TSLS with multiple instruments using Card data")

test_that("TSLS estimate without exogeneous covariates and intercept", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,c("nearc4","nearc2")]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,intercept=FALSE)
  tsls_output = as.numeric(coef(foo)["TSLS",-1])
  expect_equal(tsls_output[1],0.467542363569)
  expect_equal(tsls_output[2],0.00188597999413)
  expect_equal(tsls_output[3],247.904201012)
  expect_equal(tsls_output[4],0)
})

test_that("TSLS estimate with intercept, but without exogeneous covariates", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,c("nearc4","nearc2")]
  
  foo = ivmodel(Y=Y,D=D,Z=Z)
  tsls_output = as.numeric(coef(foo)["TSLS",-1])
  tsls_other = as.numeric(coefOther(foo)$TSLS)
  
  expect_equal(tsls_output[1],0.198413329689)
  expect_equal(tsls_output[2],0.0265814086442)
  expect_equal(tsls_output[3],7.4643647500)
  expect_equal(tsls_output[4],1.08988024378*10^(-13))
  
  expect_equal(tsls_other[1],3.630185655870)
  expect_equal(tsls_other[2],0.3527172577269)
  expect_equal(tsls_other[3],10.2920556801)
  expect_equal(tsls_other[4],1.93651427275*10^(-24))
})

test_that("TSLS estimate with exogeneous covariates, but no intercept", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,c("nearc4","nearc2")]
  Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
          "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
          "reg668", "smsa66")
  X=card.data[,Xname]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,X=X,intercept=FALSE)
  tsls_output = as.numeric(coef(foo)["TSLS",-1])
  tsls_other = coefOther(foo)$TSLS
  
  expect_equal(tsls_output[1],0.33111664127322)
  expect_equal(tsls_output[2],0.01196472190663)
  expect_equal(tsls_output[3],27.674411813098)
  expect_equal(tsls_output[4],2.61916034523*10^(-150))
  
  other_est = c(0.23479475610025,-0.00443128725739,0.03229427192418,-0.12068553496214,0.06199162862063,
                0.05485549147834 , 0.17152179874178,0.21891623853196, 0.10781815206282,0.26886047791751,
                0.31021214301357, 0.24045655835388, -0.11682019529151, 0.01738214356183)
  other_se = c(0.020088221,
               0.001096627,
               0.033011287,
               0.04400757,
               0.041235719,
               0.088908707,
               0.076211643,
               0.076328156,
               0.093886546,
               0.083977586,
               0.09170969,
               0.084744296,
               0.103934491,
               0.035850625)
  other_t = c(11.68818094,
              -4.040833265,
              0.978279703,
              -2.742381265,
              1.503347845,
              0.61698672,
              2.250598362,
              2.868092834,
              1.148387679,
              3.201574278,
              3.382544887,
              2.83743649,
              -1.123979101,
              0.484849102)
  other_pval = c(6.794225484920*10^(-31),
                 5.460107020420*10^(-5),
                 3.280150171910*10^(-1),
                 6.135823223650*10^(-3),
                 1.328548008960*10^(-1),
                 5.372903838150*10^(-1),
                 2.448325000290*10^(-2),
                 4.158423483270*10^(-3),
                 2.509002692110*10^(-1),
                 1.381094889550*10^(-3),
                 7.273967643340*10^(-4),
                 4.578327834730*10^(-3),
                 2.611119708790*10^(-1),
                 6.278188472010*10^(-1))
  expect_equal(as.numeric(tsls_other[,1]),other_est)
  expect_equal(as.numeric(tsls_other[,2]),other_se)
  expect_equal(as.numeric(tsls_other[,3]),other_t)
  expect_equal(as.numeric(tsls_other[,4]),other_pval)
})

test_that("TSLS estimate with exogeneous covariates and intercept", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,c("nearc4","nearc2")]
  Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
          "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
          "reg668", "smsa66")
  X=card.data[,Xname]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,X=X)
  tsls_output = as.numeric(coef(foo)["TSLS",-1])
  tsls_other = coefOther(foo)$TSLS
  
  expect_equal(tsls_output[1],0.157059370024375)
  expect_equal(tsls_output[2],0.05257824168154)
  expect_equal(tsls_output[3],2.98715523762946)
  expect_equal(tsls_output[4],2.83871433855*10^(-3))
  
  other_est = c(0.118814881,
                -0.002356484,
                -0.123277795,
                -0.143194462,
                0.100753,
                -0.102975996,
                -0.000228649,
                0.046955624,
                -0.055408388,
                0.051504145,
                0.069996805,
                0.03905956,
                -0.198037081,
                0.015062582,
                3.339686812)
  other_se = c(0.022806069,
               0.000347517,
               0.052150037,
               0.028444785,
               0.031519343,
               0.043422369,
               0.033794272,
               0.032649036,
               0.039182839,
               0.047567785,
               0.053304927,
               0.049749874,
               0.052534992,
               0.022335974,
               0.894537747)
  other_t = c(5.209792324,
              -6.78090708,
              -2.363906182,
              -5.034120037,
              3.196545079,
              -2.371496505,
              -0.006765914,
              1.438193294,
              -1.414098367,
              1.082752641,
              1.313139484,
              0.785118781,
              -3.769622376,
              0.67436422,
              3.733421896)
  other_pval = c(2.018681860470*10^(-7),
                 1.433054519650*10^(-11),
                 1.814690049260*10^(-2),
                 5.084491004100*10^(-7),
                 1.405319681830*10^(-3),
                 1.777913469090*10^(-2),
                 9.946020734200*10^(-1),
                 1.504837320980*10^(-1),
                 1.574369962530*10^(-1), 
                 2.790054237960*10^(-1),
                 1.892365812910*10^(-1),
                 4.324460493030*10^(-1),
                 1.666512531760*10^(-4),
                 5.001318485850*10^(-1),
                 1.924118039560*10^(-4))
  expect_equal(as.numeric(tsls_other[,1]),other_est)
  expect_equal(as.numeric(tsls_other[,2]),other_se)
  expect_equal(as.numeric(tsls_other[,3]),other_t)
  expect_equal(as.numeric(tsls_other[,4]),other_pval)
})
