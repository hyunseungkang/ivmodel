context("Check accuracy of TSLS with single instrument using Card data")

test_that("TSLS estimate without exogeneous covariates and intercept", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,"nearc4"]

  foo = ivmodel(Y=Y,D=D,Z=Z,intercept=FALSE)
  tsls_output = as.numeric(coef(foo)["TSLS",-1])
  expect_equal(tsls_output[1],0.4665769)
  expect_equal(tsls_output[2],0.001940093)
  expect_equal(tsls_output[3],240.492081152)
  expect_equal(tsls_output[4],0)
})

test_that("TSLS estimate with intercept, but without exogeneous covariates", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,"nearc4"]
  
  foo = ivmodel(Y=Y,D=D,Z=Z)
  tsls_output = as.numeric(coef(foo)["TSLS",-1])
  tsls_other = as.numeric(coefOther(foo)$TSLS)
  
  expect_equal(tsls_output[1],0.188062632758)
  expect_equal(tsls_output[2],0.026291343964)
  expect_equal(tsls_output[3],7.15302469953)
  expect_equal(tsls_output[4],1.061464*10^(-12))
  
  expect_equal(tsls_other[1],3.767471660374)
  expect_equal(tsls_other[2],0.348861744657)
  expect_equal(tsls_other[3],10.79932585924)
  expect_equal(tsls_other[4],1.063760*10^(-26))
})

test_that("TSLS estimate with exogeneous covariates, but no intercept", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,"nearc4"]
  Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
          "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
          "reg668", "smsa66")
  X=card.data[,Xname]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,X=X,intercept=FALSE)
  tsls_output = as.numeric(coef(foo)["TSLS",-1])
  tsls_other = coefOther(foo)$TSLS
  
  expect_equal(tsls_output[1],0.311872823011726)
  expect_equal(tsls_output[2],0.01637624982246)
  expect_equal(tsls_output[3],19.04421503047432)
  expect_equal(tsls_output[4],2.021866*10^(-76))
  
  other_est = c(0.264680781815391,-0.006066303123084,0.010770651971337,-0.111136194883946,0.101606674342650,
                0.158161632945826 , 0.273710902514624,0.323305448915530, 0.231086891116473,0.369085137097003,
                0.417753000654896 , 0.331919063677470, -0.000200505471685, 0.042212714810450)
  other_se = c(0.02648192,0.001447105,0.034450887,0.042937923,0.046465209,
               0.106121595,0.096012842,0.09695341,0.117274666,0.101179265,
               0.109812044,0.098736306,0.122609095,0.03775272)
  other_t = c(9.994773149,-4.192026522,0.312637873,-2.588299244,
              2.186725879,1.490381222,2.850773886,3.33464753,1.970475804,
              3.647833736,3.804254861,3.361671871,-0.001635323,1.118137034)
  other_pval = c(3.68888545695*10^(-23),2.84479425579*10^(-5),7.54577566779*10^(-1),
                 9.69174313767*10^(-3),2.88396560010*10^(-2),1.36229329392*10^(-1),
                 4.39113473619*10^(-3),8.64474417991*10^(-4),4.88758353229*10^(-2),
                 2.68968590078*10^(-4),1.45073675929*10^(-4),7.84441650405*10^(-4),
                 9.98695310526*10^(-1),2.63598148986*10^(-1))
  expect_equal(as.numeric(tsls_other[,1]),other_est)
  expect_equal(as.numeric(tsls_other[,2]),other_se)
  expect_equal(as.numeric(tsls_other[,3]),other_t)
  expect_equal(as.numeric(tsls_other[,4]),other_pval)
})

test_that("TSLS estimate with exogeneous covariates and intercept", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,"nearc4"]
  Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
          "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
          "reg668", "smsa66")
  X=card.data[,Xname]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,X=X)
  tsls_output = as.numeric(coef(foo)["TSLS",-1])
  tsls_other = coefOther(foo)$TSLS
  
  expect_equal(tsls_output[1],0.131503836)
  expect_equal(tsls_output[2],0.0549636726)
  expect_equal(tsls_output[3],2.3925591)
  expect_equal(tsls_output[4],1.679262*10^(-2))
  
  other_est = c(0.108271106,-0.002334938,-0.146775747,-0.144671501,
                0.111808309,-0.107814233,-0.007046452,0.040444546,
                -0.057917154,0.03845768,0.055088709,0.026757977,
                -0.190891226,0.018531104,3.773965141)
  other_se = c(0.023658571,0.000333497,0.053899859,0.027284623,
               0.031661988,0.041813679,0.03290727,0.031780563,
               0.037605896,0.046938701,0.052659734,0.048828702,
               0.050711333,0.021608589,0.934947017)
  other_t = c(4.5764009,-7.0013725,-2.7231193,-5.302309,
              3.5313104,-2.5784441,-0.2141305,1.2726189,
              -1.5401083,0.8193171,1.0461258,0.5479969,
              -3.7642715,0.8575805,4.0365551)
  other_pval = c(4.92290320848*10^(-6),3.11612469442*10^(-12),6.50437654555*10^(-3),
                 1.22662588576*10^(-7),4.19745945432*10^(-4),9.97198409145*10^(-3),
                 8.30459836493*10^(-1),2.03252114106*10^(-1),1.23639625283*10^(-1),
                 4.12670726440*10^(-1),2.95587380542*10^(-1),5.83734887817*10^(-1),
                 1.70242896565*10^(-4),3.91192788096*10^(-1),5.56008327986*10^(-5))
  expect_equal(as.numeric(tsls_other[,1]),other_est)
  expect_equal(as.numeric(tsls_other[,2]),other_se)
  expect_equal(as.numeric(tsls_other[,3]),other_t)
  expect_equal(as.numeric(tsls_other[,4]),other_pval)
})

  

  


