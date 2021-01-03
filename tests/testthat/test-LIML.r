context("Check accuracy of LIML using Card data")

test_that("TSLS = LIML when you have single IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,"nearc4"]
  
  # No covar, no intercept
  foo = ivmodel(Y=Y,D=D,Z=Z,intercept=FALSE)
  tsls_output = as.numeric(coef(foo)["TSLS",])
  LIML_output = as.numeric(coef(foo)["LIML",])
  expect_equal(tsls_output,LIML_output)
  
  # No covar, intercept
  foo = ivmodel(Y=Y,D=D,Z=Z)
  tsls_output = as.numeric(coef(foo)["TSLS",])
  LIML_output = as.numeric(coef(foo)["LIML",])
  expect_equal(tsls_output,LIML_output)
  expect_equal(coefOther(foo)$TSLS,coefOther(foo)$LIML)
  
  # Covar, no intercept
  Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
          "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
          "reg668", "smsa66")
  X=card.data[,Xname]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,X=X,intercept=FALSE)
  tsls_output = as.numeric(coef(foo)["TSLS",])
  LIML_output = as.numeric(coef(foo)["LIML",])
  expect_equal(tsls_output,LIML_output)
  expect_equal(coefOther(foo)$TSLS,coefOther(foo)$LIML)
  
  # Covar, intercept
  foo = ivmodel(Y=Y,D=D,Z=Z,X=X)
  tsls_output = as.numeric(coef(foo)["TSLS",])
  LIML_output = as.numeric(coef(foo)["LIML",])
  expect_equal(tsls_output,LIML_output)
  expect_equal(coefOther(foo)$TSLS,coefOther(foo)$LIML)
})

test_that("LIML estimate without exogeneous covariates, mult-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,c("nearc4","nearc2")]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,intercept=FALSE)
  est_out = as.numeric(coef(foo)["LIML",])
  expect_equal(est_out,c(1.00142462444e+00,4.67565741669e-01,1.88656148938e-03,2.47840181357e+02,0.00000000000e+00))
  
  foo = ivmodel(Y=Y,D=D,Z=Z)
  est_out = as.numeric(coef(foo)["LIML",])
  expect_equal(est_out,c(1.00110965158e+00, 2.06278670996e-01, 2.79605430592e-02, 7.37749158014e+00, 2.07611705605e-13))
  est_other_out = as.numeric(coefOther(foo)$LIML)
  expect_equal(est_other_out,c(3.525864054210,0.371009242718,9.503439936918,0))
  
})

test_that("LIML estimate with exogenous covariates, multi-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,c("nearc4","nearc2")]
  Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
          "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
          "reg668", "smsa66")
  X=card.data[,Xname]
  
  foo = ivmodel(Y=Y,D=D,Z=Z,X=X,intercept=FALSE)
  est_out = as.numeric(coef(foo)["LIML",])
  expect_equal(est_out,c(1.000859911203,0.331859012784,0.012046871349,27.547319396786, 0))
  est_other_out = coefOther(foo)$LIML
  resulting_val = c(0.23364183864109, 0.02021057771855, 11.560374072169, 0,
                   -0.00436821301809, 0.00110332841567, -3.959123100652, 7.69751447831e-05,
                    0.03312459166842, 0.03308960347837,  1.001057377133, 3.16879939962e-01,
                   -0.12105392121076, 0.04407730255489, -2.746400396440, 6.06133942841e-03,
                    0.06046339337471, 0.04137629128184,  1.461305291063, 1.44036581582e-01,
                    0.05087023578479, 0.08929079929188,  0.569714194388, 5.68914289109e-01,
                    0.16757963517009, 0.07661009933952,  2.187435293974, 2.87878235965e-02,
                    0.21488920113501, 0.07673872120887,  2.800270811786, 5.13871753283e-03,
                    0.10306279640681, 0.09436296666945,  1.092195381773, 2.74835079485e-01,
                    0.26499409700390, 0.08435145473479,  3.141547443812, 1.69697892460e-03,
                    0.30606352407280, 0.09210857765334,  3.322855828094, 9.01715695532e-04,
                    0.23692819629421, 0.08507584157808,  2.784905701776, 5.38798821139e-03,
                   -0.12131904963092, 0.10436146853516, -1.162488908347, 2.45129548098e-01,
                    0.01642425108771, 0.03593975566316,  0.456993955152, 6.47708533398e-01)
  resulting_val = matrix(resulting_val,length(Xname),4,byrow=TRUE)
  colnames(resulting_val) = c("Estimate","Std. Error","t value","Pr(>|t|)")
  rownames(resulting_val) = Xname
  expect_equal(est_other_out,resulting_val)
  
  foo = ivmodel(Y=Y,D=D,Z=Z,X=X)
  est_out = as.numeric(coef(foo)["LIML",])
  expect_equal(est_out,c(1.0004094273165,0.1640277561000,0.0554950702135,2.9557176064271,0.0031437760692))
  est_other_out = coefOther(foo)$LIML
  resulting_val = c(0.12168991721364, 0.023982112015946,  5.0741951806718, 4.12895895474e-07,
                   -0.00236235860783, 0.000352513445745, -6.7014709264197, 2.45576892155e-11,
                   -0.11687046279985, 0.054741901014118, -2.1349361391325, 3.28471930607e-02,
                   -0.14279170805679, 0.028847823460314, -4.9498260502471, 7.83906274560e-07,
                    0.09773848045461, 0.032643225029197,  2.9941429000103, 2.77469235253e-03,
                   -0.10165672445578, 0.044113721055252, -2.3044241570203, 2.12669732909e-02,
                    0.00163040341363, 0.034504112804948,  0.0472524369151, 9.62315187055e-01,
                    0.04873104057901, 0.033329378391302,  1.4621046935496, 1.43817473876e-01,
                   -0.05472430780333, 0.039747947181537, -1.3767832475320, 1.68682209993e-01,
                    0.05506160552795, 0.048860435931131,  1.1269159695088, 2.59868319889e-01,
                    0.07406188766618, 0.054781913834350,  1.3519404942683, 1.76496469376e-01,
                    0.04241390943707, 0.050976830315037,  0.8320232775352, 4.05462161154e-01,
                   -0.19998558533462, 0.053429008023303, -3.7430151285495, 1.85242733723e-04,
                    0.01411679795615, 0.022738627996378,  0.6208289241724, 5.34759446929e-01,
                    3.22126944354386, 0.944077217813185,  3.4120825953257, 6.53239525561e-04)
  resulting_val = matrix(resulting_val,length(Xname)+1,4,byrow=TRUE)
  colnames(resulting_val) = c("Estimate","Std. Error","t value","Pr(>|t|)")
  rownames(resulting_val) = c(Xname,"intercept")
  expect_equal(est_other_out,resulting_val)

})






