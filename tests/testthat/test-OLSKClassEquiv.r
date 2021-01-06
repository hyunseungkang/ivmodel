context("Check against OLS for k-class")


test_that("OLS estimates from lm should be identical to k-class with k = 0, single-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,"nearc4"]
  
  olsfit = lm(Y ~ D - 1)
  estcoef = as.numeric(coef(olsfit)); estvar = as.numeric(diag(vcov(olsfit)))
  estfit = as.numeric(predict(olsfit)); estresid = as.numeric(residuals(olsfit))
  ivfit = ivmodel(Y=Y,D=D,Z=Z,intercept=FALSE)
  estcoef_iv = as.numeric(as.numeric(coef(ivfit)["OLS","Estimate"])); 
  estvar_iv = as.numeric(vcov(ivfit)["OLS",2])
  estfit_iv = as.numeric(fitted(ivfit)[,"OLS"]); estresid_iv = as.numeric(residuals(ivfit)[,"OLS"])
  
  expect_equal(estcoef,estcoef_iv)
  expect_equal(estvar,estvar_iv)
  expect_equal(estfit,estfit_iv)
  expect_equal(estresid,estresid_iv)
  
  olsfit = lm(Y ~ D)
  estcoef = as.numeric(coef(olsfit)); estvar = as.numeric(diag(vcov(olsfit)))
  estfit = as.numeric(predict(olsfit)); estresid = as.numeric(residuals(olsfit))
  ivfit = ivmodel(Y=Y,D=D,Z=Z)
  estcoef_iv = c(as.numeric(coefOther(ivfit)$OLS[,"Estimate"]),as.numeric(coef(ivfit)["OLS","Estimate"]))
  estvar_iv = c(as.numeric(vcovOther(ivfit)["OLS",-1]),as.numeric(vcov(ivfit)["OLS",2]))
  estfit_iv = as.numeric(fitted(ivfit)[,"OLS"]); estresid_iv = as.numeric(residuals(ivfit)[,"OLS"])
  
  expect_equal(estcoef,estcoef_iv)
  expect_equal(estvar,estvar_iv)
  expect_equal(estfit,estfit_iv)
  expect_equal(estresid,estresid_iv)
  
  
  Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
          "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
          "reg668", "smsa66")
  X= as.matrix(card.data[,Xname])
  
  olsfit = lm(Y ~ D + X)
  estcoef = as.numeric(coef(olsfit)); estvar = as.numeric(diag(vcov(olsfit)))
  estfit = as.numeric(predict(olsfit)); estresid = as.numeric(residuals(olsfit))
  ivfit = ivmodel(Y=Y,D=D,Z=Z,X=X)
  estcoef_iv = c(as.numeric(coefOther(ivfit)$OLS[nrow(coefOther(ivfit)$OLS),"Estimate"]),
                 as.numeric(coef(ivfit)["OLS","Estimate"]),
                 as.numeric(coefOther(ivfit)$OLS[-nrow(coefOther(ivfit)$OLS),"Estimate"]))
  estvar_iv = c(as.numeric(vcovOther(ivfit)["OLS","intercept"]),
                as.numeric(vcov(ivfit)["OLS",2]),
                as.numeric(vcovOther(ivfit)["OLS",-c(1,ncol(vcovOther(ivfit)))]))
  estfit_iv = as.numeric(fitted(ivfit)[,"OLS"]); estresid_iv = as.numeric(residuals(ivfit)[,"OLS"])
  
  expect_equal(estcoef,estcoef_iv)
  expect_equal(estvar,estvar_iv)
  expect_equal(estfit,estfit_iv)
  expect_equal(estresid,estresid_iv)
})

test_that("AR test without exogeneous covariates and intercept, mult-IV", {
  Y=card.data[,"lwage"]
  D=card.data[,"educ"] 
  Z=card.data[,c("nearc4","nearc2")]
  
  olsfit = lm(Y ~ D - 1)
  estcoef = as.numeric(coef(olsfit)); estvar = as.numeric(diag(vcov(olsfit)))
  estfit = as.numeric(predict(olsfit)); estresid = as.numeric(residuals(olsfit))
  ivfit = ivmodel(Y=Y,D=D,Z=Z,intercept=FALSE)
  estcoef_iv = as.numeric(as.numeric(coef(ivfit)["OLS","Estimate"])); 
  estvar_iv = as.numeric(vcov(ivfit)["OLS",2])
  estfit_iv = as.numeric(fitted(ivfit)[,"OLS"]); estresid_iv = as.numeric(residuals(ivfit)[,"OLS"])
  
  expect_equal(estcoef,estcoef_iv)
  expect_equal(estvar,estvar_iv)
  expect_equal(estfit,estfit_iv)
  expect_equal(estresid,estresid_iv)
  
  olsfit = lm(Y ~ D)
  estcoef = as.numeric(coef(olsfit)); estvar = as.numeric(diag(vcov(olsfit)))
  estfit = as.numeric(predict(olsfit)); estresid = as.numeric(residuals(olsfit))
  ivfit = ivmodel(Y=Y,D=D,Z=Z)
  estcoef_iv = c(as.numeric(coefOther(ivfit)$OLS[,"Estimate"]),as.numeric(coef(ivfit)["OLS","Estimate"]))
  estvar_iv = c(as.numeric(vcovOther(ivfit)["OLS",-1]),as.numeric(vcov(ivfit)["OLS",2]))
  estfit_iv = as.numeric(fitted(ivfit)[,"OLS"]); estresid_iv = as.numeric(residuals(ivfit)[,"OLS"])
  
  expect_equal(estcoef,estcoef_iv)
  expect_equal(estvar,estvar_iv)
  expect_equal(estfit,estfit_iv)
  expect_equal(estresid,estresid_iv)
  
  
  Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
          "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
          "reg668", "smsa66")
  X= as.matrix(card.data[,Xname])
  
  olsfit = lm(Y ~ D + X)
  estcoef = as.numeric(coef(olsfit)); estvar = as.numeric(diag(vcov(olsfit)))
  estfit = as.numeric(predict(olsfit)); estresid = as.numeric(residuals(olsfit))
  ivfit = ivmodel(Y=Y,D=D,Z=Z,X=X)
  estcoef_iv = c(as.numeric(coefOther(ivfit)$OLS[nrow(coefOther(ivfit)$OLS),"Estimate"]),
                 as.numeric(coef(ivfit)["OLS","Estimate"]),
                 as.numeric(coefOther(ivfit)$OLS[-nrow(coefOther(ivfit)$OLS),"Estimate"]))
  estvar_iv = c(as.numeric(vcovOther(ivfit)["OLS","intercept"]),
                as.numeric(vcov(ivfit)["OLS",2]),
                as.numeric(vcovOther(ivfit)["OLS",-c(1,ncol(vcovOther(ivfit)))]))
  estfit_iv = as.numeric(fitted(ivfit)[,"OLS"]); estresid_iv = as.numeric(residuals(ivfit)[,"OLS"])
  
  expect_equal(estcoef,estcoef_iv)
  expect_equal(estvar,estvar_iv)
  expect_equal(estfit,estfit_iv)
  expect_equal(estresid,estresid_iv)
})





