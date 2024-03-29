### Fuller: Generates point estimates and standard errors
###            Must run ivmodel before you run this code
### INPUT: ivmodel, an object from ivmodel() function
###        b, adjustment by Fuller
###        beta0, a vector (or scalar) of null values for testing and p-values
###        alpha, significance level for confidence intervals
### OUTPUT: a list of point estimate, standard error, test statistic, and p-value
Fuller = function(ivmodel,
                  beta0=0,alpha=0.05,b=1,
                  manyweakSE=FALSE,heteroSE=FALSE,clusterID=NULL) {
  # Error checking
  if(!inherits(ivmodel,"ivmodel")) {
    print("Fuller: You must supply an ivmodel class. See ivmodel function for details")
	return(NULL)
  }
  if(b <= 0) {
    print("Fuller: Fuller adjustment cannot be less than 0")
	return(NULL)
  }

  # Extract objects from ivmodel
  Yadj = ivmodel$Yadj; Dadj = ivmodel$Dadj; Zadj = ivmodel$Zadj; ZadjQR = ivmodel$ZadjQR

  # Value of k for LIML
  LIMLMatrix1 = matrix(0,2,2)
  LIMLMatrix1[1,1] = sum(Yadj^2)
  LIMLMatrix1[1,2] = LIMLMatrix1[2,1] = sum(Dadj * Yadj)
  LIMLMatrix1[2,2] = sum(Dadj^2)

  LIMLMatrix2 = matrix(0,2,2); projYadj = qr.resid(ZadjQR,Yadj); projDadj = qr.resid(ZadjQR,Dadj)
  LIMLMatrix2[1,1] = sum(projYadj^2)
  LIMLMatrix2[1,2] = LIMLMatrix2[2,1] = sum(projDadj * projYadj)
  LIMLMatrix2[2,2] = sum(projDadj^2)

  kLIML = eigen(LIMLMatrix1 %*% invTwobyTwoSymMatrix(LIMLMatrix2))$values[2]
  kFuller = kLIML - b/(ivmodel$n - ivmodel$L - ivmodel$p)
  output = KClass(ivmodel,k=kFuller,beta0=beta0,alpha=alpha,manyweakSE=manyweakSE,heteroSE=heteroSE,clusterID=clusterID)

  # Package output
  rownames(output$point.est) = rownames(output$std.err) = rownames(output$test.stat) = rownames(output$ci) = rownames(output$p.value) = NULL
  return(c(output,list(k=kFuller)))
}
