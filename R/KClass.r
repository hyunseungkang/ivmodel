### kClass: Generates k-class estimators (point est, se, p-value, and confint)
###         Must run ivmodel before you run this code
### INPUT: ivmodel, an object from ivmodel() function
###        k, a vector (or scalar) of k values for estimation
###        beta0, a scalar null value for testing and p-values
###        alpha, significance level for confidence intervals
###        heteroSE, use heteroscedastic robust standard errors?
###        clusterID, if cluster-robust standard errors are desired, please provide a vector of
###                   length that's identical to the sample size.
###                   For example, if n = 6 and clusterID = c(1,1,1,2,2,2),
###                   there would be two clusters where the first cluster
###                   is formed by the first three observations and the
###                   second cluster is formed by the last three observations
###                   clusterID can be numeric, character, or factor.
### OUTPUT: a list of point estimate, standard error, test statistic, and p-value
KClass = function(ivmodel,
                  beta0=0,alpha=0.05,k=c(0,1),
                  manyweakSE=FALSE,heteroSE=FALSE,clusterID=NULL) {
  # Error checking
  if(class(ivmodel) != "ivmodel") {
    print("You must supply an ivmodel class. Run ivmodel() and see ivmodel() function for details")
    return(NULL)
  }
  if(missing(k)) {
    print("You must specify a value for k")
    return(NULL)
  }

  # Extract objects from ivmodel
  Y= ivmodel$Y; D = ivmodel$D; Z = ivmodel$Z; ZXQR = ivmodel$ZXQR; p = ivmodel$p
  if(p == 0) W = D
  if(p > 0) W = cbind(D,ivmodel$X)
  
  #Yadj = ivmodel$Yadj; Dadj =ivmodel$Dadj; Zadj = ivmodel$Zadj; ZadjQR = ivmodel$ZadjQR
  degF = ivmodel$n - ivmodel$p - 1 # this looks strange, but it's from Wooldridge and Stock (2002)

  # Compute k-class estimator
  #denom = (sum(Dadj^2) - k * sum(qr.resid(ZadjQR,Dadj)^2))
  #kPointEst = (sum(Yadj * Dadj) - k * sum(Yadj * qr.resid(ZadjQR,Dadj))) / denom

  # Compute std error, testStat, and confidence intervals
  kPointEst = matrix(0,length(k),1+p) # First column is treatment effect, rest are exogeneous covariates estimates.
  kVarPointEst = matrix(0,length(k),1+p)
  kCILower = matrix(0,length(k),1+p)
  kCIUpper = matrix(0,length(k),1+p)
  
  for(i in 1:length(k)) {
    rankInverse = qrrank(qr(t(W) %*% W - k[i] * t(W) %*% qr.resid(ZXQR,W)))
    if(rankInverse < ncol(W)) stop("k is out-of-range; the matrix inversion in the k-class estimator will fail.")
    inverseMat = solve(t(W) %*% W - k[i] * t(W) %*% qr.resid(ZXQR,W))
    kPointEst[i,] = as.numeric(inverseMat %*% ( t(W) %*% Y - k[i] *t(W) %*% qr.resid(ZXQR,Y)))
    if(1 <= 0) { #k exceeds the maximum possible value
	    kVarPointEst[i,] = NA; 
	    kCILower[i,] = c(NA,NA); kCIUpper[i,] = c(NA,NA);
    } else {
        if (manyweakSE && k[i] != 0) { ## From Hansen et al. (2008)
            alpha <- - (1 - k[i])/k[i]
            H <- sum((Dadj - qr.resid(ZadjQR,Dadj))^2) - alpha * sum(Dadj^2)
            u <- Yadj - Dadj * kPointEst[i]
            D.tilde <- Dadj - u * sum(u  * Dadj) / sum(u^2)
            sigmaB <- sum(u^2) / degF * ( (1 - alpha)^2 * (sum(D.tilde^2) - sum(qr.resid(ZadjQR, D.tilde)^2)) + alpha^2 * sum(qr.resid(ZadjQR, D.tilde)^2) )
            Q <- qr.Q(ZadjQR)[, 1:ivmodel$L]
            P.diag <- apply(Q * Q, 1, sum)
            Dfitted <- Dadj - qr.resid(ZadjQR, Dadj)
            kappa <- sum(P.diag^2)/ivmodel$L
            tau <- ivmodel$L/ivmodel$n
            A <- sum((P.diag - tau) * Dfitted) * mean(u^2 * qr.resid(ZadjQR, D.tilde))
            B <- ivmodel$L * (kappa - tau) * sum((u^2 - sum(u^2) / degF) * qr.resid(ZadjQR, D.tilde)^2) / (ivmodel$n * (1 - 2 * tau + kappa * tau))
            kVarPointEst[i] = (sigmaB + 2 * A + B) / H^2
        } else if (heteroSE || (manyweakSE &&k[i] == 0)) {
	    kVarPointEst[i] = ivmodel$n/degF * sum( (Yadj - Dadj * kPointEst[i])^2 * (Dadj - k[i]*qr.resid(ZadjQR,Dadj))^2 ) / (denom[i])^2
	  } else if(!is.null(clusterID)){
	  	if(length(clusterID) != ivmodel$n) {
	  		### Missing problem here needs to be taken care of ###
	  		print("Cluster ID vector is not the same length as the sample size")
	  		return(NULL)
	  	}
	  	if(!is.character(clusterID) && !is.factor(clusterID) && !is.numeric(clusterID)) {
	  		print("Cluster ID must be either a character vector, a factor vector, or a numeric vector")
	  		return(NULL)
	  	}
	  	clusterID = as.factor(clusterID)
	  	clusterID = as.numeric(clusterID)
	  	nCluster <- length(unique(clusterID))
        dfc <- nCluster/(nCluster-1)*(ivmodel$n-1)/degF
        residFit = Yadj - Dadj * kPointEst[i]
        umat = tapply(residFit*(Dadj - k[i]*qr.resid(ZadjQR,Dadj)),clusterID, sum)
        kVarPointEst[i] = dfc * t(umat) %*% umat / denom[i]^2
	  }
	  else {
        #kVarPointEst[i] = 1/(degF) *sum((Yadj - Dadj * kPointEst[i])^2) / denom[i]
	   kVarPointEst[i,] = 1/degF * sum( (Y - W %*% kPointEst[i,])^2) * diag(inverseMat)
	  }
	  # Note that R's AER package uses the z-scores instead of the t distribution. They are basically the same as n gets large.
	    kCILower[i,] = kPointEst[i,] - qt(1-alpha/2,degF) * sqrt(kVarPointEst[i,])
	    kCIUpper[i,] = kPointEst[i,] + qt(1-alpha/2,degF) * sqrt(kVarPointEst[i,])
	    
      #kCI[i,] =  c(kPointEst[i] - qt(1-alpha/2,degF) * sqrt(kVarPointEst[i]),
	    #            kPointEst[i] + qt(1-alpha/2,degF) * sqrt(kVarPointEst[i]))
    }
  }

  # Compute test statistics
  # This operation takes a vector of k estimates (pointEst, k*1 dim) and
  # another vector of null values (beta0,1 * length(beta0))
  # and subtracts pointEst from a vector of null values.
  # The last operation of dividing by sqrt() is done column-wise
  # E.g. [1,2,3 ; 2,3,4] / (2,3,4) --> [0.5,2/3,0.75;1,1,1]
  kTestStat = (kPointEst - beta0)/sqrt(kVarPointEst) #= outer(kPointEst,beta0,FUN="-") / sqrt(kVarPointEst)

  # Compute p-value
  kPValue = 2*(1 - pt(abs(kTestStat),degF))

  # Package output
  kPointEst_trt = matrix(kPointEst[,1],length(k),1)
  kVarPointEst_trt = matrix(sqrt(kVarPointEst)[,1],length(k),1)
  kCI_trt = matrix(NA,length(k),2)
  kCI_trt[,1] = kPointEst[,1] - qt(1-alpha/2,degF) * sqrt(kVarPointEst[,1])
  kCI_trt[,2] = kPointEst[,1] + qt(1-alpha/2,degF) * sqrt(kVarPointEst[,1])
  kPValue_trt = kPValue[,1,drop=FALSE]
  kTestStat_trt = kTestStat[,1,drop=FALSE]
  
  rownames(kPointEst_trt) = rownames(kVarPointEst_trt) = rownames(kCI_trt) = rownames(kPValue_trt) = rownames(kTestStat_trt) = k
  colnames(kPValue_trt) = colnames(kTestStat_trt) = beta0
  colnames(kPointEst_trt) = "Estimate"
  colnames(kVarPointEst_trt) = "Std. Error"
  colnames(kCI_trt) = c(paste(as.character(round(alpha/2 * 100,1)),"%"),paste(as.character( round((1-alpha/2) * 100,1)),"%"))
  
  return(list(point.est = kPointEst_trt,std.err = kVarPointEst_trt,test.stat = kTestStat_trt,p.value = kPValue_trt,ci = kCI_trt))
}
