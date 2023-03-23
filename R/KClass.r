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
  if(!inherits(ivmodel,"ivmodel")) {
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
  degF = ivmodel$n - ivmodel$p - 1 # this looks strange, but it's from Wooldridge and Stock (2002)

  # Compute k-class estimator
  kPointEst = matrix(0,length(k),1+p) # First column is treatment effect, rest are exogeneous covariates estimates.
  kVarPointEst = matrix(0,length(k),1+p)
  kCILower = matrix(0,length(k),1+p)
  kCIUpper = matrix(0,length(k),1+p)
  
  for(i in 1:length(k)) {
    rankInverse = qrrank(qr(t(W) %*% W - k[i] * t(W) %*% qr.resid(ZXQR,W)))
    #if(rankInverse < ncol(W)) stop("k is out-of-range; the matrix inversion in the k-class estimator will fail.")
    inverseMat = solve(t(W) %*% W - k[i] * t(W) %*% qr.resid(ZXQR,W))
    kPointEst[i,] = as.numeric(inverseMat %*% ( t(W) %*% Y - k[i] *t(W) %*% qr.resid(ZXQR,Y)))

    if (manyweakSE && k[i] != 0) { ## From Hansen et al. (2008)
        udelta = as.numeric(Y - W %*% kPointEst[i,])
        sigmahat.sq = sum(udelta^2)/(ivmodel$n-ncol(W))
        alphatilde = sum( qr.fitted(ZXQR,udelta)^2)/sum(udelta^2)
        iotahat = qr.fitted(ZXQR,W)
        Xtilde = W - t(t(udelta)) %*% (t(udelta) %*% W) / sum(udelta^2)
        Vhat = qr.resid(ZXQR,Xtilde)
        Ptt = diag(W %*% solve(t(W) %*% W) %*% t(W))
        kappa = sum(Ptt^2)/ivmodel$L 
        tauT = ivmodel$L / ivmodel$n
            
        H = t(W) %*% qr.fitted(ZXQR,W) - alphatilde * t(W) %*% W
        sigmaB = sigmahat.sq * ( (1-alphatilde)^2 * t(Xtilde) %*% qr.fitted(ZXQR,Xtilde) + alphatilde^2 *  t(Xtilde) %*% qr.resid(ZXQR,Xtilde))
        A_other = matrix(0,1,ncol(sigmaB))
        for(j in 1:ivmodel$n) {
          A_other = A_other + udelta[j]^2 * Vhat[j,,drop=FALSE] / ivmodel$n
        }
        A = matrix(0,nrow(sigmaB),ncol(sigmaB))
        B = matrix(0,nrow(sigmaB),ncol(sigmaB))
        for(j in 1:ivmodel$n) {
          A = A + (Ptt[j] - tauT) * (t(iotahat[j,,drop=FALSE]) %*% A_other)
          B = B + (udelta[j]^2 - sigmahat.sq) * (t(Vhat[j,,drop=FALSE]) %*% Vhat[j,,drop=FALSE])
        }
        B = ivmodel$L * (kappa - tauT) * B / (ivmodel$n*(1-2*tauT + kappa * tauT))
        Sigmahat = sigmaB + A + t(A) + B
            
        Hinv = solve(H)
        Lambdahat = Hinv  %*% Sigmahat %*% Hinv 
        kVarPointEst[i,] = diag(Lambdahat)
    } else if (heteroSE || (manyweakSE &&k[i] == 0)) {
        inner = matrix(0,ncol(W),ncol(W))
        for(j in 1:length(Y)) {
          adjustVec = as.numeric(W[j,]) - as.numeric(k[i]*qr.resid(ZXQR,W)[j,])
          inner = inner + (Y[j] - sum(W[j,] * kPointEst[i,]))^2 * (t(t(adjustVec)) %*% t(adjustVec) )
        }
        kVarPointEst[i,] = diag(inverseMat %*% inner %*% inverseMat)
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
	  	  nCluster <- length(unique(clusterID)); uniqueclusterID = unique(clusterID)
	  	  inner = matrix(0,ncol(W),ncol(W))
	  	  for(j in uniqueclusterID) {
	  	    clusterSame = which(clusterID == j)
	  	    innerC = rep(0,ncol(W))
	  	    for(t in clusterSame) {
	  	      innerC = innerC + (Y[t] - sum(W[t,] * kPointEst[i,])) * 
	  	                (as.numeric(W[t,]) - as.numeric(k[i]*qr.resid(ZXQR,W)[t,]))
	  	    }
	  	    inner = inner + t(t(innerC)) %*% t(innerC)
	  	  }
	  	    kVarPointEst[i,] = diag(inverseMat %*% inner %*% inverseMat)
	    }
	    else {
	      kVarPointEst[i,] = 1/degF * sum( (Y - W %*% kPointEst[i,])^2) * diag(inverseMat)
	    }
	  # Note that R's AER package uses the z-scores instead of the t distribution. They are basically the same as n gets large.
	   kCILower[i,] = kPointEst[i,] - qt(1-alpha/2,degF) * sqrt(kVarPointEst[i,])
	   kCIUpper[i,] = kPointEst[i,] + qt(1-alpha/2,degF) * sqrt(kVarPointEst[i,])
	    
  }

  # Compute test statistics
  kTestStat = (kPointEst - beta0)/sqrt(kVarPointEst) 

  # Compute p-value
  kPValue = 2*(1 - pt(abs(kTestStat),degF))

  # Package output
  kPointEst_trt = matrix(kPointEst[,1],length(k),1)
  kVarPointEst_trt = matrix(sqrt(kVarPointEst)[,1],length(k),1)
  kCI_trt = matrix(c(kCILower[,1],kCIUpper[,1]),length(k),2)
  kPValue_trt = kPValue[,1,drop=FALSE]
  kTestStat_trt = kTestStat[,1,drop=FALSE]
  
  rownames(kPointEst_trt) = rownames(kVarPointEst_trt) = rownames(kCI_trt) = rownames(kPValue_trt) = rownames(kTestStat_trt) = k
  colnames(kPValue_trt) = colnames(kTestStat_trt) = beta0
  colnames(kPointEst_trt) = "Estimate"
  colnames(kVarPointEst_trt) = "Std. Error"
  colnames(kCI_trt) = c(paste(as.character(round(alpha/2 * 100,1)),"%"),paste(as.character( round((1-alpha/2) * 100,1)),"%"))
  
  return(list(k = k,
              point.est = kPointEst_trt,std.err = kVarPointEst_trt,test.stat = kTestStat_trt,p.value = kPValue_trt,ci = kCI_trt,
              point.est.other = kPointEst[,-1,drop=FALSE],std.err.other = sqrt(kVarPointEst[,-1,drop=FALSE]),
              test.stat.other = kTestStat[,-1,drop=FALSE],p.value.other = kPValue[,-1,drop=FALSE],
              ci.other.lower = kCILower[,-1,drop=FALSE],ci.other.upper=kCIUpper[,-1,drop=FALSE]))
}
