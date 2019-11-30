print.ivmodel <- function(x, ...){
  ivmodel<-x
### print formula
  cat("\nCall:\n", paste(deparse(ivmodel$call), sep = "\n", collapse = "\n"), 
      "\n", sep = "")
  cat("sample size: ", ivmodel$n, "\n", sep="")
  
### first stage regression result
  cat(rep("_", 30), "\n")
  cat("\nFirst Stage Regression Result:\n\n")

  SSM=c(sum(qr.fitted(ivmodel$ZadjQR, ivmodel$Dadj)^2))
  SST=sum(ivmodel$Dadj^2)
  SSE=SST-SSM
  DM=ivmodel$L
  DE=ivmodel$n-ivmodel$p-ivmodel$L
  DT=DM+DE
  
  Fstat=SSM/SSE*DE/DM
  pval=1-pf(Fstat, df1=DM, df2=DE)
  RSquare=SSM/SST
  adjRSquare=1-(1-RSquare)*DT/DE
  RMSE=sqrt(SSE/DE)
  
  cat("F=", Fstat, ", df1=", DM, ", df2=", DE, ", p-value is ", format.pval(pval), "\n", sep="")
  cat("R-squared=", RSquare, ",   Adjusted R-squared=", adjRSquare, "\n", sep="")
  cat("Residual standard error: ", RMSE, " on ", DT, " degrees of freedom\n", sep="")
	  
### Sargan test
  if(ivmodel$L>1){

  cat(rep("_", 30), "\n")
  cat("\nSargan Test Result:\n\n")

  TSLS=sum(ivmodel$Dadj*qr.fitted(ivmodel$ZadjQR, ivmodel$Yadj))/
       sum(ivmodel$Dadj*qr.fitted(ivmodel$ZadjQR, ivmodel$Dadj))
  e=ivmodel$Yadj-ivmodel$Dadj*TSLS
  Sargan=sum(qr.fitted(ivmodel$ZadjQR, e)^2)/(sum(e^2)/length(e))
  pval=1-pchisq(Sargan, df=ivmodel$L-1)
  
  cat("Sargan Test Statistics=", Sargan, ", df=", ivmodel$L-1, ", p-value is ", format.pval(pval), "\n", sep="")
  
  }
  
### print TSLS, kClass, LIML, Fuller
  cat(rep("_", 30), "\n")
  cat("\nCoefficients of k-Class Estimators:\n\n")
  printCoefmat(coef(ivmodel), digits = max(3L, getOption("digits") - 3L))
	  
### print AR, ARsens, CLR
  cat(rep("_", 30), "\n")
  cat("\nAlternative tests for the treatment effect under H_0: beta=", ivmodel$beta0, ".\n",sep = "")
  if(!is.null(ivmodel$AR)){
    cat("\nAnderson-Rubin test (under F distribution):\n")
	cat("F=", ivmodel$AR$Fstat, ", df1=", ivmodel$AR$df[1], ", df2=", 
	    ivmodel$AR$df[2], ", p-value is ", format.pval(ivmodel$AR$p.value), "\n", sep="")
	cat(round((1-ivmodel$alpha)*100, digits=1), "percent confidence interval:\n", 
	    ivmodel$AR$ci.info)
  }
  if(!is.null(ivmodel$ARsens)){
    cat("\n\nSensitivity analysis with deltarange [", ivmodel$ARsens$deltarange[1], 
	    ", ", ivmodel$ARsens$deltarange[2], "]:\n")
	cat("non-central F=", ivmodel$ARsens$ncFstat, ", df1=", ivmodel$ARsens$df[1], 
	    ", df2=", ivmodel$ARsens$df[2], ", ncp=", ivmodel$ARsens$ncp, ", p-value is ", format.pval(ivmodel$ARsens$p.value), "\n", sep="")
	cat(round((1-ivmodel$alpha)*100, digits=1), "percent confidence interval:\n", 
	    ivmodel$ARsens$ci.info)
  }
  if(!is.null(ivmodel$CLR)){
    cat("\n\nConditional Likelihood Ratio test (under Normal approximation):\n")
	cat("Test Stat=", ivmodel$CLR$test.stat, ", p-value is ", format.pval(ivmodel$CLR$    p.value), "\n", sep="")
	cat(round((1-ivmodel$alpha)*100, digits=1), "percent confidence interval:\n", 
	    ivmodel$CLR$ci.info)
  }
  cat("\n")
}