vcov.ivmodel<-function(object, ...){
  ivmodel<-object
  vcovmat <- matrix(NA, ncol=2, nrow=0)
  colnames(vcovmat) <- c("k","Endogenous variable")
  if(!is.null(ivmodel$kClass)){
    temp<-cbind(as.numeric(rownames(ivmodel$kClass$ci)),ivmodel$kClass$std.err)
	rownames(temp) <- rep("k-class", nrow(temp))
    rownames(temp)[temp[,1]==0] <- "OLS"
    rownames(temp)[temp[,1]==1] <- "TSLS"
    vcovmat <- rbind(vcovmat, temp)
  }
  if(!is.null(ivmodel$LIML)){
    temp<-cbind(ivmodel$LIML$k, ivmodel$LIML$std.err)
	rownames(temp) <- "LIML"
    vcovmat <- rbind(vcovmat, temp)
  }
  if(!is.null(ivmodel$Fuller)){
    temp<-cbind(ivmodel$Fuller$k, ivmodel$Fuller$std.err)
	rownames(temp) <- "Fuller"
    vcovmat <- rbind(vcovmat, temp)
  }
  
  vcovmat<-vcovmat[sort(vcovmat[,1], index.return=T)$ix,]  
  return(vcovmat)
}