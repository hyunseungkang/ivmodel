vcovOther <-function(ivmodel){
  if(ivmodel$p == 0) stop("There are no additional exogenous covariates!")
  
  varMatrix = matrix(0,0,1+ivmodel$p)
  if(!is.null(ivmodel$kClass)){
    for(i in 1:length(ivmodel$kClass$k)) {
      temp = as.numeric(c(ivmodel$kClass$k[i],ivmodel$kClass$std.err.other[i,]))
      if(ivmodel$kClass$k[i] == 0) {
        varMatrix = rbind(varMatrix,OLS = temp)
      } else if(ivmodel$kClass$k[i] == 1) {
        varMatrix = rbind(varMatrix,TSLS = temp)
      } else {
        varMatrix = rbind(varMatrix,temp)
        print("hello")
        rownames(varMatrix) = c(rownames(varMatrix)[-nrow(varMatrix)],ivmodel$kClass$k[i])
      }
    }
  }
  if(!is.null(ivmodel$LIML)){
    temp = c(as.numeric(ivmodel$LIML$k), as.numeric(ivmodel$LIML$std.err.other))
    varMatrix = rbind(varMatrix,LIML = temp)
  }
  if(!is.null(ivmodel$Fuller)){
    temp = c(as.numeric(ivmodel$Fuller$k),as.numeric(ivmodel$Fuller$std.err.other))
    varMatrix = rbind(varMatrix,Fuller = temp)
  }
  colnames(varMatrix) = c("k",colnames(ivmodel$X))
  varMatrix = varMatrix^2
  return(varMatrix)
}