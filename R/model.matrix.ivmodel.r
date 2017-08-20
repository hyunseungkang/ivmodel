model.matrix.ivmodel <- function(object, ...){
  ivmodel <- object
  
  if(is.null(dim(ivmodel$X))) {
    return(ivmodel$D)
  } else {
    intVector = rep(1,nrow(ivmodel$X))
    whichIntIndex = NA
    for(i in ncol(ivmodel$X)) {
      if(all(ivmodel$X[,i] == intVector)) {
        whichIntIndex = i
        break
      }
    }
    if(is.na(whichIntIndex)) {
      return(cbind(ivmodel$D,ivmodel$X))
    } else {
      if(ncol(ivmodel$X) == 1) {
        matOut = cbind(1,ivmodel$D)
        colnames(matOut)[1] = "(Intercept)"
      } else {
        matOut = cbind(1,ivmodel$D,ivmodel$X[,-whichIntIndex])
        colnames(matOut)[1] = "(Intercept)"
      }
      return(matOut)
    }
  }
}