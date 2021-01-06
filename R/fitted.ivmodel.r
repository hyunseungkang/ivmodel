fitted.ivmodel <- function(object, ...){
  ivmodel <- object
  residual_est =  resid(ivmodel)
  stopifnot("Number of estimated residuals do not equal number of outcome"= nrow(residual_est) == nrow(ivmodel$Y))
  
  result = matrix(ivmodel$Y,nrow(residual_est),ncol(residual_est)) - residual_est
  result = result[which(ivmodel$naindex),]
  return(result)
}
