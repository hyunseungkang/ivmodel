#Returns the Mahalanobis distance (MD)
#according to an indicator (instrument or exposure).
#Specifically, the MD for treatment-versus-control (T or C) is defined as:
#((Nt*(N-Nt))/N)*(\bar{X}_T - \bar{X}_C)%*%[cov(X)]^{-1}%*%(\bar{X}_T - \bar{X}_C)
#Here, Nt is the number of treated units (or indicator = 1 units),
#N is the number of units,
#\bar{X}_T is the vector of covariate means for treated units,
#\bar{X}_C is the vector of covariate means for control units,
#and cov(X) is the covariate-covariance matrix.
getMD = function(X, indicator, covX.inv = NULL){

  data = data.frame(X, indicator = indicator)
  #treatment group
  treatmentData = subset(data, indicator == 1)
  #control group
  controlData = subset(data, indicator == 0)

  #now we can get rid of the indicator variable
  treatmentData = subset(treatmentData, select = -c(indicator))
  controlData = subset(controlData, select = -c(indicator))

  #covariate mean difference
  covMeanDiffs = colMeans(treatmentData) - colMeans(controlData)
  if(is.null(covX.inv)){
    covX.inv = solve(as.matrix(stats::cov(X)))
  }

  n = as.numeric(nrow(data))
  n.t = as.numeric(sum(indicator))
  md = ((n.t*(n-n.t))/n)*t(covMeanDiffs)%*%covX.inv%*%covMeanDiffs
  return(md)
}