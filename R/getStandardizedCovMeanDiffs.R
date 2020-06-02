#returns the standardized covariate mean differences between two groups (1 or 0), i.e.,
#returns (\bar{X}_T - \bar{X}_C)/sd(X), where
#\bar{X}_T is the vector of covariate means for treated units,
#\bar{X}_C is the vector of covariate means for control units,
#and sd(X) is the square root of the pooled variance between groups.
getStandardizedCovMeanDiffs = function(X, indicator){

  data = data.frame(X, indicator = indicator)
  #treatment group
  treatmentData = subset(data, indicator == 1)
  #control group
  controlData = subset(data, indicator == 0)

  #now we can get rid of the indicator variable
  treatmentData = subset(treatmentData, select = -c(indicator))
  controlData = subset(controlData, select = -c(indicator))

  #variance of each covariate within treatment/control groups
    cov.variances.treatment = apply(treatmentData, MARGIN = 2, FUN = var)
    cov.variances.control = apply(controlData, MARGIN = 2, FUN = var)

    pooled.cov.variance = (cov.variances.treatment + cov.variances.control)/2

    #covariate mean difference
    covMeanDiffs = colMeans(treatmentData) - colMeans(controlData)

    covMeanDiffs.standardized = covMeanDiffs/sqrt(pooled.cov.variance)
    return(covMeanDiffs.standardized)
}