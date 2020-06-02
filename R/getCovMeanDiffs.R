#returns the covariate mean differences between two groups (1 or 0), i.e.,
#returns (\bar{X}_T - \bar{X}_C), where
#\bar{X}_T is the vector of covariate means for treated units,
#\bar{X}_C is the vector of covariate means for control units.
getCovMeanDiffs = function(X, indicator){

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

  return(covMeanDiffs)
}