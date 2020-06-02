#create a Love plot of the standardized covariate mean differences
#across the exposure (D) and the instrument (Z).
#Can also display the permutation quantiles for these quantities.
#This function is used to create Figure 3a in Branson and Keele (2020).
balanceLovePlot = function(X, D, Z, permQuantiles = FALSE, alpha = 0.05, perms = 1000){
  #compute standardized covariate mean differences for D and Z
  covMeanDiff.d = getStandardizedCovMeanDiffs(X = X, indicator = D)
  covMeanDiff.z = getStandardizedCovMeanDiffs(X = X, indicator = Z)

  #number of covariates
  K = ncol(X)

  #get range of covMeanDiffs for plot limits
  plot.min = min( c(covMeanDiff.d, covMeanDiff.z) )
  plot.max = max( c(covMeanDiff.d, covMeanDiff.z) )

  #make love plot
  graphics::plot(covMeanDiff.d, 1:K,
    xlim = c(plot.min, plot.max), xlab = "Balance", ylab = "", main = "", yaxt = "n",
    pch = 16, col = "red")
  graphics::points(covMeanDiff.z, 1:K, pch = 17, col = "blue")
  graphics::abline(v = 0, col = "gray")
  graphics::axis(side=2, at=1:K, labels = names(covMeanDiff.z), las = 1, cex.axis = 0.5)
  if(permQuantiles == TRUE){
    #the covMeanDiffs (across permutations) are
    permutations.covMeanDiffs = getCompletePerms.balance(X = X, indicator = Z, perms = perms)
    #the quantiles are
    permutations.covMeanDiffs.lowerQuantile = apply(permutations.covMeanDiffs, MARGIN = 2, FUN = stats::quantile, probs = alpha/2)
    permutations.covMeanDiffs.upperQuantile = apply(permutations.covMeanDiffs, MARGIN = 2, FUN = stats::quantile, probs = 1 - alpha/2)
    #add the quantile lines to the plot
    graphics::lines(permutations.covMeanDiffs.lowerQuantile, 1:K, lty = 2)
    graphics::lines(permutations.covMeanDiffs.upperQuantile, 1:K, lty = 2)
  }

}