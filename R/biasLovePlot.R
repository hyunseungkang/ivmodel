#create a Love plot of the bias across the exposure (D) and the instrument (Z).
#Can also display the permutation quantiles for these quantities.
#Note that the bias is different for the exposure (D) than for the instrument (Z),
#as discussed in Equation (3) of Branson and Keele (2020).
#This function is used to create Figure 3b in Branson and Keele (2020).
biasLovePlot = function(X, D, Z, permQuantiles = FALSE, alpha = 0.05, perms = 1000){
  #compute standardized covariate mean differences for exposure and instrument
  covMeanDiff.d = getStandardizedCovMeanDiffs(X = X, indicator = D)
  covMeanDiff.z = getStandardizedCovMeanDiffs(X = X, indicator = Z)

  #number of covariates
  K = ncol(X)

  #mean difference across the instrument
  DMeanDiff = mean(D[Z == 1]) - mean(D[Z == 0])
  #thus the IV bias is
  bias.z = covMeanDiff.z/DMeanDiff

  #get range of covMeanDiffs for plot limits
  plot.min = min( c(covMeanDiff.d, bias.z) )
  plot.max = max( c(covMeanDiff.d, bias.z) )

  #make love plot
  graphics::plot(covMeanDiff.d, 1:K,
    xlim = c(plot.min, plot.max), xlab = "Bias", ylab = "", main = "", yaxt = "n",
    pch = 16, col = "red")
  graphics::points(bias.z, 1:K, pch = 17, col = "blue")
  graphics::abline(v = 0, col = "gray")
  graphics::axis(side=2, at=1:K, labels = names(covMeanDiff.z), las = 1, cex.axis = 0.5)
  if(permQuantiles == TRUE){
    #the covMeanDiffs (across permutations) are
    permutations.covMeanDiffs = getCompletePerms.balance(X = X, indicator = Z, perms = perms)
    #then, the bias is simply the covariate mean differences divided by the exposure mean difference
    permutations.bias = permutations.covMeanDiffs/DMeanDiff
    #the quantiles are
    permutations.bias.lowerQuantile = apply(permutations.bias, MARGIN = 2, FUN = stats::quantile, probs = alpha/2)
    permutations.bias.upperQuantile = apply(permutations.bias, MARGIN = 2, FUN = stats::quantile, probs = 1 - alpha/2)
    #add the quantile lines to the plot
    graphics::lines(permutations.bias.lowerQuantile, 1:K, lty = 2)
    graphics::lines(permutations.bias.upperQuantile, 1:K, lty = 2)
  }

}
