#Internal function that permutes the treatment indicator according to biased-coin randomization.
#This simply flips N-many biased coins to get a new treatment indicator;
#it ensures that treatment is not all 1s or 0s
permuteData.biasedCoin = function(N, probs){

  permutation = stats::rbinom(n = N, size = 1, prob = probs)

  #ensure that the treatment indicator is not all 1s or 0s
  while( sum(permutation) == 0 | sum(permutation) == N ){
    permutation = stats::rbinom(n = N, size = 1, prob = probs )
  }

  return(permutation)
}

#Internal function for permuting an indicator
#(instrument or exposure) within a subclass.
#This function takes a table(subclass, indicator) object.
getBlockPerm = function(subclassIndicatorTable){
  #number of subclasses
  Ns = nrow(subclassIndicatorTable)
  #total number of units
  N = sum(subclassIndicatorTable)
  #for each subclass, create a permuted indicator,
  #according to the number of 1s and 0s in that subclass.
  permutedIndicator = vector()
  for(s in 1:Ns){
    permutedIndicator.s = sample(c(rep(0, subclassIndicatorTable[s,"0"]), rep(1, subclassIndicatorTable[s,"1"])))
    permutedIndicator = append(permutedIndicator, permutedIndicator.s)
  }
  return(permutedIndicator)
}

#Internal function that returns the covariate mean differences
#across many permutations of an indicator.
getCompletePerms.meanDiffs = function(X, indicator, perms = 1000){
  #number of covariates
  K = ncol(X)

  #observed standardized covariate mean differences
  covMeanDiff.obs = getCovMeanDiffs(X, indicator)

  #permutations for the randomization test
  indicator.permutations = replicate(perms, sample(indicator), simplify = FALSE)
  #compute the vector of covariate mean differences for each permutation
  permutations.covMeanDiffs = matrix(nrow = perms, ncol = K)
  for(i in 1:perms){
    permutations.covMeanDiffs[i,] = getCovMeanDiffs(X, indicator.permutations[[i]])
  }
  return(permutations.covMeanDiffs)
}

#Internal function that returns the standardized covariate mean differences
#across many permutations of an indicator.
getCompletePerms.balance = function(X, indicator, perms = 1000){
  #number of covariates
  K = ncol(X)

  #observed standardized covariate mean differences
  standardizedCovMeanDiff.obs = getStandardizedCovMeanDiffs(X, indicator)

  #permutations for the randomization test
  indicator.permutations = replicate(perms, sample(indicator), simplify = FALSE)
  #compute the vector of covariate mean differences for each permutation
  permutations.standardizedCovMeanDiffs = matrix(nrow = perms, ncol = K)
  for(i in 1:perms){
    permutations.standardizedCovMeanDiffs[i,] = getStandardizedCovMeanDiffs(X, indicator.permutations[[i]])
  }
  return(permutations.standardizedCovMeanDiffs)
}

#Internal function that returns the Mahalanobis distance (MD)
#across many permutations of an indicator.
getCompletePerms.md = function(X, indicator, perms = 1000){
  #to efficiently compute the MD across permutations, it'll be helpful
  #to compute the inverse of the covariate covariance matrix
  #(which doesn't change across permutations)
  covX.inv = solve(as.matrix(stats::cov(X)))

  #permutations for the randomization test
  indicator.permutations = replicate(perms, sample(indicator), simplify = FALSE)
  #compute the vector of covariate mean differences for each permutation
  permutations.md = vector(length = perms)
  for(i in 1:perms){
    permutations.md[i] = getMD(X, indicator.permutations[[i]], covX.inv)
  }
  return(permutations.md)
}

#Internal function that returns the sum of absolute biases
#across many permutations of an indicator.
#Note that the bias is different for the exposure (D) than for the instrument (Z),
#as discussed in Equation (3) of Branson and Keele (2020).
getCompletePerms.absBias = function(X, D = NULL, Z = NULL, perms = 1000){
  if(is.null(D) & is.null(Z)){
    print("Error: Need to provide D or D and Z indicators.")
  }
  #if the exposure is provided, the absolute bias is just the absolute standardized covariate mean differences
  if(!is.null(D) & is.null(Z)){
    absCovMeanDiffs = abs(getCompletePerms.meanDiffs(X = X, indicator = D, perms = perms))
    #take the sum across covariates
    absCovMeanDiffs.sum = rowSums(absCovMeanDiffs)
    return(absCovMeanDiffs.sum)
  }
  #if the instrument is provided, the absolute bias is the standardized covariate mean difference
  #divided by the mean exposure difference
  #Thus, we also need the exposure to be provided.
  if( !is.null(Z) & is.null(D) ){
    return( "Error: to compute instrument bias, also need D indicator.")
  }
  else{
    #the exposure mean difference is
    DMeanDiff = mean(D[Z == 1]) - mean(D[Z == 0])
    #then, the absolute biases are
    absBias = abs(getCompletePerms.meanDiffs(X = X, indicator = Z, perms = perms)/DMeanDiff)
    absBias.sum = rowSums(absBias)
    return(absBias.sum)
  }
}

#Internal function that returns the Mahalanobis distance (MD)
#across many block permutations of an indicator within a subclass.
getBlockPerms.md = function(X, indicator, subclass, perms = 1000){
  #for efficiency purposes, it'll be helpful to order the data by subclass

  #collect covariates, indicator, and subclass into a dataframe
  data = data.frame(X, indicator = indicator, subclass = subclass)

  #order the data by subclass
  data = data[order(subclass),]

  #Now, the new X, indicator, and subclass are
  X = as.matrix(subset(data, select = -c(indicator, subclass)))
  covX.inv = solve(as.matrix(stats::cov(X)))
  indicator = data$indicator
  subclass = data$subclass

  #To efficiently get block permutations,
  #we just need a table of the indicator and subclass
  subclassIndicatorTable = table(subclass, indicator)

  #Then, the set of block permutations is
  indicator.permutations = t(replicate(perms,
    getBlockPerm(subclassIndicatorTable), simplify = TRUE))
  #compute the vector of covariate mean differences for each permutation
  permutations.md = vector(length = perms)
  for(i in 1:perms){
    permutations.md[i] = getMD(X, indicator.permutations[i,], covX.inv)
  }
  return(permutations.md)
}

#Internal function that returns the sum of absolute biases
#across many block permutations of an indicator within a subclass.
#Note that the bias is different for the exposure (D) than for the instrument (Z),
#as discussed in Equation (3) of Branson and Keele (2020).
getBlockPerms.absBias = function(X, D = NULL, Z = NULL, subclass = NULL, perms = 1000){
  if(is.null(subclass)){
    print("Error: Need to provide subclass vector for block randomization.")
  }
  if(is.null(D) & is.null(Z)){
    print("Error: Need to provide D or D and Z indicators.")
  }
  
  #permutations for the randomization test
  if(!is.null(D) & is.null(Z)){
    indicator = D
  }
  else{
    indicator = Z
  }
  #for efficiency purposes, it'll be helpful to order the data by subclass

  #collect covariates, indicator, and subclass into a dataframe
  data = data.frame(X, indicator = indicator, subclass = subclass)

  #order the data by subclass
  data = data[order(subclass),]

  #Now, the new X, indicator, and subclass are
  X = as.matrix(subset(data, select = -c(indicator, subclass)))
  covX.inv = solve(as.matrix(stats::cov(X)))
  indicator = data$indicator
  subclass = data$subclass

  #To efficiently get block permutations,
  #we just need a table of the indicator and subclass
  subclassIndicatorTable = table(subclass, indicator)

  #Then, the set of block permutations is
  indicator.permutations = t(replicate(perms,
    getBlockPerm(subclassIndicatorTable), simplify = TRUE))

  #if the exposure is provided, the absolute bias is just the absolute standardized covariate mean differences
  if(!is.null(D) & is.null(Z)){
      #compute the vector of covariate mean differences for each permutation
      permutations.covMeanDiffs = matrix(nrow = perms, ncol = ncol(X))
      for(i in 1:perms){
        permutations.covMeanDiffs[i,] = getCovMeanDiffs(X, indicator.permutations[i,])
      }
    absCovMeanDiffs = abs(permutations.covMeanDiffs)
    #take the sum across covariates
    absCovMeanDiffs.sum = rowSums(absCovMeanDiffs)
    return(absCovMeanDiffs.sum)
  }
  #if the instrument is provided, the absolute bias is the standardized covariate mean difference
  #divided by the mean exposure difference
  #Thus, we also need the exposure to be provided.
  if( !is.null(Z) & is.null(D) ){
    return( "Error: to compute instrument bias, also need D indicator.")
  }
  else{
    #the exposure mean difference is
    DMeanDiff = mean(D[Z == 1]) - mean(D[Z == 0])
    #then, the absolute biases are
      permutations.covMeanDiffs = matrix(nrow = perms, ncol = ncol(X))
      for(i in 1:perms){
        permutations.covMeanDiffs[i,] = getCovMeanDiffs(X, indicator.permutations[i,])
      }
    absBias = abs(permutations.covMeanDiffs/DMeanDiff)
    absBias.sum = rowSums(absBias)
    return(absBias.sum)
  }
}

#Internal function that returns the Mahalanobis distance (MD)
#across many Bernoulli-trial permutations of an indicator.
getBernoulliPerms.md = function(X, indicator, perms = 1000){
  #to efficiently compute the MD across permutations, it'll be helpful
  #to compute the inverse of the covariate covariance matrix
  #(which doesn't change across permutations)
  covX.inv = solve(as.matrix(stats::cov(X)))

  #observed Mahalanobis distance is
  md.obs = getMD(X, indicator, covX.inv)

  #need to ensure that X is a matrix to compute the propensity scores
  X = as.matrix(X)
  #model for the propensity scores
  psModel = glm(indicator~X, family = "binomial")
  #the propensity scores
  ps = stats::predict(psModel, type = "response")

  #permutations for the randomization test
  indicator.permutations = replicate(perms,
    permuteData.biasedCoin(N = length(indicator), probs = ps),
    simplify = FALSE)
  #compute the vector of covariate mean differences for each permutation
  permutations.md = vector(length = perms)
  for(i in 1:perms){
    permutations.md[i] = getMD(X, indicator.permutations[[i]], covX.inv)
  }
  return(permutations.md)
}

#Internal function that returns the sum of absolute biases
#across many Bernoulli-trial permutations of an indicator.
#Note that the bias is different for the exposure (D) than for the instrument (Z),
#as discussed in Equation (3) of Branson and Keele (2020).
getBernoulliPerms.absBias = function(X, D = NULL, Z = NULL, perms = 1000){
  if(is.null(D) & is.null(Z)){
    print("Error: Need to provide D or D and Z indicators.")
  }
  #need to ensure that X is a matrix to compute the propensity scores
  X = as.matrix(X)
  #model for the propensity scores
  if(!is.null(D) & is.null(Z)){
    psModel = glm(D~X, family = "binomial")
  }
  else{
    psModel = glm(Z~X, family = "binomial")
  }
  #the propensity scores
  ps = stats::predict(psModel, type = "response")

  #permutations for the randomization test
  if(!is.null(D) & is.null(Z)){
    indicator.permutations = replicate(perms,
    permuteData.biasedCoin(N = length(D), probs = ps),
    simplify = FALSE)
  }
  else{
    indicator.permutations = replicate(perms,
    permuteData.biasedCoin(N = length(Z), probs = ps),
    simplify = FALSE)
  }
  #if the exposure is provided, the absolute bias is just the absolute standardized covariate mean differences
  if(!is.null(D) & is.null(Z)){
      #compute the vector of covariate mean differences for each permutation
      permutations.covMeanDiffs = matrix(nrow = perms, ncol = ncol(X))
      for(i in 1:perms){
        permutations.covMeanDiffs[i,] = getCovMeanDiffs(X, indicator.permutations[[i]])
      }
    absCovMeanDiffs = abs(permutations.covMeanDiffs)
    #take the sum across covariates
    absCovMeanDiffs.sum = rowSums(absCovMeanDiffs)
    return(absCovMeanDiffs.sum)
  }
  #if the instrument is provided, the absolute bias is the standardized covariate mean difference
  #divided by the mean exposure difference
  #Thus, we also need the exposure to be provided.
  if( !is.null(Z) & is.null(D) ){
    return( "Error: to compute instrument bias, also need D indicator.")
  }
  else{
    #the exposure mean difference is
    DMeanDiff = mean(D[Z == 1]) - mean(D[Z == 0])
    #then, the absolute biases are
      permutations.covMeanDiffs = matrix(nrow = perms, ncol = ncol(X))
      for(i in 1:perms){
        permutations.covMeanDiffs[i,] = getCovMeanDiffs(X, indicator.permutations[[i]])
      }
    absBias = abs(permutations.covMeanDiffs/DMeanDiff)
    absBias.sum = rowSums(absBias)
    return(absBias.sum)
  }
}