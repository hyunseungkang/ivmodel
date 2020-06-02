#Perform a permutation test for complete randomization
#using the sum of absolute biases as a test statistic.
#Returns a p-value testing whether or not an indicator (exposure or instrument)
#is as-if randomized under complete randomization (i.e., random permutations),
#block randomization (i.e., random permutations within subclasses), or Bernoulli trials.
#Note that the bias is different for the exposure (D) than for the instrument (Z),
#as discussed in Equation (3) of Branson and Keele (2020).
permTest.absBias = function(X, D = NULL, Z = NULL, assignment = "complete", perms = 1000, subclass = NULL){
  if(assignment != "complete" & assignment != "block" & assignment != "bernoulli"){
    return( "Error: assignment must be ''complete'', ''block'', or ''bernoulli''." )
  }
  if(is.null(D) & is.null(Z)){
    return("Error: must enter D or Z indicator.")
  }
  #if only instrument is provided
  if(is.null(D) & !is.null(Z)){
    return("Error: to compute instrument bias, also need D indicator.")
  }
  #if only D is provided
  if(!is.null(D) & is.null(Z)){
    stat.obs = sum(abs(getCovMeanDiffs(X = X, indicator = D)))
  }
  #if D and Z are provided
  if(!is.null(D) & !is.null(Z)){
    #the exposure mean difference is
    DMeanDiff = mean(D[Z == 1]) - mean(D[Z == 0])
    #the observed sum of absolute instrument biases is
    stat.obs = sum(abs(getCovMeanDiffs(X = X, indicator = Z)/DMeanDiff))
  }

  #Now, different tests will be performed depending on the assignment mechanism provided.

  if(assignment == "complete"){
    #permutations for sum of absolute biases across permutations
    permutations.absBias = getCompletePerms.absBias(X, D = D, Z = Z, perms)
  }
  if(assignment == "block"){
    if(is.null(subclass)){ return("Error: Need to provide subclass vector for block randomization.") }
    #permutations for sum of absolute biases across permutations
    permutations.absBias = getBlockPerms.absBias(X, D = D, Z = Z, subclass, perms)
  }
  if(assignment == "bernoulli"){
    #permutations for sum of absolute biases across permutations
    permutations.absBias = getBernoulliPerms.absBias(X, D = D, Z = Z, perms)
  }

  #randomization test pvalue:
  pvalue = mean( permutations.absBias >= stat.obs ) 
  return(pvalue)
}