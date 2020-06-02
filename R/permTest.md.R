#Perform a permutation test for complete randomization, block randomization, or Bernoulli trials
#using the Mahalanobis distance as a test statistic.
#Returns a p-value testing whether or not an indicator (exposure or instrument)
#is as-if randomized under complete randomization (i.e., random permutations).
permTest.md = function(X, indicator, assignment = "complete", perms = 1000, subclass = NULL){
  if(assignment != "complete" & assignment != "block" & assignment != "bernoulli"){
    return( "Error: assignment must be ''complete'', ''block'', or ''bernoulli''." )
  }
  #observed Mahalanobis distance
  md.obs = as.numeric(getMD(X, indicator))

  if(assignment == "complete"){
    #permutations for MD across permutations
    permutations.md = getCompletePerms.md(X, indicator, perms)
  }
  if(assignment == "block"){
    if(is.null(subclass)){ return("Error: Need to input subclass vector for block randomization.") }
    #permutations for MD across permutations
    permutations.md = getBlockPerms.md(X, indicator, subclass, perms)
  }
  if(assignment == "bernoulli"){
    #permutations for MD across permutations
    permutations.md = getBernoulliPerms.md(X, indicator, perms)
  }

  #randomization test pvalue:
  pvalue = mean( permutations.md >= md.obs ) 

  return(pvalue)
}