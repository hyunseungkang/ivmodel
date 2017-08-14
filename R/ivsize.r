### calculate the sample size needed for achieving a certain power in three ways: TSLS, AR, ARsens
### some parameter values are inferred from the ivmodel

IVsize=function(ivmodel, power, alpha=0.05, beta=NULL, type="TSLS", deltarange=NULL, delta=NULL){
  esti = para(ivmodel)
  if(is.null(beta)) 
    beta = esti$beta
  if(is.null(deltarange)) 
    deltarange = ivmodel$deltarange

  if(type=="TSLS"){
     return(TSLS.size(power, beta, cor(ivmodel$Zadj, ivmodel$Dadj), esti$sigmau, var(ivmodel$Dadj), alpha))
  }
  if(type=="AR"){
    return(AR.size(power, ivmodel$p, ivmodel$L, beta, esti$gamma, var(ivmodel$Zadj), esti$sigmau, esti$sigmav, esti$rho, alpha))  
  }

  if(type=="ARsens"){
    return(ARsens.size(power, ivmodel$p, beta, esti$gamma, var(ivmodel$Zadj), esti$sigmau, esti$sigmav, esti$rho, alpha, deltarange, delta))  
  }
  print("Input error.")
  return(NULL)		  
}
