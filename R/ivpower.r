### calculate the power for testing H0 when under H_beta in three ways: TSLS, AR, ARsens
### some parameter values are inferred from the ivmodel

IVpower=function(ivmodel, n=NULL, alpha=0.05, beta=NULL, type="TSLS", deltarange=NULL, delta=NULL){
  esti = para(ivmodel)
  if(is.null(n)) 
    n = ivmodel$n
  if(is.null(beta)) 
    beta = esti$beta
  if(is.null(deltarange)) 
    deltarange = ivmodel$deltarange

  if(type=="TSLS"){
     return(TSLS.power(n, beta, cor(ivmodel$Zadj, ivmodel$Dadj), esti$sigmau, var(ivmodel$Dadj), alpha))
  }
  if(type=="AR"){
    return(AR.power(n, ivmodel$p, ivmodel$L, beta, esti$gamma, var(ivmodel$Zadj), esti$sigmau, esti$sigmav, esti$rho, alpha))  
  }

  if(type=="ARsens"){
    return(ARsens.power(n, ivmodel$p, beta, esti$gamma, var(ivmodel$Zadj), esti$sigmau, esti$sigmav, esti$rho, alpha, deltarange, delta))  
  }
  print("Input error.")
  return(NULL)
}