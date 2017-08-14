##### calculate the TSLS esitmators and covariance matrix estimator

para=function(ivmodel){
  output=list()
  fit1=lm(ivmodel$Dadj~ivmodel$Zadj-1)
  output$gamma=coef(fit1)
  names(output$gamma)=NULL
  output$beta=c(KClass(ivmodel, k=1)$point.est)
  hat_eps=ivmodel$Yadj-c(output$beta)*ivmodel$Dadj
  hat_eta=resid(fit1)
  output$sigmau=sqrt(sum(hat_eps^2)/(ivmodel$n-ivmodel$p))
  output$sigmav=sqrt(sum(hat_eta^2)/(ivmodel$n-ivmodel$p))
  output$rho=sum(hat_eta*hat_eps)/(ivmodel$n-ivmodel$p)/output$sigmau/output$sigmav
  return(output)
}
