###  power and sample size function for TSLS estimator

TSLS.power=function(n, beta, rho_ZD, sigmau, sigmaDsq, alpha=0.05){
  return(c(1+pnorm(-qnorm(1-alpha/2)-beta*rho_ZD*sqrt(n*sigmaDsq)/sigmau)-pnorm(qnorm(1-alpha/2)-beta*rho_ZD*sqrt(n*sigmaDsq)/sigmau)))
}

TSLS.size=function(power, beta, rho_ZD, sigmau, sigmaDsq, alpha=0.05){
  return(ceiling(c((qnorm(1-alpha/2)+qnorm(power))^2*sigmau^2/beta^2/rho_ZD^2/sigmaDsq)))
}
