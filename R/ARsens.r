### Sensitivity analysis and its corresponding power and sample size functions


ARsens.test=function(ivmodel, beta0=0, alpha=0.05, deltarange=NULL){
  if(is.null(deltarange))
    return(NULL)

  Yadj = ivmodel$Yadj; Dadj = ivmodel$Dadj; Zadj = ivmodel$Zadj; 
  n=ivmodel$n;  k=ivmodel$p;  l=ivmodel$L

  if(ncol(Yadj)>1){
    print("The outcome variable should be in one dimension!")
	return(NULL)
  }
  if(ncol(Dadj)>1){
    print("The treatment variable should be in one dimension!")
	return(NULL)
  }
  if(l+k>=n){
    print("Too many IVs, AR can't handle!")
    return(NULL)
  }
  
  if(l!=1){
    print("Please input exact one IV for AR sensitivity analysis")
	return(NULL)
  }

  if(is.numeric(deltarange) & length(deltarange)==2){
    ncp=max(deltarange^2)*sum(Zadj^2)
  }else{
    print("Wrong input of the sensitivity range.")
    return(NULL)
  }
  
### test  
  temp=Yadj-beta0*Dadj
  ncFstat=c(sum(qr.fitted(ivmodel$ZadjQR, temp)^2))/c(sum(temp^2)-sum(qr.fitted(ivmodel$ZadjQR, temp)^2))*(n-k-l)/l
  p.value=1-pf(ncFstat, df1=l, df2=n-k-l, ncp=ncp)
  
### confidence interval  
  cval=qf(1-alpha, df1=l, df2=n-k-l, ncp=ncp)*l/(n-k-l)
  coef.beta0sq=cval*sum(Dadj^2)-(cval+1)*sum(qr.fitted(ivmodel$ZadjQR, Dadj)^2)
  coef.beta0=-2*cval*sum(Dadj*Yadj)+2*(cval+1)*sum(Dadj*qr.fitted(ivmodel$ZadjQR, Yadj))
  coef.constant=cval*sum(Yadj^2)-(cval+1)*sum(qr.fitted(ivmodel$ZadjQR, Yadj)^2)
  Delta=coef.beta0^2-4*coef.constant*coef.beta0sq

  ci=matrix(NA, ncol=2)
  colnames(ci)<-c("lower", "upper")

  if(coef.beta0sq==0){
    if(coef.beta0>0){
      info=c("[",-coef.constant/coef.beta0,",Infinity)")
      ci[1,]=c(-coef.constant/coef.beta0, Inf)
    }
    if(coef.beta0<0){
      info=c("(-Infinity,",-coef.constant/coef.beta0,"]");
      ci[1,]=c(-Inf, -coef.constant/coef.beta0)
    }
    if(coef.beta0==0){
      if(coef.constant>=0){
        info="Whole Real Line"
        ci[1,]=c(-Inf, Inf)
      }
      if(coef.constant<0){
        info="Empty Set"
      }
    }
  }
  
  if(coef.beta0sq!=0){
    if(Delta<=0){
      if(coef.beta0sq>0){
        info="Whole Real Line"
        ci[1,]=c(-Inf, Inf)
      }
      if(coef.beta0sq<0){
        info="Empty Set"
      }
    }
    if(Delta>0){
      # Roots of quadratic equation
      root1=(-coef.beta0+sqrt(Delta))/(2*coef.beta0sq)
      root2=(-coef.beta0-sqrt(Delta))/(2*coef.beta0sq)
      upper.root=max(root1,root2)
      lower.root=min(root1,root2)
      if(coef.beta0sq<0){
        info=paste("[",lower.root,",",upper.root,"]")
        ci[1, ]=c(lower.root, upper.root)
      }
      if(coef.beta0sq>0){
        info= paste("(-Infinity,",lower.root,"] union [",upper.root,",Infinity)")
        ci[1, ]=c(-Inf, lower.root)
        ci<-rbind(ci, c(upper.root, Inf))
      }
    }
  }  

  return(list(ncFstat=ncFstat, df=c(l, n-k-l), ncp=ncp, 
              p.value=p.value, ci.info=info, ci=ci, deltarange=deltarange))
}

ARsens.size=function(power, k, beta, gamma, Zadj_sq, sigmau, sigmav, 
                     rho, alpha=0.05, deltarange=deltarange, delta=NULL){
  if(!is.numeric(gamma) | length(gamma)!=1){
	print("Wrong input of gamma or gamma is not one dimension")
	stop()
  }
  
  if(!is.numeric(Zadj_sq) | length(Zadj_sq)!=1){
	print("Wrong input of Zadj_sq or Zadj_sq is not one dimension")
	stop()
  }

  if(is.numeric(deltarange) & length(deltarange)==2){
    deltarange<-sort(deltarange)
    ncp2=max(deltarange^2)*Zadj_sq
	if(is.null(delta)){
	  if((deltarange[1] < -gamma*beta/sigmau) & 
	     (deltarange[2] > -gamma*beta/sigmau)){
	    ncp1=0
	  }else{
	    ncp1=min((beta*gamma+deltarange[1]*sigmau)^2,
		         (beta*gamma+deltarange[2]*sigmau)^2)*
		     Zadj_sq/(sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)
	  }	
	}else{
	  ncp1=(beta*gamma+delta*sigmau)^2*Zadj_sq/
	       (sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)
	}
  }else{
    print("Wrong input of the sensitivity range.")
    stop()
  }

  if(ncp1<=ncp2){
    print("Sensitivity range too large and the power is always going to be less than the significance level; see Wang, Jiang, Small, and Zhang (2018) `Sensitivity analysis...', Proposition 1.")
	stop()
  }else{
    oldn<-k+2
    state<-1
    while(state){
      temp = qf(1-alpha, df1=1, df2=oldn-k-1, ncp=ncp2*oldn)  
      temppower = 1-pf(temp, df1=1, df2=oldn-k-1, ncp=ncp1*oldn)  
      if(temppower < power){
	    oldn <- oldn*2
	  }else{
	    state <- 0
	  }
    } 

    lower <- oldn%/%2
    upper <- oldn
    while((upper-lower)>2){
      new <- (upper+lower)%/%2
      temp = qf(1-alpha, df1=1, df2=oldn-k-1, ncp=ncp2*new)  
      temppower = 1-pf(temp, df1=1, df2=oldn-k-1, ncp=ncp1*new)  
      if(temppower < power){
	    lower <- new
	  }else{
	    upper <- new
	  }
    }
  }
  return(upper)
}

#####  calculate the power of AR sensitivity analysis

ARsens.power=function(n, k, beta, gamma, Zadj_sq, sigmau, sigmav, 
                      rho, alpha=0.05, deltarange=deltarange, delta=NULL){
  if(!is.numeric(gamma) | length(gamma)!=1){
	print("Wrong input of gamma or gamma is not one dimension")
	stop()
  }
  
  if(!is.numeric(Zadj_sq) | length(Zadj_sq)!=1){
	print("Wrong input of Zadj_sq or Zadj_sq is not one dimension")
	stop()
  }

  if(is.numeric(deltarange) & length(deltarange)==2){
    deltarange<-sort(deltarange)
    ncp2=max(deltarange^2)*n*Zadj_sq
	if(is.null(delta)){
	  if((deltarange[1] < -gamma*beta/sigmau) & 
	     (deltarange[2] > -gamma*beta/sigmau)){
	    ncp1=0
	  }else{
	    ncp1=min((beta*gamma+deltarange[1]*sigmau)^2,
		         (beta*gamma+deltarange[2]*sigmau)^2)*
		     n*Zadj_sq/(sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)
	  }	
	}else{
	  ncp1=(beta*gamma+delta*sigmau)^2*n*Zadj_sq/
	       (sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)
	}
  }else{
    print("Wrong input of the sensitivity range.")
    stop()
  }

  temp = qf(1-alpha, df1=1, df2=n-k-1, ncp=ncp2)
  power = 1-pf(temp, df1=1, df2=n-k-1, ncp=ncp1)

  return(power)
}
