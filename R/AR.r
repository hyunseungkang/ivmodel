### AR test and its corresponding power and sample size functions

AR.test=function(ivmodel, beta0=0, alpha=0.05){
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

### test  
  temp=Yadj-beta0*Dadj
  Fstat=c(sum(qr.fitted(ivmodel$ZadjQR, temp)^2))/c(sum(temp^2)-sum(qr.fitted(ivmodel$ZadjQR, temp)^2))*(n-k-l)/l
  p.value=1-pf(Fstat, df1=l, df2=n-k-l)
  
### confidence interval  
  cval=qf(1-alpha, df1=l, df2=n-k-l)*l/(n-k-l)
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

  return(list(Fstat=Fstat, df=c(l, n-k-l), p.value=p.value,
              ci.info=info, ci=ci))
}


AR.power=function(n, k, l, beta, gamma, Zadj_sq, 
                  sigmau, sigmav, rho, alpha=0.05){
  if(length(gamma)!=ncol(as.matrix(Zadj_sq))|length(gamma)!=nrow(as.matrix(Zadj_sq))){
	print("The dimension of Zadj_sq doesn't match gamma")
	stop()
  }
  ncp=beta^2*n*c(t(as.matrix(gamma))%*%as.matrix(Zadj_sq)%*%as.matrix(gamma))/
	    (sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)  

  temp = qf(1-alpha, df1=l, df2=n-k-l)
  power = 1-pf(temp, df1=l, df2=n-k-l, ncp=ncp)

  return(power)
}

AR.size=function(power, k, l, beta, gamma, Zadj_sq, 
                 sigmau, sigmav, rho, alpha=0.05){
  if(length(gamma)!=ncol(as.matrix(Zadj_sq))|length(gamma)!=nrow(as.matrix(Zadj_sq))){
    print("The dimension of Zadj_sq doesn't match gamma")
	stop()
  }
  ncp=beta^2*c(t(as.matrix(gamma))%*%as.matrix(Zadj_sq)%*%as.matrix(gamma))/
	  (sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)  

  oldn<-k+l+1
  state<-1
  while(state){
    temp = qf(1-alpha, df1=l, df2=oldn-k-l)  
    temppower = 1-pf(temp, df1=l, df2=oldn-k-l, ncp=ncp*oldn)  
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
    temp = qf(1-alpha, df1=l, df2=oldn-k-l)  
    temppower = 1-pf(temp, df1=l, df2=oldn-k-l, ncp=ncp*new)  
    if(temppower < power){
	  lower <- new
	}else{
	  upper <- new
	}
  }

  return(upper)
}
