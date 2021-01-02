coefOther <-function(ivmodel){
  if(ivmodel$p == 0) stop("There are no additional exogenous covariates!")

  coefList = list(); nameList = c()
  if(!is.null(ivmodel$kClass)){
    for(i in 1:length(ivmodel$kClass$k)) {
      temp = matrix(0,length(ivmodel$kClass$point.est.other[i,]),4)
      colnames(temp) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
      temp[,1] = ivmodel$kClass$point.est.other[i,]
      temp[,2] = ivmodel$kClass$std.err.other[i,] 
      temp[,3] = ivmodel$kClass$test.stat.other[i,] 
      temp[,4] = ivmodel$kClass$p.value.other[i,]
      rownames(temp) = colnames(ivmodel$X)
      if(ivmodel$kClass$k[i] == 0) {
        nameList = c(nameList,"OLS")
      } else if(ivmodel$kClass$k[i] == 1) {
        nameList = (c(nameList,"TSLS"))
      } else {
        nameList = (c(nameList,ivmodel$kClass$k[i]))
      }
      coefList[[length(coefList) + 1]] = temp
    }
  }
  if(!is.null(ivmodel$LIML)){
    temp = matrix(0,length(ivmodel$LIML$point.est.other),4)
    colnames(temp) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    rownames(temp) = colnames(ivmodel$X)
    temp[,1] = ivmodel$LIML$point.est.other
    temp[,2] = ivmodel$LIML$std.err.other
    temp[,3] = ivmodel$LIML$test.stat.other 
    temp[,4] = ivmodel$LIML$p.value.other
    coefList[[length(coefList) + 1]] = temp
    nameList = c(nameList,"LIML")
  }
  if(!is.null(ivmodel$Fuller)){
    temp = matrix(0,length(ivmodel$Fuller$point.est.other),4)
    colnames(temp) = c( "Estimate", "Std. Error", "t value", "Pr(>|t|)")
    rownames(temp) = colnames(ivmodel$X)
    temp[,1] = ivmodel$Fuller$point.est.other
    temp[,2] = ivmodel$Fuller$std.err.other
    temp[,3] = ivmodel$Fuller$test.stat.other 
    temp[,4] = ivmodel$Fuller$p.value.other
    coefList[[length(coefList) + 1]] = temp
    nameList = c(nameList,"Fuller")
  }
  names(coefList) = nameList
  return(coefList)
}