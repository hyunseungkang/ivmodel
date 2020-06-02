#Display the randomization distribution of the Mahalanobis distance
#across the exposure and/or instrument for different assignment mechanisms.
#Currently, this function supports the following assignment mechanisms:
#Complete randomization (displayed in black),
#block randomization (displayed in green),
#and Bernoulli trials for exposure (displayed in red) and instrument (displayed in blue)
#This function is used to create Figure 4 of Branson and Keele (2020).
distributionBalancePlot = function(X,
  D = NULL, Z = NULL, subclass = NULL,
  complete = FALSE, blocked = FALSE, bernoulli = FALSE, perms = 1000){
  if(is.null(D) | is.null(Z)){
    return("Error: Must provide D and Z indicators.")
  }

  #vector of MDs across different assignment mechanisms
  compPerms = vector()
  blockPerms = vector()
  bernPerms.exposure = vector()
  bernPerms.instrument = vector()
  #get the xlim and ylim for the plot
  plot.x.vec = vector()
  plot.y.vec = vector()
  if(complete){
    compPerms = getCompletePerms.md(X = X, indicator = Z, perms = perms)
    plot.x.vec = append(plot.x.vec, sqrt(compPerms))
    plot.y.vec = append(plot.y.vec, stats::density(sqrt(compPerms))$y)
  }
  if(blocked){
    if(is.null(subclass)){ return("Error: Need to provide subclass vector for block randomization.") }
    blockPerms = getBlockPerms.md(X = X, indicator = Z, subclass = subclass, perms = perms)
    plot.x.vec = append(plot.x.vec, sqrt(blockPerms))
    plot.y.vec = append(plot.y.vec, stats::density(sqrt(blockPerms))$y)
  }
  if(bernoulli){
    bernPerms.instrument = getBernoulliPerms.md(X = X, indicator = Z, perms = perms)
    bernPerms.exposure = getBernoulliPerms.md(X = X, indicator = D, perms = perms)
    plot.x.vec = append(plot.x.vec, sqrt(bernPerms.instrument))
    plot.x.vec = append(plot.x.vec, sqrt(bernPerms.exposure))
    plot.y.vec = append(plot.y.vec, stats::density(sqrt(bernPerms.instrument))$y)
    plot.y.vec = append(plot.y.vec, stats::density(sqrt(bernPerms.exposure))$y)
  }
  #get the xlim and ylim for the plot
  #also consider the observed MDs
  md.obs.instrument = getMD(X, indicator = Z)
  md.obs.exposure = getMD(X, indicator = D)
  plot.x.max = max(c(plot.x.vec, sqrt(md.obs.exposure), sqrt(md.obs.instrument)))
  plot.y.max = max(plot.y.vec)
  if(plot.y.max == -Inf){plot.y.max = 1}
  #make a blank plot
  graphics::plot(0,0, col = "white", xlim = c(0, plot.x.max), ylim = c(0,plot.y.max),
    xlab = "Square-root of Mahalanobis Distance", ylab = "Density", main = "")
  #add various densities
  if(complete){
    graphics::lines(stats::density(sqrt(compPerms)))
  }
  if(blocked){
    graphics::lines(stats::density(sqrt(blockPerms)), col = "green", lty = 2)
  }
  if(bernoulli){
    graphics::lines(stats::density(sqrt(bernPerms.instrument)), col = "blue", lty = 3)
    graphics::lines(stats::density(sqrt(bernPerms.exposure)), col = "red", lty = 4)
  }
  #add the observed MDs
  graphics::abline(v = sqrt(md.obs.instrument), col = "gray", lty = 2)
  graphics::mtext(expression(sqrt({MD[Z]^{obs}})), side = 1, line = 1, at = sqrt(md.obs.instrument))

  graphics::abline(v = sqrt(md.obs.exposure), col = "gray", lty = 2)
  graphics::mtext(expression(sqrt({MD[D]^{obs}})), side = 1, line = 1, at = sqrt(md.obs.exposure))

}