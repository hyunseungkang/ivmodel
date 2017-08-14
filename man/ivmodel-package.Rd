\name{ivmodel-package}
\alias{ivmodel-package}
\docType{package}
\title{
\packageTitle{ivmodel}
}
\description{
The package fits an instrumental variables (IV) model of the following type. Let \eqn{Y}, \eqn{D}, \eqn{X}, and \eqn{Z} represent the outcome, endogenous variable, \emph{p} dimensional exogenous covariates, and \emph{L} dimensional instruments, respectively; note that the intercept can be considered as a vector of ones and a part of the exogenous covariates \eqn{X}.The package assumes the following IV model
\deqn{ Y = X \alpha + D \beta + \epsilon, E(\epsilon | X, Z) = 0}
It carries out several IV regressions, diagnostics, and tests associated with the parameter \eqn{\beta} in the IV model. Also, if there is only one instrument, the package runs a sensitivity analysis discussed in Jiang et al. (2015).

The package is robust to most data formats, including factor and character data, and can handle very large IV models efficiently using a sparse QR decomposition.
}
\details{
Supply the outcome \eqn{Y}, the endogenous variable \eqn{D}, and a data frame and/or matrix of instruments \eqn{Z}, and a data frame and/or matrix of exogenous covariates \eqn{X} (optional) and run \code{ivmodel}. Alternatively, one can supply a formula. \code{ivmodel} will generate all the relevant statistics for the parameter \eqn{\beta}.

The DESCRIPTION file:
\packageDESCRIPTION{ivmodel}
\packageIndices{ivmodel}
}
\author{
\packageAuthor{ivmodel}

Maintainer: \packageMaintainer{ivmodel}
}
\references{
Anderson, T. W. and Rubin, H. (1949). Estimation of the parameters of a single equation in a complete system of stochastic equations. \emph{Annals of Mathematical Statistics} 20, 46-63.\cr
  
  Andrews, D. W. K., Moreira, M. J., and Stock, J. H. (2006). Optimal two-side invariant similar tests for instrumental variables regression. \emph{Econometrica} 74, 715-752.\cr
  
  Card, D.  Using Geographic Variation in College Proximity to Estimate the Return to Schooling. In Aspects of Labor Market Behavior: Essays in Honor of John Vanderkamp, eds. L.N. Christophides, E.K. Grant and R. Swidinsky.  201-222. 
National Longitudinal Survey of Young Men: https://www.nlsinfo.org/investigator/pages/login.jsp
  
  Fuller, W. (1977). Some properties of a modification of the limited information estimator. \emph{Econometrica}, 45, 939-953.  \cr
  
  Moreira, M. J. (2003). A conditional likelihood ratio test for structural models. \emph{Econometrica} 71, 1027-1048.   \cr
  
  Wang, X., Jiang, Y., Small, D. and Zhang, N (2017), Sensitivity analysis and power for instrumental variable studies, (under review of Biometrics). \cr
  
  Sargan, J. D. (1958). The estimation of economic relationships using instrumental variables. \emph{Econometrica} , 393-415. 
}

\keyword{ package }
\keyword{ Instrumental variables}
\examples{
data(card.data)
# One instrument #
Y=card.data[,"lwage"]
D=card.data[,"educ"]
Z=card.data[,"nearc4"]
Xname=c("exper", "expersq", "black", "south", "smsa", "reg661", 
        "reg662", "reg663", "reg664", "reg665", "reg666", "reg667", 
		"reg668", "smsa66")
X=card.data[,Xname]
card.model1IV = ivmodel(Y=Y,D=D,Z=Z,X=X)
card.model1IV

# Multiple instruments
Z = card.data[,c("nearc4","nearc2")]
card.model2IV = ivmodel(Y=Y,D=D,Z=Z,X=X)
card.model2IV
}
