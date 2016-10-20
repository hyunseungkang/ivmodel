\name{ivmodelFormula}
\alias{ivmodelFormula}
\title{Fitting Instrumental Variables (IV) Models}
\description{
\code{ivmodelFormula} fits an instrumental variables (IV) model with one endogenous variable and a continuous outcome. It carries out several IV regressions, diagnostics, and tests associated this IV model. It is robust to most data formats, including factor and character data, and can handle very large IV models efficiently.
}
\usage{
ivmodelFormula(formula, data, subset, na.action,
        beta0=0,alpha=0.05,k=c(0,1), 
        heteroSE = FALSE, clusterID = NULL, 
        deltarange=NULL)
 }                      	
\arguments{
  \item{formula}{a formula describing the model to be fitted. For example, the formula \code{Y ~ D + X1 + X2 | Z1 + Z2 + X1 + X2} describes the mode where
  	             \deqn{ Y = \alpha_0 + D \beta +  X_{1} \alpha_1 + X_{2} \alpha_2 + \epsilon}
  	             and 
  	             \deqn{ D = \gamma_0 + Z_{1} \gamma_1 + Z_{2} \gamma2 + X_{1} \kappa_1 + X_{2} \kappa_2 + \xi}
  	             The outcome is \code{Y}, the endogenous variable is \code{D}. The exogenous covariates are \code{X1} and \code{X2}.
  	             The instruments are \code{Z1} and \code{Z2}. The formula environment follows the formula environment 
  	             in the \code{ivreg} function in the AER package.}
  \item{data}{an optional data frame containing the variables in the model.
              By default the variables are taken from the environment which \code{ivmodel} is called from}	             
  \item{subset}{an index vector indicating which rows should be used.}
  \item{na.action}{A function to specify the action to be taken if NAsare found. }
  \item{beta0}{Null value \eqn{\beta_0} for testing null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}. Default is $0$.}
  \item{alpha}{The significance level for hypothesis testing. Default is 0.05.}
  \item{k}{A numeric vector of k values for k-class estimation. Default is 0 (OLS) and 1 (TSLS).}
  \item{heteroSE}{Should heteroscedastic-robust standard errors be used? Default is FALSE.}
  \item{clusterID}{If cluster-robust standard errors are desired, provide a vector of 
  	               length that's identical to the sample size. For example, if n = 6 and clusterID = c(1,1,1,2,2,2), 
  	               there would be two clusters where the first cluster is formed by the first three observations and 
  	               the second cluster is formed by the last three observations. clusterID can be numeric, character, or factor.}
  \item{deltarange}{Range of \eqn{\delta} for sensitivity analysis with the Anderson-Rubin (1949) test.}
}
\details{
Let \eqn{Y}, \eqn{D}, \eqn{X}, and \eqn{Z} represent the outcome, endogenous variable, \emph{p} dimensional exogenous covariates, and \emph{L} dimensional instruments, respectively; note that the intercept can be considered as a vector of ones and a part of the exogenous covariates \eqn{X}. \code{ivmodel} assumes the following IV model
\deqn{ Y = X \alpha + D \beta + \epsilon, E(\epsilon | X, Z) = 0}
and produces statistics for \eqn{\beta}. In particular, \code{ivmodel} computes the OLS, TSLS, k-class, limited information maximum likelihood (LIML), and Fuller-k  (Fuller 1977) estimates of \eqn{\beta} using \code{KClass}, \code{LIML}, and code{Fuller}. Also, \code{ivmodel} computes confidence intervals and hypothesis tests of the type \eqn{H_0: \beta = \beta_0} versus \eqn{H_0: \beta \neq \beta_0} for the said estimators as well as two weak-IV confidence intervals, Anderson and Rubin (Anderson and Rubin 1949) confidence interval (Anderson and Rubin 1949) and the conditional likelihood ratio confidence interval (Moreira 2003). Finally, the code also conducts a sensitivity analysis if \eqn{Z} is one-dimensional (i.e. there is only one instrument) using the method in Jiang et al. (2015).

Some procedures (e.g. conditional likelihood ratio test, sensitivity analysis with Anderson-Rubin) assume an additional linear model
\deqn{D = Z \gamma + X \kappa + \xi, E(\xi | X, Z) = 0}
}
\value{
  \code{ivmodel} returns an object of class "ivmodel".
  
  An object class "ivmodel" is a list containing the following components
  \item{alpha}{Significance level for the hypothesis tests.}
  \item{beta0}{Null value of the hypothesis tests.}
  \item{kClass}{A list from \code{KClass} function.}
  \item{LIML}{A list from \code{LIML} function.}
  \item{Fuller}{A list from \code{Fuller} function.}
  \item{AR}{A list from \code{AR.test}.}
  \item{CLR}{A list from \code{CLR}.}
  In addition, if there is only one instrument, \code{ivreg} will generate an "ARsens" list within "ivmodel" object.
}
\references{
  Anderson, T. W. and Rubin, H. (1949). Estimation of the parameters of a single equation in a complete system of stochastic equations. \emph{Annals of Mathematical Statistics} 20, 46-63.\cr
  
  Freeman G, Cowling BJ, Schooling CM (2013). Power and Sample Size Calculations for Mendelian Randomization Studies Using One Genetic Instrument. International journal of epidemiology, 42(4), 1157-1163. \cr
  
  Fuller, W. (1977). Some properties of a modification of the limited information estimator. \emph{Econometrica}, 45, 939-953.  \cr
  
  Jiang, Y., Zhang, N and Small, D. (2015), Sensitivity analysis and power for instrumental variable studies, (under review of Biometrics). \cr
  
  Moreira, M. J. (2003). A conditional likelihood ratio test for structural models. \emph{Econometrica} 71, 1027-1048. \cr

  Sargan, J. D. (1958). The estimation of economic relationships using instrumental variables. Econometrica , 393-415. 
  
}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small}
\seealso{
See also \code{\link{KClass}}, \code{\link{LIML}}, \code{\link{Fuller}}, \code{\link{AR.test}}, and \code{\link{CLR}} for individual methods associated with \code{ivmodel}. For sensitivity analysis with the AR test, 
see \code{\link{ARsens.test}}. \code{ivmodel} has \code{\link{summary.ivmodel}}, \code{\link{confint.ivmodel}}, \code{\link{fitted.ivmodel}}, 
          \code{\link{residuals.ivmodel}} and \code{\link{coef.ivmodel}} methods associated with it.
}
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
card.model1IV = ivmodelFormula(lwage ~ educ + exper + expersq + black + 
                                south + smsa + reg661 + 
                                reg662 + reg663 + reg664 + 
                                reg665 + reg666 + reg667 + 
                                reg668 + smsa66 | nearc4 + 
                                exper + expersq + black + 
                                south + smsa + reg661 + 
                                reg662 + reg663 + reg664 + 
                                reg665 + reg666 + reg667 + 
                                reg668 + smsa66,data=card.data)
card.model1IV

# Multiple instruments
Z = card.data[,c("nearc4","nearc2")]
card.model2IV = ivmodelFormula(lwage ~ educ + exper + expersq + black + 
                                south + smsa + reg661 + 
                                reg662 + reg663 + reg664 + 
                                reg665 + reg666 + reg667 + 
                                reg668 + smsa66 | nearc4 + nearc2 +
                                exper + expersq + black + 
                                south + smsa + reg661 + 
                                reg662 + reg663 + reg664 + 
                                reg665 + reg666 + reg667 + 
                                reg668 + smsa66,data=card.data)
card.model2IV
}
\keyword{Instrumental variables}
\keyword{kClass estimation}
\keyword{Anderson-Rubin test}
\keyword{Sensitivity analysis with Anderson-Rubin test}
\keyword{Conditional likelihood ratio test}
