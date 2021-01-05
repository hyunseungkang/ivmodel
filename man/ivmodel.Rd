\name{ivmodel}
\alias{ivmodel}
\title{Fitting Instrumental Variables (IV) Models}
\description{
\code{ivmodel} fits an instrumental variables (IV) model with one endogenous variable and a continuous outcome. It carries out several IV regressions, diagnostics, and tests associated this IV model. It is robust to most data formats, including factor and character data, and can handle very large IV models efficiently.
}
\usage{
ivmodel(Y, D, Z, X, intercept = TRUE,
        beta0 = 0, alpha = 0.05, k = c(0, 1),
        manyweakSE = FALSE, heteroSE = FALSE, clusterID = NULL,
        deltarange = NULL, na.action = na.omit)
 }
\arguments{
  \item{Y}{A numeric vector of outcomes.}
  \item{D}{A vector of endogenous variables.}
  \item{Z}{A matrix or data frame of instruments.}
  \item{X}{A matrix or data frame of (exogenous) covariates.}
  \item{intercept}{Should the intercept be included? Default is TRUE and if so, you do not need to add a column of 1s in X.}
  \item{beta0}{Null value \eqn{\beta_0} for testing null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}. Default is $0$.}
  \item{alpha}{The significance level for hypothesis testing. Default is 0.05.}
  \item{k}{A numeric vector of k values for k-class estimation. Default
    is 0 (OLS) and 1 (TSLS).}
  \item{manyweakSE}{Should many weak instrument (and
    heteroscedastic-robust) asymptotics in Hansen, Hausman and Newey
    (2008) be used to compute standard errors? (Not supported for k ==0)}
  \item{heteroSE}{Should heteroscedastic-robust standard errors be used? Default is FALSE.}
  \item{clusterID}{If cluster-robust standard errors are desired, provide a vector of
  	               length that's identical to the sample size. For example, if n = 6 and clusterID = c(1,1,1,2,2,2),
  	               there would be two clusters where the first cluster is formed by the first three observations and
  	               the second cluster is formed by the last three observations. clusterID can be numeric, character, or factor.}
  \item{deltarange}{Range of \eqn{\delta} for sensitivity analysis with the Anderson-Rubin (1949) test.}
  \item{na.action}{NA handling. There are \code{na.fail}, \code{na.omit}, \code{na.exclude}, \code{na.pass} available. Default is na.omit.}
}
\details{
Let \eqn{Y}, \eqn{D}, \eqn{X}, and \eqn{Z} represent the outcome, endogenous variable, \emph{p} dimensional exogenous covariates, and \emph{L} dimensional instruments, respectively. Note that the intercept is a type of exogenous covariate and can be added to \eqn{X} by specifying \code{intercept} as TRUE (the default behavior); the user does not have to manually add an intercept column in \eqn{X}. \code{ivmodel} assumes the following IV model
\deqn{ Y = X \alpha + D \beta + \epsilon, E(\epsilon | X, Z) = 0}
and produces statistics for \eqn{\beta}. In particular, \code{ivmodel} computes the OLS, TSLS, k-class, limited information maximum likelihood (LIML), and Fuller-k  (Fuller 1977) estimates of \eqn{\beta} using \code{KClass}, \code{LIML}, and code{Fuller}. Also, \code{ivmodel} computes confidence intervals and hypothesis tests of the type \eqn{H_0: \beta = \beta_0} versus \eqn{H_0: \beta \neq \beta_0} for the said estimators as well as two weak-IV confidence intervals, Anderson and Rubin (Anderson and Rubin 1949) confidence interval (Anderson and Rubin 1949) and the conditional likelihood ratio confidence interval (Moreira 2003). Finally, the code also conducts a sensitivity analysis if \eqn{Z} is one-dimensional (i.e. there is only one instrument) using the method in Jiang et al. (2015).

Some procedures (e.g. conditional likelihood ratio test, sensitivity analysis with Anderson-Rubin) assume an additional linear model
\deqn{D = Z \gamma + X \kappa + \xi, E(\xi | X, Z) = 0}
}
\value{
  \code{ivmodel} returns an object of class "ivmodel".

  An object class "ivmodel" is a list containing the following components
  \item{n}{Sample size.}
  \item{L}{Number of instruments.}
  \item{p}{Number of exogenous covariates (including intercept).}
  \item{Y}{Outcome, cleaned for use in future methods.}
  \item{D}{Treatment, cleaned for use in future methods.}
  \item{Z}{Instrument(s), cleaned for use in future methods.}
  \item{X}{Exogenous covariates (if provided), cleaned for use in future methods.}
  \item{Yadj}{Adjusted outcome, projecting out X.}
  \item{Dadj}{Adjusted treatment, projecting out X.}
  \item{Zadj}{Adjusted instrument(s), projecting out X.}
  \item{ZadjQR}{QR decomposition for adjusted instrument(s).}
  \item{ZXQR}{QR decomposition for concatenated matrix of Z and X.}
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

  Freeman G., Cowling B. J., Schooling C. M. (2013). Power and Sample Size Calculations for Mendelian Randomization Studies Using One Genetic Instrument. \emph{International Journal of Epidemiology} 42(4), 1157-1163. \cr

  Fuller, W. (1977). Some properties of a modification of the limited information estimator. \emph{Econometrica}, 45, 939-953.  \cr
  
  Hansen, C., Hausman, J., and Newey, W. (2008) Estimation with many instrumental variables. \emph{Journal of Business \& Economic Statistics} 26(4), 398-422. \cr
  
  Moreira, M. J. (2003). A conditional likelihood ratio test for structural models. \emph{Econometrica} 71, 1027-1048. \cr
  
  Sargan, J. D. (1958). The estimation of economic relationships using instrumental variables. Econometrica , 393-415. \cr

  Wang, X., Jiang, Y., Small, D. and Zhang, N. (2017), Sensitivity analysis and power for instrumental variable studies. \emph{Biometrics} 74(4), 1150-1160.


}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small}
\seealso{
See also \code{\link{KClass}}, \code{\link{LIML}}, \code{\link{Fuller}}, \code{\link{AR.test}}, and \code{\link{CLR}} for individual methods associated with \code{ivmodel}. For extracting the estimated effect of the exogenous covariates on the outcome, see \code{\link{coefOther}}. For sensitivity analysis with the AR test,
see \code{\link{ARsens.test}}. \code{ivmodel} has \code{\link{vcov.ivmodel}},\code{\link{model.matrix.ivmodel}},\code{\link{summary.ivmodel}}, \code{\link{confint.ivmodel}}, \code{\link{fitted.ivmodel}},
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
card.model1IV = ivmodel(Y=Y,D=D,Z=Z,X=X)
card.model1IV

# Multiple instruments
Z = card.data[,c("nearc4","nearc2")]
card.model2IV = ivmodel(Y=Y,D=D,Z=Z,X=X)
card.model2IV
}
\keyword{Instrumental variables}
\keyword{kClass estimation}
\keyword{Anderson-Rubin test}
\keyword{Sensitivity analysis with Anderson-Rubin test}
\keyword{Conditional likelihood ratio test}
