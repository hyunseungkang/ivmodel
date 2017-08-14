\name{ARsens.test}
\alias{ARsens.test}
\title{Sensitivity Analysis for the Anderson-Rubin (1949) Test}
\description{\code{ARsens.test} computes sensitivity analysis with possibly invalid instruments, which is an extension of the Anderson-Rubin (1949) test. The formula for sensitivity analysis is derived in Jiang, Small and Zhang (2015).}
\usage{
ARsens.test(ivmodel, beta0 = 0, alpha = 0.05, deltarange = NULL)
}
\arguments{
  \item{ivmodel}{\code{ivmodel} object.}
  \item{beta0}{Null value \eqn{\beta_0} for testing null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}}
  \item{alpha}{The significance level for hypothesis testing. Default is 0.05.}
  \item{deltarange}{Range of sensitivity allowance. A numeric vector of length 2.}
}
\value{
   \code{ARsens.test} returns a list containing the following components
   \item{ncFstat}{The value of the test statistic for testing the null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}}
   \item{df}{degree of freedom for the test statistic}
   \item{ncp}{non-central parameter for the test statistic}
   \item{p.value}{The p value of the test under the null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}}
   \item{ci}{A matrix of two columns, each row contains an interval associated with the confidence interval}
   \item{ci.info}{A human-readable string describing the confidence interval} 
   \item{deltarange}{The inputted range of sensitivity allowance.}
}
\references{Anderson, T.W. and Rubin, H. (1949), Estimation of the parameters of a single equation in a complete system of stochastic equations, Annals of Mathematical Statistics, 20, 46-63.\cr
Wang, X., Jiang, Y., Small, D. and Zhang, N. (2017), Sensitivity analysis and power for instrumental variable studies, (under review of Biometrics).
}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small}
\seealso{See also \code{\link{ivmodel}} for details on the instrumental variables model.}
\examples{
data(card.data)
Y=card.data[,"lwage"]
D=card.data[,"educ"]
Z=card.data[,"nearc4"]
Xname=c("exper", "expersq", "black", "south", "smsa", "reg661", 
        "reg662", "reg663", "reg664", "reg665", "reg666", "reg667", 
		"reg668", "smsa66")
X=card.data[,Xname]
foo = ivmodel(Y=Y,D=D,Z=Z,X=X)
ARsens.test(foo, deltarange=c(-0.03, 0.03))
}
\keyword{Anderson-Rubin (1949) test}
\keyword{Sensitivity analysis}
