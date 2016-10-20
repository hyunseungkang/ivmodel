\name{AR.test}
\alias{AR.test}
\title{Anderson-Rubin (1949) Test}
\description{
\code{AR.test} computes the Anderson-Rubin (1949) test for the \code{ivmodel} object as well as the associated confidence interval.
}
\usage{
AR.test(ivmodel, beta0 = 0, alpha = 0.05)
}
\arguments{
  \item{ivmodel}{\code{ivmodel} object}
  \item{beta0}{Null value \eqn{\beta_0} for testing null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}. Default is 0.}
  \item{alpha}{The significance level for hypothesis testing. Default is 0.05.}
}
\value{
  \code{AR.test} returns a list containing the following components
   \item{Fstat}{The value of the test statistic for testing the null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}}
   \item{df}{degree of freedom for the test statistic}
   \item{p.value}{The p value of the test under the null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}}
   \item{ci}{A matrix of two columns, each row contains an interval associated with the confidence interval}
   \item{ci.info}{A human-readable string describing the confidence interval} 
}
\references{Anderson, T.W. and Rubin, H. (1949), Estimation of the parameters of a single equation in a complete system of stochastic equations, Annals of Mathematical Statistics, 20, 46-63.}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small}
\seealso{See also \code{\link{ivmodel}} for details on the instrumental variables model.
}
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
AR.test(foo)
}
\keyword{Anderson-Rubin (1949) test}
