\name{CLR}
\alias{CLR}
\title{Conditional Likelihood Ratio Test}
\description{\code{CLR} computes the conditional likelihood ratio test (Moreira, 2003) for the \code{ivmodel} object as well as the associated confidence interval.}
\usage{
CLR(ivmodel, beta0 = 0, alpha = 0.05)
}
\arguments{
  \item{ivmodel}{\code{ivmodel} object}
  \item{beta0}{Null value \eqn{\beta_0} for testing null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}. Default is 0}
  \item{alpha}{The significance level for hypothesis testing. Default is 0.05}
}
\details{
\code{CLR.test} computes the conditional likelihood ratio test for the instrumental variables model in \code{ivmodel} object, specifically for the parameter \eqn{\beta}. It also computes the \eqn{1 -\alpha} confidence interval associated with it by inverting the test. The test is fully robust to weak instruments (Moreira 2003). We use the approximation suggested in Andrews et al. (2006) to evaluate the p value and the confidence interval.}
\value{
   \code{CLR} returns a list containing the following components
   \item{test.stat}{The value of the test statistic for testing the null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}}
   \item{p.value}{The p value of the test under the null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}}
   \item{ci}{A matrix of two columns, each row contains an interval associated with the confidence interval}
   \item{ci.info}{A human-readable string describing the confidence interval} 
}
\references{
  Andrews, D. W. K., Moreira, M. J., and Stock, J. H. (2006). Optimal two-side invariant similar tests for instrumental variables regression. \emph{Econometrica} 74, 715-752.\cr
  Moreira, M. J. (2003). A conditional likelihood ratio test for structural models. \emph{Econometrica} 71, 1027-1048.
}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small}
\seealso{See also \code{\link{ivmodel}} for details on the instrumental variables model.}
\examples{
data(card.data)
Y=card.data[,"lwage"]
D=card.data[,"educ"]
Z=card.data[,c("nearc4","nearc2")]
Xname=c("exper", "expersq", "black", "south", "smsa", "reg661", 
        "reg662", "reg663", "reg664", "reg665", "reg666", "reg667", 
		"reg668", "smsa66")
X=card.data[,Xname]
card.model2IV = ivmodel(Y=Y,D=D,Z=Z,X=X)
CLR(card.model2IV,alpha=0.01)
}
\keyword{Instrumental variables}
\keyword{Conditional likelihood ratio test}
