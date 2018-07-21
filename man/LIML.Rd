\name{LIML}
\alias{LIML}
\title{Limited Information Maximum Likelihood Ratio (LIML) Estimator}
\description{\code{LIML} computes the LIML estimate for the \code{ivmodel} object.
}
\usage{
LIML(ivmodel,
     beta0 = 0, alpha = 0.05,
     manyweakSE = FALSE, heteroSE = FALSE,clusterID = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ivmodel}{\code{ivmodel} object.}
  \item{beta0}{Null value \eqn{\beta_0} for testing null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}. Default is 0.}
  \item{alpha}{The significance level for hypothesis testing. Default is
    0.05.}
    \item{manyweakSE}{Should many weak instrument (and
    heteroscedastic-robust) asymptotics in Hansen, Hausman and Newey
    (2008) be used to compute standard errors?}
  \item{heteroSE}{Should heteroscedastic-robust standard errors be used? Default is FALSE.}
  \item{clusterID}{If cluster-robust standard errors are desired, provide a vector of length that's identical to the sample size. For example, if n = 6 and clusterID = c(1,1,1,2,2,2), there would be two clusters where the first cluster is formed by the first three observations and the second cluster is formed by the last three observations. clusterID can be numeric, character, or factor.
}
}
\details{
\code{LIML} computes the LIML estimate for the instrumental variables model in \code{ivmodel}, specifically for the parameter \eqn{beta}. The computation uses \code{KClass} with the value of \eqn{k = k_{LIML}}, which is the smallest root of the equation
\deqn{det(L^T L - k L^T R_Z L) = 0}
where \eqn{L} is a matrix of two columns, the first column consisting of the outcome vector, \eqn{Y}, and the second column consisting of the endogenous variable, \eqn{D}, and \eqn{R_Z = I - Z (Z^T Z)^{-1} Z^T} with \eqn{Z} being the matrix of instruments. \code{LIML} generates a point estimate, a standard error associated with the point estimate, a test statistic and a p value under the null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel} along with a \eqn{1-\alpha} confidence interval.
}
\value{
    \code{LIML} returns a list containing the following components

	\item{k}{The k value for LIML.}
	\item{point.est}{Point estimate of \eqn{\beta}.}
    \item{std.err}{Standard error of the estimate.}
	\item{test.stat}{The value of the test statistic for testing the null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}.}
	\item{p.value}{The p value of the test under the null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}.}
	\item{ci}{A matrix of one row by two columns specifying the confidence interval associated with the Fuller estimator.}
}
\author{Yang Jiang, Hyunseung Kang, Dylan Small}
\seealso{
See also \code{\link{ivmodel}} for details on the instrumental variables model. See also \code{\link{KClass}} for more information about the k-Class estimator.
}
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
LIML(card.model2IV,alpha=0.01)
}
\keyword{Limited information maximum likelihood (LIML) estimator}
\keyword{Instrumental variables}
