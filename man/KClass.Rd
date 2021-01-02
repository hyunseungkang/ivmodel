\name{KClass}
\alias{KClass}
\title{k-Class Estimator}
\description{\code{KClass} computes the k-Class estimate for the \code{ivmodel} object.}
\usage{
KClass(ivmodel,
       beta0 = 0, alpha = 0.05, k = c(0, 1),
       manyweakSE = FALSE, heteroSE = FALSE,clusterID = NULL)
}
\arguments{
  \item{ivmodel}{\code{ivmodel} object.}
  \item{beta0}{Null value \eqn{\beta_0} for testing null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}. Default is 0.}
  \item{alpha}{The significance level for hypothesis testing. Default is 0.05.}
  \item{k}{A vector of \eqn{k} values for the k-Class estimator. Default
    is 0 (OLS) and 1 (TSLS).}
  \item{manyweakSE}{Should many weak instrument (and
    heteroscedastic-robust) asymptotics in Hansen, Hausman and Newey
    (2008) be used to compute standard errors? (Not supported for k=0)}
  \item{heteroSE}{Should heteroscedastic-robust standard errors be used? Default is FALSE.}
  \item{clusterID}{If cluster-robust standard errors are desired, provide a vector of length that's identical to the sample size. For example, if n = 6 and clusterID = c(1,1,1,2,2,2), there would be two clusters where the first cluster is formed by the first three observations and the second cluster is formed by the last three observations. clusterID can be numeric, character, or factor.
}
}
\details{
\code{KClass} computes the k-Class estimate for the instrumental variables model in \code{ivmodel}, specifically for the parameter \eqn{\beta}. It generates a point estimate, a standard error associated with the point estimate, a test statistic and a p value under the null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel} along with a \eqn{1-\alpha} confidence interval.
}
\value{
    \code{KClass} returns a list containing the following components

	\item{k}{A row matrix of k values supplied to \code{KClass}.}
	\item{point.est}{A row matrix of point estimates of \eqn{\beta}, with each row corresponding to the k values supplied.}
    \item{std.err}{A row matrix of standard errors of the estimates, with each row corresponding to the k values supplied.}
	\item{test.stat}{A row matrix of test statistics for testing the null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}, with each row corresponding to the k values supplied.}
	\item{p.value}{A row matrix of p value of the test under the null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}, with each row corresponding to the k values supplied.}
	\item{ci}{A matrix of two columns specifying the confidence interval, with each row corresponding to the k values supplied.}
}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small}
\seealso{
See also \code{\link{ivmodel}} for details on the instrumental variables model.
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
KClass(card.model2IV,
          k=c(0,1,0.5))

\dontrun{
## The following code tests the mank weak IV standard error for LIML and Fuller.
example <- function(q = 10, rho1 = 0.5, n1 = 10000,
sigma.uv = 0.5, beta = 1, gamma = rep(1/sqrt(q), q)) {

    Sigma1 <- outer(1:q, 1:q, function(i, j) rho1^abs(i - j))

    library(MASS)
    Z1 <- mvrnorm(n1, rep(1, q), Sigma1)
    Z1 <- matrix(2 * as.numeric(Z1 > 0) - 1, nrow = n1)
    UV1 <- mvrnorm(n1, rep(0, 2), matrix(c(1, sigma.uv, sigma.uv, 1), 2))
    X1 <- Z1 %*% gamma + runif(n1, 0.5, 2) * UV1[, 1] + runif(n1, 0.1, 0.5) * rexp(n1)
    Y1 <- X1 %*% beta + runif(n1, 0.5, 2) * UV1[, 2] + runif(n1, 0.1, 0.5) * rexp(n1)

    list(Z1 = Z1, X1 = X1, Y1 = Y1)

}

one.sim <- function(manyweakSE) {
    data <- example(q = 100, n1 = 200)
    fit <- ivmodel(data$Y1, data$X1, data$Z1, manyweakSE = manyweakSE)
    1 > coef(fit)[, 2] - 1.96 * coef(fit)[, 3] & 1 < coef(fit)[, 2] + 1.96 * coef(fit)[, 3]
}

res <- replicate(200, one.sim(TRUE))
apply(res, 1, mean)

res <- replicate(200, one.sim(FALSE))
apply(res, 1, mean)
}

}
\keyword{k-Class estimator}
\keyword{Instrumental variables}
