\name{AR.size}
\alias{AR.size}
\title{Sample Size Calculator for the Power of the Anderson-Rubin (1949) Test}
\description{\code{AR.size} computes the minimum sample size required for achieving certain power of Anderson-Rubin (1949) test.
}
\usage{
AR.size(power, k, l, beta, gamma, Zadj_sq, 
        sigmau, sigmav, rho, alpha = 0.05)
}
\arguments{
  \item{power}{The desired power over a constant.}
  \item{k}{Number of exogenous variables.}
  \item{l}{Number of instrumental variables.}
  \item{beta}{True causal effect minus null hypothesis causal effect.}
  \item{gamma}{Regression coefficient for the effect of instrument on treatment.}
  \item{Zadj_sq}{Variance of instruments after regressed on the observed variables.}
  \item{sigmau}{Standard deviation of potential outcome under control (structural error for y).}
  \item{sigmav}{Standard deviation of error from regressing treatment on instruments}
  \item{rho}{Correlation between u (potential outcome under control) and v (error from regressing treatment on instrument).}
  \item{alpha}{Significance level.}
}
\value{Minimum sample size required for achieving certain power of Anderson-Rubin (1949) test.
}
\references{Anderson, T.W. and Rubin, H. (1949), Estimation of the parameters of a single equation in a complete system of stochastic equations, Annals of Mathematical Statistics, 20, 46-63.
}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small}
\seealso{
See also \code{\link{ivmodel}} for details on the instrumental variables model.
}
\examples{
# Assume we performed an AR test in a study with one IV (l=1) and the 
# only one exogenous variable is the intercept (k=1). We want to know 
# the minimum sample size for this test to have an at least 0.8 power.

# Suppose the difference between the null hypothesis and true causal 
# effect is 1 (beta=1).
# The IV variance is .25 (Zadj_sq =.25).
# The standard deviation of potential outcome is 1(sigmau= 1). 
# The coefficient of regressing IV upon exposure is .5 (gamma= .5).
# The correlation between u and v is assumed to be .5 (rho=.5). 
# The standard deviation of first stage error is .4 (sigmav=.4). 
# The significance level for the study is alpha = .05.

# minimum sample size required for Anderson-Rubin test:
AR.size(power=0.8, k=1, l=1, beta=1, gamma=.5, Zadj_sq=.25, 
        sigmau=1, sigmav=.4, rho=.5, alpha = 0.05)

}
\keyword{Anderson-Rubin (1949) test}
\keyword{Power}
\keyword{Sample size}
