\name{AR.power}
\alias{AR.power}
\title{Power of the Anderson-Rubin (1949) Test}
\description{
\code{AR.power} computes the power of Anderson-Rubin (1949) test based on the given values of parameters.
}
\usage{
AR.power(n, k, l, beta, gamma, Zadj_sq, 
         sigmau, sigmav, rho, alpha = 0.05)
}
\arguments{
  \item{n}{Sample size.}
  \item{k}{Number of exogenous variables.}
  \item{l}{Number of instrumental variables.}
  \item{beta}{True causal effect minus null hypothesis causal effect.}
  \item{gamma}{Regression coefficient for effect of instruments on treatment.}
  \item{Zadj_sq}{Variance of instruments after regressed on the observed variables.}
  \item{sigmau}{Standard deviation of potential outcome under control. (structural error for y)}
  \item{sigmav}{Standard deviation of error from regressing treatment on instruments.}
  \item{rho}{Correlation between u (potential outcome under control) and v (error from regressing treatment on instrument).}
  \item{alpha}{Significance level.}
}
\value{
Power of the Anderson-Rubin test based on the given values of parameters.
}
\references{
Anderson, T.W. and Rubin, H. (1949).  Estimation of the parameters of a single equation in a complete system of stochastic equations.  Annals of Mathematical Statistics, 20, 46-63.
}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small
}

\seealso{
See also \code{\link{ivmodel}} for details on the instrumental variables model.
}
\examples{
# Assume we calculate the power of AR test in a study with one IV (l=1) 
# and the only one exogenous variable is the intercept (k=1). 

# Suppose the difference between the null hypothesis and true causal 
# effect is 1 (beta=1).
# The sample size is 250 (n=250), the IV variance is .25 (Zadj_sq =.25).
# The standard deviation of potential outcome is 1(sigmau= 1). 
# The coefficient of regressing IV upon exposure is .5 (gamma= .5).
# The correlation between u and v is assumed to be .5 (rho=.5). 
# The standard deviation of first stage error is .4 (sigmav=.4). 
# The significance level for the study is alpha = .05.

# power of Anderson-Rubin test:
AR.power(n=250, k=1, l=1, beta=1, gamma=.5, Zadj_sq=.25, 
         sigmau=1, sigmav=.4, rho=.5, alpha = 0.05)

}
\keyword{Anderson-Rubin (1949) test}
\keyword{Power}
