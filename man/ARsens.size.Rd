\name{ARsens.size}
\alias{ARsens.size}
\title{Sample Size Calculator for the Power of the Anderson-Rubin (1949) Test with Sensitivity Analysis}
\description{\code{ARsens.size} computes the minimum sample size required for achieving certain power of sensitivity analysis, which is based on an extension of Anderson-Rubin (1949) test and allows IV be possibly invalid within a certain range.}
\usage{
ARsens.size(power, k, beta, gamma, Zadj_sq, sigmau, sigmav, rho, 
            alpha = 0.05, deltarange = deltarange, delta = NULL)
}
\arguments{
  \item{power}{The desired power over a constant.}
  \item{k}{Number of exogenous variables. =}
  \item{beta}{True causal effect minus null hypothesis causal effect.}
  \item{gamma}{Regression coefficient for effect of instruments on treatment.}
  \item{Zadj_sq}{Variance of instruments after regressed on the observed covariates.}
  \item{sigmau}{Standard deviation of potential outcome under control (structural error for y).}
  \item{sigmav}{Standard deviation of error from regressing treatment on instruments.}
  \item{rho}{Correlation between u (potential outcome under control) and v (error from regressing treatment on instruments).}
  \item{alpha}{Significance level.}
  \item{deltarange}{Range of sensitivity allowance. A numeric vector of length 2.}
  \item{delta}{True value of sensitivity parameter when calculating power. Usually take delta = 0 for the favorable situation or delta = NULL for unknown delta.}
}
\value{Minimum sample size required for achieving certain power of sensitivity analysis for the proposed study, which extends the Anderson-Rubin (1949) test with possibly invalid IV. The power formula is derived in Jiang, Small and Zhang (2015).}
\references{
Anderson, T.W. and Rubin, H. (1949), Estimation of the parameters of a single equation in a complete system of stochastic equations, Annals of Mathematical Statistics, 20, 46-63.\cr
Wang, X., Jiang, Y., Small, D. and Zhang, N (2017), Sensitivity analysis and power for instrumental variable studies, (under review of Biometrics).
}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{ivmodel}} for details on the instrumental variables model.
}
\examples{
# Assume we performed a sensitivity analysis in a study with one 
# IV (l=1) and the only exogenous variable is the intercept (k=1). 
# We want to calculate the minimum sample size needed for this 
# sensitivity analysis to have an at least 0.8 power.

# Suppose the difference between the null hypothesis and true causal 
# effect is 1 (beta=1).
# The IV variance is .25 (Zadj_sq =.25).
# The standard deviation of potential outcome is 1(sigmau= 1). 
# The coefficient of regressing IV upon exposure is .5 (gamma= .5).
# The correlation between u and v is assumed to be .5 (rho=.5). 
# The standard deviation of first stage error is .4 (sigmav=.4). 
# The significance level for the study is alpha = .05.

# minimum sample size for sensitivity analysis under the favorable 
# situation, assuming the range of sensitivity allowance is (-0.1, 0.1)
ARsens.size(power=0.8, k=1, beta=1, gamma=.5, Zadj_sq=.25, sigmau=1, 
    sigmav=.4, rho=.5, alpha = 0.05, deltarange=c(-0.1, 0.1), delta=0)

# minimum sample size for sensitivity analysis with unknown delta, 
# assuming the range of sensitivity allowance is (-0.1, 0.1)
ARsens.size(power=0.8, k=1, beta=1, gamma=.5, Zadj_sq=.25, sigmau=1, 
    sigmav=.4, rho=.5, alpha = 0.05, deltarange=c(-0.1, 0.1))

}
\keyword{Anderson-Rubin (1949) test}
\keyword{Sample size}
\keyword{Sensitivity analysis}
\keyword{Power}
