\name{TSLS.size}
\alias{TSLS.size}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sample Size Calculator for the Power of Asymptotic T-test
%%  ~~function to do ... ~~
}
\description{\code{TSLS.size} computes the minimum sample size required for achieving certain power of asymptotic t-test of TSLS estimator.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TSLS.size(power, beta, rho_ZD, sigmau, sigmaDsq, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{power}{The desired power over a constant.
%%     ~~Describe \code{power} here~~
}
  \item{beta}{True causal effect minus null hypothesis causal effect.
%%     ~~Describe \code{beta} here~~
}
  \item{rho_ZD}{Correlation between the IV Z and the exposure D.
%%     ~~Describe \code{rho_ZD} here~~
}
  \item{sigmau}{Standard deviation of potential outcome under control. (structural error for y)
%%     ~~Describe \code{sigmau} here~~
}
  \item{sigmaDsq}{The variance of the exposure D.
%%     ~~Describe \code{sigmaDsq} here~~
}
  \item{alpha}{Significance level.
%%     ~~Describe \code{alpha} here~~
}
}
\details{
The calculation is based on inverting the power formula given in Freeman (2013). 
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Minimum sample size required for achieving certain power of asymptotic t-test of TSLS estimator.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Freeman G, Cowling BJ, Schooling CM (2013). Power and Sample Size Calculations for Mendelian Randomization Studies Using One Genetic Instrument. International journal of epidemiology, 42(4), 1157-1163.
%% ~put references to the literature/web site here ~
}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small
%%  ~~who you are~~
}

\seealso{
See also \code{\link{ivmodel}} for details on the instrumental variables model.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Assume we performed an asymptotic t-test of TSLS estimator in a study 
# with one IV (l=1) and the only one exogenous variable is the intercept 
# (k=1). We want to calculate the minimum sample size needed for this 
# test to have an at least 0.8 power.

# Suppose the null hypothesis causal effect is 0 and the true causal 
# effect is 1 (beta=1-0=1).
# The correlation between the IV and exposure is .5 (rho_ZD= .5).
# The standard deviation of potential outcome is 1(sigmau= 1). 
# The variance of the exposure is 1 (sigmaDsq=1).
# The significance level for the study is alpha = .05.

### minimum sample size required for aysmptotic t-test
TSLS.size(power=.8, beta=1, rho_ZD=.5, sigmau=1, sigmaDsq=1, alpha =.05)
}

