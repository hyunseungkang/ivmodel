\name{TSLS.power}
\alias{TSLS.power}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Power of TSLS Estimator
%%  ~~function to do ... ~~
}
\description{
\code{TSLS.power} computes the power of the asymptotic t-test of TSLS estimator.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TSLS.power(n, beta, rho_ZD, sigmau, sigmaDsq, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Sample size.
%%     ~~Describe \code{n} here~~
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
The power formula is given in Freeman (2013). 
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Power of the asymptotic t-test of TSLS estimator basd on given values of parameters.
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
# Assume we calculate the power of asymptotic t-test of TSLS estimator
# in a study with one IV (l=1) and the only one exogenous variable is 
# the intercept (k=1). 

# Suppose the difference between the null hypothesis and true causal
# effect is 1 (beta=1).
# The sample size is 250 (n=250). 
# The correlation between the IV and exposure is .5 (rho_ZD= .5).
# The standard deviation of potential outcome is 1(sigmau= 1). 
# The variance of the exposure is 1 (sigmaDsq=1).
# The significance level for the study is alpha = .05.

# power of asymptotic t-test of TSLS estimator
TSLS.power(n=250, beta=1, rho_ZD=.5, sigmau=1, sigmaDsq=1, alpha = 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Power }

