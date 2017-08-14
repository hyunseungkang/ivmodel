\name{IVpower}
\alias{IVpower}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Power calculation for IV models
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{IVpower} computes the power for one of the following tests: two stage least square estimates; Anderson-Rubin (1949) test; Sensitivity analysis.
}
\usage{
IVpower(ivmodel, n = NULL, alpha = 0.05, beta = NULL, type = "TSLS", 
        deltarange = NULL, delta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ivmodel}{\code{ivmodel} object.
%%     ~~Describe \code{ivmodel} here~~
}
  \item{n}{number of sample size, if missing, will use the sample size from the input \code{ivmodel} object.
%%     ~~Describe \code{n} here~~
}
  \item{alpha}{The significance level for hypothesis testing. Default is 0.05.
%%     ~~Describe \code{alpha} here~~
}
  \item{beta}{True causal effect minus null hypothesis causal effect. If missing, will use the beta calculated from the input \code{ivmodel} object.
%%     ~~Describe \code{beta} here~~
}
  \item{type}{Determines which test will be used for power calculation. "TSLS" for two stage least square estimates; "AR" for Anderson-Rubin test; "ARsens" for sensitivity analysis.
%%     ~~Describe \code{type} here~~
}
  \item{deltarange}{Range of sensitivity allowance. A numeric vector of length 2. If missing, will use the deltarange from the input \code{ivmodel} object.
%%     ~~Describe \code{deltarange} here~~
}
  \item{delta}{True value of sensitivity parameter when calculating the power. Usually take delta = 0 for the favorable situation or delta = NULL for unknown delta.
%%     ~~Describe \code{delta} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{IVpower} computes the power for one of the following tests: two stage least square estimates; Anderson-Rubin (1949) test; Sensitivity analysis. The related value of parameters will be inferred from the input of \code{ivmodel} object.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
a power value for the specified type of test.
}
\references{
%% ~put references to the literature/web site here ~
Freeman G, Cowling BJ, Schooling CM (2013). Power and Sample Size Calculations for Mendelian Randomization Studies Using One Genetic Instrument. International journal of epidemiology, 42(4), 1157-1163. \cr
Anderson, T.W. and Rubin, H. (1949).  Estimation of the parameters of a single equation in a complete system of stochastic equations.  Annals of Mathematical Statistics, 20, 46-63. \cr
ang, X., Jiang, Y., Small, D. and Zhang, N (2017), Sensitivity analysis and power for instrumental variable studies, (under review of Biometrics).
}
\author{
%%  ~~who you are~~
Yang Jiang, Hyunseung Kang, Dylan Small
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
See also \code{\link{ivmodel}} for details on the instrumental variables model. See also \code{\link{TSLS.power}}, \code{\link{AR.power}}, \code{\link{ARsens.power}} for details on the power calculation.
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
card.model = ivmodel(Y=Y,D=D,Z=Z,X=X)

IVpower(card.model)
IVpower(card.model, n=10^4, type="AR")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{power}

