\name{summary.ivmodel}
\alias{summary.ivmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{summary.ivmodel
%%  ~~function to do ... ~~
}
\description{
Summarise different IV approaches for the \code{ivmodel} object.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{summary}{ivmodel}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{object}{\code{ivmodel} object
%%     ~~Describe \code{ivmodel} here~~
}
  \item{\dots}{additional summary arguments}
}
\details{
The summary has five parts. The 1st part recalls the ivmodel expression and reports the total sample size. The 2nd part reports the first stage regression result to check the goodness of instruments. The 3rd part is the Sargan test (Sargan 1958), which checks for possible existence of invalid IVs. If there is only one IV, the 3rd part will not be displayed. The 4th part reports the result of k-Class estimations. The 5th part provides confidence interval from AR (Anderson and Rubin 1949), CLR (Moreira 2003) test, or sensitivity analysis if specified.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Display summary information in plain text.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  Andrews, D. W. K., Moreira, M. J., and Stock, J. H. (2006). Optimal two-side invariant similar tests for instrumental variables regression. \emph{Econometrica} 74, 715-752.\cr
  Moreira, M. J. (2003). A conditional likelihood ratio test for structural models. \emph{Econometrica} 71, 1027-1048. \cr
  Fuller, W. (1977). Some properties of a modification of the limited information estimator. \emph{Econometrica}, 45, 939-953.  \cr
  Sargan, J. D. (1958). The estimation of economic relationships using instrumental variables. \emph{Econometrica} , 393-415. \cr
  Jiang, Y., Zhang, N and Small, D. (2015), Sensitivity analysis and power for instrumental variable studies, (under review of Biometrics). \cr
  Anderson, T.W. and Rubin, H. (1949), Estimation of the parameters of a single equation in a complete system of stochastic equations, Annals of Mathematical Statistics, 20, 46-63.
%% ~put references to the literature/web site here ~
}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{ivmodel}} for details on the instrumental variables model.
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
summary(foo)
}

