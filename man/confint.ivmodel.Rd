\name{confint.ivmodel}
\alias{confint.ivmodel}
\title{Confidence Intervals for the Fitted Model in \code{ivmodel} Object}
\description{This confint methods returns a matrix of two columns, each row represents a confident interval for different IV approaches, which include k-Class, AR (Anderson and Rubin 1949) and CLR (Moreira 2003) estimations. 
}
\usage{
\method{confint}{ivmodel}(object,parm,level=NULL,...)
}
\arguments{
  \item{object}{\code{ivmodel} object.}
  \item{parm}{Ignored for our code.}
  \item{level}{The confidence level.}
  \item{\dots}{Additional argument(s) for methods.}
}
\value{A matrix, each row represents a confidence interval for different IV approaches.
}
\references{
  Andrews, D. W. K., Moreira, M. J., and Stock, J. H. (2006). Optimal two-side invariant similar tests for instrumental variables regression. \emph{Econometrica} 74, 715-752.\cr
  Moreira, M. J. (2003). A conditional likelihood ratio test for structural models. \emph{Econometrica} 71, 1027-1048.\cr
  Fuller, W. (1977). Some properties of a modification of the limited information estimator. \emph{Econometrica}, 45, 939-953. \cr
  Anderson, T.W. and Rubin, H. (1949), Estimation of the parameters of a single equation in a complete system of stochastic equations, Annals of Mathematical Statistics, 20, 46-63.
}
\author{Yag Jiang, Hyunseung Kang, and Dylan Small}
\seealso{See also \code{\link{ivmodel}} for details on the instrumental variables model.}
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
confint(foo)
}
