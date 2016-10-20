\name{residuals.ivmodel}
\alias{residuals.ivmodel}
\alias{resid.ivmodel}
\title{Residuals from the Fitted Model in the \code{ivmodel} Object}
\description{This function returns the residuals from the k-Class estimators inside the \code{ivmodel} object.}
\usage{
\method{residuals}{ivmodel}(object,...)
\method{resid}{ivmodel}(object,...)
}
\arguments{
  \item{object}{\code{ivmodel} object.}
  \item{\dots}{Additional arguments to \code{residuals} or \code{resid}.}
}
\value{
A matrix of residuals for each k-Class estimator.
}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small}
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
resid(foo)
residuals(foo)
}
