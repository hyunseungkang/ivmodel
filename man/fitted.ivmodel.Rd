\name{fitted.ivmodel}
\alias{fitted.ivmodel}
\title{Extract Model Fitted values in the \code{ivmodel} Object}
\description{This fitted method returns the fitted values from k-Class estimators inside \code{ivmodel}.}
\usage{
\method{fitted}{ivmodel}(object,...)
}
\arguments{
  \item{object}{\code{ivmodel} object.}
  \item{\dots}{Additional arguments to \code{fitted}.}
}
\value{
A matrix of fitted values from the k-Class estimations.
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
fitted(foo)
}
