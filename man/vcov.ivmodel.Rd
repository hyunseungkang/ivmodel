\name{vcov.ivmodel}
\alias{vcov.ivmodel}
\title{Calculate Variance-Covariance Matrix (i.e. Standard Error) for k-Class Estimators in the \code{ivmodel} Object}
\description{This \code{vcov} method returns the variance-covariance matrix, or equivalently the standard errors, for all specified k-Class estimation from an \code{ivmodel} object.}
\usage{
\method{vcov}{ivmodel}(object,...)
}
\arguments{
  \item{object}{\code{ivmodel} object.}
  \item{\dots}{Additional arguments to \code{vcov}.}
}
\value{
A matrix of standard error estimates for each k-Class estimator.
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
vcov(foo)
}