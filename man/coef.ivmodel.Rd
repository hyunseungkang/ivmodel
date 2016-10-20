\name{coef.ivmodel}
\alias{coef.ivmodel}
\title{Coefficients of the Fitted Model in the \code{ivmodel} Object}
\description{This coef methods returns the point estimation, standard error, test statistic and p value for all specified k-Class estimation from an \code{ivmodel} object.}
\usage{
\method{coef}{ivmodel}(object,...)
}
\arguments{
  \item{object}{\code{ivmodel} object.}
  \item{\dots}{Additional arguments to \code{coef}.}
}
\value{
A matrix summarizes all the k-Class estimations.
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
coef(foo)
}
