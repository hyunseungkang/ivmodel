\name{coefOther}
\alias{coefOther}
\title{Exogenous Coefficients of the Fitted Model in the \code{ivmodel} Object}
\description{This \code{coefOther} returns the point estimates, standard errors, test statistics and p values for the exogenous covariates associated with the outcome. It returns a list of matrices where each matrix is one of the k-Class estimates from an \code{ivmodel} object.}
\usage{
coefOther(ivmodel)
}
\arguments{
  \item{ivmodel}{\code{ivmodel} object.}
}
\value{
A list of matrices swhere each matrix summarizes the estimated coefficients from one of hte k-Class estimates.
}
\author{Hyunseung Kang}
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
coefOther(foo)
}