\name{vcovOther}
\alias{vcovOther}
\title{Variance of Exogenous Coefficients of the Fitted Model in the \code{ivmodel} Object}
\description{This \code{vcovOther} returns the estimated variances of the estimated coefficients for the exogenous covariates associated with the outcome. All the estimation is based on k-Class estimators.}
\usage{
  vcovOther(ivmodel)
}
\arguments{
  \item{ivmodel}{\code{ivmodel} object.}
}
\value{
  A matrix where each row represents a k-class estimator and each column represents one of the exogenous covariates. Each element is the estimated variance of the estimated coefficients.
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
  vcovOther(foo)
}