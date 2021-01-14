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
A matrix of fitted values from the k-Class estimations. Specifically, each column of the matrix represents predicted values of the outcome for each individual based on different estimates of the treatment effect from k-Class estimators. By default, one of the columns of the matrix is the predicted outcome when the treatment effect is estimated by ordinary least squares (OLS). Because OLS is generally biased in instrumental variables settings, the predictions will likely be biased. For consistent estimates, the predictions are estimates of E[Y | D,X]. In other words, they marginalize over the unmeasured confounder U
and estimate the mean outcomes among all individuals with measured confounders X if
they were to be assigned treatment value D. For example, in the Card study, if U represents
the income of the study unit's parents which were not measured and X represents experience in years, the value of fitted for E[Y |
D = 16, X = 4] is what the average log income among individuals who
had 4 years of experience would be if they were assigned 16 years of education.
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
