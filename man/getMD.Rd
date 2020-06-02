\name{getMD}
\alias{getMD}
\title{Get Mahalanobis Distance}
\description{
\code{getMD} returns the Mahalanobis distance between two groups.
}
\usage{
getMD(X, indicator, covX.inv = NULL)
}
\arguments{
  \item{X}{Covariate matrix (with units as rows and covariates as columns).}
  \item{indicator}{Binary indicator vector (must contain 1 or 0 for each unit). For example, could be a binary treatment or instrument.}
  \item{covX.inv}{Inverse of the covariate covariance matrix. Usually this is left as \code{NULL}, because \code{getMD()} will compute \code{covX.inv} for you. However, if \code{getMD()} is used many times (e.g., as in a permutation test), it can be computationally efficient to specify \code{covX.inv} beforehand.}
}
\value{
 Mahalanobis distance between two groups.
}
\references{
Branson, Z. and Keele, L. (2020). Evaluating a Key Instrumental Variable Assumption Using Randomization Tests. American Journal of Epidemiology. To appear.
}
\author{Zach Branson and Luke Keele
}
\examples{
	#load the data
	data(icu.data)
	#the covariate matrix is
	X = as.matrix(subset(icu.data, select = -c(open_bin, icu_bed)))
	#mahalanobis distance across the treatment
	getMD(X = X, indicator = icu.data$icu_bed)
	#mahalanobis distance across the instrument
	getMD(X = X, indicator = icu.data$open_bin)
}
\keyword{Covariate Mean Differences}
