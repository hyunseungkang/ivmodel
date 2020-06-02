\name{getCovMeanDiffs}
\alias{getCovMeanDiffs}
\title{Get Covariate Mean Differences}
\description{
\code{getCovMeanDiffs} returns the covariate mean differences between two groups.
}
\usage{
getCovMeanDiffs(X, indicator)
}
\arguments{
  \item{X}{Covariate matrix (with units as rows and covariates as columns).}
  \item{indicator}{Binary indicator vector (must contain 1 or 0 for each unit). For example, could be a binary treatment or instrument.}
}
\value{
 Covariate mean differences between two groups.
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
	#covariate mean differences across the treatment
	getCovMeanDiffs(X = X, indicator = icu.data$icu_bed)
	#covariate mean differences across the instrument
	getCovMeanDiffs(X = X, indicator = icu.data$open_bin)
}
\keyword{Covariate Mean Differences}
