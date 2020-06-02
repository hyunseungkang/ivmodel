\name{balanceLovePlot}
\alias{balanceLovePlot}
\title{Create Love plot of standardized covariate mean differences}
\description{
\code{balanceLovePlot} creates a Love plot of the standardized covariate mean differences across the treatment and the instrument. Can also display the permutation quantiles for these quantities. This function is used to create Figure 3a in Branson and Keele (2020).
}
\usage{
balanceLovePlot(X, D, Z, permQuantiles = FALSE, alpha = 0.05, perms = 1000)
}
\arguments{
  \item{X}{Covariate matrix (with units as rows and covariates as columns).}
  \item{D}{Indicator vector for a binary treatment (must contain 1 or 0 for each unit).}
  \item{Z}{Indicator vector for a binary instrument (must contain 1 or 0 for each unit).}
  \item{permQuantiles}{If \code{TRUE}, displays the permutation quantiles for the standardized covariate mean differences.}
  \item{alpha}{The significance level used for the permutation quantiles. For example, if \code{alpha = 0.05}, then the 2.5\% and 97.5\% permutation quantiles are displayed.}
  \item{perms}{Number of permutations used to approximate the permutation quantiles.}
}
\value{
Plot of the standardized covariate mean differences across the treatment and the instrument.
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
	#the treatment
	D = icu.data$icu_bed
	#the instrument
	Z = icu.data$open_bin
	#make the Love plot with permutation quantiles
	\dontrun{balanceLovePlot(X = X, D = D, Z = Z, permQuantiles = TRUE, perms = 500)}
}
\keyword{Balance Plots}
