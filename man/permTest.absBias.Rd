\name{permTest.absBias}
\alias{permTest.absBias}
\title{Perform a permutation test using the sum of absolute biases}
\description{
\code{permTest.absBias} performs a permutation test for complete randomization using the sum of absolute biases as a test statistic.
}
\usage{
permTest.absBias(X, D = NULL, Z = NULL,
assignment = "complete",
perms = 1000, subclass = NULL)
}
\arguments{
  \item{X}{Covariate matrix (with units as rows and covariates as columns).}
  \item{D}{Indicator vector for a binary treatment (must contain 1 or 0 for each unit).}
  \item{Z}{Indicator vector for a binary instrument (must contain 1 or 0 for each unit).}
  \item{assignment}{Must be "complete", "block", or "bernoulli". Designates whether to test for complete randomization, block randomization, or Bernoulli trials.}
  \item{subclass}{Vector of subclasses (one for each unit). Subclasses can be numbers or characters, as long as there is one specified for each unit. Only needed if \code{assignment = "block"}.}
  \item{perms}{Number of permutations used to approximate the permutation test.}
}
\value{
p-value testing whether or not an indicator (treatment or instrument) is as-if randomized under complete randomization (i.e., random permutations), block randomization (i.e., random permutations within subclasses), or Bernoulli trials.
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
	#the subclass
	subclass = icu.data$site
	
	#can uncomment the following code for examples
	
	#permutation test for complete randomization (for the treatment)
	#permTest.absBias(X = X, D = D,
	#assignment = "complete", perms = 500)
	#permutation test for complete randomization (for the instrument)
	#permTest.absBias(X = X, D = D, Z = Z,
	#assignment = "complete", perms = 500)
	#permutation test for block randomization (for the treatment)
	#permTest.absBias(X = X, D = D,
	#assignment = "block", subclass = subclass, perms = 500)
	#permutation test for block randomization (for the instrument)
	#permTest.absBias(X = X, D = D, Z = Z,
	#assignment = "block",
	#subclass = subclass, perms = 500)
	#permutation test for bernoulli trials (for the treatment)
	#permTest.absBias(X = X, D = D,
	#assignment = "bernoulli", perms = 500)
	#permutation test for bernoulli randomization (for the instrument)
	#permTest.absBias(X = X, D = D, Z = Z,
	#assignment = "bernoulli", perms = 500)
}
\keyword{Permutation Tests}
