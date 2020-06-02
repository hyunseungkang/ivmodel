\name{distributionBalancePlot}
\alias{distributionBalancePlot}
\title{Plot randomization distributions of the Mahalanobis distance}
\description{
\code{distributionBalancePlot} displays the randomization distribution of the square root of the Mahalanobis distance across the treatment and/or instrument for different assignment mechanisms. This function supports complete randomization (displayed in black), block randomization (displayed in green), and Bernoulli trials for exposure (displayed in red) and instrument (displayed in blue). This function is used to create Figure 4 of Branson and Keele (2020).
}
\usage{
distributionBalancePlot(X, D = NULL, Z = NULL, subclass = NULL,
complete = FALSE, blocked = FALSE, bernoulli = FALSE, perms = 1000)
}
\arguments{
  \item{X}{Covariate matrix (with units as rows and covariates as columns).}
  \item{D}{Indicator vector for a binary treatment (must contain 1 or 0 for each unit).}
  \item{Z}{Indicator vector for a binary instrument (must contain 1 or 0 for each unit).}
  \item{subclass}{Vector of subclasses (one for each unit). Subclasses can be numbers or characters, as long as there is one specified for each unit. Only needed if \code{blocked = TRUE}.}
  \item{complete}{If \code{TRUE}, displays the randomization distribution of the Mahalanobis distance under complete randomization.}
  \item{blocked}{If \code{TRUE}, displays the randomization distribution of the Mahalanobis distance under block randomization. Needs \code{subclass} specified.}
  \item{bernoulli}{If \code{TRUE}, displays the randomization distribution of the Mahalanobis distance under Bernoulli trials for the treatment and for the instrument.}
  \item{perms}{Number of permutations used to approximate the randomization distributions.}
}
\value{
Plot of randomization distributions of the square root of the Mahalanobis distance across the treatment and/or instrument for different assignment mechanisms.
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
  #make distribution plot of sqrt(MD) for
  #complete randomization, block randomization, and bernoulli trials
  #(just uncomment the code below)
  #distributionBalancePlot(X = X, D = D, Z = Z, subclass = subclass,
  #complete = TRUE, blocked = TRUE, bernoulli = TRUE, perms = 500)
}
\keyword{Balance Plots}
