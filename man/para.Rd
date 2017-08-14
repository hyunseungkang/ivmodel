\name{para}
\alias{para}
\title{Parameter Estimation from Ivmodel
}
\description{\code{para} computes the estimation of several parameters for the \code{ivmodel} object.
}
\usage{
para(ivmodel)
}
\arguments{
  \item{ivmodel}{\code{ivmodel} object.
}
}
\details{\code{para} computes the coefficients of 1st and 2nd stage regression (gamma and beta). It also computes the covariance matrix of the error term of 1st and 2nd stage. (sigmau, sigmav, and rho)
}
\value{
    \code{para} returns a list containing the following components
    
	\item{gamma}{The coefficient of IV in first stage, calculated by linear regression}
	\item{beta}{The TSLS estimator of the exposure effect}
	\item{sigmau}{Standard deviation of potential outcome under control (structural error for y).}
	\item{sigmav}{Standard deviation of error from regressing treatment on instruments}
	\item{rho}{Correlation between u (potential outcome under control) and v (error from regressing treatment on instrument).}

}

\author{Yang Jiang, Hyunseung Kang, Dylan Small
}

\seealso{
See also \code{\link{ivmodel}} for details on the instrumental variables model.
}
\examples{
data(card.data)
Y=card.data[,"lwage"]
D=card.data[,"educ"]
Z=card.data[,"nearc4"]
Xname=c("exper", "expersq", "black", "south", "smsa", "reg661", 
        "reg662", "reg663", "reg664", "reg665", "reg666", "reg667", 
		"reg668", "smsa66")
X=card.data[,Xname]
cardfit=ivmodel(Y=Y, D=D, Z=Z, X=X)
para(cardfit)
}

