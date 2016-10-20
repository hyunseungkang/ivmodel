\name{ivmodel-internal}
\title{Internal \code{ivmodel} functions}
\alias{invTwobyTwoSymMatrix}
\alias{condPvalue}
\alias{qrrank}
\alias{qrRM}
\alias{quadSolver}
\alias{print.ivmodel}
\description{Internal \code{ivmodel} functions}
\usage{
  invTwobyTwoSymMatrix(M)
  quadSolver(a,b,c)
  condPvalue(m,qT,k,df2,eps = 0.02)
  qrrank(QR,tol=1e-8)
  qrRM(QR)
  \method{print}{ivmodel}(x,...)
}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small}
\details{
  \code{invTwobyTwoSymMatrix} is a helper function to invert two-by-two symmetric matrices. 
  \code{quadSolver} is a helper function to invert to solve quadratic equations.
  \code{condPvalue} is a helper function for the conditional likelihood ratio test in \code{CLR.test}.
  \code{qrrank} and \code{qrRM} are helper functions for the sparse QR package.
  \code{print} is a method to print the ivmodel object.
}
\keyword{internal}