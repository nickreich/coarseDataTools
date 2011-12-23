\name{dic.getSE}
\alias{dic.getSE}
\title{
Calculates the standard errors for estimates from dic.fit() using the delta method.}
\description{
This function is used at the end of the \code{\link{dic.fit}} routine to calculate the standard errors of point estimates for the desired characteristics of the log-normal distribution.  It is not meant to be run by end users.
}
\usage{
dic.getSE(mu, log.s, Sig, ptiles)
}
\arguments{
  \item{mu}{
the log of the median of the assumed underlying distribution, i.e. the mean of the underlying normal distribution.
}
  \item{log.s}{
the log-log of the dispersion of the assumed underlying distribution, i.e. the log of the standard deviation of the underlying normal distribution.
} 
  \item{Sig }{
The inverse of the hessian matrix for the likelihood surface at the MLE. 
}
  \item{ptiles}{
percentiles of the distribution to be estimated in this routine.  Note that the median and dispersion are estimated by default.  See the ests matrix in the returned list below for more information.
}
}
\value{
A vector of standard errors for (in this order) the 50th percentile, the dispersion, and any percentile specified in ptiles.
}
\author{
Nicholas G. Reich
}


\seealso{
\code{\link{dic.fit}}
}
