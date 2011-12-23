\name{fwFuncs}
\alias{fwFuncs}
\alias{fw1}
\alias{fw3}
\title{
Functions to calculate the product of the density and the likelihood weight functions.   
}
\description{
The fw1() and fw3() functions are not meant to be called by the end user but are used by \code{\link{dic.fit}} when calculating the likelihood of a given observation.  For a given time, t, they calculate the likelihood, given a set of parameters (mu and sigma) and the observed endpoints of exposure and symptom onset windows.  Details on the likelihood functions can be found in Reich et al. (2009) -- see reference below.
}
\usage{
fw1(t, EL, ER, SL, SR, mu, sigma)
fw3(t, EL, ER, SL, SR, mu, sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
a real number, representing time.
}
  \item{EL}{
the left endpoint of the exposure window.
}
  \item{ER}{
the right endpoint of the exposure window.
}
  \item{SL}{
the left endpoint of the symptom onset window.
}
  \item{SR}{
the right endpoint of the symptom onset window.

}
  \item{mu}{
the log of the median of the assumed underlying distribution, i.e. the mean of the underlying normal distribution.
}
  \item{sigma}{
the log-log of the dispersion of the assumed underlying distribution, i.e. the log of the standard deviation of the underlying normal distribution.
}
}
\value{
fw1() returns the value for the first third of the integral for the doubly interval-censored likelihood calculation.  fw3() returns the value for the last third of the integral for the doubly interval-censored likelihood calculation.  }
\references{
Reich NG et al.  Statistics in Medicine.  Estimating incubation periods with coarse data. 2009.  \url{http://www3.interscience.wiley.com/journal/122507367/abstract}
}
\author{
Nicholas G. Reich
}