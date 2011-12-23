\name{LikelihoodFuncs}
\alias{LikelihoodFuncs}
\alias{lik}
\alias{diclik}
\alias{diclik2}
\alias{diclik2.helper1}
\alias{diclik2.helper2}
\alias{siclik}
\alias{exactlik}
\title{
Functions to calculate the likelihood for a given, possibly coarse, survival-time observation.
}
\description{
The functions documented here calculate the likelihood for a given survival-time observation.  Observations may be doubly interval-censored, single interval-censored or exact.  Several helper functions are included.  None of these functions are meant to be run by end users, as they are all used by the function \code{\link{dic.fit}}.  
}
\usage{
lik(mu, sigma, EL, ER, SL, SR, type)
diclik(mu, sigma, EL, ER, SL, SR)
diclik2(mu, sigma, EL, ER, SL, SR)
diclik2.helper1(x, SL, SR, mu, sigma)
diclik2.helper2(x, SR, mu, sigma)
siclik(mu, sigma, EL, ER, SL, SR)
exactlik(mu, sigma, EL, ER, SL, SR)
}
\arguments{
  \item{mu}{
the log of the median of the assumed underlying distribution, i.e. the mean of the underlying normal distribution.
}
  \item{sigma}{
the log of the dispersion of the assumed underlying distribution, i.e. the standard deviation of the underlying normal distribution.
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
\item{type}{
a value indicating one of three data types: 0 = doubly interval-censored; 1 = single interval-censored; 2 = exact.
}
\item{x}{
a real number, used as an argument over which to integrate.
}
}

\value{
Each of these functions returns a likelihood (or a portion of a likelihood calculation).  Again, none of these are meant to be called by the end user. 
}
\author{
Nicholas G. Reich
}

\seealso{
\code{\link{dic.fit}}
}
