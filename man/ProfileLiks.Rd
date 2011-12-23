\name{ProfileLiks}
\alias{ProfileLiks}
\alias{pl.mu}
\alias{pl.sigma}
\title{Profile likelihood helper functions used by dic.fit().}
\description{
Calculates the likeliood of a dataset given parameters of a log-normal distribution.  These functions are not meant to be called directly by the end user but are used by the \code{\link{dic.fit}} function to determine starting values for the log-scale median and log-log-scale dispersion of the distribution.}
\usage{
pl.mu(mu, log.sigma, dat)
pl.sigma(log.sigma, mu, dat)
}
\arguments{
  \item{mu}{
the log of the median of the log-normal distribution.
}
  \item{log.sigma}{
the log-log (i.e. the log taken twice) of the dispersion of the log-normal distribution.
}
  \item{dat}{
a dataset, as described in the help file of \code{\link{dic.fit}}.
}
}
\value{
the negative log-likelihood for the dataset, given the parameters mu and log.sigma.
}
\author{
Nicholas G. Reich
}
\note{Two versions of the function are needed because the function optimize only can optimize on the first parameter specified by the function.  So pl.mu() is used to find the value of mu that maximizes the likelihood and pl.sigma() is used to find the value of the dispersion (on the log-log scale) that maximizes the likelihood. }

\seealso{
\code{\link{dic.fit}}
}
