\name{loglik}
\alias{loglik}
\title{
General likelihood function used in the maximization performed in dic.fit(). 
}
\description{
This function calculates the negative log-likelihood of a dataset of coarse observations and parameters for a log-normal distribution.  This function is not meant to be called by the end user but is rather to be used by \code{\link{dic.fit}} when fitting a log-normal distribution to coarse data.
}
\usage{
loglik(pars, dat)
}
\arguments{
  \item{pars}{
A vector with parameters for a log-normal distribution.  pars[1] is the mean of the underlying normal distribution assumed for the data.  pars[2] is the log of the standard deviation of the underlying normal distribution.
}
  \item{dat}{
A dataset as described in the documentation for \code{\link{dic.fit}}.
}
}
\value{
The negative log likelihood of the entire dataset.  The negative log-likelihood is returned so that the default settings of optim() (i.e. to minimize and not to maximize) can be used to find the maximum likelihood estimate. 
}
\author{
Nicholas G Reich
}

\seealso{
\code{\link{dic.fit}}
}