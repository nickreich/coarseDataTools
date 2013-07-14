\name{cd.fit.mcmc}
\alias{cd.fit.mcmc}
\alias{cd.fit.mcmc-class}
\title{An S4 Class that stores a MCMC fit coarse data object}
\description{
  This is the output from \code{dic.fit.mcmc()}, which
  contains the important bits of information about the
  model fit and key options used.
}
\section{Slots}{
  \describe{ \item{\code{ests}:}{Matrix of class
  \code{"numeric"}. This matrix summarizes the results of
  fitting the model. Rows correspond to the first parameter
  , the second paramaeter and then percentiles specified by
  the ptiles argument. Columns correspond to the point
  estimate, the lower and upper bounds on the 95\% credible
  interval and the standard error of the point estimate.}
  \item{\code{conv}:}{Object of class \code{"numeric"}. Not
  used in with \code{dic.fit.mcmc}}
  \item{\code{MSG}:}{Object of class \code{"character"}.
  The error message returned from optim() if the routine
  fails to converge.} \item{\code{loglik}:}{Object of class
  \code{"numeric"}.  Not used in with \code{dic.fit.mcmc}.}
  \item{\code{samples}:}{Object of class
  \code{"data.frame"}. Data frame of posterior draws of
  parameters.} \item{\code{data}:}{Object of class
  \code{"data.frame"}. Original data used to fit model.}
  \item{\code{dist}:}{Object of class \code{"character"}.
  Failure time distribution fit to data. "L" for
  log-normal, "G" for gamma, "W" for Weibull, and "E" for
  Erlang.} \item{\code{inv.hessian}:}{Object of class
  \code{"matrix"}. Not used in with \code{dic.fit.mcmc}.}
  \item{\code{est.method}:}{Object of class
  \code{"character"}. Method used for estimation.}
  \item{\code{ci.method}:}{Object of class
  \code{"character"}. Method used for estimation of
  confidence/credible intervals.} }
}

