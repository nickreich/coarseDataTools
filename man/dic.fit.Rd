\name{dic.fit}
\alias{dic.fit}
\title{Fits a log-normal, Gamma, Erlang, or Weibull model to doubly interval
censored survival data}
\usage{
  dic.fit(dat, start.par2 = log(2),
    opt.method = "L-BFGS-B",
    par1.int = c(log(0.5), log(13)),
    par2.int = c(log(1.01), log(log(5))),
    ptiles = c(0.05, 0.95, 0.99), dist = "L", n.boots = 0,
    ...)
}
\arguments{
  \item{dat}{a matrix with columns named "EL", "ER", "SL",
  "SR", corresponding to the left (L) and right (R)
  endpoints of the windows of possible exposure (E) and
  symptom onset (S). Also, a "type" column must be
  specified and have entries with 0, 1, or 2, corresponding
  to doubly interval-censored, single interval-censored or
  exact observations, respsectively.}

  \item{start.par2}{starting value for 2nd parameter of
  desired distribtution}

  \item{opt.method}{method used by optim}

  \item{par1.int}{the log-scale interval of possible median
  values (in the same units as the observations in dat).
  Narrowing this interval can help speed up convergence of
  the algorithm, but care must be taken so that possible
  values are not excluded or that the maximization does not
  return a value at an endpoint of this interval.}

  \item{par2.int}{the log-scale interval of possible
  dispersion values}

  \item{ptiles}{percentiles of interest}

  \item{dist}{what distribution to use to fit the data.
  Default "L" for log-normal. "G" for gamma, and "W" for
  Weibull. Note: If dist is Gamma (G) or Weibull (W), the
  mu refers to the shape and sigma refers to the scale
  param.}

  \item{n.boots}{number of bootstrap resamples if non-log
  normal model}

  \item{...}{additional options passed to optim}

  \item{start.sigma}{the log-scale starting value for the
  dispersion}
}
\value{
  a cd.fit S4 object.
}
\description{
  \code{dic.fit} fits a parametric accelerated failure time
  model to survival data.  It was developed with the
  application to incubation periods in mind. The data can
  be a mixture of doubly interval-censored, single
  interval-censored or exact observations from a single
  univariate distribution. Currently, three distributions
  are supported: log-normal, gamma, and Weibull. (The
  Erlang distribution is supported in the
  \code{dic.fit.mcmc} function, which implements an MCMC
  version of this code.) We use a consistent (par1, par2)
  notation for each distribution, they map in the following
  manner: Log-normal(meanlog=par1, sdlog=par2)
  Gamma(shape=par1, scale=par2) Weibull(shape=par1,
  scale=par2)
}
\examples{
data(fluA.inc.per)
dic.fit(fluA.inc.per, dist="L")
}
\references{
  Reich NG et al.  Statistics in Medicine.  Estimating
  incubation periods with coarse data. 2009.
  \url{http://www3.interscience.wiley.com/journal/122507367/abstract}
}
\seealso{
  \code{\link{cd.fit}}
}

