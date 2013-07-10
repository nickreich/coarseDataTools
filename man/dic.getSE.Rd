\name{dic.getSE}
\alias{dic.getSE}
\usage{
  dic.getSE(par1, log.par2, Sig, ptiles, dist, dat,
    opt.method)
}
\arguments{
  \item{par1}{}

  \item{log.par2}{log.scale param 2 (log-log dispersion for
  log-normal)}

  \item{Sig}{- var-cov matrix from hessian}

  \item{ptiles}{- percentiles of interest}

  \item{dist}{- failure time distribtion}

  \item{dat}{- data}

  \item{opt.method}{- optimization method for optim (see
  ?optim for options)}
}
\value{
  asymptotic standard errors for Log-Normal or Weibull
  models
}

