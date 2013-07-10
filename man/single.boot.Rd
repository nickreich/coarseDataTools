\name{single.boot}
\alias{single.boot}
\usage{
  single.boot(par1.s, par2.s, opt.method, dat.tmp, dist,
    ...)
}
\arguments{
  \item{par1.s}{starting value for first param}

  \item{par2.s}{starting value for second param}

  \item{opt.method}{optimization method}

  \item{dat.tmp}{one relaization of resampled data}

  \item{dist}{distribution}

  \item{...}{}
}
\value{
  returns optim list object with estimates for the
  untransformed two parameters of the specified dist
}

