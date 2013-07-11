\name{dic.fit.mcmc}
\alias{dic.fit.mcmc}
\title{Fits the distribution to the passed-in data using MCMC
as implemented in MCMCpack.}
\usage{
  dic.fit.mcmc(dat, prior.par1 = c(0, 0.001),
    prior.par2 = if (dist == "L") {     c(1000, 1000)
    } else {     c(1000, 0.001) }, init.pars = c(1, 1),
    ptiles = c(0.05, 0.95, 0.99), verbose = 1000,
    burnin = 3000, n.samples = 5000, dist = "L", ...)
}
\arguments{
  \item{dat}{the data}

  \item{prior.par1}{vector of first prior parameters}

  \item{prior.par2}{vector of second prior parameters}

  \item{init.pars}{the initial parameters, defaults to
  par.prior.par1}

  \item{ptiles}{what percentiles of the incubation period
  to return estimates for}

  \item{verbose}{how often do you want a print out from
  MCMCpack on iteration number and MH acceptance rate}

  \item{burnin}{number of burnin samples}

  \item{n.samples}{number of samples to draw from the
  posterior}

  \item{dist}{distribution to be used (L for log-normal,W
  for weibull, G for Gamma, and E for erlang)}

  \item{...}{additional parameters to MCMCmetrop1R}
}
\value{
  list with (1) ests - a matrix of estimates with columns
  est (e.g., the median estimate), (2) CIlow (0.025
  quantile) and CIhigh (0.975 quantile), and (3) an mcmc
  object as defined in MCMC pack containing the posterior
  samples
}
\description{
  Similar to \code{dic.fit} but uses MCMC instead of a
  direct likelihood optimization routine to fit the model.
  Currently, four distributions are supported: log-normal,
  gamma, Weibull, and Erlang. See Details for prior
  specification.
}
\details{
  The following priors are used: \deqn{Log-normal: (par1,
  par2) ~ Gamma()} \deqn{Weibull: par1 ~ Gamma(), par2 ~
  Normal()} \deqn{Gamma: (par1,par2) ~ 1/beta}
  \deqn{Erlang: p(par1,par2) \propto 1}
}

