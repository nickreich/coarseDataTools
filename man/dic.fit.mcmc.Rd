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
  a cd.fit.mcmc S4 object
}
\description{
  Similar to \code{dic.fit} but uses MCMC instead of a
  direct likelihood optimization routine to fit the model.
  Currently, four distributions are supported: log-normal,
  gamma, Weibull, and Erlang. See Details for prior
  specification.
}
\details{
  The following models are used: \deqn{Log-normal model:
  f(x) = \frac{1}{x*\sigma \sqrt{2 * \pi}}
  exp\{-\frac{(\log x - \mu)^2}{2 * \sigma^2}\}}
  \deqn{Log-normal Default Prior: \mu ~ N( 0, 1000), \sigma
  ~ N(0,1000)} \deqn{Weibull model: f(x) =
  \frac{\alpha}{\beta}(\frac{x}{\beta})^{\alpha-1}
  exp\{-(\frac{x}{\beta})^{\alpha}\}} \deqn{Weibull Default
  Prior Specification: \alpha ~ N( 0, 1000), \beta ~
  Gamma(0.001,0.001)} \deqn{Gamma model: f(x) =
  \frac{1}{\theta^k \Gamma(k)} x^{k-1}
  exp\{-\frac{x}{\theta}\}} \deqn{Gamma Default Prior
  Specification: p(k,\theta) \propto \frac{1}{\theta}}
  \deqn{Erlang model: f(x) = \frac{1}{\theta^k (k-1)!}
  x^{k-1} exp\{-\frac{x}{\theta}\}} \deqn{Erlang Default
  Prior Specification: p(k,\theta) \propto 1}
}

