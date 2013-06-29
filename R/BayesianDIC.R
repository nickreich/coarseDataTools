##' Bayesian version of Log-normal model
##' Fits the distribution to the passed in data using MCMC
##' as implemented in MCMCpack asssuming lognormal prior for the median
##' and log-log normal prior for the dispersion.
##' @param dat the data
##' @param par.prior.mu the mean for the prior distribution for the parameters a vector of [log median, log log dispersion]
##' @param par.prior.sd vector of standard deviations for the prior on the log scale
##' @param init.pars the initial parameters, defaults to par.prior.mu
##' @param ptiles what percentiles of the incubation period to return estimates for
##' @param verbose how often do you want a print out from MCMCpack on iteration number and acceptance rate
##' @param ... additional parameters to MCMCmetrop1R
##' @return list with (1) ests - a matrix of estimates with columns est (e.g., the median estimate), (2) CIlow (0.025 quantile) and CIhigh (0.975 quantile), and (3) an mcmc object as defined in MCMC pack containing the posterior samples
dic.fit.mcmc <- function(dat,
                         par.prior.mu = c(0,0),
                         par.prior.sd = c(1000,1000),
                         init.pars = par.prior.mu,
                         ptiles = c(0.05,0.95,0.99),
                         verbose=1000,#how often to print update
                         burnin = 30000,
                         n.samples = 50000
                         ...) {

    require(MCMCpack)

    if (n.samples <= burnin) stop("Number of samples needs to be more than the number of burnin samples")

    #log liklihood function to pass to MCMCpack sampler
    local.ll <- function(pars,
                         dat,
                         par.prior.mu,
                         par.prior.sd,dist) {
        if (dist == "L"){
            ll <- tryCatch(-loglik(pars,dat,dist)+sum(dnorm(pars, par.prior.mu, par.prior.sd, log=T)),
                           error=function(e) {
                               warning("Loglik failure, returning -Inf")
                               return(-Inf)
                           })
        } else if (dist == "W"){
            ll <- tryCatch(-loglik(pars,dat,dist) +
                           dnorm(pars[1], par.prior.mu[1], par.prior.sd[1],log=T) +
                           dgamma(pars[2],par.prior.mu[2],par.prior.sd[2],log=T),
                           error=function(e) {
                               warning("Loglik failure, returning -Inf")
                               return(-Inf)
                           })
            ##use normal prior on shape and gamma on scale
        } else if (dist == "G"){
            #
        } else {
            stop("Unknown Distribution")
        }

        return(ll)
    }

    #run the MCMC chains
    mcmc.run <- MCMCmetrop1R(local.ll,
                             init.pars,
                             dat = dat,
                             par.prior.mu = par.prior.mu,
                             par.prior.sd = par.prior.sd,
                             verbose=verbose,
                             dist="L",
                             burnin = burnin,
                             ...)


    #make the return matrix
    est.pars <- matrix(nrow=length(ptiles)+2,
                       ncol=3)
    colnames(est.pars) <- c("est","CIlow", "CIhigh")
    rownames(est.pars) <- c("m", "disp", sprintf("p%d", 100*ptiles))

    est.pars[1,] <- quantile(exp(mcmc.run[,1]), c(0.5,0.025,0.975))
    est.pars[2,] <- quantile(exp(exp(mcmc.run[,2])), c(0.5,0.025,0.975))

    #make sure at least one ptile was requested
    if (length(ptiles)>0) {
        for (i in 1:length(ptiles)) {
            est.pars[i+2,] <-
                quantile(qlnorm(ptiles[i],mcmc.run[,1], exp(mcmc.run[,2])),
                         c(0.5, 0.025, 0.975))
        }
    }

    return(list(ests=est.pars, mcmc=mcmc.run))

}





#Fits the distribution to the passed in data using MCMC
#as implemented in MCMCpack asssuming lognormal prior for the median
#and log-log normal prior for the dispersion.
##
#Parameters -
#   dat - the data
#   par.prior.1 - prior parameters for the log-median. First entry is the mean, decond is sd
#   par.prior.2 - prior paramter for the log-log-dispersion, prior is gamma(a,b)
#   init.pars - the initial parameters, defaults to par.prior.mu
#   ptiles - what percentiles of the incubation period to return estimates for
#   ... - additional parameters to MCMCmetrop1R
#
#Returns -
#  A list with the following elements:
#  ests - a matrix of estimates with columns est (e.g., the median estimate),
#    CIlow (0.025 quantile) and CIhigh (0.975 quantile)
#  mcmc - an mcmc object as defined in MCMC pack containing the posterior
#    samples
dic.fit.mcmc.test<- function(dat,
                       #  par.prior.1 = c(0,100),
                       #  par.prior.2 = c(0.0001,0.0001),
                         init.pars = c(0,1),
                         ptiles = c(0.05,0.95,0.99),
                         verbose=1000, #how often to print update
                         ...) {

    require(MCMCpack)

    #log liklihood function to pass to MCMCpack sampler
    local.ll <- function(pars,
                         dat,
                       #  par.prior.1,
                       #  par.prior.2,
                         dist) {
        if (dist == "L"){
            ll <- tryCatch(-loglik(pars,dat,dist)+log(1/pars[2]),
                           error=function(e) {
                               warning("Loglik failure, returning -Inf")
                               return(-Inf)
                           })
        } else if (dist == "W"){
        } else if (dist == "G"){
        } else {
            stop("Unknown Distribution")
        }

        return(ll)
    }


    #run the MCMC chains
    mcmc.run <- MCMCmetrop1R(local.ll,
                             init.pars,
                             dat = dat,
                            # par.prior.1 = par.prior.1,
                            # par.prior.2 = par.prior.2,
                             verbose=verbose,
                             dist="L",
                             burnin = 50000,
                             ...)


    #make the return matrix
    est.pars <- matrix(nrow=length(ptiles)+2,
                       ncol=3)
    colnames(est.pars) <- c("est","CIlow", "CIhigh")
    rownames(est.pars) <- c("m", "disp", sprintf("p%d", 100*ptiles))

    est.pars[1,] <- quantile(exp(mcmc.run[,1]), c(0.5,0.025,0.975))
    est.pars[2,] <- quantile(exp(exp(mcmc.run[,2])), c(0.5,0.025,0.975))

    #make sure at least one ptile was requested
    if (length(ptiles)>0) {
        for (i in 1:length(ptiles)) {
            est.pars[i+2,] <-
                quantile(qlnorm(ptiles[i],mcmc.run[,1], exp(mcmc.run[,2])),
                         c(0.5, 0.025, 0.975))
        }
    }

    return(list(ests=est.pars, mcmc=mcmc.run))

}


