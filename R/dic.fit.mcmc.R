##' Fits the distribution to the passed-in data using MCMC
##' as implemented in MCMCpack.
##'
##' Similar to \code{dic.fit} but uses MCMC instead of a direct likelihood optimization routine to fit the model. Currently, four distributions are supported: log-normal, gamma, Weibull, and Erlang
##'
##'   The following priors are used:
##'   Survival Model = Log-normal --> $(par1,par2) ~ Gamma()$
##'   Survival Model = Weibull --> $par1 ~ Gamma()$, $par2 ~ Normal()$
##'   Survival Model = Gamma --> $(par1,par2) ~ 1/beta$
##'   Survival Model = Erlang --> $p(par1,par2) proportionalto 1$
##' @param dat the data
##' @param prior.par1 vector of first prior parameters
##' @param prior.par2 vector of second prior parameters
##' @param init.pars the initial parameters, defaults to par.prior.par1
##' @param ptiles what percentiles of the incubation period to return estimates for
##' @param verbose how often do you want a print out from MCMCpack on iteration number and MH acceptance rate
##' @param burnin number of burnin samples
##' @param n.samples number of samples to draw from the posterior
##' @param dist distribution to be used (L for log-normal,W for weibull, G for Gamma, and E for erlang)
##' @param ... additional parameters to MCMCmetrop1R
##' @return list with (1) ests - a matrix of estimates with columns est (e.g., the median estimate), (2) CIlow (0.025 quantile) and CIhigh (0.975 quantile), and (3) an mcmc object as defined in MCMC pack containing the posterior samples
dic.fit.mcmc <- function(dat,
                         prior.par1 = c(0,0.001),
                         prior.par2 = c(1000,0.001),
                         init.pars = c(1,1),
                         ptiles = c(0.05,0.95,0.99),
                         verbose=1000,#how often to print update
                         burnin = 3000,
                         n.samples = 5000,
                         dist = "L",
                         ...){

        require(MCMCpack)

        ## check to make sure data is well formed for CDT use:
        check.data.structure(dat)

        ## check to make sure distribution is supported
        if(!dist %in% c("G","W","L","E")) stop("Please use one of the following distributions Log-Normal (L) , Weibull (W), Gamma (G), or Erlang (E)")

        ## log liklihood function to pass to MCMCpack sampler
        local.ll <- function(pars,
                             dat,
                             prior.par1,
                             prior.par2,
                             dist) {

                ## get parameters on untransformed scale
                pars.untrans <- dist.optim.untransform(dist,pars)

                if (dist == "L"){
                        ## default gamma on scale param and (inproper) uniform on location
                        ll <- tryCatch(-loglikhd(pars,dat,dist) +
                                               ## dgamma(pars.untrans[2],shape=par.prior.param1[2],
                                               ## rate=par.prior.param2[2],log=T),
                                               sum(dnorm(pars.untrans,
                                                         prior.par1,
                                                         prior.par2,log=T)),
                                       error=function(e) {
                                               warning("Loglik failure, returning -Inf")
                                               return(-Inf)
                                       })

                } else if (dist == "W"){
                        ## using normal prior on the first param and gamma on second
                        ll <- tryCatch(
                                -loglikhd(pars,dat,dist) +
                                        dnorm(pars.untrans[1],
                                              prior.par1[1],
                                              prior.par2[1],log=T) +
                                        dgamma(pars.untrans[2],
                                               shape=prior.par1[2],
                                               rate=prior.par2[2],log=T),
                                error=function(e) {
                                        warning("Loglik failure, returning -Inf")
                                        return(-Inf)
                                })
                } else if (dist == "G"){
                        ## using "non-informative" prior 1/scale for the joint prior \pi(a,b) \propto \frac{1}{\beta}
                        ll <- tryCatch(-loglikhd(pars,dat,dist) + log(1/pars.untrans[2]),
                                       error=function(e) {
                                               warning("Loglik failure, returning -Inf")
                                               return(-Inf)
                                       })

                } else if (dist == "E"){ # for Erlang
                        ## Erlang is just a gamma so we are going to use this trick
                        ll <- tryCatch(-loglikhd(pars,dat,dist="G"),
                                       # no priors for now will add later
                                       error=function(e) {
                                               warning("Loglik failure, returning -Inf")
                                               return(-Inf)
                                       })

                } else {
                        stop("Sorry, unknown distribution type. Check the 'dist' option.")
                }
                return(ll)
        }

        cat(sprintf("Running %.0f MCMC iterations \n",n.samples+burnin))

        msg <- NULL
        fail <- FALSE

        ## run the MCMC chains
        tryCatch(mcmc.run <- MCMCmetrop1R(local.ll,
                                          init.pars,
                                          dat = dat,
                                          prior.par1 = prior.par1,
                                          prior.par2 = prior.par2,
                                          verbose=verbose,
                                          dist=dist,
                                          burnin = burnin,
                                          mcmc=n.samples,
                                          logfun=TRUE,
                                          ...),
                 error=function(e){
                         msg <<- e$message
                         fail <<- TRUE
                 },
                 warning = function(w){
                         msg <<- w$message
                         fail <<- TRUE
                 })

        if (!fail){
                ## untransform MCMC parameter draws to natural scale
                untrans.mcmcs <- t(apply(mcmc.run[,1:2],1,function(x) dist.optim.untransform(dist=dist,pars=x)))

                ## append median to the percentiles in case it isn't there
                ptiles.appended <- sort(union(0.5,ptiles))
                est.pars <- matrix(nrow=length(ptiles.appended)+2,ncol=3)

                if (dist == "L"){
                        par1.name <- "meanlog"
                        par2.name <- "sdlog"
                        mcmc.quantiles <- apply(untrans.mcmcs,1,function(x) qlnorm(ptiles.appended,meanlog=x[1],sdlog=x[2]))
                } else if (dist == "G"){
                        par1.name <- "shape"
                        par2.name <- "scale"
                        mcmc.quantiles <- apply(untrans.mcmcs,1,function(x) qgamma(ptiles.appended,shape=x[1],scale=x[2]))
                } else if (dist == "W"){
                        par1.name <- "shape"
                        par2.name <- "scale"
                        mcmc.quantiles <- apply(untrans.mcmcs,1,function(x) qweibull(ptiles.appended,shape=x[1],scale=x[2]))
                } else if (dist == "E"){
                        par1.name <- "shape"
                        par2.name <- "scale"
                        mcmc.quantiles <- apply(untrans.mcmcs,1,function(x) qgamma(ptiles.appended,shape=x[1],scale=x[2]))
                } else {
                        stop("Sorry, unknown distribution type. Check the 'dist' option.")
                        ## not actually needed but just in case
                }

                ## make the return matrix
                colnames(est.pars) <- c("est","CIlow", "CIhigh")
                rownames(est.pars) <- c(par1.name,par2.name,paste0("p", 100*ptiles.appended))

                est.pars[1,] <- quantile(untrans.mcmcs[,1], c(0.5,0.025,0.975))
                est.pars[2,] <- quantile(untrans.mcmcs[,2], c(0.5,0.025,0.975))
                cis.ptiles <- t(apply(mcmc.quantiles,1,function(x) quantile(x,c(0.5,.025,.975))))
                est.pars[3:nrow(est.pars),1:3] <- cis.ptiles

                rc <- new("cd.fit.mcmc",
                          ests=round(est.pars,3),
                          conv = numeric(),
                          MSG = "",
                          loglik=numeric(),
                          samples = data.frame(untrans.mcmcs),
                          data=data.frame(dat),
                          dist=dist,
                          inv.hessian = matrix(),
                          est.method = "MCMC",
                          ci.method = "MCMC"
                )

                return(rc)

        } else {
                rc <- new("cd.fit.mcmc",
                          ests=matrix(NA, nrow=5, ncol=3),
                          conv = numeric(),
                          MSG = msg,
                          loglik=numeric(),
                          samples = data.frame(),
                          data=data.frame(dat),
                          dist=dist,
                          inv.hessian = matrix(),
                          est.method = "MCMC",
                          ci.method = "MCMC"
                )

                print("Try adjusting the starting parameters init.pars")
                return(rc)
        }
}
