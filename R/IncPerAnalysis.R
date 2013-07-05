##' Fits a log-normal, Gamma, or Weibull model to doubly interval censored survival data
##' @param dat a matrix with columns named "EL", "ER", "SL", "SR", type = either "dic" or "sic" [because optim() can be sensitive to starting parameters, some of them may be changed in the options of dic.fit().]
##' @param start.sigma the log-scale starting value for the dispersion
##' @param opt.method method used by optim
##' @param mu.int the log-scale interval of possible median values (in days)
##' @param sigma.int the log-scale interval of possible dispersion values
##' @param ptiles percentiles of interest
##' @param dist what distribution to use. Default "L" for log-normal. "G" for gamma, and "W" for Weibull. Note: If dist is Gamma (G) or Weibull (W), the mu refers to the shape and sigma refers to the scale param.
##' @param n.boots number of bootstrap resamples if non-log normal model
##' @param ... additional options passed to optim
##' @return
dic.fit <- function(dat,
		    start.sigma=log(2),
		    opt.method="L-BFGS-B",
		    mu.int=c(log(.5), log(13)),
		    sigma.int=c(log(1.01), log(log(5))),
		    ptiles=c(.05, .95, .99),
                    dist="L",
                    n.boots=0,
                    ...) {

    ## check format of dat
    cnames <- colnames(dat)
    if(!("EL" %in% cnames)) stop("dat must have column named EL")
    if(!("ER" %in% cnames)) stop("dat must have column named ER")
    if(!("SL" %in% cnames)) stop("dat must have column named SL")
    if(!("SR" %in% cnames)) stop("dat must have column named SR")

    if(!("type" %in% cnames)) stop("dat must have column named type")
    if(!all(dat[,"type"] %in% c(0,1,2)))
        stop("values in type column must be either 0, 1 or 2.")

    ## check to make sure distribution is supported
    if(!dist %in% c("G","W","L")) stop("Please use one of the following distributions Log-Normal (L) , Weibull (W), or Gamma (G)")

    ## no asymptotic results for disributions other than gamma at the moment so will need bootstrap to be larger tha 0 if dist != "L"
    if(dist %in% c("G","W") & n.boots <=0) stop("You must use bootstraping with this distrbution at the moment.  Please increase n.boots to something larger than 0")

    ## fix sample size
    n <- nrow(dat)

    ## make sure dat is a matrix
    dat <- as.matrix(dat[,c("EL", "ER", "SL", "SR", "type")])
    if(class(dat)=="data.frame") stop("dat should be a matrix.")

    ## find starting values for DIC analysis using profile likelihoods
    start.mu <- optimize(f=pl.mu, interval=mu.int,
                         sigma=start.sigma, dat=dat,dist=dist)$min
    start.sigma <- optimize(f=pl.sigma, interval=sigma.int, mu=start.mu,
                                dat=dat,dist=dist)$min

    ## find MLEs for doubly censored data using optim
    tmp <- list(convergence=1)
    msg <- NULL
    fail <- FALSE
    tryCatch(tmp <- optim(par=c(start.mu, start.sigma),
                          method=opt.method, hessian=TRUE,
                          lower=c(log(0.5), log(log(1.04))),
                          fn=loglik, dat=dat,dist=dist, ...),
             error = function(e) {
                 msg <<- e$message
                 fail <<- TRUE
             },
             warning = function(w){
                 msg <<- w$message
                 fail <<- TRUE
             })

    ## also, to catch a few more errors
    if(tmp$convergence!=0 | all(tmp$hessian==0) ){
        msg <- tmp$message
        if(all(tmp$hessian==0)) msg <- paste(msg, "& hessian is singular")
        fail <- TRUE
    }

    ## back transform optim fit
    untransformed.fit.params <- dist.optim.untransform(dist,tmp$par)

    ## check if optimaization went well
    if(!fail){

        ## get asymtotic CIs and SEs
        if (dist == "L" & n.boots<=0 ){

            med <- exp(untransformed.fit.params[1])
            disp <- exp(untransformed.fit.params[2])

            norm.quants <- qnorm(ptiles)
            ests <- c(med,
                      disp,
                      med*disp^norm.quants)

            Sig <- solve(tmp$hessian)
            ses <- dic.getSE(dat=dat,mu=log(med),log.s=log(log(disp)),Sig=Sig,ptiles=ptiles,dist=dist,opt.method=opt.method,n.boots=0)

            cil <- ests - qt(.975, n-1)*ses
            cih <- ests + qt(.975, n-1)*ses
            ## save the quantile estimates
            quant.matrix <- matrix(c(ests, cil, cih, ses),
                                   nrow=2+length(ptiles), byrow=FALSE)

            ptiles.names <- paste("p", 100*ptiles, sep="")

            rownames(quant.matrix) <- c("p50", "disp", ptiles.names)
            colnames(quant.matrix) <- c("est", "CIlow", "CIhigh", "StdErr")

        } else { ## for other distributions

            Sig <- solve(tmp$hessian)

            ##get estimates and cis for shape and scale
            boot.params <- dic.getSE(dat=dat,
                                     mu=untransformed.fit.params[1], # param 1
                                     log.s=log(untransformed.fit.params[2]), # log param 2, keeping it logged to stay consistent with previous function
                                     Sig=NULL,
                                     ptiles=ptiles,
                                     dist=dist,
                                     opt.method=opt.method,
                                     n.boots=n.boots)


            na.rows <- is.na(rowSums(boot.params))
            ## if we have any bootstraps that we couldn't get the MLE for:
            if (sum(na.rows) > 0) {
                warning(sprintf("Could not estimate the MLEs for %.0f of %.0f bootstrap replications. Excluding these from the calculation of confidence intervals and standard errors so interpret with caution. \n",sum(na.rows),n.boots))
            boot.params <- boot.params[-which(na.rows),]
            }

            cis.params <- apply(boot.params,2,function(x) quantile(x,c(.025,0.975)))

            ## adding median to  below since the exp(shape) paramter no longer has the nice interpretration
            ## of the log-normal model
            ptiles.appended <- union(0.5,ptiles)

            if (dist == "L"){
                boot.funcs <- apply(boot.params,1,function(x) qlnorm(ptiles.appended,meanlog=x[1],sdlog=x[2]))
                ests <- qlnorm(ptiles.appended,untransformed.fit.params[1],untransformed.fit.params[2])
                param1.name <- "meanlog"
                param2.name <- "sdlog"

            } else if (dist == "W"){
                boot.funcs <- apply(boot.params,1,function(x) qweibull(ptiles.appended,shape=x[1],scale=x[2]))
                ests <- qweibull(ptiles.appended,shape=untransformed.fit.params[1],scale=untransformed.fit.params[2])

                param1.name <- "shape"
                param2.name <- "scale"

            } else if (dist == "G"){
                boot.funcs <- apply(boot.params,1,function(x) qgamma(ptiles.appended,shape=x[1],scale=x[2]))

                ests <- qgamma(ptiles.appended,shape=untransformed.fit.params[1],scale=untransformed.fit.params[2])

                param1.name <- "shape"
                param2.name <- "scale"
            }

            sds.params <- apply(boot.params,2,sd)

            ## get percentile estimates including the median
            ## ests.ptiles <- apply(boot.funcs,1,mean)

            cis.ptiles <- apply(boot.funcs,1,function(x) quantile(x,c(.025,.975)))
            sds.ptiles <- apply(boot.funcs,1,sd)
            quant.matrix <- matrix(c(untransformed.fit.params,ests,cis.params[1,],cis.ptiles[1,],cis.params[2,],cis.ptiles[2,],sds.params,sds.ptiles), nrow=2+1+length(ptiles), byrow=FALSE)
            ptiles.names <- paste0("p", 100*ptiles.appended)
            rownames(quant.matrix) <- c(param1.name, param2.name, ptiles.names)
            colnames(quant.matrix) <- c("est", "CIlow", "CIhigh", "SD")
        }

        if ("boot.params" %in% ls()) {
            bp <- data.frame(boot.params)
            ci.method <- "Bootstrap"
        } else {
            bp <- data.frame()
            ci.method <- "Asymtotic"
        }

        return(
            new("cd.fit",
                ests=round(quant.matrix,3),
                conv = 1,
                MSG = "",
                loglik=-tmp$value,
                samples = bp,
                data=data.frame(dat),
                dist=dist,
                inv.hessian = Sig,
                est.method = "Maximum Likihood - optim",
                ci.method = ci.method
                )
            )

        ## return(list(ests=round(quant.matrix, 3),
        ##             conv=1,
        ##             MSG=NULL,
        ##             Sig.log.scale=Sig,
        ##             loglik=-tmp$value,
        ##             dist=dist))

    } else { ## if optimization fails:
        ## return(list(ests=matrix(NA, nrow=5, ncol=4),
        ##             conv=0,
        ##             MSG=msg,
        ##             Sig.log.scale=NULL,
        ##             loglik=NULL,
        ##             dist=NULL))
                return(
                    new("cd.fit",
                        ests=matrix(NA, nrow=5, ncol=4),
                        conv = 0,
                        MSG = msg,
                        loglik=-tmp$value,
                        samples = data.frame(),
                        data=data.frame(dat),
                        dist=dist,
                        inv.hessian = NULL,
                        est.method = "Maximum Likihood - optim",
                        ci.method = ci.method
                        )
                    )
    }
}


## profile likelihood for mu -- used by dic.fit() to get starting values
pl.mu <- function(mu, sigma, dat, dist){
    loglik(pars=c(mu, sigma),dist=dist,dat=dat)
}


## profile likelihood for sigma -- used by dic.fit() to get starting values
pl.sigma <- function(sigma, mu, dat, dist){
    loglik(pars=c(mu, sigma), dist=dist, dat=dat)
}

## functions that manipulate/calculate the likelihood for the censored data

## the functions coded here are taken directly from the
## doubly interval censored likelihood notes.

fw1 <- function(t, EL, ER, SL, SR, mu, sigma, dist){
    ## function that calculates the first function for the DIC integral
    if (dist=="W"){
        (ER-SL+t) * dweibull(x=t,shape=mu,scale=sigma)
    } else if (dist=="G") {
        (ER-SL+t) * dgamma(x=t, shape=mu, scale=sigma)
    } else {
        (ER-SL+t) * dlnorm(x=t, meanlog=mu, sdlog=sigma)
    }
}

fw3 <- function(t, EL, ER, SL, SR, mu, sigma, dist){
    ## function that calculates the third function for the DIC integral
    if (dist == "W"){
	(SR-EL-t) * dweibull(x=t, shape=mu, scale=sigma)
    } else if (dist == "G"){
    	(SR-EL-t) * dgamma(x=t, shape=mu, scale=sigma)
    } else {
        (SR-EL-t) * dlnorm(x=t, meanlog=mu, sdlog=sigma)
    }
}


lik <- function(mu, sigma, EL, ER, SL, SR, type, dist){
    ## returns the right likelihood for the type of data
    ## 0 = DIC, 1=SIC, 2=exact
    if(type==0) return(diclik2(mu, sigma, EL, ER, SL, SR, dist))
    if(type==1) return(siclik(mu, sigma, EL, ER, SL, SR, dist))
    if(type==2) return(exactlik(mu, sigma, EL, ER, SL, SR, dist))
}


diclik <- function(mu, sigma, EL, ER, SL, SR, dist){
    ## calculates the DIC likelihood by integration

    ## if symptom window is bigger than exposure window
    if(SR-SL>ER-EL){
        dic1 <- integrate(fw1, lower=SL-ER, upper=SL-EL,
                          subdivisions=10,
                          mu=mu, sigma=sigma,
                          EL=EL, ER=ER, SL=SL, SR=SR,
                          dist=dist)$value
        if (dist == "W"){
            dic2 <- (ER-EL)*
                (pweibull(SR-ER, shape=mu, scale=sigma) - pweibull(SL-EL, shape=mu, scale=sigma))
        } else if (dist == "G"){
            dic2 <- (ER-EL)*
                (pgamma(SR-ER, shape=mu, scale=sigma) - pgamma(SL-EL, shape=mu, scale=sigma))
        } else {
            dic2 <- (ER-EL)*
                (plnorm(SR-ER, mu, sigma) - plnorm(SL-EL, mu, sigma))
        }
        dic3 <- integrate(fw3, lower=SR-ER, upper=SR-EL,
                          subdivisions=10,
                          mu=mu, sigma=sigma,
                          EL=EL, ER=ER, SL=SL, SR=SR,
                          dist=dist)$value
        return(dic1 + dic2 + dic3)
    }

    ## if exposure window is bigger than symptom window
    else{
        dic1 <- integrate(fw1, lower=SL-ER, upper=SR-ER,                          subdivisions=10,
                          mu=mu, sigma=sigma,
                          EL=EL, ER=ER, SL=SL, SR=SR,
                          dist=dist)$value
        if (dist == "W"){
            dic2 <- (SR-SL)*
                (pweibull(SL-EL, shape=mu, scale=sigma) - pweibull(SR-ER, shape=mu, scale=sigma))
        } else if (dist == "G"){
            dic2 <- (SR-SL)*
                (pgamma(SL-EL, shape=mu, scale=sigma) - pgamma(SR-ER, shape=mu, scale=sigma))
        } else {
            dic2 <- (SR-SL)*
                (plnorm(SL-EL, mu, sigma) - plnorm(SR-ER, mu, sigma))
        }
        dic3 <- integrate(fw3, lower=SL-EL, upper=SR-EL,
                          subdivisions=10,
                          mu=mu, sigma=sigma,
                          EL=EL, ER=ER, SL=SL, SR=SR,
                          dist=dist)$value
        return(dic1 + dic2 + dic3)
    }
}

## this dic likelihood is designed for data that has overlapping intervals
diclik2 <- function(mu, sigma, EL, ER, SL, SR, dist){
    if(SL>ER) {
        return(diclik(mu, sigma, EL, ER, SL, SR, dist))
    } else {
        lik1 <- integrate(diclik2.helper1, lower=EL, upper=SL,
                          SL=SL, SR=SR, mu=mu, sigma=sigma, dist=dist)$value
        lik2 <- integrate(diclik2.helper2, lower=SL, upper=ER,
                          SR=SR, mu=mu, sigma=sigma, dist=dist)$value
        return(lik1+lik2)
    }
}

## likelihood functions for diclik2
diclik2.helper1 <- function(x, SL, SR, mu, sigma, dist){
    if (dist =="W"){
        pweibull(SR-x, shape=mu, scale=sigma) - pweibull(SL-x, shape=mu, scale=sigma)
    } else if (dist =="G") {
        pgamma(SR-x, shape=mu, scale=sigma) - pgamma(SL-x, shape=mu, scale=sigma)
    } else {
        plnorm(SR-x, mu, sigma) - plnorm(SL-x, mu, sigma)
    }
}

diclik2.helper2 <- function(x, SR, mu, sigma, dist){
    if (dist =="W"){
        pweibull(SR-x, shape=mu, scale=sigma)
    } else if (dist =="G") {
        pgamma(SR-x, shape=mu, scale=sigma)
    } else {
	plnorm(SR-x, mu, sigma)
    }
}

siclik <- function(mu, sigma, EL, ER, SL, SR, dist){
    ## calculates the SIC likelihood as the difference in CDFs
    if (dist =="W"){
        pweibull(SR-EL, shape=mu, scale=sigma) - pweibull(SL-ER, shape=mu, scale=sigma)
    } else if (dist =="G") {
        pgamma(SR-EL, shape=mu, scale=sigma) - pgamma(SL-ER, shape=mu, scale=sigma)
    } else {
        plnorm(SR-EL, mu, sigma) - plnorm(SL-ER, mu, sigma)
    }
}

exactlik <- function(mu, sigma, EL, ER, SL, SR, dist){
    ## calculates the likelihood for an exact observation

    ## NB: the two Ss should be equal and the two Es should be equal
    ##     so it doesn't matter which pair we use in the formula below.
    if (dist =="W"){
        dweibull(SR-EL, shape=mu, scale=sigma)
    } else if (dist =="G") {
        dgamma(SR-EL, shape=mu, scale=sigma)
    } else {
        dlnorm(SR-EL, mu, sigma)
    }
}


##' negative log likelihood for a set of parameters, data, and a distribution
##' @param pars transformed parameters \in (-\infty,\infty)
##' @param dat data
##' @param dist distribution
##' @return negative log likelihood
loglik <- function(pars, dat, dist) {
    ## calculates the log-likelihood of DIC data
    ## dat must have EL, ER, SL, SR and type columns

    ## expecting transformed params from optimiztion
    ## e.g. for log-normal expecting c(mu,log.sigma)
    pars <- dist.optim.untransform(dist,pars)
    mu <- pars[1]
    sigma <- pars[2]

    ## cat(sprintf("mu = %.2f, sigma = %.2f \n",mu, sigma))  ## for debugging
    n <- nrow(dat)
    totlik <- 0
    for(i in 1:n){
        totlik <- totlik +
            log(lik(mu, sigma, type=dat[i,"type"],
                    EL=dat[i,"EL"], ER=dat[i,"ER"],
                    SL=dat[i,"SL"], SR=dat[i,"SR"],
                    dist=dist))
    }
    return(-totlik) ## NB: NEEDS TO BE -totlik IF WE ARE MAXIMIZING USING OPTIM!
    ## May want to change this name later to reflect that is it negative log lik
}


## calculates the standard errors for estimates from dic.fit() using delta method or bootstrap
##' @param mu - shape param
##' @param log.s - log.scale param
##' @param Sig
##' @param ptiles
##' @param dist
##' @param dat
##' @param opt.method
##' @param n.boots number of bootstraps
##' @return
dic.getSE <- function(mu, log.s, Sig, ptiles, dist, dat, opt.method, n.boots=100){
    boots <- vector("list",n.boots)

    if (dist == "L" & n.boots <= 0) {
        cat(sprintf("Computing Asymtotic Confidence Intervals for Log Normal Model \n"))
        s <- exp(log.s)
        qnorms <- qnorm(ptiles)
        df <- matrix(c(exp(mu), 0, exp(mu+qnorms*s),
                       0, exp(s+log.s), qnorms * exp(mu + qnorms*s + log.s)),
                     nrow=2, ncol=2+length(ptiles), byrow=TRUE)
        ses <- sqrt(diag(t(df)%*%Sig%*%df))
        return(ses)

    } else {

        cat(sprintf("Bootstrapping (n=%i) Standard Errors for %s \n",n.boots,dist))
        ## sample line numbers from the data
        line.nums <- matrix(sample(1:nrow(dat),nrow(dat)*n.boots,replace=T),nrow=nrow(dat),ncol=n.boots)
        ## set up progress bar
        pb <- txtProgressBar(min = 0, max = n.boots, style = 3)
        for (i in 1:n.boots){
            boots[[i]] <-
                single.boot(mu.s=mu,sigma.s=exp(log.s),opt.method=opt.method,dat.tmp=dat[line.nums[,i],],dist=dist)
            setTxtProgressBar(pb, i)
        }

        close(pb)

        ## grab the params from each
        ## remember if any failed there will be NAs here
        mus <- sapply(boots,function(x) x$par[1])
        sigmas <- sapply(boots,function(x) x$par[2])

        return(cbind(mus=mus,sigmas=sigmas))
    }
}

## estimates one set of parameters for a single bootstrap resample
##' @param mu.s starting value for first param
##' @param s.s starting value for second param
##' @param opt.method optimization method
##' @param dat.tmp one relaization of resampled data
##' @param dist distribution
##' @param ...
##' @return returns optim list object with estiamtes for the untransformed two parameters of the specified dist
single.boot <- function(mu.s,sigma.s,opt.method,dat.tmp,dist,...){
    tmp <- list(convergence=1)
    msg <- NULL
    fail <- FALSE
    pars.transformed <- dist.optim.transform(dist,c(mu.s,sigma.s))
    tryCatch(tmp <- optim(par=pars.transformed,
                          method=opt.method, hessian=FALSE,
                          lower=c(-10,-10),
                          fn=loglik, dat=dat.tmp,dist=dist,...),
             error = function(e) {
                 msg <- e$message
                 fail <- TRUE
             },
             warning = function(w){
                 msg <- w$message
                 fail <- TRUE
             },
             if(tmp$convergence!=0 | all(tmp$hessian==0) ){
                 msg <- tmp$message
                 if(all(tmp$hessian==0)) msg <- paste(msg, "& hessian is singular")
                 fail <- TRUE
             })

    ## transform back to original scale
    ## return NAs if we can't find the min for this param set
    if(is.null(tmp$par)){
        tmp$par <- c(NA,NA)
        } else {
            tmp$par <- dist.optim.untransform(dist,tmp$par)
        }

    return(tmp)
}


get.obs.type <- function(dat) {
    type <- rep(0, nrow(dat))
    ## get the single interval censored
    type[dat[,"EL"]==dat[,"ER"]]<-1
    type[dat[,"SL"]==dat[,"SR"]]<-1
    type[dat[,"ER"]>=dat[,"SL"]]<-1

    ## some of those are actually exact!
    type[(dat[,"EL"]==dat[,"ER"]) & (dat[,"SL"]==dat[,"SR"])]<- 2
    return(type)
}


##' Fits the distribution to the passed in data using MCMC
##' as implemented in MCMCpack. The following priors are used:
##' Survival Model = Log-normal --> $(\mu,\sigma) \sim Gamma()$
##' Survival Model = Weibull --> $\alpha \sim Gamma()$, $\beta \sim Normal()$
##' Survival Model = Gamma --> $\alpha,\beta) \sim \frac{1}{\beta}$
##' @param dat the data
##' @param par.prior.param1 vector of first prior parameters
##' @param par.prior.param2 vector of second prior parameters
##' @param init.pars the initial parameters, defaults to par.prior.mu
##' @param ptiles what percentiles of the incubation period to return estimates for
##' @param verbose how often do you want a print out from MCMCpack on iteration number and MH acceptance rate
##' @param burnin number of burnin samples
##' @param n.samples number of samples to draw from the posterior
##' @param dist distribution to be used
##' @param ... additional parameters to MCMCmetrop1R
##' @param par.prior.mu the mean for the prior distribution for the parameters a vector of [log median, log log dispersion]
##' @param par.prior.sd vector of standard deviations for the prior on the log scale
##' @return list with (1) ests - a matrix of estimates with columns est (e.g., the median estimate), (2) CIlow (0.025 quantile) and CIhigh (0.975 quantile), and (3) an mcmc object as defined in MCMC pack containing the posterior samples
dic.fit.mcmc <- function(dat,
                         par.prior.param1 = c(0,0.001),
                         par.prior.param2 = c(1000,0.001),
                         init.pars = c(1,1),
                         ptiles = c(0.05,0.95,0.99),
                         verbose=1000,#how often to print update
                         burnin = 3000,
                         n.samples = 5000,
                         dist = "L",
                         ...){

    require(MCMCpack)

    ## log liklihood function to pass to MCMCpack sampler
    local.ll <- function(pars,
                         dat,
                         par.prior.1,
                         par.prior.2,
                         dist) {

        ## get parameters on untransformed scale
        pars.untrans <- dist.optim.untransform(dist,pars)

        if (dist == "L"){
            ## default gamma on scale param and (inproper) uniform on location
            ll <- tryCatch(-loglik(pars,dat,dist) +
                           dgamma(pars.untrans[2],shape=par.prior.param1[2],
                                  rate=par.prior.param2[2],log=T),
                           ## sum(dnorm(pars.untrans,
                           ##           par.prior.param1,
                           ##           par.prior.param2,log=T)),
                           error=function(e) {
                               warning("Loglik failure, returning -Inf")
                               return(-Inf)
                           })

        } else if (dist == "W"){
            ## using normal prior on the first param and gamma on second
            ll <- tryCatch(-loglik(pars,dat,dist) +
                           dnorm(pars.untrans[1],
                                 par.prior.param1[1],
                                 par.prior.param2[1],log=T) +
                           dgamma(pars.untrans[2],
                                  shape=par.prior.param1[2],
                                  rate=par.prior.param2[2],log=T),
                           error=function(e) {
                               recover()
                               warning("Loglik failure, returning -Inf")
                               return(-Inf)
                           })
        } else if (dist == "G"){
            ## using "non-informative" prior 1/scale for the joint prior \pi(a,b) \propto \frac{1}{\beta}
            ll <- tryCatch(-loglik(pars,dat,dist) + log(1/pars.untrans[2]),
                           error=function(e) {
                               warning("Loglik failure, returning -Inf")
                               return(-Inf)
                           })

        } else {
            stop("Sorry, unknown distribution type. Check the 'dist' option.")
        }
        return(ll)
    }


    cat(sprintf("Running %.0f MCMC iterations to get Bayesian estimates \n",n.samples+burnin))

    msg <- NULL
    fail <- FALSE

    ## run the MCMC chains
    tryCatch(mcmc.run <- MCMCmetrop1R(local.ll,
                                      init.pars,
                                      dat = dat,
                                      par.prior.1 = par.prior.param1,
                                      par.prior.2 = par.prior.param2,
                                      verbose=verbose,
                                      dist=dist,
                                      burnin = burnin,
                                      mcmc=n.samples,
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
        ptiles.appended <- union(0.5,ptiles)
        est.pars <- matrix(nrow=length(ptiles.appended)+2,ncol=3)

        if (dist == "L"){
            param1.name <- "meanlog"
            param2.name <- "sdlog"

            mcmc.quantiles <- apply(untrans.mcmcs,1,function(x) qlnorm(ptiles.appended,meanlog=x[1],sdlog=x[2]))


        } else if (dist == "G"){
            param1.name <- "shape"
            param2.name <- "scale"
            mcmc.quantiles <- apply(untrans.mcmcs,1,function(x) qgamma(ptiles.appended,shape=x[1],scale=x[2]))

        } else if (dist == "W"){
            param1.name <- "shape"
            param2.name <- "scale"
            mcmc.quantiles <- apply(untrans.mcmcs,1,function(x) qweibull(ptiles.appended,shape=x[1],scale=x[2]))

        } else {
            stop("Sorry, unknown distribution type. Check the 'dist' option.")
            ## not actually needed but just in case
        }

        ## make the return matrix
        colnames(est.pars) <- c("est","CIlow", "CIhigh")
        rownames(est.pars) <- c(param1.name,param2.name,paste0("p", 100*ptiles.appended))

        est.pars[1,] <- quantile(untrans.mcmcs[,1], c(0.5,0.025,0.975))
        est.pars[2,] <- quantile(untrans.mcmcs[,2], c(0.5,0.025,0.975))
        cis.ptiles <- t(apply(mcmc.quantiles,1,function(x) quantile(x,c(0.5,.025,.975))))
        est.pars[3:nrow(est.pars),1:3] <- cis.ptiles

        ## return(list(ests=round(est.pars,3),
        ##             MSG=NULL,
        ##             mcmc=mcmc.run,
        ##             dist=dist))

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


##' Transforms parameters of a specific distriution for unbounded optimization
##' @param dist string representing distirbution
##' @param pars vector of parameters
##' @return vector of transformed parameters
dist.optim.transform <- function(dist,pars){
    if (dist == "G"){
        log(pars) # for shape and scale
    } else if (dist == "W"){
        log(pars) # for shape and scale
    } else if (dist == "Erlang"){
        stop("Not yet complete")
    } else if (dist == "L"){
        c(pars[1],log(pars[2])) # for meanlog, sdlog
    } else {
        stop(sprintf("Distribtion (%s) not supported",dist))
    }
}

##' Untransforms parameters
##' @param dist
##' @param pars
##' @return vector of untransformed parameters
dist.optim.untransform <- function(dist,pars){
    if (dist == "G"){
        exp(pars) # for shape and scale
    } else if (dist == "W"){
        exp(pars) # for shape and scale
    } else if (dist == "Erlang"){
        stop("Not yet complete")
    } else if (dist == "L"){
        c(pars[1],exp(pars[2])) # for meanlog, sdlog
    } else {
        stop(sprintf("Distribtion (%s) not supported",dist))
    }
}
