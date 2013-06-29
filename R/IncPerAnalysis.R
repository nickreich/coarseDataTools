##' dic.fit() fits a log-normal model to doubly interval censored incubation period data
##' @param dat a matrix with columns named "EL", "ER", "SL", "SR", type = either "dic" or "sic" [because optim() can be sensitive to starting parameters, some of them may be changed in the options of dic.fit().]
##' @param start.log.sigma the log-log-scale starting value for the dispersion
##' @param opt.method method used by optim
##' @param mu.int the log-scale interval of possible median values (in days)
##' @param log.sigma.int the log-log-scale interval of possible dispersion values
##' @param ptiles percentiles of interest
##' @param dist what distribution to use. Default "L" for log-normal. "G" for gamma, and "W" for Weibull. Note: If dist is Gamma (G) or Weibull (W), the mu refers to the shape and sigma refers to the scale param.
##' @param n.boots number of bootstrap resamples if non-log normal model
##' @param bayesian if yes, estimates paramters with a Bayesian model using non-informative priors (for log-normal model only)
##' @param ... additional options passed to optim
##' @return
dic.fit <- function(dat,
		    start.log.sigma=log(log(2)),
		    opt.method="L-BFGS-B",
		    mu.int=c(log(.5), log(13)),
		    log.sigma.int=c(log(log(1.01)), log(log(5))),
		    ptiles=c(.05, .95, .99),
                    dist="L",
                    n.boots=100,
                    bayesian=FALSE,
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

#    if(bayesian == TRUE & dist != "L") stop("Bayesian analysis only available for Log-Normal model at the moment")

    ## check to make sure disitribution is supported
    if(!dist %in% c("G","W","L")) stop("Please use one of the following distributions Log-Normal (L) , Weibull (W), or Gamma (G)")
    ## fix sample size
    n <- nrow(dat)

    ## make sure dat is a matrix
    dat <- as.matrix(dat[,c("EL", "ER", "SL", "SR", "type")])
    if(class(dat)=="data.frame") stop("dat should be a matrix.")

    ## find starting values for DIC analysis using profile likelihoods
    start.mu <- optimize(f=pl.mu, interval=mu.int,
                         log.sigma=start.log.sigma, dat=dat,dist=dist)$min
    start.log.sigma <- optimize(f=pl.sigma, interval=log.sigma.int, mu=start.mu,
                                dat=dat,dist=dist)$min

    ## find MLEs for doubly censored data using optim
    tmp <- list(convergence=1)
    msg <- NULL
    fail <- FALSE
    if (bayesian == TRUE){
        ## source("BayesianDIC.R")
        cat(sprintf("Running MCMCs to get Bayesian estimates \n"))
        tryCatch(tmp <- dic.fit.mcmc(dat=dat,
                                     ptiles=ptiles,
                                     ...),
                 error = function(e) {
                     msg <<- e$message
                     fail <<- TRUE
                 },
                 warning = function(w){
                         msg <<- w$message
                         fail <<- TRUE
                     })
        if (!fail){
            ## return list with results mcmc draws
            ## TO DO: may want to consider wrapping all of this into an S4 object and setting it so it doesn't spit out mcmc by default
            return(list(ests=round(tmp$est, 3),
                        conv=1,
                        MSG=NULL,
                        Sig.log.scale=NULL,
                        loglik=NULL,
                        dist=dist,
                        mcmc=tmp$mcmc))
        } else {
            return(list(ests=matrix(NA, nrow=5, ncol=4),
                        conv=0,
                        MSG=msg,
                        Sig.log.scale=NULL,
                        loglik=NULL,
                        dist=NULL))

        }} else {
            tryCatch(tmp <- optim(par=c(start.mu, start.log.sigma),
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

            if(!fail){

                if (dist == "L"){
                    med <- exp(tmp$par[1])
                    disp <- exp(exp(tmp$par[2]))
                    norm.quants <- qnorm(ptiles)
                    ests <- c(med,
                              disp,
                              med*disp^norm.quants)
                    Sig <- solve(tmp$hessian)
                    ses <- dic.getSE(log(med), log(log(disp)), Sig, ptiles,dist=dist)
                    cil <- ests - qt(.975, n-1)*ses
                    cih <- ests + qt(.975, n-1)*ses
                    ## save the quantile estimates
                    quant.matrix <- matrix(c(ests, cil, cih, ses),
                                           nrow=2+length(ptiles), byrow=FALSE)
                    ptiles.names <- paste("p", 100*ptiles, sep="")
                    rownames(quant.matrix) <- c("p50", "disp", ptiles.names)
                    colnames(quant.matrix) <- c("est", "CIlow", "CIhigh", "StdErr")

                } else {
                    shape <- exp(tmp$par[1])
                    scale <- exp(exp(tmp$par[2]))
                    Sig <- solve(tmp$hessian)

                    ##get estimates and cis for shape and scale
                    boot.params <- dic.getSE(dat=dat,mu=tmp$par[1],log.s=tmp$par[2], Sig=NULL, ptiles=ptiles,dist=dist,opt.method=opt.method,n.boots=n.boots)
                    cis.params <- apply(boot.params,2,function(x) quantile(x,c(.025,0.975)))

                    ##adding 0.5 to percentiles below since the exp(shape) paramter no longer has the nice interpretration of the log-normal model
                    if (dist == "W"){
                        boot.funcs <- apply(boot.params,1,function(x) qweibull(c(0.5,ptiles),shape=x[1],scale=x[2]))
                    } else if (dist == "G"){
                        boot.funcs <- apply(boot.params,1,function(x) qgamma(c(0.5,ptiles),shape=x[1],scale=x[2]))
                    }

                    sds.params <- apply(boot.params,2,sd)
                                        #get percentile estimates including the median
                    ests.ptiles <- apply(boot.funcs,1,mean)
                    cis.ptiles <- apply(boot.funcs,1,function(x) quantile(x,c(.025,.975)))
                    sds.ptiles <- apply(boot.funcs,1,sd)
                    quant.matrix <- matrix(c(shape,scale,ests.ptiles,cis.params[1,],cis.ptiles[1,],cis.params[2,],cis.ptiles[2,],sds.params,sds.ptiles),
                                           nrow=2+1+length(ptiles), byrow=FALSE)
                    ptiles.names <- paste("p", 100*c(.5,ptiles), sep="")
                    rownames(quant.matrix) <- c("shape", "scale", ptiles.names)
                    colnames(quant.matrix) <- c("est", "CIlow", "CIhigh", "SD")
                }

                return(list(ests=round(quant.matrix, 3),
                            conv=1,
                            MSG=NULL,
                            Sig.log.scale=Sig,
                            loglik=-tmp$value,
                            dist=dist))
            }
            else {
		return(list(ests=matrix(NA, nrow=5, ncol=4),
			    conv=0,
			    MSG=msg,
			    Sig.log.scale=NULL,
                            loglik=NULL,
                            dist=NULL))
            }
        }
}


## profile likelihood for mu -- used by dic.fit() to get starting values
pl.mu <- function(mu, log.sigma, dat, dist){
    loglik(pars=c(mu, log.sigma),dist=dist,dat=dat)
}


## profile likelihood for sigma -- used by dic.fit() to get starting values
pl.sigma <- function(log.sigma, mu, dat, dist){
    loglik(pars=c(mu, log.sigma), dist=dist, dat=dat)
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

loglik <- function(pars, dat, dist) {
    ## calculates the log-likelihood of DIC data
    ## dat must have EL, ER, SL, SR and type columns
    mu <- pars[1]
    sigma <- exp(pars[2])
    sprintf("mu = %.2f, sigma = %.2f",mu, sigma)  ## for debugging
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
}


## calculates the standard errors for estimates from dic.fit() using delta method or bootstrap
dic.getSE <- function(mu, log.s, Sig, ptiles, dist,dat=dat,opt.method,n.boots=500){
    boots <- vector("list",n.boots)
    if (dist == "L") {
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

        line.nums <- matrix(sample(1:nrow(dat),nrow(dat)*n.boots,replace=T),nrow=nrow(dat),ncol=n.boots)
        pb <- txtProgressBar(min = 0, max = n.boots, style = 3)
        for (i in 1:n.boots){
            boots[[i]] <-
                single.boot(mu.s=mu,log.s.s=log.s,opt.method=opt.method,dat.tmp=dat[line.nums[,i],],dist=dist)
            setTxtProgressBar(pb, i)
        }
        close(pb)
        mus <- sapply(boots,function(x) x$par[1])
        sigmas <- sapply(boots,function(x) exp(x$par[2]))
        return(cbind(shape=mus,scale=sigmas))
    }
}

## estimates one set of parameters for bootstraps
single.boot <- function(mu.s,log.s.s,opt.method,dat.tmp,dist,...){
    tmp <- list(convergence=1)
    msg <- NULL
    fail <- FALSE
    tryCatch(tmp <- optim(par=c(mu.s,log.s.s),
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
