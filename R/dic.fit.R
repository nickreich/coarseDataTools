##' Fits a log-normal, Gamma, Erlang, or Weibull model to doubly interval 
##' censored survival data
##' 
##' \code{dic.fit} fits a parametric accelerated failure time model to survival 
##' data.   It was developed with the application to incubation periods in mind.
##' The data can be a mixture of doubly interval-censored, single
##' interval-censored or exact observations from a single univariate
##' distribution. Currently, three distributions are supported: log-normal,
##' gamma, and Weibull. (The Erlang distribution is supported in the
##' \code{dic.fit.mcmc} function, which implements an MCMC version of this
##' code.) We use a consistent (par1, par2) notation for each distribution, they
##' map in the following manner: 
##' Log-normal(meanlog=par1, sdlog=par2) 
##' Gamma(shape=par1, scale=par2) 
##' Weibull(shape=par1, scale=par2)
##' 
##' 
##' 
##' @param dat a matrix with columns named "EL", "ER", "SL", "SR", corresponding
##'   to the left (L) and right (R) endpoints of the windows of possible 
##'   exposure (E) and symptom onset (S). Also, a "type" column must be 
##'   specified and have entries with 0, 1, or 2, corresponding to doubly 
##'   interval-censored, single interval-censored or exact observations, 
##'   respsectively.
##' @param start.par2 starting value for 2nd parameter of desired distribtution
##' @param opt.method method used by optim
##' @param par1.int the log-scale interval of possible median values (in the 
##'   same units as the observations in dat).  Narrowing this interval can help 
##'   speed up convergence of the algorithm, but care must be taken so that 
##'   possible values are not excluded or that the maximization does not return 
##'   a value at an endpoint of this interval.
##' @param par2.int the log-scale interval of possible dispersion values
##' @param ptiles percentiles of interest
##' @param dist what distribution to use to fit the data. Default "L" for 
##'   log-normal. "G" for gamma, and "W" for Weibull. Note: If dist is Gamma (G)
##'   or Weibull (W), the mu refers to the shape and sigma refers to the scale 
##'   param.
##' @param n.boots number of bootstrap resamples if non-log normal model
##' @param ... additional options passed to optim
##' @param start.sigma the log-scale starting value for the dispersion
##' @return a cd.fit S4 object.
##' @seealso \code{\link{cd.fit}}
##' @export
##' @examples
##' data(fluA.inc.per)
##' dic.fit(fluA.inc.per, dist="L")
##' @references Reich NG et al.  Statistics in Medicine.  Estimating incubation periods with coarse data. 2009.  \url{http://www3.interscience.wiley.com/journal/122507367/abstract}

dic.fit <- function(dat,
		    start.par2=log(2),
		    opt.method="L-BFGS-B",
		    par1.int=c(log(.5), log(13)),
		    par2.int=c(log(1.01), log(log(5))),
		    ptiles=c(.05, .95, .99),
                    dist="L",
                    n.boots=0,
                    ...) {

    ## check format of dat
    check.data.structure(dat)

    ## check to make sure distribution is supported
    if(!dist %in% c("G","W","L")) stop("Please use one of the following distributions Log-Normal (L) , Weibull (W), or Gamma (G)")

    ## no asymptotic results for disributions other than gamma at the moment so will need bootstrap to be larger tha 0 if dist != "L"
    if(dist %in% c("G") & n.boots <=0) stop("You must use bootstraping with this distrbution at the moment.  Please increase n.boots to something larger than 0")

    ## fix sample size
    n <- nrow(dat)

    ## make sure dat is a matrix
    dat <- as.matrix(dat[,c("EL", "ER", "SL", "SR", "type")])
    if(class(dat)=="data.frame") stop("dat should be a matrix.")

    ## find starting values for DIC analysis using profile likelihoods
    start.par1 <- optimize(f=pl.par1, interval=par1.int,
                           par2=start.par2, dat=dat,dist=dist)$min
    start.par2 <- optimize(f=pl.par2, interval=par2.int, par1=start.par1,
                           dat=dat,dist=dist)$min

    ## find MLEs for doubly censored data using optim
    tmp <- list(convergence=1)
    msg <- NULL
    fail <- FALSE
    tryCatch(tmp <- optim(par=c(start.par1, start.par2),
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

        ## always going to report median even if not requested
        ptiles.appended <- union(0.5,ptiles)

        ## get asymtotic CIs and SEs
        if (dist == "L" & n.boots<=0 ){

            med <- exp(untransformed.fit.params[1])
            disp <- exp(untransformed.fit.params[2])

            norm.quants <- qnorm(ptiles.appended)
            ests <- c(untransformed.fit.params[1],
                      untransformed.fit.params[2],
                      med*disp^norm.quants)
            Sig <- solve(tmp$hessian)
            ses <- dic.getSE(dat=dat,par1=log(med),log.par2=log(log(disp)),Sig=Sig,ptiles=ptiles.appended,dist=dist,opt.method=opt.method)
            ## get cis
            cil <- ests - qt(.975, n-1)*ses
            cih <- ests + qt(.975, n-1)*ses
            ## save the quantile estimates
            quant.matrix <- matrix(c(ests, cil, cih, ses),
                                   nrow=2+length(ptiles.appended), byrow=FALSE)
            ptiles.names <- paste0("p", 100*ptiles.appended)

            rownames(quant.matrix) <- c("meanlog", "sdlog", ptiles.names)
            colnames(quant.matrix) <- c("est", "CIlow", "CIhigh", "StdErr")

        } else if (dist == "W" & n.boots <=0){
            shape <- untransformed.fit.params[1]
            scale <- untransformed.fit.params[2]

            ests <- c(shape,
                      scale,
                      scale*(-log(1-ptiles.appended))^(1/shape))

            Sig <- solve(tmp$hessian)
            ses <- dic.getSE(dat=dat,
                             par1=shape,
                             log.par2=log(scale),
                             Sig=Sig,
                             ptiles=ptiles.appended,
                             dist=dist,
                             opt.method=opt.method)
            ## get cis
            cil <- ests - qt(.975, n-1)*ses
            cih <- ests + qt(.975, n-1)*ses

            ## save the quantile estimates
            quant.matrix <- matrix(c(ests, cil, cih, ses),
                                   nrow=2+length(ptiles.appended), byrow=FALSE)

            ptiles.names <- paste0("p", 100*ptiles.appended)

            rownames(quant.matrix) <- c("shape", "scale", ptiles.names)
            colnames(quant.matrix) <- c("est", "CIlow", "CIhigh", "StdErr")

        } else { ## for other distributions

            Sig <- solve(tmp$hessian)

            ##get estimates and cis for shape and scale
            boot.params <- dic.get.boots(dat=dat,
                                         par1=untransformed.fit.params[1],
                                         par2=untransformed.fit.params[2], # keeping it logged to stay consistent with previous function
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

            ## std deviations of bootstraps for parameters
            sds.params <- apply(boot.params,2,sd)

            ## get percentile estimates
            cis.ptiles <- apply(boot.funcs,1,function(x) quantile(x,c(.025,.975)))
            sds.ptiles <- apply(boot.funcs,1,sd)
            quant.matrix <- matrix(c(untransformed.fit.params,ests,cis.params[1,],cis.ptiles[1,],cis.params[2,],cis.ptiles[2,],sds.params,sds.ptiles), nrow=2+1+length(ptiles), byrow=FALSE)
            ## deal with row and column names
            ptiles.names <- paste0("p", 100*ptiles.appended)
            rownames(quant.matrix) <- c(param1.name, param2.name, ptiles.names)
            colnames(quant.matrix) <- c("est", "CIlow", "CIhigh", "SD")
        }

        if ("boot.params" %in% ls()) {
            bp <- data.frame(boot.params)
            ci.method <- "Bootstrap"
        } else {
            bp <- data.frame()
            ci.method <- "Asymptotic"
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

    } else { ## if optimization fails:

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
##' @param par1
##' @param par2
##' @param dat
##' @param dist
pl.par1 <- function(par1, par2, dat, dist){
    loglik(pars=c(par1, par2),dist=dist, dat=dat)
}


## profile likelihood for sigma -- used by dic.fit() to get starting values
##' @param par2
##' @param par1
##' @param dat
##' @param dist
pl.par2 <- function(par2, par1, dat, dist){
    loglik(pars=c(par1, par2), dist=dist, dat=dat)
}

## functions that manipulate/calculate the likelihood for the censored data

##' the functions coded here are taken directly from the
##' doubly interval censored likelihood notes.
##' @param t
##' @param EL
##' @param ER
##' @param SL
##' @param SR
##' @param par1
##' @param par2
##' @param dist
fw1 <- function(t, EL, ER, SL, SR, par1, par2, dist){
    ## function that calculates the first function for the DIC integral
    if (dist=="W"){
        (ER-SL+t) * dweibull(x=t,shape=par1,scale=par2)
    } else if (dist=="G") {
        (ER-SL+t) * dgamma(x=t, shape=par1, scale=par2)
    } else {
        (ER-SL+t) * dlnorm(x=t, meanlog=par1, sdlog=par2)
    }
}

##' @param t
##' @param EL
##' @param ER
##' @param SL
##' @param SR
##' @param par1
##' @param par2
##' @param dist
fw3 <- function(t, EL, ER, SL, SR, par1, par2, dist){
    ## function that calculates the third function for the DIC integral
    if (dist == "W"){
	(SR-EL-t) * dweibull(x=t, shape=par1, scale=par2)
    } else if (dist == "G"){
    	(SR-EL-t) * dgamma(x=t, shape=par1, scale=par2)
    } else {
        (SR-EL-t) * dlnorm(x=t, meanlog=par1, sdlog=par2)
    }
}


##' @param par1
##' @param par2
##' @param EL
##' @param ER
##' @param SL
##' @param SR
##' @param type
##' @param dist
lik <- function(par1, par2, EL, ER, SL, SR, type, dist){
    ## returns the right likelihood for the type of data
    ## 0 = DIC, 1=SIC, 2=exact
    if(type==0) return(diclik2(par1, par2, EL, ER, SL, SR, dist))
    if(type==1) return(siclik(par1, par2, EL, ER, SL, SR, dist))
    if(type==2) return(exactlik(par1, par2, EL, ER, SL, SR, dist))
}


##' calculates the DIC likelihood by integration
##' @param par1
##' @param par2
##' @param EL
##' @param ER
##' @param SL
##' @param SR
##' @param dist
diclik <- function(par1, par2, EL, ER, SL, SR, dist){

    ## if symptom window is bigger than exposure window
    if(SR-SL>ER-EL){
        dic1 <- integrate(fw1, lower=SL-ER, upper=SL-EL,
                          subdivisions=10,
                          par1=par1, par2=par2,
                          EL=EL, ER=ER, SL=SL, SR=SR,
                          dist=dist)$value
        if (dist == "W"){
            dic2 <- (ER-EL)*
                (pweibull(SR-ER, shape=par1, scale=par2) - pweibull(SL-EL, shape=par1, scale=par2))
        } else if (dist == "G"){
            dic2 <- (ER-EL)*
                (pgamma(SR-ER, shape=par1, scale=par2) - pgamma(SL-EL, shape=par1, scale=par2))
        } else {
            dic2 <- (ER-EL)*
                (plnorm(SR-ER, par1, par2) - plnorm(SL-EL, par1, par2))
        }
        dic3 <- integrate(fw3, lower=SR-ER, upper=SR-EL,
                          subdivisions=10,
                          par1=par1, par2=par2,
                          EL=EL, ER=ER, SL=SL, SR=SR,
                          dist=dist)$value
        return(dic1 + dic2 + dic3)
    }

    ## if exposure window is bigger than symptom window
    else{
        dic1 <- integrate(fw1, lower=SL-ER, upper=SR-ER,                          subdivisions=10,
                          par1=par1, par2=par2,
                          EL=EL, ER=ER, SL=SL, SR=SR,
                          dist=dist)$value
        if (dist == "W"){
            dic2 <- (SR-SL)*
                (pweibull(SL-EL, shape=par1, scale=par2) - pweibull(SR-ER, shape=par1, scale=par2))
        } else if (dist == "G"){
            dic2 <- (SR-SL)*
                (pgamma(SL-EL, shape=par1, scale=par2) - pgamma(SR-ER, shape=par1, scale=par2))
        } else {
            dic2 <- (SR-SL)*
                (plnorm(SL-EL, par1, par2) - plnorm(SR-ER, par1, par2))
        }
        dic3 <- integrate(fw3, lower=SL-EL, upper=SR-EL,
                          subdivisions=10,
                          par1=par1, par2=par2,
                          EL=EL, ER=ER, SL=SL, SR=SR,
                          dist=dist)$value
        return(dic1 + dic2 + dic3)
    }
}

##' this dic likelihood is designed for data that has overlapping intervals
##' @param par1
##' @param par2
##' @param EL
##' @param ER
##' @param SL
##' @param SR
##' @param dist
diclik2 <- function(par1, par2, EL, ER, SL, SR, dist){
    if(SL>ER) {
        return(diclik(par1, par2, EL, ER, SL, SR, dist))
    } else {
        lik1 <- integrate(diclik2.helper1, lower=EL, upper=SL,
                          SL=SL, SR=SR, par1=par1, par2=par2, dist=dist)$value
        lik2 <- integrate(diclik2.helper2, lower=SL, upper=ER,
                          SR=SR, par1=par1, par2=par2, dist=dist)$value
        return(lik1+lik2)
    }
}

##' likelihood functions for diclik2
##' @param x
##' @param SL
##' @param SR
##' @param par1
##' @param par2
##' @param dist
diclik2.helper1 <- function(x, SL, SR, par1, par2, dist){
    if (dist =="W"){
        pweibull(SR-x, shape=par1, scale=par2) - pweibull(SL-x, shape=par1, scale=par2)
    } else if (dist =="G") {
        pgamma(SR-x, shape=par1, scale=par2) - pgamma(SL-x, shape=par1, scale=par2)
    } else {
        plnorm(SR-x, par1, par2) - plnorm(SL-x, par1, par2)
    }
}

##' @param x
##' @param SR
##' @param par1
##' @param par2
##' @param dist
diclik2.helper2 <- function(x, SR, par1, par2, dist){
    if (dist =="W"){
        pweibull(SR-x, shape=par1, scale=par2)
    } else if (dist =="G") {
        pgamma(SR-x, shape=par1, scale=par2)
    } else {
	plnorm(SR-x, par1, par2)
    }
}


##' @param par1
##' @param par2
##' @param EL
##' @param ER
##' @param SL
##' @param SR
##' @param dist
siclik <- function(par1, par2, EL, ER, SL, SR, dist){
    ## calculates the SIC likelihood as the difference in CDFs
    if (dist =="W"){
        pweibull(SR-EL, shape=par1, scale=par2) - pweibull(SL-ER, shape=par1, scale=par2)
    } else if (dist =="G") {
        pgamma(SR-EL, shape=par1, scale=par2) - pgamma(SL-ER, shape=par1, scale=par2)
    } else {
        plnorm(SR-EL, par1, par2) - plnorm(SL-ER, par1, par2)
    }
}

##' @param par1
##' @param par2
##' @param EL
##' @param ER
##' @param SL
##' @param SR
##' @param dist
exactlik <- function(par1, par2, EL, ER, SL, SR, dist){
    ## calculates the likelihood for an exact observation

    ## NB: the two Ss should be equal and the two Es should be equal
    ##     so it doesn't matter which pair we use in the forpar1la below.
    if (dist =="W"){
        dweibull(SR-EL, shape=par1, scale=par2)
    } else if (dist =="G") {
        dgamma(SR-EL, shape=par1, scale=par2)
    } else {
        dlnorm(SR-EL, par1, par2)
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
    ## e.g. for log-normal expecting c(meanlog,log(sdlog))
    pars <- dist.optim.untransform(dist,pars)
    par1 <- pars[1]
    par2 <- pars[2]

    ## cat(sprintf("par1 = %.2f, par2 = %.2f \n",par1, par2))  ## for debugging
    n <- nrow(dat)
    totlik <- 0
    for(i in 1:n){
        totlik <- totlik +
            log(lik(par1, par2, type=dat[i,"type"],
                    EL=dat[i,"EL"], ER=dat[i,"ER"],
                    SL=dat[i,"SL"], SR=dat[i,"SR"],
                    dist=dist))
    }
    return(-totlik) ## NB: NEEDS TO BE -totlik IF WE ARE MAXIMIZING USING OPTIM!
    ## May want to change this name later to reflect that is it negative log lik
}


## calculates the standard errors for estimates from dic.fit() using delta method (NOTE: only works for Log Normal and Weibull Models at the moment)
##' @param par1
##' @param log.par2 log.scale param 2 (log-log dispersion for log-normal)
##' @param Sig - var-cov matrix from hessian
##' @param ptiles - percentiles of interest
##' @param dist - failure time distribtion
##' @param dat - data
##' @param opt.method - optimization method for optim (see ?optim for options)
##' @return asymptotic standard errors for Log-Normal or Weibull models
dic.getSE <- function(par1, log.par2, Sig, ptiles, dist, dat, opt.method){

        cat(sprintf("Computing Asymtotic Confidence Intervals \n"))
        par2 <- exp(log.par2) # log.par2 input historically so I kept it as is

        if (dist == "L"){
            qnorms <- qnorm(ptiles)
            df <- matrix(c(1, 0,exp(par1+qnorms*par2),
                           0, par2, qnorms * exp(par1 + qnorms*par2 + log.par2)),
                         nrow=2, ncol=2+length(ptiles), byrow=TRUE)
        } else if (dist == "W"){
            df <- matrix(c(par1,0,
                           -(par2*(-log(1-ptiles))^(1/par1)*log(-log(1-ptiles)))/par1^2, #d/dmu
                           0,par2,
                           (-log(1-ptiles))^(1/par1)), #d/ds
                         nrow=2, ncol=2+length(ptiles), byrow=TRUE)
        }
        ses <- sqrt(diag(t(df)%*%Sig%*%df))
        return(ses)
    }

##' @param par1
##' @param par2
##' @param dist distribution
##' @param dat data
##' @param opt.method optim method
##' @param n.boots number of bootstraps
##' @return matrix of bootstrap estimates of untransformed parameters for distrbution
dic.get.boots <- function(par1, par2, dist, dat, opt.method, n.boots=100){
    cat(sprintf("Bootstrapping (n=%i) Standard Errors for %s \n",n.boots,dist))
    boots <- vector("list",n.boots)

    ## sample line numbers from the data
    line.nums <- matrix(sample(1:nrow(dat),nrow(dat)*n.boots,replace=T),nrow=nrow(dat),ncol=n.boots)
    ## set up progress bar
    pb <- txtProgressBar(min = 0, max = n.boots, style = 3)
    for (i in 1:n.boots){
        boots[[i]] <-
            single.boot(par1.s=par1,par2.s=par2,opt.method=opt.method,dat.tmp=dat[line.nums[,i],],dist=dist)
        setTxtProgressBar(pb, i)
    }
    close(pb)

    ## grab the params from each
    ## remember if any failed there will be NAs here
    par1s <- sapply(boots,function(x) x$par[1])
    par2s <- sapply(boots,function(x) x$par[2])

    return(cbind(par1=par1s,par2=par2s))
}

## estimates one set of parameters for a single bootstrap resample
##' @param par1.s starting value for first param
##' @param par2.s starting value for second param
##' @param opt.method optimization method
##' @param dat.tmp one relaization of resampled data
##' @param dist distribution
##' @param ...
##' @return returns optim list object with estimates for the untransformed two parameters of the specified dist
single.boot <- function(par1.s,par2.s,opt.method,dat.tmp,dist,...){

    tmp <- list(convergence=1)
    msg <- NULL
    fail <- FALSE
    pars.transformed <- dist.optim.transform(dist,c(par1.s,par2.s))
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


##' Transforms parameters of a specific distriution for unbounded optimization
##' @param dist string representing distirbution
##' @param pars vector of parameters
##' @return vector of transformed parameters
dist.optim.transform <- function(dist,pars){
    if (dist == "G"){
        log(pars) # for shape and scale
    } else if (dist == "W"){
        log(pars) # for shape and scale
    } else if (dist == "E"){
        log(pars)
    } else if (dist == "L"){
        c(pars[1],log(pars[2])) # for meanlog, sdlog
    } else {
        stop(sprintf("Distribtion (%s) not supported",dist))
    }
}

##' Untransforms parameters before entering likelihood
##' @param dist
##' @param pars
##' @return vector of untransformed parameters
dist.optim.untransform <- function(dist,pars){
    if (dist == "G"){
        exp(pars) # for shape and scale
    } else if (dist == "W"){
        exp(pars) # for shape and scale
    } else if (dist == "E"){
        ## we want shape to be restricted to integers
        tmp <- exp(pars)
        tmp[1] <- round(tmp[1],0)
        tmp
    } else if (dist == "L"){
        c(pars[1],exp(pars[2])) # for meanlog, sdlog
    } else {
        stop(sprintf("Distribtion (%s) not supported",dist))
    }
}

##' Issues a stop if the data does not conform with the expected structure
##' @param dat data for use in CDT
##' @return NULL
##' @author Andrew Azman
check.data.structure <- function(dat){
    ## check format of dat
    cnames <- colnames(dat)
    if(!("EL" %in% cnames)) stop("dat must have column named EL")
    if(!("ER" %in% cnames)) stop("dat must have column named ER")
    if(!("SL" %in% cnames)) stop("dat must have column named SL")
    if(!("SR" %in% cnames)) stop("dat must have column named SR")

    if(!("type" %in% cnames)) stop("dat must have column named type")
    if(!all(dat[,"type"] %in% c(0,1,2)))
        stop("values in type column must be either 0, 1 or 2.")

    if (any(is.na(dat[,c("EL","ER","SL","SR","type")]))) stop("Missing (NA) values not permitted")
    NULL
}
