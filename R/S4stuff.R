#' An S4 Class that stores a fitted coarse data object
##' @rdname cd.fit
##' @export
setClass("cd.fit",
         representation(ests = "matrix",
                        conv = "numeric",
                        MSG="character",
                        loglik="numeric",
                        samples="data.frame",
                        data="data.frame",
                        dist="character",
                        inv.hessian="matrix",
                        est.method="character",
                        ci.method="character"))

setClass("cd.fit.mcmc",
         contains="cd.fit")

## we don't want boots and data printing out all the time
setMethod("show","cd.fit",function(object){
    cat(sprintf("Coarse Data Model Parameter and Quantile Estimates: \n"))
    print(object@ests)
    cat(sprintf("\n-2*Log Likelihood = %.1f \n",-2*object@loglik))
    if (object@conv == 0) warning("This model did not converge. Try different starting values or increase the number of iterations")
    if (object@dist == "L" & nrow(object@samples) > 0) {
        se.dispersion <- sqrt(object@inv.hessian[2,2] * (exp(object@ests[2,1] + log(object@ests[2,1])))^2) #using delta method
        cat(sprintf("\nNote: dispersion parameter is exp(sdlog). In this case it is %.3f (95%% CI %.3f-%.3f). \n",
                    exp(object@ests[2,1]),
                    exp(object@ests[2,1])-qt(.975, nrow(object@data)-1)*se.dispersion,
                    exp(object@ests[2,1])+qt(.975, nrow(object@data)-1)*se.dispersion
                    ))
    }
})

## we don't want boots and data printing out all the time
setMethod("show","cd.fit.mcmc",function(object){
    cat(sprintf("Coarse Data Model Parameter and Quantile Estimates: \n"))
    print(object@ests)
    if (object@dist == "L") cat(sprintf("\n Note: dispersion parameter is exp(sdlog). In this case it is %.3f. \n",exp(object@ests[2,1])))
    cat(sprintf("Note: please check that the MCMC converged on the target distribution by running multiple chains. MCMC samples are available in the mcmc slot (e.g. my.fit@mcmc) \n"))
})


##' log likelihood of a cd.fit object
##' 
##' returns the log-likelihood of a cd.fit object
##' 
##' @export
##' @param object cd.fit object
setMethod("logLik",
          "cd.fit",
          function(object){
              object@loglik
          })

##' For now this is going to plot the estimated survival
##' function and if bootstrap or mcmc samples are present, it will plot these samples with alpha
##' @export
##' @param x cd.fit object
##' @param col.main - color for plotting the main estimate
##' @param col.samps - color for the samples (should include some alpha transparency)
##' @param plot.n.samps - how many posterior or boostrap samples do you want to plot?
##' @param ... - other options to pass to plot
setMethod("plot",
          "cd.fit",
          function(x,
                   col.main=rgb(230, 85, 13,maxColorValue=255),
                   col.samps=rgb(99,99,99,10,maxColorValue=255),
                   plot.n.samps=200,
                   ...) {

              par1 <- x@ests[1,1]
              par2 <- x@ests[2,1]

              ## which samples do we want to plot
              if (nrow(x@samples) > 0){
                  plot.samples <- x@samples[sample(nrow(x@samples),min(plot.n.samps,nrow(x@samples))),]
              }

              xlims <- range(x@data)
              xs <-seq(xlims[1],xlims[2],length=100)
              plot(-100,-100,xlim=xlims,ylim=c(0,1),xlab="Time",ylab="Proportion with Symptoms",...)
              if (x@dist == "L"){
                  if (plot.n.samps> 0 & nrow(x@samples) > 0){
                      apply(plot.samples,1,function(y) lines(xs,plnorm(xs,meanlog=y[1],sdlog=y[2]),col=col.samps))
                  }
                  lines(xs,plnorm(xs,meanlog=par1,sdlog=par2),col=col.main,lwd=2)
              } else if (x@dist == "G" | x@dist == "E"){
                  if (plot.n.samps > 0 & nrow(x@samples) > 0){
                      apply(plot.samples,1,function(y) lines(xs,pgamma(xs,shape=y[1],scale=y[2]),col=col.samps))
                  }
                  lines(xs,pgamma(xs,shape=par1,scale=par2),col=col.main,lwd=2)
              } else if (x@dist == "W"){
                  if (plot.n.samps > 0 & nrow(x@samples) > 0){
                      apply(plot.samples,1,function(y) lines(xs,pweibull(xs,y[1],y[2]),col=col.samps))
                  }
                  lines(xs,pweibull(xs,par1,par2),col=col.main,lwd=2)
              }
          })
