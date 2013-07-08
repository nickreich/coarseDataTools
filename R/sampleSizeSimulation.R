######################################
## code for sample size simulations ##
## Nicholas Reich
## August 2010
######################################

## this is the central user function
## N = overall sample size
## med = median of log normal distribution
## disp = dispersion of log normal distribution
## percentile = what percentile(s) should be simulated
## nsim = how many simulations to run
## exact data = T/F should simulation be done with exact data
## pct.type.A = percent type A data, will be rounded up to nearest integer
## exp.win.dat = a vector of exposure window lengths to sample from
## verb = whether to print 10 iteration counts during the course of the sim
precision.simulation <- function(N,
			       med=2,
			       disp=1.3,
			       percentile=.5,
			       nsim=100,
			       exact.data=FALSE,
			       pct.type.A=.5,
 			       exp.win.dat=NULL,
			       verb=FALSE) {
	## logic check
	if(percentile <= 0 | percentile >=1)
		stop("percentile must be between 0 and 1.")
	if(pct.type.A < 0 | pct.type.A >1)
		stop("% of data that is type A must be between 0 and 1.")
	if(is.null(exp.win.dat)){
		if(verb) message("NYC exposure window data used")
		#data("exp.win.lengths")
		exp.win.dat <- get(data("exp.win.lengths"))
	}

	## TODO: add default exp.win.dat to package and load if NULL

	if(exact.data) {
		out <- precision.simulation.exact(N=N,
						med=med,
						disp=disp,
						percentile=percentile,
						nsim=nsim,
						verb=verb)
	} else {
		out <- precision.simulation.coarse(N=N,
						 med=med,
						 disp=disp,
						 percentile=percentile,
						 nsim=nsim,
						 pct.type.A=pct.type.A,
						 exp.win.dat=exp.win.dat,
						 verb=verb)
	}

	target <- qlnorm(percentile, log(med), log(disp))
	bias <- out[,"ests"]-target
	out <- cbind(out, bias)

	return(out)
}


precision.simulation.exact <- function(N,
                                       med,
                                       disp,
                                       percentile,
                                       nsim,
                                       verb) {
    storage <- matrix(NA, ncol=3, nrow=nsim)
    colnames(storage) <- c("ests", "SE", "conv")

    data <- matrix(0, ncol=6, nrow=nsim*N)
    colnames(data) <- c("dataset.id", "EL", "ER", "SL", "SR", "type")
    data[,"dataset.id"] <- rep(1:nsim, each=N)
    data[,"SL"] <- rlnorm(N*nsim, meanlog=log(med), sdlog=log(disp))
    data[,"SR"] <- data[,"SL"]
    data[,"type"] <- 2

    for(i in 1:nsim){
        tmp.dat <- data[which(data[,"dataset.id"]==i),]
        tmp.fit <- dic.fit(tmp.dat, ptiles=percentile)
        if(tmp.fit$conv==1){
            row.name <- paste("p", round(percentile*100), sep="")
            which.row <-
                which(rownames(tmp.fit$ests)==row.name)[1]
            storage[i,c("ests", "SE")] <-
                tmp.fit$ests[which.row, c("est", "StdErr")]
        } else {
            storage[i,c("ests", "SE")] <- NA
        }
        storage[i,"conv"] <- tmp.fit$conv
        if(verb & i%%(round(nsim/10))==0)
			print(paste("iteration",i,"complete ::", Sys.time()))

    }
    return(storage)
}


precision.simulation.coarse <- function(N,
				      med,
				      disp,
				      percentile,
				      nsim,
				      pct.type.A,
				      exp.win.dat,
				      verb) {
	## create storage
	storage <- matrix(NA, ncol=3, nrow=nsim)
	colnames(storage) <- c("ests", "SE", "conv")

	for(i in 1:nsim){
		tmp.dat <- generate.coarse.data(N=N,
						med=med,
						disp=disp,
						pct.type.A=pct.type.A,
						exp.win.dat=exp.win.dat)
		tmp.fit <- dic.fit(tmp.dat, ptiles=percentile)
		if(tmp.fit$conv==1){
			row.name <- paste("p", round(percentile*100), sep="")
			which.row <-
                            which(rownames(tmp.fit$ests)==row.name)[1]
			storage[i,c("ests", "SE")] <-
				tmp.fit$ests[which.row, c("est", "StdErr")]
		} else {
			storage[i,c("ests", "SE")] <- NA
		}
		storage[i,"conv"] <- tmp.fit$conv
		if(verb & i%%(round(nsim/10))==0)
			print(paste("iteration",i,"complete ::", Sys.time()))
	}
	return(storage)
}



generate.coarse.data <- function(N, med, disp, pct.type.A, exp.win.dat) {

 	n.type.A <- round(N*pct.type.A)
	n.type.B <- N-n.type.A


	E <- runif(N, 10, 11)
	T <- rlnorm(N, log(med), log(disp))
	S <- T + E
	SR <- ceiling(S)
	SL <- floor(S)

	## generate window types
	##   0 = short with bounded ER = type A
	##   1 = "long" with no ER = type B
	win.type <- rep(0:1, times=c(n.type.A, n.type.B))

	## generate window lengths,
	potential.lengths <- sample(exp.win.dat, size=N, replace=TRUE)
	win.length <- 1*(win.type==0) + potential.lengths*(win.type>0)

	## fix data
	ER <- ceiling(E)*(win.type==0) + SR*(win.type==1)
	EL <- pmin(ER-win.length, floor(E))

	## return the data
	cbind(EL, E, ER, SL, S, SR, win.length, win.type, type=0)
}
