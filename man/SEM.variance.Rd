\name{SEM.variance}
\alias{SEM.variance}
\title{Implementation of the Supplemented EM algorithm}
\description{
This function is meant to be run only through the function EMforCFR() and is used to calculate the variance via the Supplemented EM algorithm (see Meng and Rubin, 1991)
}
\usage{
SEM.variance(full.data, 
	     dat, 
	     phi, 
	     max.iter, 
	     tol, 
	     nlag, 
	     alpha.start.values, 
	     assumed.nu)
}
\arguments{
  \item{full.data}{ A matrix of observed data. See description in EMforCFR helpfile.}

  \item{dat}{A data frame.}
  \item{phi}{ A vector of fitted parameters from the final EM iteration.}
  \item{max.iter}{ The maximum number of iterations for SEM algorithm. }
  \item{tol}{ A tolerance to use to test for convergence of the EM algorithm. }
  \item{nlag}{ The number of time units for lagged data.  Corresponds to length(assumed.nu).}  
  \item{alpha.start.values}{ a vector starting values for the reporting rate parameter of the GLM model. This must have length which corresponds to one less than the number of unique integer values of full.dat[,"new.times"].}
 \item{assumed.nu}{ a vector of probabilities corresponding to the survival distribution, i.e. nu[i]=Pr(surviving i days | fatal case) }

 
}

\value{ A list with the following components
  \item{DM}{The estimate of the variance-covariance matrix for the model parameters.  Only converged rows are returned.}
  \item{DMiter}{A vector whose ith entry is the number of iterations needed for convergence of the ith row of the DM matrix.}
  \item{loop.idx}{If not NULL, the values correspond to the original indices of DM which have been omitted because of lack of convergence.}

}
\references{
Meng, X.L. and Rubin, D.B.  JASA. 1991: 86 (416), 899-909.}
\author{ Nicholas G. Reich }
\seealso{\code{\link{EMforCFR}}}
\keyword{ supplemented EM algorithm }
\keyword{ EM algorithm }
