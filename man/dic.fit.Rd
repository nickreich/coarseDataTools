\name{dic.fit}
\alias{dic.fit}
\title{
Fitting log-normal models to doubly interval-censored data
}
\description{
This function fits a log-normal model to doubly interval-censored survival data.  It was developed with the application to incubation periods in mind.  The data can be a mixture of doubly interval-censored, single interval-censored or exact observations from a single univariate distribution.
}
\usage{
dic.fit(dat, 
        start.log.sigma = log(log(2)), 
        opt.method = "L-BFGS-B", 
        mu.int = c(log(0.5), log(13)), 
        log.sigma.int = c(log(log(1.01)), log(log(5))), 
        ptiles = c(0.05, 0.95, 0.99), ...)
}

\arguments{
  \item{dat}{
a matrix with columns named "EL", "ER", "SL", "SR", corresponding to the left (L) and right (R) endpoints of the windows of possible exposure (E) and symptom onset (S). Also, a "type" column must be specified and have entries with 0, 1, or 2, corresponding to doubly interval-censored, single interval-censored or exact observations, respsectively.
}
  \item{start.log.sigma}{
the log-log-scale starting value for the dispersion, used in a profile likelihood maximization which determines the starting values for the full likelihood maximization routine.
}
  \item{opt.method}{
the optimization method to use in the optim() function 
}
  \item{mu.int}{
the log-scale interval of possible median values (in the same units as the observations in dat).  Narrowing this interval can help speed up convergence of the algorithm, but care must be taken so that possible values are not excluded or that the maximization does not return a value at an endpoint of this interval. 
}
  \item{log.sigma.int}{
the log-log-scale interval of possible dispersion values.  See notes for mu.int above.
}
  \item{ptiles}{
percentiles of the distribution to be estimated in this routine.  Note that the median and dispersion are estimated by default.  See the ests matrix in the returned list below for more information.
}
  \item{\dots}{
other parameters to be passed to optim() which runs the full likelihood maximization. 
}
}
\details{
See instructions in the arguments above.
}
\value{
A list is returned with elements as follows:
  \item{ests }{a matrix summarizing the results of fitting the model.  Rows correspond to the median, dispersion and then percentiles specified by the ptiles argument.  Columns correspond to the point estimate, the lower and upper bounds on the 95\% confidence interval and the standard error of the point estimate. If the maximization does not converge, this matrixis filled with NAs.}
  \item{conv }{A value of 1 indicates successful convergence; 0 indicates unsuccessful convergence. }
  \item{MSG }{The error message returned from optim() if the routine fails to converge. }
  \item{Sig.log.scale }{The inverse of the hessian matrix for the likelihood surface at the MLE.  Used to determine the standard errors for the percentiles.}
}
\references{
Reich NG et al.  Statistics in Medicine.  Estimating incubation periods with coarse data. 2009.  \url{http://www3.interscience.wiley.com/journal/122507367/abstract}

Lessler J et al.  Lancet Infectious Diseases.  Incubation periods of acute respiratory viral infections: a systematic review.  \url{http://www.thelancet.com/journals/laninf/article/PIIS1473309909700696/abstract}
}
\author{
Nicholas G. Reich
}

\examples{
data(fluA.inc.per)
dic.fit(fluA.inc.per)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ models }
