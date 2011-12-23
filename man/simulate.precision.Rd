\name{precision.simulation}
\alias{precision.simulation}
\alias{precision.simulation.coarse}
\alias{precision.simulation.exact}
\alias{generate.coarse.data}

\title{
Simulate incubation period analyses with coarse data
}
\description{
These functions simulate coarse incubation period data sets and analyze them.  The goal is for these simulations to provide evidence for how much information a given dataset contains about a characteristic of the incubation period distribution. 
}
\usage{
precision.simulation(N, 	   
		     med = 2, 
		     disp = 1.3, 
		     percentile = 0.5, 
		     nsim = 100, 
		     exact.data = FALSE, 
		     pct.type.A = 0.5, 
		     exp.win.dat = NULL, 
		     verb = FALSE)
precision.simulation.coarse(N, 
			    med, 
			    disp, 
			    percentile, 
			    nsim, 
			    pct.type.A, 
			    exp.win.dat, 
			    verb)
precision.simulation.exact(N,
			   med, 
			   disp, 
			   percentile, 
			   nsim, 
			   verb)
generate.coarse.data(N,
		     med, 
		     disp, 
		     pct.type.A, 
		     exp.win.dat)

}

\arguments{
  \item{N}{
Overall sample size for the datasets to be simulated.
}
  \item{med}{
Median for the assumed log normal distribution of the incubation periods.
}
  \item{disp}{
Dispersion for the assumed log normal distribution of the incubation periods.
}
  \item{percentile}{
Percentile of the incubation period distribution which we want to estimate. 
}
  \item{nsim}{
Number of datasets to analyze in the simulation.
}
  \item{exact.data}{
Either TRUE/FALSE.  Incidates whether the data generated should be coarsened at all.  If TRUE, pct.type.A and exp.win.dat are ignored.
}
  \item{pct.type.A}{
Percent of the N observations that are assumed to be type A data.  If N*pct.type.A is not an integer, it will be rounded to the nearest integer.
}
  \item{exp.win.dat}{
A vector of exposure window lengths.  Defaults to the observed window lengths from Lessler et al. (see below).
}
  \item{verb}{
If TRUE, a message with the system time and iteration number will be printed ten times during the simulation run.
}
}

\value{
The simulation.precision functions return a matrix with four columns and nsim rows.  The "ests" column gives the estimated percentiles for the incubation period distribution.  The "SE" column gives the standard error for the estimate.  The "conv" column is 1 if the doubly interval-censored likelihood maximization converged.  Otherwise, it is 0. The "bias" column gives the estimated percentile - true percentile.

The generate.coarse.data function returns a matrix with data suitable for analysis by the dic.fit() function.
}
\references{
Lessler J et al.  New England Journal of Medicine. Outbreak of 2009 Pandemic Influenza A (H1N1) at a New York City School. 2009. 361(27):2628-2636. \url{http://content.nejm.org/cgi/content/full/361/27/2628}

Reich NG et al.  Statistics in Medicine.  Estimating incubation periods with coarse data. 2009.  \url{http://www3.interscience.wiley.com/journal/122507367/abstract}
}
\author{
Nicholas G. Reich
}

\seealso{
\code{\link{dic.fit}}
}


\examples{
a <- precision.simulation(N=100, nsim=10, exp.win.dat=c(1,2), verb=TRUE)
head(a)
summary(a)
}
