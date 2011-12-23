\name{simulated.outbreak.deaths}
\alias{simulated.outbreak.deaths}
\docType{data}
\title{
Simulated case and death reports from a fictional outbreak
}
\description{
This dataset provides reported counts of cases and deaths occuring at different time points across a simulated outbreak. Details of the data simulation algorithm are provided in the manuscript "Estimating case fatality ratios from infectious disease surveillance data" (Reich et al., under review, available upon request).
}
\usage{data(simulated.outbreak.deaths)}
\format{

\describe{
	\item{\code{time}}{time, t,  after start of outbreak}
		
	\item{\code{grp}}{an categorical variable indicating membership in one of two groups of a covariate, j}

	\item{\code{R}}{number of recovered cases reported at the given t and j}

	\item{\code{D}}{number of deaths reported at the given t and j}

	\item{\code{N}}{total number of cases and deaths reported at t and j, or D+R}
	}

}
\references{
Reich NG, Lessler J, Brookmeyer R.  Estimating case fatality ratios from infectious disease surveillance data.  [currently under review but availble from the authors upon request]
}
\examples{
data(simulated.outbreak.deaths)
head(simulated.outbreak.deaths)
plot(simulated.outbreak.deaths[simulated.outbreak.deaths[,"grp"]==1,"D"], type="l")
}
