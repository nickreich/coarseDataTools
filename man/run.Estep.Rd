\name{run.Estep}
\alias{run.Estep}
\title{
The expectation step for a specialized EM algorithm}
\description{
The Expectation- or E-step for the EM algorithm implemented in EMforCFR().  Not meant to be called by users.}
\usage{
run.Estep(alpha, full.data, nlag, assumed.nu)
}
\arguments{
  \item{alpha}{Current estimates of the alpha parameters from the GLM model.}
  \item{full.data}{ A matrix of observed data. See description in EMforCFR helpfile.}
  \item{nlag}{ The number of time units for lagged data.  Corresponds to length(assumed.nu).}  
 \item{assumed.nu}{ a vector of probabilities corresponding to the survival distribution, i.e. nu[i]=Pr(surviving i days | fatal case) }
}
\value{ A data matrix with the same format as full.data from the EMforCFR() documentation.
}
\author{
Nicholas G. Reich
}

\seealso{
\code{\link{EMforCFR}}
}
