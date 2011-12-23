\name{run.Mstep}
\alias{run.Mstep}
\title{
The maximization step for a specialized EM algorithm}
\description{
The Maximization- or M-step for the EM algorithm implemented in EMforCFR().  Not meant to be called by users.}
\usage{
run.Mstep(dat)
}
\arguments{
  \item{dat}{
data matrix passed from EMforCFR().
}
}
\value{ A list with two components
  \item{phi }{fitted vector of parameters}
  \item{Var }{variance-covariance matrix from the fitted model}
}
\author{
Nicholas G. Reich}


\seealso{
\code{\link{EMforCFR}}
}