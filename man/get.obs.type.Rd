\name{get.obs.type}
\alias{get.obs.type}

\title{
Generation of data types 
}
\description{
This function makes an educated guess at what type of data (doubly interval-censored, single interval-censored or exact) each observation in a doubly interval-censored dataset is.
}
\usage{
get.obs.type(dat)
}

\arguments{
  \item{dat}{
a matrix with columns named "EL", "ER", "SL", "SR", corresponding to the left (L) and right (R) endpoints of the windows of possible exposure (E) and symptom onset (S). 
}
}
\value{
  \item{type }{A numeric vector with entries corresponding to each row of 'dat'. The values 0, 1, or 2, correspond to doubly interval-censored, single interval-censored or exact observations, respsectively.}
}
\author{
Justin Lessler
}

\examples{
data(fluA.inc.per)
get.obs.type(fluA.inc.per)
}


\keyword{ survival }
