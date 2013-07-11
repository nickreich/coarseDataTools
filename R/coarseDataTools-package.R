##' @name exp.win.lengths.rda 
##' @title Exposure window lengths from an influenza outbreak at a NYC school
##' @description A numeric vector of exposure window lengths taken from a dataset of doubly interval-censored incubation period observations.  All observations came from a NYC public school.  The outbreak has been described in full in Lessler et al. (see citation below). 
##' @docType data
##' @format A numeric vector with 134 positive values.  Each value represents an exposure window length from an observation of the incubation period for that individual.  The exposure window length is the length of time during which exposure could have occured.  For example, if an individual could have been exposed anytime between 6am on Monday to 6am on Wednesday, her exposure window length would be 2 days. 
##' 
##' @usage data(exp.win.lengths)
##' @source Lessler J et al.  New England Journal of Medicine. Outbreak of 2009 Pandemic Influenza A (H1N1) at a New York City School. 2009. 361(27):2628-2636. \url{http://content.nejm.org/cgi/content/full/361/27/2628}
##' 
##' @examples
##' data(exp.win.lengthsab)
##' summary(exp.win.lengths)
##' hist(exp.win.lengths)
##' @keywords datasets
NULL

##' @name fluA.inc.per
##' @docType data
##' @title Coarse incubation period data for influenza A
##' @description These observations on the incubation period of influenza A come from a variety of sources, and were gathered for a literature review.  They report doubly interval-censored, single interval-censored or exact observations for the incubation period. 
##' @usage data(fluA.inc.per)
##' @format A data frame with 151 observations on the following 7 variables.
##' \item{\code{author}}{the name of the primary author for the source of the observation}
##' \item{\code{year}}{the year of the study which is the source of the observation}
##' \item{\code{EL}}{the earliest possible time of infection}
##' \item{\code{ER}}{the latest possible time of infection}
##' \item{\code{SL}}{the earliest possible time of symptom onset}
##' \item{\code{SR}}{the latest possible time of symptom onset}
##' \item{\code{type}}{an indicator of the type of observation: 0 for doubly interval-censored, 1 for single-interval censored, 2 for exact}
##'@source Lessler J, Reich NG, Brookmeyer R, Perl TM, Nelson KE, Cummings DAT. (2009) A systematic review of the incubation periods of acute respiratory viral infections. Lancet Infectious Diseases. 9(5):291-300.
##' @examples
##' data(fluA.inc.per)
##' head(fluA.inc.per)
##' @keywords datasets
NULL