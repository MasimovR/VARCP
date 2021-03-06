% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{changetest}
\alias{changetest}
\title{Computes P value of estimated test statistic}
\usage{
changetest(VARCP)
}
\arguments{
\item{VARCP}{A list of class "VARCP"}
}
\value{
changetest returns a list of class "htest" with following items
\item{method}{Character string which states type of the input test statistic}
\item{data.name}{Character string which states name of the input argument}
\item{p.value}{P value estimated using asymptotic distribution of the test statistic}
\item{estimate}{Maximal value of the test statistic}
\item{statistic}{Time index of estimated changepoint}
}
\description{
This function allows testing the presence of the change point by computing P value using corresponding asymptotic distribution of test statistic.
}
\examples{
## Simulation of time series of length 200
phi <- matrix(c(0.2,-0.6,0.3,1),2,2)
sig <- matrix(c(4,0.8,0.8,1),2,2)
data <- VARMAsim(200, arlags = 1, phi = phi, sigma = sig)
model <- VAR(data[["series"]], p = 1)

## Tests estimation
CUSUM <- statcomp(model, type = "CUSUM")
DET <- statcomp(model, type = "DET")

## Computing test P values
changetest(CUSUM)
changetest(DET)

}
