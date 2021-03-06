% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rand-mechanisms.R
\name{summarize.exposures}
\alias{summarize.exposures}
\title{Counts the number of units receiving each exposure}
\usage{
summarize.exposures(W.mat)
}
\arguments{
\item{W.mat}{A matrix. Each row represents an exposure vector, so if
the matrix is K x N, then K is the number of exposure vectors, and N
the number of units in the population.}
}
\value{
A matrix with the same number of rows as `W.mat', and whose
number of columns is equal to the number of unique exposures in `W.mat'.
}
\description{
For each row of the input matrix `W.mat', computes the number of units
receiving each exposure.
}
\examples{
W.mat <- rbind(c(0, 1, 1, 2, 0, 1), c(1, 1, 0, 0, 0, 0))
summarize.exposures(W.mat)
#      0 1 2
# [1,] 2 3 1
# [2,] 4 2 0
}
