% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rand-mechanisms.R
\name{rsL}
\alias{rsL}
\title{Generates a random room assignment from a template L0}
\usage{
rsL(L0, A = NULL, n.draws = 10)
}
\arguments{
\item{L0}{A room assignment vector.}

\item{A}{A vector of attributes. Should be of the same length as `L0`.}

\item{n.draws}{An integer. The number of permutations of `L0` that should be
drawn (defaults to 10).}
}
\value{
If `n.draws` = 1, then a vector is returned, representing a random
permutation of `L0`. If `n.draws > 1`, returns a matrix whose rows are
random permutations of `L0`.
}
\description{
`rsL` draws a random permutation of `L0`. If `A` is null, then the
permutation will be drawn uniformly at random from the symmetric group `S`.
If A is not null, then the permutation will be drawn from `S_A`, the
stabilizer of `A` in the symmetric group.
}
