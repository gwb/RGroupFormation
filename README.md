
# RGroupFormation

<!-- badges: start -->
<!-- badges: end -->

This R package provides utilities for running randomization tests of peer-effects 
in group-formation experiments. It also includes functions for designing such 
experiments.

## Installation

Install the package from GitHub with:
``` r
# install.packages("devtools")
devtools::install_github("gwb/RGroupFormation")
```

## Example

The following code describes a simulation computing the power of a conditional 
randomization test for different alternatives of the 
form `H1: Y(1) = Y(0) + tau`. It uses the functions `find.good.L0`, 
`summarize.exposures`, `test.FRT.strat`, and `T.dim` from the package.

```r
library(RGroupFormation)

set.seed(123) # for reproducibility

# initialization
N <- 90 # total number of units in the experiment
M <- 3 # numher of units in each (equal-sized) group
A <- rep(c(0, 1), each=N/2) # each unit has a binary attribute
L0 <- find.good.L0(N, M, 1000) # draws 1000 seed room assignments.

# Displays the exposure for each seed assignment.
summarize.exposures(W.star.fn(L0, A))

# Initializes the potential outcomes
# Each entry in tau.10.ls is a value for the effect of having a single neighbor 
# with attribute `A=1`. These values are used to define the alternative hypotheses. 
# tau.20 is the effect of having 2 neighbors with attribute `A=1`.
Y0 <- rnorm(N)
tau.10.ls <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1) 
tau.20 <- 1


# Initializes result vector
res.ls <- vector('list', length=length(tau.10.ls)) 

# Iterates over values of tau.10.ls and computes a p-value for each
for(i in seq(length(tau.10.ls))){
    print(paste0("tau spill = ", tau.10.ls[i], collapse=""))
	
	# constructs a science table
    science.i <- cbind(Y0, Y0 + tau.10.ls[i], Y0 + tau.20) 
    colnames(science.i) <- NULL
	
	# Computes the p-value
    res.ls[[i]] <- test.FRT.strat(L0, science.i, T.dim, A, 1, 0, n.rand=300, n.draws=1000)
    print(paste0("power = ", mean(res.ls[[i]] <= 0.05), collapse=""))
}

# Computes the power!
power.ls <- sapply(res.ls, function(res.i) mean(res.i <= 0.05))

```


