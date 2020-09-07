
#' Runs a Conditional Fisher Randomization Test
#'
#' This function tests the null that H_0: Y(k) = Y(l) where `k` and `l` are the
#' exposures of interest. The test statistic `T` should be inputable (see
#' description of argument).
#'
#' @param Zobs A list of numeric vectors. A group-assignment structure.
#' @param Yobs A numeric vector of observed outcomes.
#' @param T A function representing the test statistic. This function
#' is called as followed in the code: `T(W, Yobs.u, k, l, A.u)` where
#' W, Yobs.u and A are the subsets of the exposure vector, observed
#' outcomes vector and attribute vectors (respectively) for the focal
#' units. The statistics `T.dim` and `T.studentized` can be used out
#' of the box.
#' @param A A numeric vector of attributes.
#' @param k An integer representing an exposure value.
#' @param l An integer representing an exposure value.
#' @param n.rand An integer. The number of draws to use to build the
#' (conditional) null distribution of the test statistic (defaults to 1000).
#' @param two.sided A boolean. Whether to report a one-sided or two-sided
#' p-value (defaults to FALSE).
#' @examples
#'
#' # setup
#' N <- 90 # number of units
#' M <- 3 # number of units in each (equal-sized) group
#' A <- rep(c(0,1), each = N/2) # vector of binary attributes
#' L0 <- find.good.L0(N, M, 1000)
#'
#' # generating a schedule of potential outcomes
#' # (for which H0^{0,1} is true, but H0^{0,2} is not).
#' Y0 <- rnorm(N)
#' science = cbind(Y0, Y0, Y0+1)
#' colnames(science) <- NULL
#'
#' # Observing the assignments and outcomes
#' L <- rsL(L0, A, 1000)
#' Zobs <- Z.fn(L)
#' Yobs <- observe.outcomes(Zobs, A, science)
#'
#' # running the tests
#' FRT(Zobs, Yobs, T.dim, A, 1, 0)
#' FRT(Zobs, Yobs, T.dim, A, 2, 0)
#' 
#' @export
FRT <- function(Zobs, Yobs, T, A, k, l, n.rand=1000, two.sided=FALSE) {
    Wobs <- W.fn(Zobs, A)
    N <- length(Wobs)
    
    U <- rep(0,N)
    U[which(Wobs %in% c(k,l))] <- 1
   
    Wobs.u <- Wobs[U==1]
    Yobs.u <- Yobs[U==1]
    A.u <- A[U==1]

    Tobs <- T(Wobs.u, Yobs.u, k, l, A.u)

    W.mat <- rsW(Wobs.u, A.u, n.draws=n.rand)
    T.ls <- apply(W.mat, 1,
                  function(W) T(W, Yobs.u, k, l, A.u))

    if(!two.sided) {
        return(mean(T.ls >= Tobs))
    } else {
        right.pval <- mean(T.ls >= Tobs)
        left.pval <- mean(T.ls <= Tobs)
        res <- min(1, 2*min(right.pval, left.pval))
        return(res)        
    }
}


#' Runs a standard Fisher Randomization Test (for the sharp null hypothesis).
#'
#' @param Zobs A list of numeric vectors representing the group-assignment
#' structure.
#' @param Yobs A numeric vector of observed outcomes.
#' @param T A function representing the test statistic. It is called as
#' `T(W, Yobs)` in the code, where `W` is a vector of exposures, and `Yobs`
#' the vector of observed outcomes. IMPORTANT NOTE `T.dim` and `T.studentized`
#' cannot be used here, but `T.F` can.
#' @param A A numeric vector of attributes.
#' @param n.rand An integer. The number of draws to use to build the
#' null distribution of the test statistic (defaults to 1000).
#' @export
FRT.sharp <- function(Zobs, Yobs, T, A, n.rand=1000) {
    Wobs <- W.fn(Zobs, A)
    Tobs <- T(Wobs, Yobs)
    W.mat <- rsW(Wobs, A, n.draws=n.rand)
    T.ls <- apply(W.mat, 1,
                  function(W) T(W, Yobs))
    return(mean(T.ls >= Tobs))        
}

generate.null.science <- function(Yobs, Wobs, k, l, tau) {
    ## generates science with Y(k) = Y(l) + tau
    science <- matrix(NA, nrow=length(Yobs), ncol=max(k+1,l+1))

    ##
    science[Wobs==k, k+1] <- Yobs[Wobs==k]
    science[Wobs==l, k+1] <- Yobs[Wobs==l] + tau
    science[Wobs==l, l+1] <- Yobs[Wobs==l]
    science[Wobs==k, l+1] <- Yobs[Wobs==k] - tau

    ## returns a science with a bunch of NAs, but the entries
    ## with NAs will never be used, if the rest of the code
    ## is correct.
    return(science)
}

partial_observe_outcomes <- function(science, W, k, l) {
    Yobs <- rep(NA, length(W))
    Yobs[W == k]  <-  science[W == k, k+1]
    Yobs[W == l]  <-  science[W == l, l+1]

    ## Yobs will contain NAs, but the NA entries will never
    ## be used if code is correct
    return(Yobs)
}

FRT.tau <- function(Zobs, Yobs, T, A, k, l, tau=0, n.rand=1000, two.sided=FALSE, return.T.ls=FALSE) {
    Wobs <- W.fn(Zobs, A)
    science.tau <- generate.null.science(Yobs, Wobs, k, l, tau)
    
    N <- length(Wobs)
    
    U <- rep(0,N)
    U[which(Wobs %in% c(k,l))] <- 1

    science.tau.u <- science.tau[U==1,]
    Wobs.u <- Wobs[U==1]
    Yobs.u <- Yobs[U==1]
    A.u <- A[U==1]

    
    Tobs <- T(Wobs.u, Yobs.u, k, l, A.u)
    
    W.mat <- rsW(Wobs.u, A.u, n.draws=n.rand)

    Yobs.u.ls <- lapply(seq(n.rand),
                        function(i) partial_observe_outcomes(science.tau.u,
                                                             W.mat[i,],
                                                             k,
                                                             l))
    T.ls <- sapply(seq(n.rand),
                   function(i) T(W.mat[i,],
                                 Yobs.u.ls[[i]],
                                 k, l,
                                 A.u))
                                 
##    T.ls <- apply(W.mat, 1,
##                  function(W) T(W, Yobs.u, k, l, A.u))

    if(return.T.ls==TRUE) {
        return(T.ls)
    }
    if(!two.sided) {
        return(mean(T.ls >= Tobs))
    } else {
        right.pval <- mean(T.ls >= Tobs)
        left.pval <- mean(T.ls <= Tobs)
        res <- min(1, 2*min(right.pval, left.pval))
        return(res)
    }
}

