##source("group.R")
##suppressMessages(library(pbapply))

## TODO: current code requires A to be integer and contiguous I think. Expand
## to arbitrary attribute
##
##


### utilities

Z.fn <- function(L) {
    res <- lapply(seq_along(L),
                  function(i) setdiff(which(L == L[i]), i))
    return(res)
}

W.fn <- function(Z, A) {
    res <- sapply(seq_along(A),
                  function(i) sum(A[Z[[i]]]))
}

.W.star.fn <- function(L, A) {
    Z.tmp <- Z.fn(L)
    return(W.fn(Z.tmp, A))
}

W.star.fn <- function(L, A) {
    if(is.null(nrow(L))){
        res <- .W.star.fn(L, A)
    } else {
        res <- t(apply(L, 1, function(L.i) .W.star.fn(L.i,A)))
    }
    return(res)
}

## Randomization in latent scale



## Finding good L0

#' Generates one (or more) random vectors of room assignments
#'
#' Generates a room assignment vector (or a room assignment matrix)
#' depending on `n.draws`. For a vector, the entry at index `i`
#' represents the room to which unit `i` is to be assigned.
#' 
#' @param N An integer. The total number of units in the population.
#' @param M An integer. The number of unit per room.
#' @param n.draws An integer. The number of room assignments vectors to
#' return (defaults to 10).
#' @return If `n.draws = 1` then returns a vector of room assignments.
#' If `n.draws > 1', returns a matrix where each row represents a room
#' assignment vector.
#' @export
rL <- function(N, M, n.draws=10) {
    L0 <- rep(seq(N / M), each = M)
    res <- t(replicate(n.draws, rp(N) %p% L0))
    if(n.draws == 1) {
        res <- as.vector(res)
    }
    return(res)
}


count.matches <- function(v.ls, elts.ls) {
    res <- sapply(elts.ls, function(i) sum(v.ls == i))
    names(res) <- elts.ls
    return(res)
}

#' Counts the number of units receiving each exposure
#'
#' For each row of the input matrix `W.mat', computes the number of units
#' receiving each exposure.
#' 
#' @param W.mat A matrix. Each row represents an exposure vector, so if
#' the matrix is K x N, then K is the number of exposure vectors, and N
#' the number of units in the population.
#' @return A matrix with the same number of rows as `W.mat', and whose
#' number of columns is equal to the number of unique exposures in `W.mat'.
#' @examples
#' W.mat <- rbind(c(0, 1, 1, 2, 0, 1), c(1, 1, 0, 0, 0, 0))
#' summarize.exposures(W.mat)
#' #      0 1 2
#' # [1,] 2 3 1
#' # [2,] 4 2 0
#' @export
summarize.exposures <- function(W.mat) {
    if(is.vector(W.mat)) {
        W.mat <- matrix(W.mat, nrow=1)
    }
    exp.ls <- sort(unique(as.vector(W.mat)))
    res <- apply(W.mat, 1,
                 function(W.mat.row) count.matches(as.vector(W.mat.row), exp.ls))
    return(t(res))
}

#' Finds an initial room-assignment vector
#'
#' Here, a `good' room-assignment vector is one whose corresponding exposure
#' vector is approximately balanced. The function generates `n.draws` potential
#' room assignment vectors and selects the one that is most balanced.
#'
#' @param N An integer. The number of units in the population.
#' @param M An integer. The number of units per room. `M` must be a divisor of
#' `N`.
#' @param A A numeric vector. The vector of attributes (must be of length `N`).
#' @param n.draws An integer. The number of potential room assignment vector to
#' draw, and from which to choose the best one.
#' @return A numeric vector, representing the `best` room assignment vector among
#' those randomly generated.
#' @export
find.good.L0 <- function(N, M, A, n.draws=10) {
    Ls <- rL(N, M, n.draws)
    Ws <- W.star.fn(Ls, A)
    exp.mat <- summarize.exposures(Ws)
    mse.ls <- apply(exp.mat, 1, var)
    idx.min <- which.min(mse.ls)
    return(Ls[idx.min,])
}

find.target.L0 <- function(N, M, A, target.frac=0.1, n.draws=10) {
    Ls <- rL(N, M, n.draws)
    Ws <- W.star.fn(Ls, A)
    exp.mat <- summarize.exposures(Ws)
    nexp <- ncol(exp.mat)
    ##    mse.ls <- apply(exp.mat, 1,
    ##                    function(x) var(x[-nexp]) + (x[nexp] - target.frac*N)^2)
    ##mse.ls <- apply(exp.mat, 1,
    ##                  function(x) 0.01* var(x[-nexp]) + (x[nexp] - target.frac*N)^2)
    mse.ls <- apply(exp.mat, 1,
                    function(x) (x[nexp] - target.frac*N)^2)
    ##res <- cbind(exp.mat, mse.ls)
    idx.min <- which.min(mse.ls)    
    return(Ls[idx.min,])
}


## T stats

#' Computes the difference-in-means test statistic
#'
#' \code{T.dim} computes the average difference in outcomes between the
#' units assigned to treatment \code{k} and those assigned to treatment
#' \code{l}.
#'
#' The value of \code{A} is ignored. The parameter is here to allow a unified
#' treatment of test statistics by other functions.
#' 
#' @param W A vector representing the assignment of units to treatments
#' (or exposures).
#' @param Y A vector representing the vector of observed outcomes. Must be
#' of the same length as \code{W}.
#' @param k An integer. One of the treatments to be contrasted.
#' @param l An integer. The other treatment to be contrasted.
#' @param A This parameter is ignored.
#' @return A scalar.
#' @export
T.dim <- function(W, Y, k, l, A=NULL) {
    return(mean(Y[W==k]) - mean(Y[W==l]))
}

#' Compute the stratified studentized test statistic.
#'
#' \code{T.studentizedd} computes the studentized difference in outcomes
#' between untis assigned to treatment \code{k} and those assigned to
#' treatment \code{l}, stratifying by values in \code{A}.
#'
#' The test statistic is described in Basse et al. (2020+).
#'
#' @param W A vector representing the assignment of units to treatments
#' (or exposures).
#' @param Y A vector representing the vector of observed outcomes. Must be
#' of the same length as \code{W}.
#' @param k An integer. One of the treatments to be contrasted.
#' @param l An integer. The other treatment to be contrasted.
#' @param A A vector representing the values to stratify on. Must be of
#' the same length as \code{A}.
#' @return A scalar.
#' @export
T.studentized <- function(W, Y, k, l, A) {
    A.all <- unique(A)
    pa.ls <- count.matches(A, A.all) / length(A)
    num.ls <- vector("numeric", length=length(pa.ls))
    denom.ls <- vector("numeric", length=length(pa.ls))
    for(i in seq_along(A.all)) {
        a <- A.all[i]
        pa <- pa.ls[as.character(a)]

        ## numerator
        Y.hat.a.k <- mean(Y[W==k & A==a])
        Y.hat.a.l <- mean(Y[W==l & A==a])
        num.ls[i] <- pa * (Y.hat.a.k - Y.hat.a.l)

        ## denominator
        n.a.k <- sum(W==k & A==a)
        n.a.l <- sum(W==l & A==a)
        S.hat.a.k <- sum( (Y[W==k & A==a] - Y.hat.a.k)^2 ) / (n.a.k - 1)
        S.hat.a.l <- sum( (Y[W==l & A==a] - Y.hat.a.l)^2 ) / (n.a.l - 1)
        denom.ls[i] <- pa^2 * (S.hat.a.k / n.a.k + S.hat.a.l / n.a.l)      
    }
    res <- sum(num.ls) / sqrt(sum(denom.ls))
    return(res)
}

## FIX: need to weight by group size. In our simulation
## it didn't matter because all groups had equal size, but
## we shouldn't assume this.
##T.F <- function(W, Y) {
##    
##    num <- sum(sapply(unique(W),
##                  function(w) var(Y[W==w])))
##
##    return(var(Y) / num)
##}
##

##T.F <- function(W, Y) {
##    K <- length(unique(W))
##    df.num <- K-1
##    df.denom <- length(W) - K
##    m <- mean(Y)
##    m.i <- tapply(Y, W, mean)
##    n.i <- tapply(Y, W, length)
##    num <- sum(n.i * (m.i - m)^2)
##    denom <- sum(tapply(Y, W, var) * (n.i-1))
##    return( (num / df.num) / (denom / df.denom))
##}

T.F <- function(W, Y) {
    K <- length(unique(W))
    df.num <- K-1
    df.denom <- length(W) - K
    uW <- sort(unique(W))
    m <- mean(Y)
    m.i <- sapply(uW, function(wi) mean(Y[W==wi]))
    n.i <- sapply(uW, function(wi) sum(W==wi))
    num <- sum(n.i * (m.i - m)^2)
    v.i <- sapply(uW, function(wi) var(Y[W==wi]))
    denom <- sum(v.i * (n.i-1))
    return( (num / df.num) / (denom / df.denom))
}

##T.F.new <- function(W, Y) {
##    F.val <- summary(aov(Y ~ W))[[1]][["F value"]][1]
##    return(F.val)
##}
##

T.F.2 <- function(W, Y, k, l, A=NULL) {
    Yk <- Y[W==k]
    Yl <- Y[W==l]
    nk <- length(Yk)
    nl <- length(Yl)
    mk <- mean(Yk)
    ml <- mean(Yl)
    m <- (sum(Yk) + sum(Yl)) / (nk + nl)
    num <- nk * (mk - m)^2 + nl * (ml - m)^2
    denom <- sum((Yk - mk)^2) + sum((Yl - ml)^2)
    res <- num / (denom / (nk + nl - 2))
    return(res)
}
##
##T.F.2 <- function(W, Y, k, l, A=NULL) {
##    num <- var(Y[W==k]) + var(Y[W==l])
##    return(var(Y) / num)
##}

## Observe outcomes

#' Reveals the observed outcomes.
#'
#' @param Z A list of numeric vectors. A group-assignment structure.
#' @param A A numeric vector. A vector of attributes.
#' @param science A matrix. A full schedule of potential outcomes. Its number
#' of rows should equal the length of `A` and `Z`. Each column corresponds to
#' a different exposure. For now, exposures should be contiguous integers
#' starting from 0.
#' @return A vector of observed outcomes, of the same length as `A` and `Z`.
#' @export
observe.outcomes <- function(Z, A, science) {
    ## !!! assumes exposures of the form 0, 1, 2, ...
    W <- W.fn(Z, A)
    N <- length(W)
    Y <- sapply(seq(N), function(i) science[i, W[i] + 1])
    return(Y)
}

## permutation distributions

#' Generates a random room assignment
#'
#' `rsL` draws a random permutation of `L0`. If `A` is null, then the
#' permutation will be drawn uniformly at random from the symmetric group `S`.
#' If A is not null, then the permutation will be drawn from `S_A`, the
#' stabilizer of `A` in the symmetric group.
#'
#' @param L0 A room assignment vector.
#' @param A A vector of attributes. Should be of the same length as `L0`.
#' @param n.draws An integer. The number of permutations of `L0` that should be
#' drawn.
#' @return If `n.draws` = 1, then a vector is returned, representing a random
#' permutation of `L0`. If `n.draws > 1`, returns a matrix whose rows are
#' random permutations of `L0`.
#' @export
rsL <- function(L0, A=NULL, n.draws=10) {
    N <- length(L0)
    if(is.null(A)) {
        res <- t(replicate(n.draws, rp(N) %p% L0))
    } else {
        res <- t(replicate(n.draws, rsp(A) %p% L0))
    }
    if(n.draws==1) {
        res <- as.vector(res)
    }
    return(res)
}


rsW <- function(W0, A=NULL, n.draws=10) {
    N <- length(W0)
    if(is.null(A)) {
        res <- t(replicate(n.draws, rp(N) %p% W0))
    } else {
        res <- t(replicate(n.draws, rsp(A) %p% W0))
    }
    return(res)
}


## randomization tests

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

FRT.sharp <- function(Zobs, Yobs, T, A, n.rand=1000) {
    Wobs <- W.fn(Zobs, A)
    Tobs <- T(Wobs, Yobs)
    W.mat <- rsW(Wobs, A, n.draws=n.rand)
    T.ls <- apply(W.mat, 1,
                  function(W) T(W, Yobs))
    return(mean(T.ls >= Tobs))        
}

###
## Hodges-Lehman
###

## A few things to fix for package:
## - currently assumes that est_HL is inside the interval
##   [lb, ub] which I'm not sure is always the case
## - for this to work properly, there are some implicit
##   "monotony" assumptions about how increasing tau
##    affects E[T]. This may need to be investigated.
##

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


hl_find_nonreject <- function(Zobs, Yobs, T, A, k, l, n.rand=300, step.size=0.1, max.iter=100) {
    pval <- FRT(Zobs, Yobs, T, A, k, l, n.rand, two.sided=TRUE)
    Wobs <- W.fn(Zobs, A)
    counter <- 0
    tau <- 0
    
    while(pval <= 0.05 && counter <= max.iter) {
        counter <- counter+1
        tau.up <- step.size * counter
        tau.down <- - step.size * counter
        pval.up <- FRT.tau(Zobs, Yobs, T, A, k, l, tau=tau.up, n.rand=n.rand, two.sided=TRUE)
        pval.down <- FRT.tau(Zobs, Yobs, T, A, k, l, tau=tau.down, n.rand=n.rand, two.sided=TRUE)

        if(pval.up > 0.05) {
            pval <- pval.up
            tau <- tau.up
        } else {
            pval <- pval.down
            tau <- tau.down
        }
    }

    if(counter >= max.iter) {
        stop("max iterations reached.. now solutions found")
    } else {
        return(tau)
    }
}

hl_find_ub <- function(Zobs, Yobs, T, A, k, l, n.rand=300, init.step.size=0.1,
                       init.max.iter=100, precision=0.1, debug=FALSE) {
    tau.init <- hl_find_nonreject(Zobs, Yobs, T, A, k, l, n.rand, step.size=init.step.size,
                                  max.iter=init.max.iter)

    if(debug) { print(paste0("Initial Tau = ", tau.init)) }
    
    last_try_wrong <- Inf
    current_ub <- tau.init
    try_ub <- tau.init + abs(tau.init)

    counter <- 1
    while(last_try_wrong - current_ub > precision & counter < 20) {
        counter <- counter + 1
        pval_try_ub <- FRT.tau(Zobs, Yobs, T, A, k, l, tau=try_ub, n.rand=n.rand, two.sided=TRUE)
        if(pval_try_ub <= 0.05) {
            last_try_wrong <- try_ub
            try_ub <- (current_ub + try_ub) / 2
        } else {
            if(is.infinite(last_try_wrong)) {
                former_ub <- current_ub
                current_ub <- try_ub
                try_ub <- 2 * former_ub
            } else {
                current_ub <- try_ub
                try_ub <- current_ub + (last_try_wrong - current_ub) / 2
            }
        }
        if(debug) { print(paste0("Current Upper Bound = ", current_ub)) }
    }
    if(counter >= 20) {
        stop("reached max iter in hl_find_ub")
    }
    return(current_ub)
}

hl_find_lb <- function(Zobs, Yobs, T, A, k, l, n.rand=300, init.step.size=0.1,
                       init.max.iter=100, precision=0.1, debug=FALSE) {

    tau.init <- hl_find_nonreject(Zobs, Yobs, T, A, k, l, n.rand, step.size=init.step.size,
                                  max.iter=init.max.iter)

    if(debug) { print(paste0("Initial Tau = ", tau.init)) }
    
    last_try_wrong <- -Inf
    current_lb <- tau.init
    try_lb <- tau.init - abs(tau.init)

    counter <- 1
    while(abs(last_try_wrong - current_lb) > precision & counter < 20) {
        counter <- counter + 1
        pval_try_lb <- FRT.tau(Zobs, Yobs, T, A, k, l, tau=try_lb, n.rand=n.rand, two.sided=TRUE)

        if(pval_try_lb <= 0.05) {
            last_try_wrong <- try_lb
            try_lb <- (current_lb + try_lb) / 2
        } else {
            if(is.infinite(last_try_wrong)) {
                former_lb <- current_lb
                current_lb <- try_lb
                try_lb <- former_lb - abs(former_lb)
            } else {
                current_lb <- try_lb
                try_lb <- current_lb - abs(last_try_wrong - current_lb) / 2
            }
        }
        if(debug) { print(paste0("Current Lower Bound = ",
                                 current_lb,
                                 " | Try Lower Bound = ",
                                 try_lb)) }
    }
    if(counter >= 20) {
        stop("reached max iter in hl_find_lb")
    }
    return(current_lb)
}

hl_find_pe <- function(Zobs, Yobs, T, A, k, l, lb, ub, 
                       n.rand=300, precision=0.1, debug=FALSE) {
    Wobs <- W.fn(Zobs, A)

    U <- rep(0,N)
    U[which(Wobs %in% c(k,l))] <- 1

    Wobs.u <- Wobs[U==1]
    Yobs.u <- Yobs[U==1]
    A.u <- A[U==1]
    
    Tobs <- T(Wobs.u, Yobs.u, k, l, A.u)
    tau.try <- (ub + lb) / 2
    T.mean <- mean(FRT.tau(Zobs, Yobs, T, A, k, l, tau=tau.try, n.rand, return.T.ls=TRUE))
    
    current.gap <- Tobs - T.mean

    if(debug) {print(paste0("Init Tau = ", tau.try, " | gap = ", current.gap))}

    previous.ub <- ub
    previous.lb <- lb
    counter <- 1
    while(abs(current.gap) > precision & counter < 20) {
        counter <- counter + 1
        if(current.gap > 0) {
            previous.lb <- tau.try
            tau.try <- tau.try + (previous.ub - tau.try) / 2
        } else {
            previous.ub <- tau.try
            tau.try <- tau.try - (tau.try-previous.lb) / 2
        }        
        T.mean <- mean(FRT.tau(Zobs, Yobs, T, A, k, l, tau=tau.try, n.rand, return.T.ls=TRUE))
        current.gap <- Tobs - T.mean
        if(debug) {print(paste0("Tau = ", tau.try, " | gap = ", current.gap))}        
    }
    if(counter >= 20) {
        stop("reached max iter in hl_find_ub")
    }
    return(tau.try)        
}

hl_inference <- function(Zobs, Yobs, T, A, k, l,
                         n.rand=300,
                         precision=0.01,
                         init.step.size=0.1,
                         init.max.iter=100,
                         debug=FALSE) {
    
    ub0 <- hl_find_ub(Zobs, Yobs, T, A, 0, 1, n.rand,
                      init.step.size, init.max.iter, precision, debug)
    lb0 <- hl_find_lb(Zobs, Yobs, T, A, 0, 1, n.rand,
                      init.step.size, init.max.iter, precision, debug)
    est0 <- hl_find_pe(Zobs, Yobs, T.dim, A, 0, 1, lb0, ub0, n.rand, precision, debug)

    return(c(lb0, est0, ub0))
}

test_hl_inference <- function(L0, science, T, A, k, l,
                              n.rand=300,
                              precision=0.01,
                              init.step.size=0.1,
                              init.max.iter=100,
                              n.draws=10) {
    Ls <- rsL(L0, A, n.draws)
    Zobs.ls <- apply(Ls, 1, Z.fn)
    Yobs.ls <- lapply(seq(n.draws),
                      function(i) observe.outcomes(Zobs.ls[[i]], A, science))

    res <- vector('list', length=n.draws)
    hl_inference_safely <- safely(hl_inference)
#    for(i in seq_len(n.draws)) {
#        res[[i]] <- hl_inference(Zobs.ls[[i]], Yobs.ls[[i]], T, A, k, l,
#                                 n.rand, precision, init.step.size, init.max.iter)
                                        #    }
    res <- pbapply::pblapply(seq(n.draws),
                             function(i) hl_inference_safely(Zobs.ls[[i]],
                                                             Yobs.ls[[i]],
                                                             T, A, k, l,
                                                             n.rand,
                                                             precision, 
                                                             init.step.size,
                                                             init.max.iter,
                                                             debug=FALSE))
    
    return(res)    
}
                              



## FRT simulation functions

test.FRT.strat <- function(L0, science, T, A, k, l, n.rand=1000, n.draws=100, pb=TRUE) {
    Ls <- rsL(L0, A, n.draws)
    Zobs.ls <- apply(Ls, 1, Z.fn)
    Yobs.ls <- lapply(seq(n.draws),
                      function(i) observe.outcomes(Zobs.ls[[i]], A, science))
    if(pb){
        ## shows progress bar
        res <- pbapply::pbsapply(seq(n.draws),
                        function(i) FRT(Zobs.ls[[i]], Yobs.ls[[i]], T, A, k, l, n.rand))
    } else {
        res <- sapply(seq(n.draws),
                        function(i) FRT(Zobs.ls[[i]], Yobs.ls[[i]], T, A, k, l, n.rand))
    }
    return(res)
}

test.FRT.sharp <- function(L0, science, T, A, n.rand=1000, n.draws=100, pb=TRUE) {
    Ls <- rsL(L0, A, n.draws)
    Zobs.ls <- apply(Ls, 1, Z.fn)
    Yobs.ls <- lapply(seq(n.draws),
                      function(i) observe.outcomes(Zobs.ls[[i]], A, science))
    if(pb){
        ## shows progress bar
        res <- pbapply::pbsapply(seq(n.draws),
                        function(i) FRT.sharp(Zobs.ls[[i]], Yobs.ls[[i]], T, A, n.rand))
    } else {
        res <- sapply(seq(n.draws),
                      function(i) FRT(Zobs.ls[[i]], Yobs.ls[[i]], T, A, n.rand))
    }
    return(res)           
}


## Rejection sampling

rejection.sample.W <- function(rfn, A, U, k, l){
    valid <- FALSE
    while(!valid) {
        L <- rfn()
        W <- W.star.fn(L, A)
        U.new <- ifelse(W %in% c(k, l), 1, 0)
        valid <- all(U == U.new)
        ##valid <- all(W[U==1] %in% c(k,l))
    }
    return(list(L,W))
}

get.rfn.1 <- function(L0, A) {
    fn <- function() {
        return(rsL(L0, A, 1))
    }
    return(fn)
}

get.rfn.2 <- function(N, M) {
    fn <- function() {
        return(rL(N, M, 1))
    }
    return(fn)
}




