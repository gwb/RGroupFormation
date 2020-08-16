
##suppressMessages(library(purrr))


#' Generates a random permutation
#'
#' If \code{x} is an integer, then returns an element of the symmetric
#' group on \code{x} elements. If \code{x} is a vector, then generates a 
#' permutation of the elements of \code{x}.
#' 
#' @param x An integer or a vector.
#' @return A vector representing a permutation.
#' @examples
#' rp(5)
#' rp(c(1, 4, 6, 10, 2))
#' @export
rp <- function(x) {
    if(length(x) == 1) {
        ## x is an integer
        res <- sample(seq(x))
    } else {
        ## x is a set
        res <- sample(x)
    }
    return(res)
}

count.uniques <- function(v.ls) {
    uniq.v <- unique(v.ls)
    res <- sapply(uniq.v, function(v) sum(v.ls == v))
    names(res) <- uniq.v
    return(res)
}


.rsp <- function(A) {
    stopifnot(all(count.uniques(A) > 1))
    idx.strata.ls <- lapply(unique(A),
                            function(a) which(A==a))

    res <- rep(NA, length(A))
    for(idx.strata in idx.strata.ls) {
        res[idx.strata] <- rp(idx.strata)
    }

    return(res)
}

#'@importFrom purrr array_branch
#' @importFrom purrr map_chr
#' @importFrom purrr %>%
format.attributes <- function(A) {
    if(is.matrix(A)) {
        new.A <- array_branch(A, 1) %>% map_chr(~ paste0(., collapse="."))
    } else if(is.numeric(A)){        
        new.A <- as.character(A)
    } else {
        stop("unsupported format")
    }    
    return(new.A)
}

#' Generates a random stratified permutation
#'
#' \code{rsp} returns an element of the stabilizer of \code{A} in the
#' symmetric group on \code{|A|} elements.
#' 
#' @param A A vector.
#' @return A vector representing a permutation
#' @examples
#' A <- c(0, 0, 0, 0, 1, 1, 1, 1)
#' rsp(A)
#' @export
rsp <- function(A) {
    res <- .rsp(format.attributes(A))
    return(res)
}

#' Group action operator
#'
#' Applies a permutation \code{s} to a vector \code{v}.
#' 
#' @param s A vector representing a permutation.
#' @param v A vector of the same length as \code{s}.
#' @return A vector of the same length as \code{s} and \code{v}.
#' @examples
#' X <- seq(8)
#' A <- c(0, 0, 0, 0, 1, 1, 1, 1)
#' rsp(A) %p% X
#' rp(8) %p% X
#' @export
`%p%` <- function(s, v) {
    stopifnot(length(s) == length(v))
    res <- v[s]
    return(res)
}

inv <- function(s){
    stopifnot(setequal(s, seq_along(s)))
    res <- sort(s, index.return=TRUE)$ix
    return(res)
}


#' Pretty print permutations
#'
#' \code{pprint} shows both the permutationa and its inverse. The column
#' labelled \code{p.i} shows the image of the column labelled \code{i} while
#' the image labelled \code{pm.i} shows the inverse image of \code{i}.
#'
#' @param s A vector representing a permutation.
#' @return A dataframe with 3 columns: \code{i}, \code{p.i}, \code{pm.i}.
#' @export
pprint <- function(s) {
    res <- data.frame(i = seq_along(s), pm.i = s, p.i = inv(s))
    return(res)
}
