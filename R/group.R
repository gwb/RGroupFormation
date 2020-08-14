
##suppressMessages(library(purrr))




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

rsp <- function(A) {
    res <- .rsp(format.attributes(A))
    return(res)
}

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

pprint <- function(s) {
    res <- data.frame(i = seq_along(s), pm.i = s, p.i = inv(s))
    return(res)
}
