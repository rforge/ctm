### Box product of two matrices
.boxprod <- function(x, y) {
    x[, rep(1:ncol(x), ncol(y)), drop = FALSE] * 
    y[, rep(1:ncol(y), rep(ncol(x), ncol(y))), drop = FALSE]
}

### box product of design matices
.box <- function(...) {
    args <- list(...)
    ret <- args[[1]]
    if (length(args) == 1) return(ret)
    for (i in 2:length(args))
        ret <- .boxprod(ret, args[[i]])
    ret
}

### of linear functions for constraints
.box_ui <- function(...) {
    args <- list(...)
    ret <- args[[1]]
    if (length(args) == 1) return(ret)
    for (i in 2:length(args))
        ret <- rBind(kronecker(Diagonal(ncol(args[[i]])), ret),
                     kronecker(args[[i]], Diagonal(ncol(ret))))
    ret
}

### of constraints
.box_ci <- function(...) {
    args <- list(...)
    ret <- args[[1]]
    if (length(args) == 1) return(ret)
    for (i in 2:length(args))
        ret <- c(rep(ret, length(args[[i]])), 
                 rep(args[[i]], rep(length(ret), length(args[[i]]))))
    ret
}

.rec2flat <- function(x) {

    ret <- c()
    foo <- function(x) {
        a <- sapply(x, is.atomic)
        if (any(a))
            ret <<- c(ret, x[a])
        if (all(a)) {
            return(NULL)
        } else {
            sapply(x[!a], foo)
        }
    }
    foo(x)
    return(ret)
}


