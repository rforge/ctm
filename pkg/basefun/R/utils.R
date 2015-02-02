
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

### Box product of two matrices
.boxprod_char <- function(x, y) {
    matrix(paste(x[rep(1:nrow(x), length = nrow(y)), rep(1:ncol(x), ncol(y)), drop = FALSE], 
                 y[rep(1:nrow(y), length = nrow(x)), rep(1:ncol(y), rep(ncol(x), ncol(y))), drop = FALSE], sep = ":"), nrow = nrow(x))
}

### box product of design matices
.box_char <- function(...) {
    args <- list(...)
    ret <- args[[1]]
    if (length(args) == 1) return(ret)
    for (i in 2:length(args))
        ret <- .boxprod_char(ret, args[[i]])
    ret
}



### of linear functions for constraints
.box_ui_ci <- function(...) {
    args <- list(...)
    ret <- args[[1]]
    if (length(args) == 1) return(ret)
    for (i in 2:length(args)) {
        nci <- ncol(args[[i]]$ui)
        ncr <- ncol(ret$ui)
        ret$ui <- rBind(kronecker(Diagonal(nci), ret$ui),
                        kronecker(args[[i]]$ui, Diagonal(ncr)))
        ret$ci <- c(as(kronecker(Diagonal(nci), 
                                 matrix(ret$ci, ncol = 1)) %*% rep(1, nci), 
                       "vector"),
                    as(kronecker(matrix(args[[i]]$ci, ncol = 1), 
                                 Diagonal(ncr)) %*% rep(1, ncr), 
                       "vector"))
    }
    ret$ci <- 
    ret
}
