
### compute predicted values of (multiple) basis objects
predict.basis <- function(object, newdata, coef, 
                          dim = !is.data.frame(newdata), ...) {

    if (isTRUE(dim))
        dim <- sapply(newdata, NROW) 
    else if (is.logical(dim)) 
        dim <- NULL
    
    X <- model.matrix(object = object, data = newdata, dim = dim, ...)
    lp <- c(X %*% coef)
    if (is.null(dim)) return(lp)
    nd <- names(dim)
    ### <FIXME> essentially handle the length(dim) == 2 case
    if (any(nd %in% variable.names(object))) {
        nd <- nd[nd %in% variable.names(object)]
    } else {
        stopifnot(length(dim) == 2)
        nd <- names(dim)[2]
    }
    ### <FIXME>
    return(.const_array(dim, nd, lp))
}

nparm <- function(object, data)
    UseMethod("nparm")

nparm.basis <- function(object, data) {
    if (!is.data.frame(data)) 
        data <- expand.grid(data[variable.names(object)])
    ncol(object(data))
}

nparm.box_bases <- function(object, data)
    prod(sapply(object, nparm, data = data))

nparm.cbind_bases <- function(object, data)
    sapply(object, nparm, data = data)

variable.names.basis <- function(object, ...)
    attr(object, "varnames")

variable.names.bases <- function(object, ...)
    unique(sapply(object, variable.names))

support <- function(x)
    UseMethod("support")

support.basis <- function(x)
    attr(x, "support")

support.bases <- function(x)
    lapply(x, support)

mkgrid <- function(object, n)
    UseMethod("mkgrid")

### <FIXME> can we generate only a subset of variables??? </FIXME>
mkgrid.basis <- function(object, n) {
    ret <- list()
    .mkgrid <- function(s, n) {
        if (!is.atomic(s)) {
            tmp <- lapply(s, .mkgrid, n = n)
            if (all(sapply(s, is.atomic)))
                ret[names(tmp)] <<- tmp
        }
        ### <FIXME> better type check as in mlt </FIXME>
        if ((length(s) == 2) && (storage.mode(s) == "double")) {
            x <- seq(from = s[1], to = s[2], length.out = n)
        } else if (is.factor(s)) {
            x <- s
        } else if (is.integer(s)) {
            x <- s
        } else {
            x <- NULL
        }
        return(x)
    }
    ret2 <- .mkgrid(support(object), n = n)
    if (length(ret) == 0) return(ret2)
    ret
}

mkgrid.bases <- mkgrid.basis

as.basis <- function(object, ...)
    UseMethod("as.basis")

intercept <- function(object, ...)
    UseMethod("intercept")

intercept.default <- function(object, ...)
    attr(object, "intercept")

intercept.box_bases <- function(object, ...)
    any(sapply(object, intercept))

intercept.cbind_bases <- function(object, ...)
    any(sapply(object, intercept))
