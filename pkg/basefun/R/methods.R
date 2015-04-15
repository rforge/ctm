
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
    return(.const_array(dim, nd[nd %in% varnames(object)], lp))
}

nparm <- function(object, data)
    UseMethod("nparm")

nparm.basis <- function(object, data) {
    if (!is.data.frame(data)) 
        data <- expand.grid(data[varnames(object)])
    ncol(object(data))
}

nparm.bases <- function(object, data)
    sapply(object, nparm, data = data)

varnames <- function(x)
    UseMethod("varnames")

varnames.basis <- function(x)
    attr(x, "varnames")

varnames.bases <- function(x)
    unique(sapply(x, varnames))

support <- function(x)
    UseMethod("support")

support.basis <- function(x)
    attr(x, "support")

support.bases <- function(x)
    lapply(x, support)

generate <- function(object, n)
    UseMethod("generate")

### <FIXME> can we generate only a subset of variables??? </FIXME>
generate.basis <- function(object, n) {
    ret <- list()
    .generate <- function(s, n) {
        if (!is.atomic(s)) {
            tmp <- lapply(s, .generate, n = n)
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
    ret2 <- .generate(support(object), n = n)
    if (length(ret) == 0) return(ret2)
    ret
}

generate.bases <- generate.basis

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
