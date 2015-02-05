
### evaluate model.matrix of multiple basis objects
model.matrix.basis <- function(object, data, ...)
    object(data, ...)

### compute predicted values of (multiple) basis objects
predict.basis <- function(object, newdata, coef, ...) {
    X <- model.matrix(object = object, data = newdata, ...)
    return(drop(X %*% coef))
}

nparm <- function(object, data)
    UseMethod("nparm")

nparm.basis <- function(object, data)
    ncol(object(data))

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

generate.basis <- function(object, n) {
    ret <- list()
    .generate <- function(s, n) {
        if (!is.atomic(s)) {
            tmp <- lapply(s, .generate, n = n)
            if (all(sapply(s, is.atomic)))
                ret[names(tmp)] <<- tmp
        }
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
