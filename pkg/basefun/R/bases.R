
### evaluate model.matrix of multiple basis objects
model.matrix.basis <- function(object, data, ...)
    object(data, ...)

model.matrix.log_basis <- function(object, data, deriv = 0L, ...)
    object(data, deriv = deriv, ...)

model.matrix.polynomial_basis <- function(object, data, deriv = 0L, ...)
    object(data, deriv = deriv, ...)

### evaluate model.matrix of Bernstein polynom
model.matrix.Bernstein_basis <- function(object, data,
                                         deriv = 0L, integrate = FALSE, ...)
    model.matrix.basis(object = object, data = data, 
                       deriv = deriv, integrate = integrate, ...)

### compute predicted values of (multiple) basis objects
predict.basis <- function(object, newdata, coef, ...) {
    X <- model.matrix(object = object, data = newdata, ...)
    X %*% coef
}

predict.box_bases <- function(object, newdata, coef, ...) {
    if (is.data.frame(newdata)) 
        return(predict.basis(object = object, newdata = newdata, 
                             coef = coef, ...))
    X <- model.matrix(object = object, data = newdata, model.matrix = FALSE, ...)
    X <- lapply(X, function(x) as(x, "matrix"))
    X$beta <- array(coef, sapply(X, NCOL))
    do.call(".cXb", X)
}

predict.cbind_bases <- function(object, newdata, coef, ...) {
    stopifnot(is.data.frame(newdata)) 
    predict.basis(object = object, newdata = newdata, 
                  coef = coef, ...)
}

varnames <- function(x)
    UseMethod("varnames")

varnames.default <- function(x)
    attr(x, "varnames")

support <- function(x)
    UseMethod("support")

support.default <- function(x)
    attr(x, "support")

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
        } else {
            x <- NULL
        }
        return(x)
    }
    ret2 <- .generate(support(object), n = n)
    if (length(ret) == 0) return(ret2)
    ret
}

as.basis <- function(object, ...)
    UseMethod("as.basis")
