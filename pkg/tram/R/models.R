
Coxph <- function(formula, data, subset, weights, offset, cluster, na.action = na.omit, order = 6, 
                  prob = c(.1, .9), ...)
{
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(tram_data)
    td <- eval(mf, parent.frame())

    stopifnot(inherits(td$response, "Surv"))

    ret <- tram(td, order = order, prob = prob, transformation = "smooth", distribution = "MinExtrVal", ...)
    if (!inherits(ret, "mlt")) return(ret)
    ret$call <- match.call(expand.dots = TRUE)
    class(ret) <- c("Coxph", class(ret))
    ret
}

Colr <- function(formula, data, subset, weights, offset, cluster, na.action = na.omit, order = 6, 
                 prob = c(.1, .9), ...)
{
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(tram_data)
    td <- eval(mf, parent.frame())

    stopifnot(is.numeric(td$response))

    ret <- tram(td, order = order, prob = prob, transformation = "smooth", distribution = "Logistic", ...)
    if (!inherits(ret, "mlt")) return(ret)
    ret$call <- match.call(expand.dots = TRUE)
    class(ret) <- c("colr", class(ret))
    ret
}

Polr <- function(formula, data, subset, weights, offset, cluster, na.action = na.omit, ...)
{
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(tram_data)
    td <- eval(mf, parent.frame())

    stopifnot(is.ordered(td$response))

    ret <- tram(td, transformation = "discrete", distribution = "Logistic", negative = TRUE, ...)
    if (!inherits(ret, "mlt")) return(ret)
    ret$call <- match.call(expand.dots = TRUE)
    class(ret) <- c("Polr", class(ret))
    ret
}

