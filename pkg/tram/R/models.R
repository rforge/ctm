
Coxph <- function(formula, data, subset, weights, offset, cluster, na.action = na.omit, order = 6, 
                  prob = c(.1, .9), ...)
{
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(tram_data)
    td <- eval(mf, parent.frame())

    stopifnot(inherits(td$response, "Surv") ||
              inherits(td$response, "response") ||
              is.numeric(td$response))

    ret <- tram(td, order = order, prob = prob, transformation = "smooth", distribution = "MinExtrVal", ...)
    if (!inherits(ret, "mlt")) return(ret)
    ret$call <- match.call(expand.dots = TRUE)
    ret$tram <- paste(ifelse(is.null(td$terms$s), "", "(Stratified)"), 
                      "Parametric Linear Cox' Regression Model")
    class(ret) <- c("Coxph", class(ret))
    ret
}

Survreg <- function(formula, data, subset, weights, offset, cluster, na.action = na.omit, 
                    dist = c("logistic", "gaussian", "weibull", "exponential", 
                             "rayleigh", "loggaussian", "lognormal", "loglogistic"), 
                    scale = 0, ...)
{
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(tram_data)
    td <- eval(mf, parent.frame())

    stopifnot(inherits(td$response, "Surv") ||
              inherits(td$response, "response") ||
              is.numeric(td$response))

    dist <- match.arg(dist)
    distribution <- switch(dist, "logistic" = "Logistic",
                                 "loglogistic" = "Logistic",
                                 "gaussian" = "Normal",
                                 "loggaussian" = "Normal",
                                 "lognormal" = "Normal",
                                 "weibull" = "MinExtrVal",
                                 "exponential" = "MinExtrVal",
                                 "rayleigh" = "MinExtrVal")

    transformation <- switch(dist, "logistic" = "linear",
                                   "loglogistic" = "logarithmic",
                                   "gaussian" = "linear",
                                   "loggaussian" = "logarithmic",
                                   "lognormal" = "logarithmic",
                                   "weibull" = "logarithmic",
                                   "exponential" = "logarithmic",
                                   "rayleigh" = "logarithmic")

    ret <- tram(td, transformation = transformation, 
                distribution = distribution, negative = TRUE, model_only = TRUE, ...)

    rname <- names(td$mf)[1]    
    cfnm <- names(coef(ret))

    scalecf <- grep(rname, cfnm, fixed = TRUE)
    if (length(scalecf) > 0) {
        if (length(grep("Intercept", cfnm[scalecf])) > 0)
            stop("strata contain explicit intercept term")
    }
    if (dist == "exponential") 
        scale <- 1
    if (dist == "rayleigh")
        scale <- 2
    if (scale > 0) {
        fixed <- rep(1 / scale, length(scalecf))
        names(fixed) <- cfnm[scalecf]
        ret <- tram(td, transformation = transformation, 
                    distribution = distribution, negative = TRUE, fixed = fixed, ...)
    } else {
        ret <- tram(td, transformation = transformation, 
                    distribution = distribution, negative = TRUE, ...)
    }
    if (!inherits(ret, "mlt")) return(ret)
    if (length(scalecf) > 0)
        ret$invscale <- 1 / coef(as.mlt(ret))[scalecf]
    ret$call <- match.call(expand.dots = TRUE)

    .simpleCap <- function(x) {
         s <- strsplit(x, " ")[[1]]
         paste(toupper(substring(s, 1, 1)), substring(s, 2),
               sep = "", collapse = " ")
    }

    ret$tram <- paste(ifelse(is.null(td$terms$s), "", "(Stratified)"),
                      .simpleCap(dist), "Linear Regression Model")
    class(ret) <- c("Survreg", class(ret))
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

    stopifnot(inherits(td$response, "Surv") ||
              inherits(td$response, "response") ||
              is.numeric(td$response))

    ret <- tram(td, order = order, prob = prob, transformation = "smooth", distribution = "Logistic", ...)
    if (!inherits(ret, "mlt")) return(ret)
    ret$call <- match.call(expand.dots = TRUE)
    ret$tram <- paste(ifelse(is.null(td$terms$s), "", "(Stratified)"),
                      "Continuous Outcome Logistic Regression")
    class(ret) <- c("colr", class(ret))
    ret
}

Polr <- function(formula, data, subset, weights, offset, cluster, na.action = na.omit, 
                 method = c("logistic", "probit", "cloglog"), ...)
{
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(tram_data)
    td <- eval(mf, parent.frame())

    stopifnot(is.ordered(td$response) || inherits(td$response, "response"))

    method <- match.arg(method)
    distribution <- c("logistic" = "Logistic", "probit" = "Normal", "cloglog" = "MinExtrVal")
    distribution <- distribution[method]

    ret <- tram(td, transformation = "discrete", distribution = distribution, negative = TRUE, ...)
    if (!inherits(ret, "mlt")) return(ret)
    ret$call <- match.call(expand.dots = TRUE)
    if (method != "probit") {
        ret$tram <- paste(ifelse(is.null(td$terms$s), "", "(Stratified)"),
                          "Proportional", ifelse(method == "logistic",
                          "Odds", "Hazards"), "Regression Model")
    } else {
        ret$tram <- paste(ifelse(is.null(td$terms$s), "", "(Stratified)"),
                          "Ordered Probit Regression Model")
    }                   
    class(ret) <- c("Polr", class(ret))
    ret
}

Lm <- function(formula, data, subset, weights, offset, cluster, na.action = na.omit, ...)
{
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(tram_data)
    td <- eval(mf, parent.frame())

    stopifnot(inherits(td$response, "Surv") ||
              inherits(td$response, "response") ||
              is.numeric(td$response))

    ret <- tram(td, transformation = "linear", distribution = "normal", ...)
    if (!inherits(ret, "mlt")) return(ret)
    ret$call <- match.call(expand.dots = TRUE)
    ret$tram <- paste(ifelse(is.null(td$terms$s), "", "(Stratified)"),
                      "Normal Linear Regression Model")
    class(ret) <- c("Lm", class(ret))
    ret
}

BoxCox <- function(formula, data, subset, weights, offset, cluster, na.action = na.omit, ...)
{
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(tram_data)
    td <- eval(mf, parent.frame())

    stopifnot(inherits(td$response, "Surv") ||
              inherits(td$response, "response") ||
              is.numeric(td$response))

    ret <- tram(td, transformation = "smooth", distribution = "normal", ...)
    if (!inherits(ret, "mlt")) return(ret)
    ret$call <- match.call(expand.dots = TRUE)
    ret$tram <- paste(ifelse(is.null(td$terms$s), "", "(Stratified)"),
                      "Transformed Normal Linear Regression Model")
    class(ret) <- c("BoxCox", class(ret))
    ret
}
