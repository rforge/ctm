
.tramFamily <- function(object, data, weights) {

    model <- mlt(object, data, fixed = c("(Intercept)" = 0), weights = weights)

    ngradient <- function(y, f, w = 1) {
        if (length(f) == 1) f <- rep(f, nrow(data))
        if (length(w) == 1) w <- rep(w, nrow(data))
        ### F(a(y) + (Intercept) + offset)
        model <<- mlt(object, data = data, weights = w, fixed = c("(Intercept)" = 0),
                      theta = coef(model, fixed = FALSE), offset = -f)
        tmp <- mlt(object, data = data, dofit = FALSE, offset = -f)
        coef(tmp) <- coef(model, fixed = TRUE)
        i <- names(coef(tmp)) == "(Intercept)"
        ### gradient wrt offset == score wrt (Intercept) = 1
        ### We need UNWEIGHTED scores
        iw <- w
        iw[iw > 0] <- 1/(iw[iw > 0])
        estfun(tmp, w = w)[,i] * iw
    }

    risk <- function(y, f, w = 1) {
        if (length(w) == 1) w <- rep(w, nrow(data))
        tmp <- mlt(object, data = data, dofit = FALSE, offset = -f)
        -logLik(tmp, w = w, parm = coef(model, fixed = TRUE))
    }

    mboost::Family(ngradient = ngradient,
           risk = risk,
           offset = function(y, w) 0,
           check_y = function(y) y,
           nuisance = function() return(coef(model, fixed = FALSE)),
           name = "Conditional Transformation Model Family",
           response = function(f) {
               data$f <- f
               predict(model, newdata = data, type = "quantile", prob = .5)
           })
}

.ctmFamily <- function(model, data, weights) {
    mf <- mlt(model, data, weights  = weights)
    theta <- coef(mf)
    offset <- c(theta[1], diff(theta))
    OM <- matrix(offset, nrow = NROW(data), ncol = length(offset),
                 byrow = TRUE)
    CS <- diag(length(coef(mf)))
    CS[upper.tri(CS)] <- 1
    ngradient <- function(y, f, w) {
        ### update to weights w if necessary
        if (!isTRUE(all.equal(w, weights))) {
            mf <<- mlt(model, data, weights  = w)
            theta <<- coef(mf)
            ### start low!
            offset <<- c(theta[1] - 2, diff(theta))
            OM <<- matrix(offset, nrow = NROW(data), ncol = length(offset),
                          byrow = TRUE)
            weights <<- w
        }
        f <- matrix(f, nrow = NROW(y), ncol = length(coef(mf))) + OM
        theta <- f %*% CS
        ret <- - estfun(mf, parm = theta, w = w)
        ### we need unweighted scores!
        iw <- w
        iw[iw > 0] <- 1/(iw[iw > 0])
        ret
    }
    risk <- function(y, f, w) {
        f <- matrix(f, nrow = NROW(y), ncol = length(coef(mf))) + OM
        theta <- f %*% CS
        -logLik(mf, parm = theta, w = w)
    }
    Family(ngradient = ngradient, risk = risk, 
           offset = function(...) return(0),
           nuisance = function() offset)
}

ctmboost <- function(model, formula, data = list(), weights = NULL, ...) {

    ### note: This defines the response and MUST match data
    basedata <- model$data

    mf <- match.call(expand.dots = TRUE)
    mf$model <- NULL
    mf$gradient <- NULL
    if(missing(data)) data <- environment(formula)

    stopifnot(is.null(model$model$bases$interacting))
    stopifnot(is.null(model$model$bases$shifting))
    myctm <- model$model
    mf$family <- .ctmFamily(myctm, basedata, weights)
    mf[[1L]] <- quote(mboost::mboost)
    ret <- eval(mf, parent.frame())
    ret$model <- mlt(myctm, data = basedata, dofit = FALSE)
    class(ret) <- c("ctmboost", "tbm", class(ret))
    return(ret) 
}

drboost <- function(model, formula, data = list(), weights = NULL, ...) {

    mf <- match.call(expand.dots = TRUE)
    mf$baselearner <- "bols"
    mf[[1L]] <- quote(ctmboost)
    ret <- eval(mf, parent.frame())
    return(ret) 
}

wildboost <- function(model, formula, data = list(), weights = NULL, ...) {

    ### note: This defines the response and MUST match data
    basedata <- model$data

    mf <- match.call(expand.dots = TRUE)
    mf$model <- NULL
    if(missing(data)) data <- environment(formula)

    stopifnot(is.null(model$model$bases$interacting))
    stopifnot(is.null(model$model$bases$shifting))
    myctm <- model$model
    mf$family <- .ctmFamily(myctm, basedata, weights)
    mf[[1L]] <- quote(mboost::blackboost)	
    ret <- eval(mf, parent.frame())
    ret$model <- mlt(myctm, data = basedata, dofit = FALSE)
    class(ret) <- c("ctmboost", "tbm", class(ret))
    return(ret) 
}

tramboost <- function(model, formula, data = list(), weights = NULL, 
                      method = quote(mboost::mboost), ...) {

    ### note: This defines the response and MUST match data
    basedata <- model$data

    mf <- match.call(expand.dots = TRUE)
    mf$model <- mf$method <- NULL
    if(missing(data)) data <- environment(formula)

    tmp <- model$model$bases
    if (!is.null(tmp$shifting)) {
        tmp$shifting <- c(int = intercept_basis(), shift = tmp$shifting)
    } else {
        tmp$shifting <- intercept_basis()
    }
    td <- model$model$todistr$name
    td <- switch(td, "normal" = "Normal", "logistic" = "Logistic",
                 "minimum extreme value" = "MinExtrVal")
    myctm <- ctm(tmp$response, tmp$interacting, tmp$shifting, 
                 todistr = td)
    mf$family <- .tramFamily(myctm, basedata, weights)
    mf[[1L]] <- method

    ret <- eval(mf, parent.frame())
    ret$model <- mlt(myctm, data = basedata, dofit = FALSE)
    class(ret) <- c("tramboost", "tbm", class(ret))
    return(ret) 
}


predict.ctmboost <- function(object, newdata = NULL, which = NULL, 
                             coef = FALSE, ...) {

    class(object) <- class(object)[-1L]
    pr <- predict(object, newdata, which = which)

    pr <- matrix(pr, ncol = length(coef(object$model)))
    pr <- t(t(pr) + nuisance(object)) ### this is the OFFSET!!!
    CS <- diag(ncol(pr))
    CS[upper.tri(CS)] <- 1
    pr <- pr %*% CS
    if (coef) return(pr)
    ret <- c()
    tmpm <- object$model$model
    for (i in 1:nrow(pr)) {
        coef(tmpm) <- pr[i,]
        ret <- cbind(ret, predict(tmpm, newdata = data.frame(1), ...))
    }
    ret
}

predict.tramboost <- function(object, newdata = NULL, which = NULL, 
                              coef = FALSE, ...) {

    class(object) <- class(object)[-1L]
    pr <- predict(object, newdata, which = which)

    if (coef) {
        cf <- nuisance(object)
        ret <- matrix(cf, nrow = NROW(pr), ncol = length(cf), byrow = TRUE)
        ret <- cbind(ret, -pr)
        return(ret)
    }
    ret <- c()
    tmpm <- object$model$model
    coef(tmpm) <- c(nuisance(object), "(Intercept)" = 0)
    ret <- c()
    nd <- newdata[, !colnames(newdata) %in% variable.names(tmpm)[1], drop = FALSE]
    for (i in 1:NROW(pr)) {
        coef(tmpm)["(Intercept)"] <- -pr[i]
        ret <- cbind(ret, predict(tmpm, newdata = nd, ...))
    }
    ret
}

coef.tbm <- function(object, newdata = NULL, ...)
    predict(object, newdata = newdata, coef = TRUE, ...)

logLik.tbm <- function(object, newdata = NULL, coef = coef(object, newdata = newdata), 
                       weights = model.weights(object), ...)
    logLik(object$model, parm = coef, newdata = newdata, weights = weights, ...)
