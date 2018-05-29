
.loglinFamily <-  function(model, data, weights) {
    mf <- mlt(model, data, weights  = weights)
    theta <- coef(mf)
    offset <- c(theta[1], log(theta[2]))
    OM <- matrix(offset, nrow = NROW(data), ncol = length(offset),
                 byrow = TRUE)
    nd <- as.data.frame(mkgrid(model, n = 50))
    mm <- model.matrix(model$model, data = nd)
    ngradient <- function(y, f, w) {
        ### update to weights w if necessary
        if (!isTRUE(all.equal(w, weights))) {
            mf <<- mlt(model, data, weights  = w, theta = coef(model))
            theta <<- coef(mf)
            ### start low!
            offset <<- c(theta[1], log(theta[2]))
            OM <<- matrix(offset, nrow = NROW(data), ncol = length(offset),
                          byrow = TRUE)
            weights <<- w
        }
        f <- matrix(f, nrow = NROW(y), ncol = length(coef(mf))) + OM
        theta <- cbind(f[,1], exp(f[,2]))
        ret <- - estfun(mf, parm = theta, w = w) * cbind(1, exp(f[,2]))
        ### we need unweighted scores!
        iw <- w
        iw[iw > 0] <- 1/(iw[iw > 0])
        ret
    }
    risk <- function(y, f, w) {
        f <- matrix(f, nrow = NROW(y), ncol = length(coef(mf))) + OM
        theta <- cbind(f[,1], exp(f[,2]))
        -logLik(mf, parm = theta, w = w)
    }
    Family(ngradient = ngradient, risk = risk, 
           offset = function(...) return(0),
           nuisance = function() offset)
}

.stmFamily <- function(object, data, weights) {

    model <- mlt(object, data, fixed = c("(Intercept)" = 0), weights = weights, 
                 theta = coef(object)[-length(coef(object))])

    tmp0 <- mlt(object, data = data, dofit = FALSE)

    ngradient <- function(y, f, w = 1) {
        if (length(f) == 1) f <- rep(f, nrow(data))
        if (length(w) == 1) w <- rep(w, nrow(data))
        ### F(a(y) + (Intercept) + offset)
        ### <FIXME> this requires a complete setup (model.matrix, etc)
        ### of the model and takes tiiiiime...
        tmp <- mlt(object, data = data, weights = w, fixed = c("(Intercept)" = 0),
                   theta = coef(model, fixed = FALSE), offset = -f)
        if (logLik(tmp) < logLik(model))
            warning("risk increase; decrease stepsize nu")
        model <<- tmp
        ### this is a very bad trick to speed things up
        ### this function generates the score function
        ofuns <- get(".ofuns", envir = environment(tmp0$score))
        ### inject the offset
        assign("offset", -f, envir = environment(ofuns))
        ### compute score wrt to an intercept term constant one
        ret <- drop(-ofuns(w)$sc(coef(model, fixed = TRUE), Xmult = FALSE))
        return(ret)
    }

    risk <- function(y, f, w = 1) {
        if (length(w) == 1) w <- rep(w, nrow(data))
        ### see comments above
        ofuns <- get(".ofuns", envir = environment(tmp0$score))
        assign("offset", -f, envir = environment(ofuns))
        return(-sum(w * ofuns(w)$ll(coef(model, fixed = TRUE))))
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
    mf <- mlt(model, data, weights  = weights, theta = coef(model))
    theta <- coef(mf)
    nd <- as.data.frame(mkgrid(model, n = 50))
    dr <- 1
    names(dr) <- names(nd)
    offset <- c(theta[1], diff(theta))
    OM <- matrix(offset, nrow = NROW(data), ncol = length(offset),
                 byrow = TRUE)
    CS <- diag(length(coef(mf)))
    CS[upper.tri(CS)] <- 1
    ### a very bad trick to speed things up
    ofuns <- get(".ofuns", envir = environment(mf$score))
    score  <- ofuns(weights)$sc
    ngradient <- function(y, f, w) {
        ### update to weights w if necessary
        if (!isTRUE(all.equal(w, weights))) {
            mf <<- mlt(model, data, weights  = w, theta = coef(model))
            theta <<- coef(mf)
            offset <- c(theta[1], diff(theta))
            OM <<- matrix(offset, nrow = NROW(data), ncol = length(offset),
                          byrow = TRUE)
            weights <<- w
            ofuns <- get(".ofuns", envir = environment(mf$score))
            score  <- ofuns(weights)$sc
        }
        f <- matrix(f, nrow = NROW(y), ncol = length(coef(mf))) + OM
        theta <- f %*% CS
        return(tcrossprod(score(theta), CS))
    }
    risk <- function(y, f, w) {
        f <- matrix(f, nrow = NROW(y), ncol = length(coef(mf))) + OM
        theta <- f %*% CS
        return(-sum(w * ofuns(w)$ll(theta)))
    }
    Family(ngradient = ngradient, risk = risk, 
           offset = function(...) return(0),
           nuisance = function() offset)
}

ctmboost <- function(model, formula, data = list(), weights = NULL, 
                     method = quote(mboost::mboost), ...) {

    ### note: This defines the response and MUST match data
    basedata <- model$data

    if (inherits(model, "tram"))
        model <- as.mlt(model)

    mf <- match.call(expand.dots = TRUE)
    mf$model <- NULL
    mf$gradient <- NULL
    mf$method <- NULL
    if(missing(data)) data <- environment(formula)

    stopifnot(is.null(model$model$bases$interacting))
    stopifnot(is.null(model$model$bases$shifting))
    myctm <- model$model
    coef(myctm) <- coef(model)
    if (inherits(myctm$bases$response, "Bernstein_basis") ||
        inherits(myctm$bases$response, "formula_basis")) {
        mf$family <- .ctmFamily(myctm, basedata, weights)
        class <- "ctmboost"
    } else if (inherits(myctm$bases$response, "log_basis") ||
               inherits(myctm$bases$response, "polynomial_basis")) {
        stopifnot(length(coef(as.mlt(model))) == 2)
        mf$family <- .loglinFamily(myctm, basedata, weights)
        class <- "loglinboost"
    } else {
       stop("Cannot deal with model.")
    } 
    mf[[1L]] <- method
    ret <- eval(mf, parent.frame())
    ret$model <- mlt(myctm, data = basedata, dofit = FALSE)
    class(ret) <- c(class, "tbm", class(ret))
    return(ret) 
}

stmboost <- function(model, formula, data = list(), weights = NULL, 
                      method = quote(mboost::mboost), ...) {

    ### note: This defines the response and MUST match data
    basedata <- model$data

    if (inherits(model, "tram"))
        model <- as.mlt(model)

    mf <- match.call(expand.dots = TRUE)
    mf$model <- mf$method <- NULL
    mf$method <- NULL
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
    cf <- coef(myctm)
    cf[] <- 0
    cf[names(coef(model))] <- coef(model)
    coef(myctm) <- cf
    mf$family <- .stmFamily(myctm, basedata, weights)
    mf[[1L]] <- method

    ret <- eval(mf, parent.frame())
    ret$model <- mlt(myctm, data = basedata, dofit = FALSE)
    class(ret) <- c("stmboost", "tbm", class(ret))
    return(ret) 
}

predict.ctmboost <- function(object, newdata = NULL, which = NULL, 
                             coef = FALSE, cleanup = TRUE, ...) {

    class(object) <- class(object)[-1L]
    pr0 <- predict(object, newdata, which = which)

    pr0 <- matrix(pr0, ncol = length(coef(object$model)))
    pr0 <- t(t(pr0) + nuisance(object)) ### this is the OFFSET!!!

    CS <- diag(ncol(pr0))
    CS[upper.tri(CS)] <- 1

    pr <- pr0 %*% CS
    nonmono <- pr0[,-1] < 0

    tmpm <- object$model$model

    ### check if transformation function is non-monotone
    ### and fix
    if (cleanup & any(nonmono)) {
        idx <- which(rowSums(nonmono) > 0)
        q <- as.data.frame(mkgrid(tmpm, n = 100)[1])
        X <- model.matrix(tmpm, data = q)
        for (i in idx) {
            cf <- coef(tmpm)
            cf[] <- pr[i,]
            coef(tmpm) <- cf
            tr <- predict(tmpm, newdata = data.frame(1), q = q[[1]], 
                          type = "trafo")
            tr[!is.finite(tr)] <- sign(tr[!is.finite(tr)]) * 10
            if (any(diff(tr) < 0)) {
                w <- c(0.01, c(.01, 10)[(diff(tr) > 0) + 1L])
                l <- coneproj::qprog(crossprod(X * sqrt(w)), crossprod(X * w, tr),
                                     as(attr(X, "constraint")$ui, "matrix"), 
                                     attr(X, "constraint")$ci, msg = FALSE)
                pr[i,] <- l$theta
            }
        }
    } 

    if (coef) return(pr)
    ret <- c()
    for (i in 1:nrow(pr)) {
        cf <- coef(tmpm)
        cf[] <- pr[i,]
        coef(tmpm) <- cf
        ret <- cbind(ret, predict(tmpm, newdata = data.frame(1), ...))
    }
    ret
}

predict.loglinboost <- function(object, newdata = NULL, which = NULL, 
                                coef = FALSE, ...) {

    class(object) <- class(object)[-1L]
    pr <- predict(object, newdata, which = which)

    pr <- matrix(pr, ncol = length(coef(object$model)))
    pr <- t(t(pr) + nuisance(object)) ### this is the OFFSET!!!
    pr <- cbind(pr[,1], exp(pr[,2]))

    if (coef) return(pr)
    ret <- c()
    tmpm <- object$model$model
    for (i in 1:nrow(pr)) {
        cf <- coef(tmpm)
        cf[] <- pr[i,]
        coef(tmpm) <- cf
        ret <- cbind(ret, predict(tmpm, newdata = data.frame(1), ...))
    }
    ret
}

predict.stmboost <- function(object, newdata = NULL, which = NULL, 
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
    cf <- coef(tmpm)
    cf[] <- c(nuisance(object), "(Intercept)" = 0)
    coef(tmpm) <- cf
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

logLik.tbm <- function(object, parm = coef(object, newdata = newdata), 
                       w = NULL, newdata = NULL, ...)
{
    parm <- parm ### evaluate here
    if (is.null(newdata)) newdata <- object$model$data
    logLik(object$model, parm = parm, w = w, newdata = newdata, ...)
}

