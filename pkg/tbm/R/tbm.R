
shiftFamily <- function(object, data, weights) {

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
        estfun(tmp, w = w)[,i]
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

ctmFamily <- function(model, data, weights) {
    mf <- mlt(model, data, weights  = weights)
    theta <- coef(mf)
    offset <- c(theta[1], log(pmax(sqrt(.Machine$double.eps), diff(theta))))
    OM <- matrix(offset, nrow = NROW(data), ncol = length(offset),
                 byrow = TRUE)
    CS <- diag(length(coef(mf)))
    CS[upper.tri(CS)] <- 1
    ngradient <- function(y, f, w) {
        f <- matrix(f, nrow = NROW(y), ncol = length(coef(mf))) + OM
        theta <- cbind(f[,1], exp(f[,-1])) %*% CS
        ret <- - estfun(mf, parm = theta, w = w)
        ### we need unweighted scores!
        ret <- ret * cbind(1, exp(f[,-1]))
        ret
    }
    risk <- function(y, f, w) {
        f <- matrix(f, nrow = NROW(y), ncol = length(coef(mf))) + OM
        theta <- cbind(f[,1], exp(f[,-1])) %*% CS
        -logLik(mf, parm = theta, w = w)
    }
    Family(ngradient = ngradient, risk = risk, 
           offset = function(...) return(0),
           nuisance = function() offset)
}

tbm <- function(model, formula, data = list(), na.action = na.omit, weights = NULL, 
                shiftonly = FALSE, control = boost_control(), oobweights = NULL, 
                baselearner = c("bbs", "bols"), ...) {

    baselearner <- match.arg(baselearner)
    if (shiftonly) {
        tmp <- model$bases
        if (!is.null(tmp$shifting)) {
            tmp$shifting <- c(int = intercept_basis(), shift = tmp$shifting)
        } else {
            tmp$shifting <- intercept_basis()
        }
        td <- model$todistr$name
        td <- switch(td, "normal" = "Normal", "logistic" = "Logistic",
                     "minimum extreme value" = "MinExtrVal")
        model <- ctm(tmp$response, tmp$interacting, tmp$shifting, 
                     todistr = td)
        family <- shiftFamily(model, data, weights)
        class <- "tbm_shift"
    } else {
        stopifnot(is.null(model$model$bases$interacting))
        stopifnot(is.null(model$model$bases$shifting))
        model <- model$model
        family <- ctmFamily(model, data, weights)
        class <- "tbm_ctm"
    }

    mf <- match.call(expand.dots = FALSE)
    mf$model <- NULL
    mf$shiftonly <- NULL
    mf$na.action <- na.action ### evaluate na.action
    if(missing(data)) data <- environment(formula)
    mf$family <- family
    mf$baselearner <- baselearner
    mf[[1L]] <- quote(mboost::mboost)
    ret <- eval(mf, parent.frame())
    ret$ctm <- model
    class(ret) <- c(class, class(ret))
    return(ret) 
}

predict.tbm_ctm <- function(object, newdata = NULL, which = NULL, 
                            coef = FALSE, ...) {

    if (is.null(newdata) && is.null(which)) {
        pr <- fitted(object)
    } else {
        class(object) <- class(object)[-1L]
        if (is.null(newdata)) newdata <- object$data
        pr <- predict(object, newdata, which = which)
    }
    pr <- matrix(pr, ncol = length(coef(object$ctm)))
    CS <- diag(ncol(pr))
    CS[upper.tri(CS)] <- 1
    pr <- cbind(pr[,1], exp(pr[,-1])) %*% CS
    if (coef) return(pr)
    ret <- c()
    tmpm <- object$ctm
    for (i in 1:nrow(pr)) {
        coef(tmpm) <- pr[i,]
        ret <- cbind(ret, predict(tmpm, newdata = data.frame(1), ...))
    }
    ret
}

predict.tbm_shift <- function(object, newdata = NULL, which = NULL, 
                              coef = FALSE, ...) {

    if (is.null(newdata) && is.null(which)) {
        pr <- fitted(object)
    } else {
        class(object) <- class(object)[-1L]
        if (is.null(newdata)) newdata <- object$data
        pr <- predict(object, newdata, which = which)
    }
    if (coef) {
        cf <- nuisance(object)
        ret <- matrix(cf, nrow = NROW(pr), ncol = length(cf), byrow = TRUE)
        return(ret)
    }
    ret <- c()
    tmpm <- object$ctm
    coef(tmpm) <- c(nuisance(object), "(Intercept)" = 0)
    ret <- c()
    for (i in 1:nrow(pr)) {
        coef(tmpm)["(Intercept)"] <- -pr[i]
        ret <- cbind(ret, predict(tmpm, newdata = newdata[i,,drop = FALSE], ...))
    }
    ret
}
