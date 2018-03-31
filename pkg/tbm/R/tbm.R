
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
        ### update to weights w if necessary
        if (!isTRUE(all.equal(w, weights))) {
            mf <<- mlt(model, data, weights  = w)
            theta <<- coef(mf)
            offset <<- c(theta[1], log(pmax(sqrt(.Machine$double.eps), diff(theta))))
            OM <<- matrix(offset, nrow = NROW(data), ncol = length(offset),
                          byrow = TRUE)
            weights <<- w
        }
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

schwarzboost <- function (formula, data = list(), weights = NULL, na.action = na.pass, 
    offset = NULL, family = Gaussian(), control = boost_control(), 
    oobweights = NULL, tree_controls = partykit::ctree_control(teststat = "quad", 
        testtype = "Teststatistic", mincriterion = 0, maxdepth = 2, 
        saveinfo = FALSE), ...) 
{
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights", "na.action"), 
        names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    response <- model.response(mf)
    weights <- model.weights(mf)
    mf <- mf[, -1, drop = FALSE]
    mf$"(weights)" <- NULL
    bl <- list(bbaum(mf, tree_controls = tree_controls))
    ret <- mboost_fit(bl, response = response, weights = weights, 
        offset = offset, family = family, control = control, 
        oobweights = oobweights, ...)
    ret$call <- cl
    ret$rownames <- rownames(mf)
    class(ret) <- c("blackboost", class(ret))
    ret
}

tbm <- function(model, formula, data = list(), weights = NULL, 
                gradient = c("ctm", "shift"), baselearner = c("bbs", "bols", "bbaum"), ...) {

    baselearner <- match.arg(baselearner)
    gradient <- match.arg(gradient)
    ### note: This defines the response and MUST match data
    basedata <- model$data

    mf <- match.call(expand.dots = FALSE)
    mf$model <- NULL
    mf$gradient <- NULL
    if(missing(data)) data <- environment(formula)

    if (gradient == "shift") {
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
        family <- shiftFamily(myctm, basedata, weights)
        class <- "tbm_shift"
        mf$family <- family
        if (baselearner == "bols") {
            mf$baselearner <- NULL
            mf[[1L]] <- quote(mboost:::glmboost.formula)
        } else if (baselearner == "bbaum") {
            mf$baselearner <- NULL
            mf[[1L]] <- quote(schwarzboost)
        } else {          
            mf$baselearner <- baselearner
            mf[[1L]] <- quote(mboost::mboost)
        }
    } else {
        stopifnot(is.null(model$model$bases$interacting))
        stopifnot(is.null(model$model$bases$shifting))
        myctm <- model$model
        family <- ctmFamily(myctm, basedata, weights)
        class <- "tbm_ctm"
        mf$family <- family
        if (baselearner == "bbaum") {
            mf$baselearner <- NULL
            mf[[1L]] <- quote(schwarzboost)	
        } else {          
            mf$baselearner <- baselearner
            mf[[1L]] <- quote(mboost::mboost)
        }
    }

    ret <- eval(mf, parent.frame())
    ret$model <- mlt(myctm, data = basedata, dofit = FALSE)
    class(ret) <- c(class, "tbm", class(ret))

    if (gradient != "shift") {
    mypredict <- function(newdata = NULL, which = NULL,
                          aggregate = c("sum")) {

        if (mstop == 0) {
            if (length(offset) == 1) {
                if (!is.null(newdata))
                    return(rep(offset, NCOL(newdata)))
                return(rep(offset, NROW(y)))
            } 
            if (!is.null(newdata)) {
                warning("User-specified offset is not a scalar, ",
                        "thus it cannot be used for predictions when ",
                        sQuote("newdata"), " is specified.")
                return(rep(0, NCOL(newdata)))
            }
            return(offset)
        }
        if (!is.null(xselect))
            indx <- ((1:length(xselect)) <= mstop)
        which <- thiswhich(which, usedonly = nw <- is.null(which))
        if (length(which) == 0) return(NULL)

        aggregate <- match.arg(aggregate)

        pfun <- function(w, agg) {
            ix <- xselect == w & indx
            if (!any(ix))
                return(0)
            if (cwlin) w <- 1
            ret <- nu * bl[[w]]$predict(ens[ix],
                   newdata = newdata, aggregate = agg)
            if (agg == "sum") return(ret)
            m <- Matrix(0, nrow = nrow(ret), ncol = sum(indx))
            m[, which(ix)] <- ret
            m
        }

        pr <- do.call("cbind", lapply(which, pfun, agg = "sum"))
        attr(pr, "offset") <- offset
        return(pr)
    }

    environment(mypredict) <- environment(ret$predict)
    ret$predict <- mypredict
    }
    return(ret) 
}

predict.tbm_ctm <- function(object, newdata = NULL, which = NULL, 
                            coef = FALSE, ...) {

    class(object) <- class(object)[-1L]
    pr <- predict(object, newdata, which = which)

    pr <- matrix(pr, ncol = length(coef(object$model)))
    pr <- t(t(pr) + nuisance(object)) ### this is the OFFSET!!!
    CS <- diag(ncol(pr))
    CS[upper.tri(CS)] <- 1
    pr <- cbind(pr[,1], exp(pr[,-1])) %*% CS
    if (coef) return(pr)
    ret <- c()
    tmpm <- object$model$model
    for (i in 1:nrow(pr)) {
        coef(tmpm) <- pr[i,]
        ret <- cbind(ret, predict(tmpm, newdata = data.frame(1), ...))
    }
    ret
}

predict.tbm_shift <- function(object, newdata = NULL, which = NULL, 
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
