

shiftFamily <- function(object, data, weights) {

    tmp <- object$bases
    if (!is.null(tmp$shifting)) {
        tmp$shifting <- c(int = intercept_basis(), shift = tmp$shifting)
    } else {
        tmp$shifting <- intercept_basis()
    }
    td <- object$todistr$name
    td <- switch(td, "normal" = "Normal", "logistic" = "Logistic",
                 "minimum extreme value" = "MinExtrVal")
    object <- ctm(tmp$response, tmp$interacting, tmp$shifting, 
                  todistr = td)

    model <- mlt(object, data, fixed = c("(Intercept)" = 0), weights = weights)

    ngradient <- function(y, f, w = 1) {
        if (length(f) == 1) f <- rep(f, nrow(data))
        if (length(w) == 1) w <- rep(w, nrow(data))
        ### F(a(y) + (Intercept) + offset)
        model <<- mlt(object, data = data, weights = w, fixed = c("(Intercept)" = 0),
                      theta = coef(model, fixed = FALSE), offset = -f)
        tmp <- mlt(object, data = data, dofit = FALSE, offset = -f)
        i <- names(coef(tmp)) == "(Intercept)"
        ### gradient wrt offset == score wrt (Intercept) = 1
        estfun(tmp, parm = coef(model, fixed = TRUE))[,i]
    }

    risk <- function(y, f, w = 1) {
        if (length(w) == 1) w <- rep(w, nrow(data))
        tmp <- mlt(object, data = data, dofit = FALSE, offset = -f)
        logLik(tmp, w = w, parm = coef(model, fixed = TRUE))
    }

    mboost::Family(ngradient = ngradient,
           risk = risk,
           offset = function(y, w) 0,
           check_y = function(y) y,
           nuisance = function() return(coef(model)),
           name = "Conditional Transformation Model Family",
           response = function(f) {
               data$f <- f
               predict(model, newdata = data, type = "quantile", prob = .5)
           })
}

ctmFamily <- function(model, data, weights) {
    mf <- mlt(model, data, weights  = weights)
    ui <- attr(X <- model.matrix(model, data = data[1:2,,drop = FALSE]), "constraint")$ui
    ci <- attr(X, "constraint")$ci
    ngradient <- function(y, fit, weights) {
        ret <- matrix(0, nrow = NROW(data), ncol = length(coef(mf)))
        for (i in 1:NROW(data)) {
            of <- get(".ofuns", environment(mf$score))(weights, subset = i)
            ### already weighted!
            ret[i, ] <- of$sc(fit[i,])[i,,drop = TRUE] / weights[i]
        }
        attr(ret, "ui") <- ui
        attr(ret, "ci") <- ci
        attr(ret, "fit") <- fit
        attr(ret, "nu") <- get("nu", envir = parent.frame())
        ret
    }
    risk <- function(y, fit, weights) {
        ret <- 0
        for (i in 1:NROW(data)) {
            of <- get(".ofuns", environment(mf$score))(weights, subset = i)
            ### already weighted!
            ret <- ret + of$ll(fit[i,])[i]
        }
        ret
    }
    offset <- function(y, weights) {
        cf <- coef(mf)
        matrix(cf, nrow = NROW(data), ncol = length(cf), byrow = TRUE)
    }
    Family(ngradient = ngradient, risk = risk, offset = offset)
}


trboost <- function(model, formula, data = list(), na.action = na.omit, weights = NULL, 
                    shiftonly = FALSE, control = boost_control(), oobweights = NULL, 
                    baselearner = c("bbs", "bols"), ...) {

    baselearner <- match.arg(baselearner)
    if (shiftonly) {
        family <- shiftFamily(model$model, data, weights)
        mf <- match.call(expand.dots = FALSE)
        mf$model <- NULL
        mf$shiftonly <- NULL
        mf$na.action <- na.action ### evaluate na.action
        if(missing(data)) data <- environment(formula)
        mf$family <- family
        mf$baselearner <- baselearner
        mf[[1L]] <- quote(mboost::mboost)
        return(eval(mf, parent.frame()))
    }

    tmpcontrol <- control
    tmpcontrol$mstop <- 2
    tmpcontrol$trace <- FALSE
    m1 <- mboost(formula = formula, data = data, na.action = na.action, weights =
                 weights, offset = NULL, family = Gaussian(), control = tmpcontrol,
                 oobweights = oobweights, baselearner = baselearner, ...)
    if (is.null(weights)) weights <- rep(1, NROW(fitted(m1)))
#     weights[!complete.cases(m1$data)] <- 0

    ### Now, this is _very_ bad but works
    K <- 1
    df2lambda <- function() NA
    index <- 1
    get_index <- function() NA
    newX <- function() NA

    dpp_constr <- function(weights) {
        w <- weights
        if (!is.null(index))
            w <- tapply(as.double(weights), as.integer(index), sum)
        XtX <- crossprod(X * w, X)
        lambdadf <- df2lambda(X, df = args$df, lambda = args$lambda,
                              dmat = K, weights = w, XtX = XtX)
        lambda <- lambdadf["lambda"]
        XtX <- XtX + lambda * K

        X <- as(X, "matrix")
        XtX <- as(forceSymmetric(XtX), "matrix")

        Xi <- X
        if (!is.null(index)) Xi <- Xi[index,,drop = FALSE]

        mysolve <- function(y) {
            si <- y
            if (!is.null(index))
                si <- apply(y, 2, function(w) tapply(w, index, sum))
            Xs <- crossprod(X, si)
            G <- solve(XtX, Xs)
            ui <- attr(y, "ui")
            lastparm <- attr(y, "fit")
            nu <- attr(y, "nu")
            ci <- attr(y, "ci")
            chk <- all(ui %*% t(nu * X %*% G + lastparm) - ci >= -.Machine$double.eps)
            if (!chk) {
                Xs2 <- as.vector(t(Xs))
                tmp <- vector(mode = "list", length = ncol(y))
                for (i in 1:ncol(y)) tmp[[i]] <- XtX
                XtX2 <- do.call("bdiag", tmp)
                i <- as.vector(t(matrix(1:(nrow(XtX) * ncol(y)), nrow = nrow(XtX))))
                ### make this baselearner monotone
                sv <- solve.QP(XtX2[i,i], Xs2, t(kronecker(Xi, ui)))
                G <- matrix(sv$solution, ncol = ncol(y), byrow = TRUE)
            }
            return(G)
        }

        fit <- function(y) {
            if (!is.null(index)) {
                y <- apply(y * weights, 2, function(w) tapply(w, index, sum))
            } else {
                y <- y * weights
            }
            coef <- mysolve(y)
            ret <- list(model = coef,
                        fitted = function() {
                            ret <- as(X %*% coef, "matrix")
                            if (NCOL(coef) == 1)
                            ret <- as.vector(coef)
                            if (is.null(index)) return(ret)
                            if (is.matrix(ret)) 
                            return(ret[index,,drop = FALSE])
                        })
            class(ret) <- c("bm_lin", "bm")
            ret
        }

        ### <FIXME> check for large n, option?
        hatvalues <- function() {
            ret <- as.matrix(tcrossprod(X %*% solve(XtX), X * w))
            if (is.null(index)) return(ret)
            return(ret[index, index])
        }
        ### </FIXME>

        ### actually used degrees of freedom (trace of hat matrix)
        df <- function() lambdadf

        ### prepare for computing predictions
        predict <- function(bm, newdata = NULL, aggregate = c("sum", "cumsum", "none")) {
            cf <- sapply(bm, coef)
            if (!is.matrix(cf))
                cf <- matrix(cf, nrow = 1)
            if(!is.null(newdata)) {
                index <- NULL
                ## Use sparse data represenation if data set is huge
                ## and a data.frame
                if (is.data.frame(newdata) && nrow(newdata) > options("mboost_indexmin")[[1]]) {
                    index <- get_index(newdata)
                    newdata <- newdata[index[[1]], , drop = FALSE]
                    index <- index[[2]]
                }
                X <- newX(newdata, prediction = TRUE)$X
            }
            aggregate <- match.arg(aggregate)
            pr <- switch(aggregate, "sum" =
                 as(X %*% rowSums(cf), "matrix"),
            "cumsum" = {
                 stop("cumsum not available")  
#                 as(X %*% .Call("R_mcumsum", as(cf, "matrix"),
#                                PACKAGE = "mboost"), "matrix")
            },
            "none" = as(X %*% cf, "matrix"))
            if (is.null(index))
                return(pr[ , , drop = FALSE])
            return(pr[index, ,drop = FALSE])
        }

        ret <- list(fit = fit, hatvalues = hatvalues,
                    predict = predict, df = df,
                    Xnames = colnames(X))
        class(ret) <- c("bl_lin", "bl")
        return(ret)
    }

    blg <- m1$baselearner
    for (i in 1:length(blg)) {
        env <- environment(blg[[i]]$dpp)
        blg[[i]]$dpp <- dpp_constr
        environment(blg[[i]]$dpp) <- env
    }

    family <- ctmFamily(model$model, data, weights)

    ret <- mboost_fit(blg, response = m1$response, weights = weights, family = family,
                      offset = NULL, control = control, oobweights = oobweights, ...)
    if (is.data.frame(data) && nrow(data) == length(m1$response)) 
        ret$rownames <- rownames(data)
    else ret$rownames <- 1:NROW(m1$response)
    ret$call <- match.call()
    class(ret) <- c("trboost", class(ret))
    ret
}

