
.is.formula <- function(x) {
    if (is.null(x)) return(FALSE)
    inherits(x, "formula")
}

model <- function(response = NULL, interacting = NULL, shifting = NULL) {

    if (.is.formula(response)) 
        response <- as.bases(response)
    if (.is.formula(interacting)) 
        interacting <- as.bases(interacting)
    if (.is.formula(shifting)) 
        shifting <- as.bases(shifting, remove_intercept = TRUE)

    if (is.null(shifting)) {
        if (!is.null(interacting)) {
            ret <- b(bresponse = response, binteracting = interacting)
        } else {
            ret <- c(bresponse = response)
        }   
    } else {
        if (!is.null(interacting)) {
            ret <- c(b(bresponse = response, binteracting = interacting), 
                    bshifting = shifting)
        } else {
            ret <- c(bresponse = response, bshifting = shifting)
        }
    }
    attr(ret, "response") <- varnames(response)
    return(ret)
}

.mlt_fit <- function(model, data, todistr = c("Normal", "Logistic", "MinExtrVal"),
                     weights = NULL, offset = NULL, trunc = c(-Inf, Inf)) {

    if (is.null(weights)) weights <- rep(1, nrow(data))
    if (is.null(offset)) offset <- rep(0, nrow(data))
    stopifnot(nrow(data) == length(weights))
    stopifnot(nrow(data) == length(offset))

    response <- attr(model, "response")
    stopifnot(length(response) == 1)
    y <- data[[response]]
    tmpdata <- data

    by <- model
    distr <- .distr(todistr)

    if (is.null(dim(y)) & (storage.mode(y) == "double")) {
        exact <- rep(TRUE, nrow(data))
        edata <- data
    } else if (is.Surv(y)) {
        sy <- .Surv2matrix(y)
        lna <- is.na(sy[, "left"])
        rna <- is.na(sy[, "right"])
        exact <- lna | rna
        if (any(exact)) {
            edata <- data[exact,]
            ct <- ifelse(all(lna[exact]), "right", "left")
            edata[[response]] <- sy[exact, ct]
        }
        ldata <- data[!exact,]
        ldata[[response]] <- sy[!exact, "left"]
        lfinite <- is.finite(ldata[[response]])
        rdata <- data[!exact,]
        rdata[[response]] <- sy[!exact, "right"]
        rfinite <- is.finite(rdata[[response]])
    } else if (is.ordered(y)) {
        exact <- rep(FALSE, nrow(data))
        ldata <- rdata <- data
        ldata[[response]] <- y
        ldata[[response]][] <- levels(y)[pmax(unclass(y) - 1, 1)]
        lfinite <- (ldata[[response]] != levels(y)[1])
        rdata[[response]] <- y
        rfinite <- rdata[[response]] != rev(levels(y))[1]
    } else if (is.integer(y)) {
        exact <- rep(FALSE, nrow(data))
        ldata <- rdata <- data
        ldata[[response]] <- y - 1
        lfinite <- (ldata[[response]] >= 0)
        rdata[[response]] <- y
        rfinite <- rep(TRUE, nrow(rdata))
    } else {
        stop("cannot deal with response class", class(y))
    }

    nbeta <- NULL
    if (any(exact)) {    
        Y <- model.matrix(by, data = edata)
        nbeta <- ncol(Y)
        Yprime <- model.matrix(by, data = edata, bresponse = list(deriv = 1))
        ### <FIXME> handle deriv better!!!
        Yprime[abs(Yprime - Y) == 0] <- 0 
        ### </FIXME>
        ui <- attr(Y, "constraint")$ui
        ci <- attr(Y, "constraint")$ci
        stopifnot(any(is.finite(ci)))
        ui <- as(ui[is.finite(ci),,drop = FALSE], "matrix")
        ci <- ci[is.finite(ci)]
        ci[ci == 0] <- sqrt(.Machine$double.eps)
    } 

    if (all(exact)) { 
        ll <- function(beta) .mlt_loglik_exact(distr, Y, Yprime, offset)(beta)
        sc <- function(beta) .mlt_score_exact(distr, Y, Yprime, offset)(beta)
    } else {
        .makeY <- function(data, finite, nc = NULL) {
            if (any(finite)) {
                tmp <- model.matrix(by, data = data[finite,,drop = FALSE])
                nc <- ncol(tmp)
            } else {
                tmp <- -Inf
            }
            ret <- matrix(-Inf, nrow = nrow(data), ncol = nc)
            ret[finite,] <- tmp
            ret
        }
        if (!is.null(nbeta)) {
            lY <- .makeY(ldata, lfinite, nbeta)
            rY <- .makeY(rdata, rfinite, nbeta)
        } else {
            if (any(lfinite)) {
                lY <- .makeY(ldata, lfinite)
                rY <- .makeY(rdata, rfinite, ncol(lY))
            } else {
                rY <- .makeY(rdata, rfinite)
                lY <- .makeY(ldata, lfinite, ncol(rY))
            }
        }
        if (all(!exact)) {
            ui <- attr(lY, "constraint")$ui
            ci <- attr(lY, "constraint")$ci
            stopifnot(any(is.finite(ci)))
            ui <- as(ui[is.finite(ci),,drop = FALSE], "matrix")
            ci <- ci[is.finite(ci)]
            ci[ci == 0] <- sqrt(.Machine$double.eps)
        }
        ll <- function(beta) {
            ret <- numeric(nrow(data))
            ret[exact] <- .mlt_loglik_exact(distr, Y, Yprime, offset[exact])(beta)
            ret[!exact] <- .mlt_loglik_interval(distr, lY, rY, offset[!exact])(beta)
            ret
        }
        sc <- function(beta) {
            ret <- matrix(0, nrow = nrow(data), ncol = ncol(Y))
            ret[exact,] <- .mlt_score_exact(distr, Y, Yprime)(beta)
            ret[!exact,] <- .mlt_score_interval(distr, lY, rY)(beta)
            ret
        }
    }


    z <- scale(edata[[response]])
    theta <- coef(lm(z ~ Y - 1))

    ret <- constrOptim(theta = theta, f = function(beta) -sum(weights * ll(beta)),
                           grad = function(beta) -colSums(weights * sc(beta)), ui = ui,
                           ci = ci, hessian = TRUE)
    ret$by <- by
    ret$response <- response
    ret$distr <- distr
    ret$loglik <- ll(ret$par)
    class(ret) <- "mlt"
    return(ret)
}

coef.mlt <- function(object, ...)
    object$par

logLik.mlt <- function(object, ...)
    object$loglik

mlt <- .mlt_fit

predict.mlt <- function(object, newdata = NULL, 
                        ...) {
    ret <- function(y, type = c("trafo", "prob"),...) {
        type <- match.arg(type)
        if (is.null(newdata)) {
            newdata <- data.frame(y)
            colnames(newdata) <- object$response
        } else {
            idx <- 1:nrow(newdata)
            tmp <- expand.grid(y = y, idx = idx)
            newdata <- newdata[tmp$idx,,drop = FALSE]
            newdata[[object$response]] <- tmp$y
        }
        Y <- model.matrix(object$by, data = newdata)
        ret <- drop(Y %*% object$par)
        ret <- switch(type, 
                      "trafo" = ret,
                      "prob" = object$distr$p(ret))
        ret
    }
    ret
}
