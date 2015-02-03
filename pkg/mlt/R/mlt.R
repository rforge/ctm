
.mlt_fit <- function(model, data, todistr = c("Normal", "Logistic", "MinExtrVal"),
                     weights = NULL, offset = NULL, fixed = NULL, trunc = NULL) {

    if (is.null(weights)) weights <- rep(1, nrow(data))
    if (is.null(offset)) offset <- rep(0, nrow(data))
    stopifnot(nrow(data) == length(weights))
    stopifnot(nrow(data) == length(offset))

    response <- attr(model, "response")
    stopifnot(length(response) == 1)
    y <- data[[response]]

    if (is.character(todistr))
        todistr <- .distr(todistr)

    if (!is.null(trunc)) {
        stopifnot(is.character(trunc))
        stopifnot(all(names(trunc) %in% c("left", "right")))
        lefttrunc <- righttrunc <- NULL
        if ("left" %in% names(trunc))
            lefttrunc <- data[trunc["left"]]
        if ("right" %in% names(trunc))
            righttrunc <- data[trunc["right"]]
        trunc <- list(left = lefttrunc, right = righttrunc)
    }

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
        if ("lefttrunc" %in% colnames(sy))
            trunc$left <- sy[, "lefttrunc"]
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

    Y <- NULL
    if (any(exact)) {    
        Y <- model.matrix(model, data = edata)
        deriv <- 1
        names(deriv) <- response
        Yprime <- model.matrix(model, data = edata, deriv = deriv)
        ui <- attr(Y, "constraint")$ui
        ci <- attr(Y, "constraint")$ci
    } 

    if (!is.null(trunc)) {
        if (!is.null(trunc$left)) {
            ltdata <- data
            ltdata[[response]] <- trunc$left
            trunc$left <- model.matrix(model, data = ltdata)
        }
        if (!is.null(trunc$right)) {
            rtdata <- data
            rtdata[[response]] <- trunc$right
            trunc$right <- model.matrix(model, data = rtdata)
        }
    }

    if (all(exact)) { 
        if (!is.null(fixed)) {
            stopifnot(all(names(fixed) %in% colnames(Y)))
            fix <- colnames(Y) %in% names(fixed)
            ui <- ui[,!fix,drop = FALSE]
            .parm <- function(beta) {
                ret <- numeric(ncol(Y))
                ret[fix] <- fixed
                ret[!fix] <- beta
                ret
            }
        } else {
            .parm <- function(beta) beta
            fix <- rep(FALSE, ncol(Y))
        } 
        ll <- function(beta)
            .mlt_loglik_exact(todistr, Y, Yprime, offset, trunc)(.parm(beta))
        sc <- function(beta)
            .mlt_score_exact(todistr, Y, Yprime, offset, trunc)(.parm(beta))[, !fix, drop = FALSE]
    } else {
        .makeY <- function(data, finite, nc = NULL) {
            if (any(finite)) {
                tmp <- model.matrix(model, data = data[finite,,drop = FALSE])
                nc <- colnames(tmp)
            } else {
                tmp <- -Inf
            }
            ret <- matrix(-Inf, nrow = nrow(data), ncol = length(nc))
            ret[finite,] <- tmp
            colnames(ret) <- nc
            ret
        }
        if (!is.null(Y)) {
            lY <- .makeY(ldata, lfinite, colnames(Y))
            rY <- .makeY(rdata, rfinite, colnames(Y))
        } else {
            if (any(lfinite)) {
                lY <- .makeY(ldata, lfinite)
                rY <- .makeY(rdata, rfinite, colnames(lY))
            } else {
                rY <- .makeY(rdata, rfinite)
                lY <- .makeY(ldata, lfinite, colnames(rY))
            }
        }
        if (all(!exact)) {
            ui <- attr(lY, "constraint")$ui
            ci <- attr(lY, "constraint")$ci
        }
        if (!is.null(fixed)) {
            stopifnot(all(names(fixed) %in% colnames(lY)))
            fix <- colnames(lY) %in% names(fixed)
            ui <- ui[,!fix,drop = FALSE]
            .parm <- function(beta) {
                ret <- numeric(ncol(lY))
                ret[fix] <- fixed
                ret[!fix] <- beta
                ret
            }
        } else {
            .parm <- function(beta) beta
            fix <- rep(FALSE, ncol(lY))
        } 
        ll <- function(beta) {
            ret <- numeric(nrow(data))
            ret[exact] <- .mlt_loglik_exact(todistr, Y, Yprime, offset[exact], trunc)(.parm(beta))
            ret[!exact] <- .mlt_loglik_interval(todistr, lY, rY, offset[!exact], trunc)(.parm(beta))
            ret
        }
        sc <- function(beta) {
            ret <- matrix(0, nrow = nrow(data), ncol = ncol(Y))
            ret[exact,] <- .mlt_score_exact(todistr, Y, Yprime, offset[exact], trunc)(.parm(beta))
            ret[!exact,] <- .mlt_score_interval(todistr, lY, rY, offset[!exact], trunc)(.parm(beta))
            ret[, !fix, drop = FALSE]
        }
    }


    z <- scale(edata[[response]])
    theta <- coef(lm(z ~ Y - 1))
    if (!is.null(fixed)) theta <- theta[!fix]

    stopifnot(any(is.finite(ci)))
    ui <- as(ui[is.finite(ci),,drop = FALSE], "matrix")
    ci <- ci[is.finite(ci)]

    loglikfct <- function(beta) -sum(weights * ll(beta))
    scorefct <- function(beta) -colSums(weights * sc(beta))
    if (sum(abs(ui)) > 0) {
        optimfct <- function(beta, hessian = FALSE) 
            constrOptim(theta = beta, f = loglikfct,
                        grad = scorefct, ui = ui,ci = ci, hessian = hessian)
    } else {
        optimfct <- function(beta, hessian = FALSE) 
            optim(par = beta, fn = loglikfct, 
                  gr = scorefct, hessian = hessian)
    } ### refit for hessian in vcov?

    ret <- optimfct(theta, hessian = FALSE)
    if (ret$convergence != 0)
        warning("algorithm did not converge")

    ret$model <- model
    ret$response <- response
    ret$todistr <- todistr
    ret$loglik <- loglikfct
    ret$score <- scorefct
    ret$optim <- optimfct
    class(ret) <- "mlt"
    return(ret)
}

mlt <- .mlt_fit

