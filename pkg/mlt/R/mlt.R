
.findstart <- function(model, data, ui = NULL, ci = NULL, fix = NULL, fixed = NULL) {

    y <- data[[response <- model$response]]
    ytype <- .type_of_response(y)
    if (is.na(ytype))
        stop("cannot deal with response class", class(y))

    if (ytype == "double") {
        z <- y
    } else if (ytype == "survival") {
        sy <- .Surv2matrix(y)
        lna <- is.na(sy[, "left"])
        rna <- is.na(sy[, "right"])
        interval <- !lna & !rna
        sy[!is.finite(sy[, "left"]), "left"] <- 
            min(sy[is.finite(sy)], na.rm = TRUE)
        sy[!is.finite(sy[, "right"]), "right"] <- 
            max(sy[is.finite(sy)], na.rm = TRUE)
        z <- numeric(ncol(sy))
        z[lna] <- sy[lna, "right"]
        z[rna] <- sy[rna, "left"]
        z[interval] <- sy[interval, "left"] + (sy[interval, "right"] - 
                       sy[interval, "left"]) / 2
        y <- z - min(z) + .1
    } else if (ytype %in% c("ordered", "unordered")) {
        z <- (1:nlevels(y))[y]
    } else if (ytype == "integer") {
        z <- y
    } else {
        stop("cannot deal with response class", class(y))
    }
    z <- scale(z)
    data[[response]] <- y
    Y <- model.matrix(model, data)
    if (any(fix)) {
        z <- z - Y[, fix, drop = FALSE] %*% fixed
        Y <- Y[, !fix, drop = FALSE]
    }
    theta <- lm.fit(y = z, x = Y)$coef
    if (any(is.na(theta))) {
        warning("cannot determine some starting values")
        theta[is.na(theta)] <- 0
    }
    if (is.null(ui)) return(theta)
    if (any(ui %*% theta - ci <= 0)) {
        citmp <- ci + .1
        theta <- MASS::ginv(ui) %*% citmp
    }
    theta
}

.mlt_fit <- function(model, data, weights = NULL, 
                     offset = NULL, fixed = NULL, trunc = NULL, 
                     theta = NULL, ...) {

    if (is.null(weights)) weights <- rep(1, nrow(data))
    if (is.null(offset)) offset <- rep(0, nrow(data))
    stopifnot(nrow(data) == length(weights))
    stopifnot(nrow(data) == length(offset))

    response <- model$response
    stopifnot(length(response) == 1)
    y <- data[[response]]

    todistr <- model$todistr

    if (!is.null(trunc)) {
        stopifnot(is.character(trunc))
        stopifnot(all(names(trunc) %in% c("left", "right")))
        lefttrunc <- righttrunc <- NULL
        if ("left" %in% names(trunc))
            lefttrunc <- data[[trunc["left"]]]
        if ("right" %in% names(trunc))
            righttrunc <- data[[trunc["right"]]]
        trunc <- list(left = lefttrunc, right = righttrunc)
    }

    ytype <- .type_of_response(y)
    if (is.na(ytype))
        stop("cannot deal with response class", class(y))

    if (ytype == "double") {
        exact <- rep(TRUE, nrow(data))
        edata <- data
    } else if (ytype == "survival") {
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
    } else if (ytype %in% c("ordered", "unordered")) {
        if (ytype == "unordered") {
            warning("results may depend on ordering of levels")
            y <- ordered(y)
        }
        exact <- rep(FALSE, nrow(data))
        ldata <- rdata <- data
        sy <- sort(unique(y))
        ldata[[response]] <- sy[pmax(1, unclass(y) - 1)]
        lfinite <- (y > levels(y)[1])
        rdata[[response]] <- y
        rfinite <- (y < levels(y)[nlevels(y)])
    } else if (ytype == "integer") {
        exact <- rep(FALSE, nrow(data))
        ldata <- rdata <- data
        ldata[[response]] <- y - 1
        lfinite <- (ldata[[response]] >= 0)
        rdata[[response]] <- y
        rfinite <- rep(TRUE, nrow(rdata))
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
        checktheta <- function(theta) {
            if (all(Yprime %*% .parm(theta) > .Machine$double.eps))
                return(theta)
            tmp <- Yprime[, !fix, drop = FALSE]
            i <- colSums(abs(tmp)) > 0
            nt <- solve(crossprod(tmp[, i, drop = FALSE]), colSums(tmp[,i]))
            ret <- numeric(ncol(tmp))
            ret[i] <- nt
            ret
        }
        ll <- function(beta)
            .mlt_loglik_exact(todistr, Y, Yprime, offset, trunc)(.parm(beta))
        sc <- function(beta)
            .mlt_score_exact(todistr, Y, Yprime, offset, trunc)(.parm(beta))[, !fix, drop = FALSE]
        he <- function(beta) {
            ret <- .mlt_hessian_exact(todistr, Y, Yprime, offset[exact], trunc, weights)(.parm(beta))
            return(ret[!fix, !fix, drop = FALSE])
        }
    } else {
        .makeY <- function(data, finite, nc = NULL, INF = -Inf) {
            if (any(finite)) {
                tmp <- model.matrix(model, data = data[finite,,drop = FALSE])
                nc <- colnames(tmp)
            } else {
                tmp <- INF
            }
            ret <- matrix(INF, nrow = nrow(data), ncol = length(nc))
            ret[finite,] <- tmp
            colnames(ret) <- nc
            if (!is.null(attr(tmp, "constraint"))) {
                attr(ret, "constraint") <- attr(tmp, "constraint") ### constraints!
                attr(ret, "Assign") <- attr(tmp, "Assign")
            }
            ret
        }
        if (!is.null(Y)) {
            lY <- .makeY(ldata, lfinite, colnames(Y), INF = -Inf)
            rY <- .makeY(rdata, rfinite, colnames(Y), INF = Inf)
        } else {
            if (any(lfinite)) {
                lY <- .makeY(ldata, lfinite, INF = -Inf)
                rY <- .makeY(rdata, rfinite, colnames(lY), INF = Inf)
            } else {
                rY <- .makeY(rdata, rfinite, INF = Inf)
                lY <- .makeY(ldata, lfinite, colnames(rY), INF = -Inf)
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
        checktheta <- function(theta) {
            tmp <- NULL
            if (any(exact))
                tmp <- Yprime
            rtmp <- rY
            if (any(!is.finite(rtmp))) 
                rtmp[!is.finite(rtmp)] <- lY[!is.finite(rtmp)] + 
                    seq(from = .1, to = 1, length = length(lY[!is.finite(rtmp)]))
            ltmp <- lY
            if (any(!is.finite(ltmp)))
               ltmp[!is.finite(ltmp)] <- rY[!is.finite(ltmp)] - 1
#                    rev(seq(from = .1, to = 1, length = length(rY[!is.finite(ltmp)])))
            tmp <- rbind(tmp, rtmp - ltmp)
            if (all(tmp %*% .parm(theta) > .Machine$double.eps))
                return(theta)
            tmp <- tmp[, !fix, drop = FALSE]
            i <- colSums(abs(tmp)) > 0
            nt <- solve(crossprod(tmp[, i, drop = FALSE]), colSums(tmp[,i]))
            ret <- numeric(ncol(tmp))
            ret[i] <- nt
            ret
        }
        ll <- function(beta) {
            ret <- numeric(nrow(data))
            if (any(exact))
                ret[exact] <- .mlt_loglik_exact(todistr, Y, Yprime, offset[exact], trunc)(.parm(beta))
            ret[!exact] <- .mlt_loglik_interval(todistr, lY, rY, offset[!exact], trunc)(.parm(beta))
            return(ret)
        }
        sc <- function(beta) {
            ret <- matrix(0, nrow = nrow(data), ncol = length(fix))
            if (any(exact))
                ret[exact,] <- .mlt_score_exact(todistr, Y, Yprime, offset[exact], trunc)(.parm(beta))
            ret[!exact,] <- .mlt_score_interval(todistr, lY, rY, offset[!exact], trunc)(.parm(beta))
            return(ret[, !fix, drop = FALSE])
        }
        he <- function(beta) {
            ret <- 0
            if (any(exact))
                ret <- ret + .mlt_hessian_exact(todistr, Y, Yprime, offset[exact], trunc, weights[exact])(.parm(beta))
            ret <- ret + .mlt_hessian_interval(todistr, lY, rY, offset[!exact], trunc, weights[!exact])(.parm(beta))
            return(ret[!fix, !fix, drop = FALSE])
        }
    }

    loglikfct <- function(beta) -sum(weights * ll(beta))
    scorefct <- function(beta) -colSums(weights * sc(beta))

    if (all(!is.finite(ci))) {
        ui <- ci <- NULL
    } else {
        ui <- as(ui[is.finite(ci),,drop = FALSE], "matrix")
        ci <- ci[is.finite(ci)]
        r0 <- rowSums(abs(ui)) == 0
        ui <- ui[!r0,,drop = FALSE]
        ci <- ci[!r0]
        if (nrow(ui) == 0) ui <- ci <- NULL
    }

    optimfct <- function(beta, usescore = TRUE, ...) {
        control <- list(...)
        if (!is.null(ui)) {
            if (usescore)
                return(spg(par = beta, fn = loglikfct, gr = scorefct, project = "projectLinear",
                           projectArgs = list(A = ui, b = ci, meq = 0), control = control))
            return(spg(par = beta, fn = loglikfct, project = "projectLinear",
                       projectArgs = list(A = ui, b = ci, meq = 0), control = control))
        }
        if (usescore)
            return(spg(par = runif(length(beta)), fn = loglikfct, gr = scorefct, control = control))
        return(spg(par = runif(length(beta)), fn = loglikfct, control = control))
    }

    if (is.null(theta))
        theta <- .findstart(model, data, ui, ci, fix, fixed)

    ### at least one serious constraint
    if (!is.null(ui))
        if(!all(ui %*% theta - ci >= 0)) 
            warning("start from inadmissible parameters")

    # theta <- checktheta(theta)
    theta <- optimfct(theta, usescore = TRUE, maxit = 200, checkGrad = FALSE)$par

    ret <- try(optimfct(theta, ...))    
    if (inherits(ret, "try-error") || ret$convergence != 0 || ret$gradient > 1)
        ret <- optimfct(ret$par, ...)

    if (ret$convergence != 0)
        warning("algorithm did not converge")

    ### check gradient and hessian
    gr <- numDeriv::grad(loglikfct, ret$par)
    s <- scorefct(ret$par)
    cat("Score:", max(abs(gr / s)), " ")

    H1 <- numDeriv::hessian(loglikfct, ret$par)
    H2 <- he(ret$par)
    cat("Hessian:", max(abs(drop(H1 - H2))), "\n")

    coef <- numeric(length(fix))
    coef[fix] <- fixed
    coef[!fix] <- ret$par
    if (is.null(Y)) {
        names(coef) <- colnames(lY)
    } else {
        names(coef) <- colnames(Y)
    }
    ret$coef <- coef
    ret$model <- model
    ret$response <- response
    ret$todistr <- todistr
    ret$loglik <- loglikfct
    ret$score <- scorefct
    ret$hessian <- he
    ret$optim <- optimfct
    ret$data <- data
    class(ret) <- "mlt"
    return(ret)
}

mlt <- .mlt_fit

