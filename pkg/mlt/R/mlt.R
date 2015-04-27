
.findstart <- function(model, y, data, ui = NULL, ci = NULL, fix = NULL, fixed = NULL) {

    response <- model$response

    r <- range(unlist(y), na.rm = TRUE, finite = TRUE)
    if (any(!is.na(y$exact))) {
        ytmp <- ifelse(!is.na(y$exact), y$exact, 
                      (ifelse(is.finite(y$cleft), y$cleft, r[1]) + 
                       ifelse(is.finite(y$cright), y$cright, r[2])) / 2)
    } else {
        if (is.factor(y$cright) || is.integer(y$right)) {
            ytmp <- y$cright
            ytmp[is.na(ytmp)] <- levels(y$cright)[nlevels(y$cright)]
        } else {
            ytmp <- (ifelse(is.finite(y$cleft), y$cleft, r[1]) +
                     ifelse(is.finite(y$cright), y$cright, r[2])) / 2
        }
    }
    if (is.factor(ytmp)) ytmps <- (1:nlevels(ytmp))[ytmp] else ytmps <- ytmp

    z <- scale(ytmps)
    data[[response]] <- ytmp
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

### <FIXME> rename fixed to coef and allow for specification of coefs, ie fitted models? </FIXME>
.mlt_setup <- function(model, data, weights = NULL, 
                       offset = NULL, fixed = NULL) {

    if (is.null(weights)) weights <- rep(1, nrow(data))
    if (is.null(offset)) offset <- rep(0, nrow(data))
    stopifnot(nrow(data) == length(weights))
    stopifnot(nrow(data) == length(offset))

    response <- model$response
    stopifnot(length(response) == 1)
    y <- data[[response]]

    todistr <- model$todistr

    if (!inherits(y, "response")) {
        ytype <- .type_of_response(y)
        if (is.na(ytype))
            stop("cannot deal with response class", class(y))

        if (ytype == "survival") {
            y <- .Surv2R(y)
        } else {
            y <- R(exact = y)
        }
    }
    
    eY <- .mm_exact(model, data = data, response = response, object = y)
    iY <- .mm_interval(model, data = data, response = response, object = y)

    if (is.null(eY)) {
        Y <- iY$Yleft
    } else {
        Y <- eY$Y
    }

    ui <- attr(Y, "constraint")$ui
    ci <- attr(Y, "constraint")$ci

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

    exact <- .exact(y)

    ll <- function(beta) {
        ret <- numeric(nrow(data))
        if (any(exact))
            ret[exact] <- .mlt_loglik_exact(todistr, eY$Y, eY$Yprime, offset[exact], eY$trunc)(.parm(beta))
        ret[!exact] <- .mlt_loglik_interval(todistr, iY$Yleft, iY$Yright, offset[!exact], iY$trunc)(.parm(beta))
        return(ret)
    }
    sc <- function(beta) {
        ret <- matrix(0, nrow = nrow(data), ncol = length(fix))
        if (any(exact))
            ret[exact,] <- .mlt_score_exact(todistr, eY$Y, eY$Yprime, offset[exact], eY$trunc)(.parm(beta))
        ret[!exact,] <- .mlt_score_interval(todistr, iY$Yleft, iY$Yright, offset[!exact], iY$trunc)(.parm(beta))
        return(ret[, !fix, drop = FALSE])
    }
    he <- function(beta) {
        ret <- 0
        if (any(exact))
            ret <- ret + .mlt_hessian_exact(todistr, eY$Y, eY$Yprime, offset[exact], eY$trunc, weights[exact])(.parm(beta))
        if (any(!exact))
            ret <- ret + .mlt_hessian_interval(todistr, iY$Yleft, iY$Yright, offset[!exact], iY$trunc, weights[!exact])(.parm(beta))
        return(ret[!fix, !fix, drop = FALSE])
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
        ci <- ci + sqrt(.Machine$double.eps) ### we need ui %*% theta > ci, not >= ci
    }

    optimfct <- function(theta, ...) {
        control <- list(...)
        if (!is.null(ui))
            return(spg(par = theta, fn = loglikfct, gr = scorefct, project = "projectLinear",
                       projectArgs = list(A = ui, b = ci, meq = 0), control = control))
        return(spg(par = theta, fn = loglikfct, gr = scorefct, control = control))
    }

    theta <- .findstart(model, y, data, ui, ci, fix, fixed)

    ### at least one serious constraint
    if (!is.null(ui))
        if(!all(ui %*% theta - ci >= 0)) 
            warning("start from inadmissible parameters")

    coef <- rep(NA, length(fix))
    coef[fix] <- fixed
    names(coef) <- colnames(Y)

    ret <- list()
    ret$theta <- theta
    ret$parm <- .parm
    ret$coef <- coef
    ret$model <- model
    ret$data <- data
    ret$weights <- weights
    ret$offset <- offset
    ret$response <- response
    ret$todistr <- todistr
    ret$loglik <- loglikfct
    ret$score <- scorefct
    ret$hessian <- he
    ret$optimfct <- optimfct
    class(ret) <- c("mlt_setup", "mlt")
    return(ret)
}

.mlt_fit <- function(object, theta = object$theta, check = TRUE, trace = FALSE, ...) {

    ret <- try(object$optimfct(theta, trace = trace, ...))    
    if (inherits(ret, "try-error") || ret$convergence != 0 || ret$gradient > 1)
        ret <- object$optimfct(ret$par, trace = trace, ...)

    if (ret$convergence != 0)
        warning("algorithm did not converge")

    cls <- class(object)
    object <- c(object, ret)
    object$coef <- object$parm(ret$par)
    class(object) <- c("mlt_fit", cls)
    
    if (check) {
        ### check gradient  and hessian
        gr <- numDeriv::grad(object$loglik, ret$par)
        s <- estfun(object)
        cat("Gradient")
        print(all.equal(gr, -s, check.attributes = FALSE))

        H1 <- numDeriv::hessian(object$loglik, ret$par)
        H2 <- Hessian(object)
        cat("Hessian:")
        print(all.equal(H1, H2, check.attributes = FALSE))
    }

    return(object)
}

mlt <- function(..., dofit = TRUE) {
    args <- list(...)
    if (sum(names(args) == "") > 1)
        stop("function mlt requires names arguments")
    sa <- c("", "model", "data", "weights", "offset", "trunc", "fixed")
    setupargs <- args[names(args) %in% sa]
    fitargs <- args[!(names(args) %in% sa)]
    s <- do.call(".mlt_setup", setupargs)
    if (!dofit) return(s)
    fitargs$object <- s
    do.call(".mlt_fit", fitargs)
}

