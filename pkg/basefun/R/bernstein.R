
.B <- function(j, n, integrate = FALSE) {
    function(x) {
        if (j < 0 || (n - j) < 0)
            return(rep(0, length(x)))
        if (integrate)
            return(pbeta(x, shape1 = j + 1, shape2 = n - j + 1) / (n + 1))
        return(dbeta(x, shape1 = j + 1, shape2 = n - j + 1) / (n + 1))
    }
}

.Bx <- function(x, j, n, deriv = 0L, integrate = FALSE) {
    if (deriv == 0L)
        return(.B(j, n, integrate = integrate)(x))
    return(n * (.Bx(x, j - 1, n - 1, deriv = deriv - 1L, integrate = integrate) -
                .Bx(x, j,     n - 1, deriv = deriv - 1L, integrate = integrate)))
}

### set-up Bernstein polynom basis functions
### http://en.wikipedia.org/wiki/Bernstein_polynomial
### http://dx.doi.org/10.1080/02664761003692423
### arXiv preprint arXiv:1404.2293
Bernstein_basis <- function(var, order = 2, 
                            ui = c("none", "increasing", "decreasing", 
                                   "cyclic", "zerointegral", "positive", "negative"),
                            extrapolate = FALSE, log_first = FALSE) {

    zeroint <- FALSE
    ui <- match.arg(ui, several.ok = TRUE)
    if (length(ui) > 2) ui <- "none"
    if (length(ui) > 1) 
        stopifnot(length(ui) == 2)
    if ("zerointegral" %in% ui) {
        zeroint <- TRUE
        ui <- ifelse(length(ui) == 2, ui[ui != "zerointegral"], "none")
    } else {
        ui <- paste(sort(ui), collapse = ".")
    }

    s <- os <- support(var)
    b <- ob <- bounds(var)
    if (log_first) {
        stopifnot(bounds(var)[[1]][1] > .Machine$double.eps)
        s <- lapply(s, log)
        b <- lapply(b, log)
    }    

    constr <- switch(ui,
        "none" = list(ui = Diagonal(order + 1), 
                      ci = rep(-Inf, order + 1)),
        "increasing" = list(ui = diff(Diagonal(order + 1), differences = 1), 
                            ci = rep(0, order)),
        "decreasing" = list(ui = diff(Diagonal(order + 1), differences = 1) * -1,
                            ci = rep(0, order)),
        "increasing.positive" = {
            tmp <- Bernstein_basis(var, order = order, extrapolate = FALSE)
            tmpdf <- as.data.frame(mkgrid(var))
            B0 <- model.matrix(tmp, data = tmpdf)[1,]
            list(ui = rbind(B0, diff(Diagonal(order + 1), differences = 1)), 
                        ci = rep(0, order + 1))
        },
        "positive" = list(ui = Diagonal(order + 1), ci = rep(0, order + 1)),
        "negative" = list(ui = -Diagonal(order + 1), ci = rep(0, order + 1)),
        "default" = stop(paste(ui, "not yet implemented")))

    ### linear extrapolation, f''(support) = 0
    if (extrapolate) {
        os[[1]] <- range(os[[1]])
        tmp <- Bernstein_basis(var, order = order, extrapolate = FALSE,
                               log_first = log_first)
        tmpdf <- as.data.frame(os)
        left <- os[[1]][1] > ob[[1]][1]
        right <- os[[1]][2] < ob[[1]][2]
        if (left || right) {
            if (left && !right) tmpdf <- tmpdf[-2,,drop = FALSE]
            if (!left && right) tmpdf <- tmpdf[-1,,drop = FALSE]
            dr <- 2
            names(dr) <- names(os)
            B0 <- model.matrix(tmp, data = tmpdf, deriv = dr)
            ### <FIXME> we don't have infrastructure for equality
            ### constraints, so use <= and >= 
            ### </FIXME>
            constr$ui <- rbind(constr$ui, B0, -B0)
            constr$ci <- c(constr$ci, rep(0, 2 * nrow(B0)))
        }
    }

    stopifnot(inherits(var, "numeric_var"))
    varname <- variable.names(var)
    support <- range(s[[varname]])
    stopifnot(all(diff(support) > 0))

    basis <- function(data, deriv = 0L, integrate = FALSE) {

        stopifnot(check(var, data))
        if (is.atomic(data)) {
            x <- data
        } else {
            if (is.null(varname)) varname <- colnames(data)[1]
            x <- data[[varname]]
        }
        ox <- x
        if (log_first) x <- log(x)

        ### applies to all basis functions
        ### deriv = -1 => 0
        ### deriv = 0 => f(x)
        ### deriv = 1 => f'(x)
        ### deriv = 2 => f''(x)
        ### ...

        if (deriv < 0) x <- rep(mean(support), NROW(data))

        stopifnot(order > max(c(0, deriv - 1L)))
        x <- (x - support[1]) / diff(support)
        stopifnot(all(x >= 0 & x <= 1))

        X <- do.call("cbind", lapply(0:order, function(j) 
                     .Bx(x, j, order, deriv = max(c(0, deriv)), integrate = integrate)))

        if (deriv < 0) {
            X[] <- 0
        } else{
            if (!log_first && deriv > 0) {
                X <- X * (1 / diff(support)^deriv)
            } else {
                if (log_first && deriv > 0) {
                    X <- switch(as.character(deriv),
                        "1" = {
                            X * (1 / diff(support)^deriv) / ox
                        },
                        "2" = {
                            X1 <- do.call("cbind", lapply(0:order, function(j) 
                               .Bx(x, j, order, deriv = 1L, integrate = integrate)))
                            (X - X1) / (ox^2)
                        },
                        stop("deriv > 2 not implemented for log_first"))
                }
                ### do nothing for deriv == 0
            }
        }
        colnames(X) <- paste("Bs", 1:ncol(X), "(", varname, ")", sep = "")
        if (zeroint) {
            X <- X[, -ncol(X), drop = FALSE] - X[, ncol(X), drop = TRUE]
            constr$ui <- constr$ui[, -ncol(constr$ui),drop = FALSE]
        }
        attr(X, "constraint") <- constr
        attr(X, "Assign") <- matrix(varname, ncol = ncol(X))
        return(X)
    }

    attr(basis, "variables") <- var
    attr(basis, "intercept") <- zeroint
    attr(basis, "log_first") <- log_first

    class(basis) <- c("Bernstein_basis", "basis", class(basis))
    return(basis)
}

### evaluate model.matrix of Bernstein polynom
model.matrix.Bernstein_basis <- function(object, data,
                                         deriv = 0L,
                                         maxderiv = 1L, 
                                         integrate = FALSE, ...) {

    varname <- variable.names(object)
    deriv <- .deriv(varname, deriv)
    if (deriv < 0)
        return(object(data = data, deriv = deriv))
    x <- data[[varname]]
    s <- range(support(as.vars(object))[[varname]])
    small <- x < s[1]
    large <- x > s[2]
    if (all(!small) && all(!large))
        return(object(data = data, deriv = deriv, integrate = integrate))

    ### extrapolate linearily for x outside s
    stopifnot(!integrate)
    ### interpolate f^d(s + x) = f^d(s) + f^d+1(s) * (x - s)
    ### note: the code allows for maxderiv <- deriv + d, d > 1
    ### Maybe: constrain second deriv in s to be zero
    order <- get("order", envir = environment(object))
    data[[varname]][small] <- s[1]
    data[[varname]][large] <- s[2]
    ret <- object(data = data, deriv = deriv, integrate = integrate)

    if (any(small)) {
        dsmall <- data.frame(x = rep(s[1], sum(small)))
        names(dsmall) <- varname
        if (!attr(object, "log_first")) {
            xdiff <- x[small] - s[1]
            ret[small,] <- switch(as.character(deriv),
                "0" = {
                    object(data = dsmall, deriv = deriv) +
                    object(data = dsmall, deriv = deriv + 1L) * xdiff
                },
                "1" = {
                    object(data = dsmall, deriv = deriv)
                },
                0)
        } else {
            xdiff <- (log(x[small]) - log(s[1]))
            ret[small,] <- switch(as.character(deriv),
                "0" = {
                    object(data = dsmall, deriv = deriv) +
                    ### note: object(data = dsmall, deriv = deriv + 1)
                    ### returns a'(log(s)) / s; we need a'(log(s))
                    object(data = dsmall, deriv = deriv + 1L) * s[1] * xdiff
                },
                "1" = {
                    object(data = dsmall, deriv = deriv) * s[1] / x[small]
                },
                "2" = {
                    -object(data = dsmall, deriv = deriv - 1L) * s[1] / (x[small]^2)
                },
                stop("deriv >= 2 not implemented"))
        } 
    }
    if (any(large)) {
        dlarge <- data.frame(x = rep(s[2], sum(large)))
        names(dlarge) <- varname
        if (!attr(object, "log_first")) {
            xdiff <- x[large] - s[2]
            ret[large,] <- switch(as.character(deriv),
                "0" = {
                    object(data = dlarge, deriv = deriv) +
                    object(data = dlarge, deriv = deriv + 1L) * xdiff
                },
                "1" = {
                    object(data = dlarge, deriv = deriv)
                },
                0)
        } else {
            xdiff <- (log(x[large]) - log(s[2]))
            ret[large,] <- switch(as.character(deriv),
                "0" = {
                    object(data = dlarge, deriv = deriv) +
                    ### note: object(data = dlarge, deriv = deriv + 1)
                    ### returns a'(log(s)) / s; we need a'(log(s))
                    object(data = dlarge, deriv = deriv + 1L) * s[2] * xdiff
                },
                "1" = {
                    object(data = dlarge, deriv = deriv) * s[2] / x[large]
                },
                "2" = {
                    -object(data = dlarge, deriv = deriv - 1L) * s[2] / (x[large]^2)
                },
                stop("deriv >= 2 not implemented"))
        } 
    }
    return(ret)
}
