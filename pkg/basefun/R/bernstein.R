
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

    s <- support(var)
    b <- bounds(var)
    if (log_first) {
        stopifnot(bounds(var)[[1]][1] >= 0)
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
        s[[1]] <- range(s[[1]])
        tmp <- Bernstein_basis(var, order = order, extrapolate = FALSE,
                               log_first = log_first)
        tmpdf <- as.data.frame(s)
        left <- s[[1]][1] > b[[1]][1]
        right <- s[[1]][2] < b[[1]][2]
        if (left || right) {
            if (left && !right) tmpdf <- tmpdf[-2,,drop = FALSE]
            if (!left && right) tmpdf <- tmpdf[-1,,drop = FALSE]
            dr <- 2
            names(dr) <- names(s)
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

        dlog <- 1
        if (log_first && deriv > 0) {
            expr <- D(quote(log(ox)), "ox")
            if (deriv > 1) {
                for (d in 2:deriv)
                    expr <- D(expr, "ox")
            }
            dlog <- eval(expr)
        } 

        if (deriv > 0)
            X <- X * (1 / diff(support)^deriv) * dlog
        if (deriv < 0)
            X[] <- 0
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

    expr <- D(quote(log(s)), "s")
    if (deriv > 1) {
        for (d in 2:deriv)
            expr <- D(expr, "ox")
        }
    dlog <- eval(expr)


    if (any(small)) {
        dsmall <- data.frame(x = rep(s[1], sum(small)))
        names(dsmall) <- varname
        if (attr(object, "log_first")) {
            ### f'(s) = b'(log(s)) / s but we need b'(log(s)) only
            xdiff <- (log(x[small]) - log(s[1]))
        } else {
            xdiff <- x[small] - s[1]
        }
        dfun <- function(deriv) {
            dlog <- rep(1L, 2)
            if (attr(object, "log_first")) {
                expr <- D(quote(log(s)), "s")
                if (deriv > 1) {
                    for (d in 2:deriv)
                        expr <- D(expr, "s")
                    }
                dlog <- eval(expr)
            }
            if (deriv >= maxderiv)
                return(object(data = dsmall, deriv = deriv))
            return(object(data = dsmall, deriv = deriv) +
                   dfun(deriv + 1L) / dlog[1] * xdiff)
        }
        ret[small,] <- dfun(deriv)
    }
    if (any(large)) {
        dlarge <- data.frame(x = rep(s[2], sum(large)))
        names(dlarge) <- varname
        X <- object(data = dlarge)
        if (attr(object, "log_first")) {
            xdiff <- (log(x[large]) - log(s[2]))
        } else {
            xdiff <- x[large] - s[2]
        }
        dfun <- function(deriv) {
            dlog <- rep(1L, 2)
            if (attr(object, "log_first")) {
                expr <- D(quote(log(s)), "s")
                if (deriv > 1) {
                    for (d in 2:deriv)
                        expr <- D(expr, "s")
                    }
                dlog <- eval(expr)
            }
            if (deriv >= maxderiv)
                return(object(data = dlarge, deriv = deriv))
            return(object(data = dlarge, deriv = deriv) +   
                   dfun(deriv + 1L) / dlog[2] * xdiff)
        }
        ret[large,] <- dfun(deriv)
    }
    return(ret)
}
