
.frmt <- function(q) {
    if (is.factor(q)) 
        return(as.character(q))
    return(formatC(q, digits = 3, width = 5))
}

### transformation function
tmlt <- function(object, newdata = object$data, q = NULL, n = 50, ...) {

    vn <- unlist(varnames(object$model$model))
    vnx <- vn[!(vn %in% object$response)]
    y <- object$response
    model <- object$model$model

    if (is.data.frame(newdata)) {

        ### in sample predictions
        ### this will _not_ work for censored responses
        if (!is.null(newdata[[y]]) & is.null(q)) {
            ret <- c(predict(model, newdata = newdata, 
                             coef = coef(object), ...))
            names(ret) <- rownames(newdata)
            ### P(Y \le y_K) = 1 but trafo can be < Inf
            ### depending on parameterisation
            if (is.factor(f <- newdata[[y]])) {
                i <- f == levels(f)[nlevels(f)]
                if (any(i)) 
                    ret[i] <- Inf
            }
            return(ret)
        }

        ### extra quantiles, compute transformation
        ### for each q and each row of newdata
        if (is.null(q)) 
            q <- generate(object, n = n)[[y]]
        stopifnot(length(unique(q)) == length(q))
        dim <- c(length(q), nrow(newdata))

        ### <FIXME> this triggers a trick in 
        ### basefun:::predict.formula_basis; better checks needed </FIXME>
        names(dim) <- c(y, vnx[1])
        newdata <- as.list(newdata)
        newdata[[y]] <- q
        ret <- predict(object$model$model, newdata = newdata, 
                       coef = coef(object), dim = dim, ...)
        dn <- vector(mode = "list", length = 2)
        names(dn) <- c(y, "newdata") ### deparse(substitute(newdata))) ?
        dn[[y]] <- .frmt(q)
        dn[["newdata"]] <- rownames(newdata)
        dimnames(ret) <- dn

        ### trafo of last level is always Inf, see above
        if (is.factor(f <- newdata[[y]])) {
            i <- f == levels(f)[nlevels(f)]
            if (any(i))
                ret[i,] <- Inf
        }
        return(ret)
    }

    ### need to generate newdata outside tmlt such that
    ### the rows of expand.grid(newdata) match the elements of
    ### the return value
    stopifnot(is.null(q))
    stopifnot(y %in% names(newdata))
    ret <- predict(object$model$model, newdata = newdata, 
                   coef = coef(object), dim = TRUE, ...)
    dn <- lapply(newdata, .frmt)
    dimnames(ret) <- dn

    ### trafo of last level is always Inf, see above
    if (is.factor(f <- newdata[[y]])) {
        i <- f == levels(f)[nlevels(f)]
        if (any(i)) {
            args <- lapply(names(dn), function(d) {
                if (d == y)
                    return(i)
                return(1:(dim(ret)[which(names(dn) == d)]))
            })
            ret <- do.call("[<-", c(list(i = ret), args, 
                                    list(value = Inf)))
        }
    }
    return(ret)
}

### distribution function
pmlt <- function(object, newdata = object$data, q = NULL, n = 50)
    object$model$todistr$p(tmlt(object = object, newdata = newdata,
                                q = q, n = n))

### survivor function
smlt <- function(object, newdata = object$data, q = NULL, n = 50)
    1 - pmlt(object = object, newdata = newdata, q = q, n = n)

### cumulative hazard function
Hmlt <- function(object, newdata = object$data, q = NULL, n = 50)
    -log(smlt(object = object, newdata = newdata, q = q, n = n))

### numerical inversion of distribution function
### to get quantile function
.p2q <- function(prob, q, p, interpolate = FALSE, fact = 1.1) {

    prob <- cbind(0, prob, 1)
    i <- rowSums(prob < p)
    if (.type_of_response(q) != "double")
        return(q[i])

    ### return interval censored quantiles
    if (!interpolate) {
        qq <- c(-Inf, q, Inf)
        return(data.frame(left = qq[i], right = qq[i + 1]))
    }

    ### interpolate linearily
    qq <- c(abs(min(q)) * fact * sign(min(q)), q, max(q) * fact)
    ptmp <- cbind(prob[cbind(1:nrow(prob), i)], 
                  prob[cbind(1:nrow(prob), i + 1)])
    beta <- (ptmp[,2] - ptmp[,1]) / (qq[i + 1] - qq[i])
    alpha <- ptmp[,1] - beta * qq[i]
    return((p - alpha) / beta) 
}

### quantile function
qmlt <- function(object, newdata = object$data, p = .5, n = 50, 
                 interpolate = TRUE) {

    y <- object$response
    ### don't accept user-generated quantiles
    newdata[[y]] <- NULL
    q <- generate(object, n = n)[[y]]
    if (is.data.frame(newdata)) {
        prob <- pmlt(object, newdata, q = q)
    } else {
        nm <- names(newdata)
        newdata[[y]] <- q
        newdata <- newdata[c(y, nm)]
        prob <- pmlt(object, newdata)
    }

    ### convert potential array-values distribition function
    ### to matrix where rows correspond to observations newdata 
    ### and columns to quantiles q
    ptmp <- t(matrix(prob, nrow = length(q)))
    nr <- nrow(ptmp)
    ptmp <- ptmp[rep(1:nr, each = length(p)),,drop = FALSE]
    pp <- rep(p, nr) ### p varies fastest
    ret <- .p2q(ptmp, q, pp, interpolate = interpolate)


    ### arrays of factors are not allowed
    if (is.factor(q)) return(ret)

    ### return list of length(p)
    if (!interpolate) {
        tmp <- vector(mode = "list", length = length(p))
        names(tmp) <- .frmt(p)
        for (i in 1:length(p)) {
            idx <- 1:nr + (i - 1) * nr
            if (is.data.frame(ret))  
                tmp[[i]] <- ret[idx,]
            else 
                tmp[[i]] <- ret[idx]
        }
        return(tmp)
    }

    dim <- dim(prob)
    dim[1] <- length(p)
    dn <- c(list(p = .frmt(p)), dimnames(prob)[-1])
    return(array(ret, dim = dim, dimnames = dn))
}
    
### simulate from model object with data newdata
simulate.mlt <- function(object, nsim = 1, seed = NULL, 
                         newdata = object$data, n = 50, 
                         interpolate = TRUE, ...) {

    ### from stats:::simulate.lm
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }

    y <- object$response
    newdata[[y]] <- NULL
    ### don't accept user-generated quantiles
    q <- generate(object, n = n)[[y]]
    if (is.data.frame(newdata)) {
        p <- runif(nsim * NROW(newdata))
        ### basically compute quantiles for p; see qmlt
        prob <- pmlt(object, newdata, q = q)
        prob <- t(prob[, rep(1:ncol(prob), nsim),drop = FALSE])
        ret <- .p2q(prob, q, p, interpolate = interpolate)
        if (nsim > 1) {
            tmp <- vector(mode = "list", length = nsim)
            for (i in 1:nsim) {
                idx <- 1:nrow(newdata) + (i - 1) * nrow(newdata)
                if (is.data.frame(ret)) 
                    tmp[[i]] <- ret[idx,]
                else 
                    tmp[[i]] <- ret[idx]
            }
            ret <- tmp
        }
    } else {
        stop("not yet implemented")
    }
    return(ret)
}

### density
dmlt <- function(object, newdata = object$data, q = NULL, n = 50, 
                 log = FALSE) {

    response <- object$data[[y <- object$response]]

    ### Lebesgue density only for double
    if (.type_of_response(response) == "double") {
        trafo <- tmlt(object, newdata = newdata, q = q, n = n)
        deriv <- 1
        names(deriv) <- y
        trafoprime <- tmlt(object, newdata = newdata, q = q, n = n, 
                           deriv = deriv)
        if (log)
            return(object$model$todistr$d(trafo, log = TRUE) + log(trafoprime))
        return(object$model$todistr$d(trafo) * trafoprime)
    }

    ### for factors and integers compute density as F(y) - F(y - 1)
    lev <- levels(response)
    if (is.data.frame(newdata)) {

        ### in sample density
        if (!is.null(newdata[[y]]) & is.null(q)) {
            q <- newdata[[y]]
            first <- q == lev[1]
            qwoK <- factor(lev[pmax(unclass(q) - 1), 1], 
                           levels = lev, labels = lev)
            p <- pmlt(model, newdata = newdata)
            newdata[[y]] <- qwoK
            pwoK <- pmlt(model, newdata = newdata)
            pwoK[first] <- 0
            ret <- p - pwoK
        } else {
            ### extra quantiles, compute density
            ### for each q and each row of newdata 
            if (is.null(q))
                q <- generate(object, n = n)[[y]]

            first <- q == lev[1]
            qfirst <- q[first]
            qwoK <- q[q != lev[length(lev)]]
            qwo1 <- q[q != lev[1]]

            pfirst <- pmlt(object, newdata = newdata, q = qfirst)
            pwo1 <- pmlt(object, newdata = newdata, q = qwo1)
            pwoK <- pmlt(object, newdata = newdata, q = qwoK)
            ret <- matrix(0, nrow = length(first), ncol = ncol(pfirst))
            ret[!first,] <- pwo1 - pwoK
            ret[first,] <- pfirst
       }
    } else {

        ### need to generate newdata outside tmlt such that
        ### the rows of expand.grid(newdata) match the elements of
        ### the return value

        stopifnot(is.null(q))
        stopifnot(y %in% names(newdata))
        dim <- sapply(newdata, NROW)
        q <- newdata[[y]]

        first <- q == lev[1]
        qfirst <- q[first]
        qwoK <- q[q != lev[length(lev)]]
        qwo1 <- q[q != lev[1]]

        newdata[[y]] <- qfirst
        pfirst <- pmlt(object, newdata = newdata)
        newdata[[y]] <- qwo1
        pwo1 <- pmlt(object, newdata = newdata)
        newdata[[y]] <- qwoK
        pwoK <- pmlt(object, newdata = newdata)

        dn <- dim(pfirst)
        names(dn) <- names(dimnames(pfirst))

        frst <- lapply(names(dn), function(d) {
            if (d != y) return(1:dn[d])
            return(first)
        })
        ntfrst <- lapply(names(dn), function(d) {
            if (d != y) return(1:dn[d])
            return(!first)
        })

        dn <- lapply(names(dn), function(d) {
            if (d != y) return(dimnames(pfirst)[[d]])
            return(as.character(q))
        })

        ret <- array(0, dim = dim, dimnames = dn)
        ret <- do.call("[<-", 
            c(list(i = ret), frst, list(value = pfirst)))
        ret <- do.call("[<-", 
            c(list(i = ret), ntfrst, list(value = pwo1 - pwoK)))
    }

    if (log) return(log(ret))
    return(ret)
}

### hazard function
hmlt <- function(object, newdata = object$data, q = NULL, n = 50, 
                 log = FALSE) {
    if (log) 
        return(dmlt(object, newdata = newdata, q = q, n = n, log = log) -
               log(smlt(object, newdata = newdata, q = q, n = n)))
    return(dmlt(object, newdata = object$data, q = q, n = n, log = log) /
           smlt(object, newdata = object$data, q = q, n = n))
}
