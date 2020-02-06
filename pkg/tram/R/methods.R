
as.mlt.tram <- function(object) {
    cls <- which(class(object) == "mlt_fit")
    class(object) <- class(object)[-(1:(cls - 1))]
    object
}    

update.tram <- function(object, ...)
    update(as.mlt(object), ...)

model.frame.tram <- function(formula, ...) {
    ret <- formula$data
    ### <FIXME>: we need different options here:
    ### model.frame for all variables, for x and z, etc.
    attr(ret, "terms") <- formula$terms$x
    ### </FIXME>
    ret
}

terms.tram <- function(x, ...)
    terms(model.frame(x))

model.matrix.tram <- function(object, data = object$data, 
                              with_baseline = FALSE, ...) 
{
    if (with_baseline) 
        return(model.matrix(as.mlt(object), data = data, ...))
    if (is.null(object$shiftcoef)) 
        return(NULL)
    ret <- model.matrix(as.mlt(object)$model$model$bshifting, 
                        data = data, ...)
    ret
}	

coef.tram <- function(object, with_baseline = FALSE, ...) 
{
    cf <- coef(as.mlt(object), ...)
    if (with_baseline) return(cf)
    if (is.null(object$shiftcoef)) return(NULL)
    return(cf[names(cf) %in% object$shiftcoef])
}

coef.Lm <- function(object, as.lm = FALSE, ...) {

    class(object) <- class(object)[-1L]
    if (!as.lm)
        return(coef(object, ...))

    if (!is.null(object$stratacoef))
        stop("cannot compute scaled coefficients with strata")

    cf <- coef(object, with_baseline = TRUE, ...)
    cfx <- coef(object, with_baseline = FALSE, ...)
    cfs <- cf[!(names(cf) %in% names(cfx))]
    sd <- 1 / cfs[names(cfs) != "(Intercept)"]

    if (is.null(object$shiftcoef)) {
        ret <- -cfs["(Intercept)"] * sd
    } else {
        ret <- c(-cfs["(Intercept)"], cfx) * sd
    }
    attr(ret, "scale") <- sd
    ret
}

coef.Survreg <- function(object, as.survreg = FALSE, ...)
    coef.Lm(object, as.lm = as.survreg, ...)
        
vcov.tram <- function(object, with_baseline = FALSE, complete = FALSE, ...) 
{

    if (complete) stop("complete not implemented")

    ### full covariance matrix
    if (with_baseline) {
        if (is.null(object$cluster)) {
            return(vcov(as.mlt(object), ...))
        } else {
            return(sandwich::vcovCL(as.mlt(object), 
                                    cluster = object$cluster))
        }
    }
    if (is.null(object$shiftcoef)) return(NULL)
   
    ### covariance matrix for shift terms only
    ### return Schur complement
    H <- Hessian(as.mlt(object), ...)
    shift <- which(colnames(H) %in% object$shiftcoef)
    Hlin <- H[shift, shift]
    Hbase <- H[-shift, -shift]
    Hoff <- H[shift, -shift]
    H <- try(Hlin - tcrossprod(Hoff %*% solve(Hbase), Hoff))
    if (inherits(H, "try-error"))
        return(vcov(as.mlt(object))[shift, shift])
    ret <- solve(H)
    if (inherits(ret, "try-error"))
        return(vcov(as.mlt(object))[shift, shift])
    colnames(ret) <- rownames(ret) <- object$shiftcoef
    return(ret)
}

nobs.tram <- function(object, ...) {
    if (!is.null(object$weights)) 
        return(sum(object$weights != 0))
    return(NROW(object$data))
}

logLik.tram <- function(object, parm = coef(as.mlt(object), fixed = FALSE), ...)
    logLik(as.mlt(object), parm = parm, ...)

Hessian.tram <- function(object, parm = coef(as.mlt(object), fixed = FALSE), ...)
    Hessian(as.mlt(object), parm = parm, ...)

Gradient.tram <- function(object, parm = coef(as.mlt(object), fixed = FALSE), ...)
   Gradient(as.mlt(object), parm = parm, ...)

estfun.tram <- function(object, parm = coef(as.mlt(object), fixed = FALSE), ...)
    estfun(as.mlt(object), parm = parm, ...)

predict.tram <- function(object, newdata = model.frame(object), 
    type = c("lp", "trafo", "distribution", "survivor", "density", 
             "logdensity", "hazard", "loghazard", "cumhazard", "quantile"), ...) {

    type <- match.arg(type)
    if (type == "lp") {
        ret <- model.matrix(object, data = newdata) %*% 
               coef(object, with_baseline = FALSE)
        if (object$negative) return(-ret)
        return(ret)
    }
    predict(as.mlt(object), newdata = newdata, type = type, ...)
}

print.tram <- function(x, ...) {
    cat("\n", x$tram, "\n")
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(coef(x, with_baseline = FALSE))
    ll <- logLik(x)
    cat("\nLog-Likelihood:\n ", ll, " (df = ", attr(ll, "df"), ")", sep = "")
    cat("\n\n")
    invisible(x)
}

summary.tram <- function(object, ...) {
    ret <- list(call = object$call,
                tram = object$tram,
                test = cftest(object, parm = names(coef(object, with_baseline = FALSE))),
                ll = logLik(object))
    if (!is.null(object$LRtest)) {
        ret$LRstat <- object$LRtest["LRstat"]
        ret$df <- floor(object$LRtest["df"])
        ret$p.value <- pchisq(object$LRtest["LRstat"], 
                              df = object$LRtest["df"], lower.tail = FALSE)
    }
    class(ret) <- "summary.tram"
    ret
}

print.summary.tram <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\n", x$tram, "\n")
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    pq <- x$test$test
    mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
    colnames(mtests) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    sig <- .Machine$double.eps
    printCoefmat(mtests, digits = digits, has.Pvalue = TRUE, 
        P.values = TRUE, eps.Pvalue = sig)
    cat("\nLog-Likelihood:\n ", x$ll, " (df = ", attr(x$ll, "df"), ")", sep = "")
    if (!is.null(x$LRstat))
        cat("\nLikelihood-ratio Test: Chisq =", x$LRstat, "on",
            x$df, "degrees of freedom; p =", format.pval(x$p.value, digits = digits, ...))
    cat("\n\n")
    invisible(x)
}

residuals.tram <- function(object, ...)
    residuals(as.mlt(object), ...)

# profile.tram was modified from
# File MASS/profiles.q copyright (C) 1996 D. M. Bates and W. N. Venables.
#
# port to R by B. D. Ripley copyright (C) 1998
#
# corrections copyright (C) 2000,3,6,7 B. D. Ripley

profile.tram <- function(fitted, which = 1:p, alpha = 0.01,
                         maxsteps = 10, del = zmax/5, trace = FALSE, ...)
{
    Pnames <- names(B0 <- coef(fitted))
    nonA <- !is.na(B0)
    pv0 <- t(as.matrix(B0))
    p <- length(Pnames)
    if(is.character(which)) which <- match(which, Pnames)

    diff <- sqrt(diag(vcov(fitted)))
    names(diff) <- Pnames

    X <- model.matrix(fitted)
    n <- NROW(X)
    O <- fitted$offset
    if(!length(O)) O <- rep(0, n)
    W <- fitted$weights
    if(length(W) == 0L) W <- rep(1, n)

    OriginalDeviance <- -2 * logLik(fitted)
    zmax <- sqrt(qchisq(1 - alpha, 1))

    prof <- vector("list", length=length(which))
    names(prof) <- Pnames[which]
    for(i in which) {
        if(!nonA[i]) next
        zi <- 0
        pvi <- pv0
        pi <- Pnames[i]
        fx <- 0
        names(fx) <- pi

        ### set-up new model with fixed parameter pi such that
        ### update can be called
        theta <- coef(as.mlt(fitted))
        theta <- theta[!names(theta) %in% pi]
        fm <- m <- mlt(fitted$model, data = fitted$data, weights = W, 
                       fixed = fx, theta = theta,
                       scale = fitted$scale, 
                       optim = fitted$optim)

        for(sgn in c(-1, 1)) {
            if(trace)
                message("\nParameter: ", pi, " ",
                        c("down", "up")[(sgn + 1)/2 + 1])
            step <- 0
            z <- 0

            while((step <- step + 1) < maxsteps && abs(z) < zmax) {
                bi <- B0[i] + sgn * step * del * diff[i]
                o <- O + X[, i] * bi

                ### compute profile likelihood
                fm <- update(m, weights = W, offset = o, 
                             theta = coef(fm, fixed = FALSE))
                pl <- logLik(fm)

                ri <- pv0
                ri[, Pnames] <- coef(fm)[Pnames]
                ri[, pi] <- bi
                pvi <- rbind(pvi, ri)
                zz <- -2 * pl - OriginalDeviance

                if(zz > - 1e-3) zz <- max(zz, 0)
                else stop("profiling has found a better solution, so original fit had not converged")

                z <- sgn * sqrt(zz)
                zi <- c(zi, z)
            }
        }
        si <- order(zi)
        prof[[pi]] <- structure(data.frame(zi[si]), names = "z")
        prof[[pi]]$par.vals <- pvi[si, ,drop=FALSE]
    }
    val <- structure(prof, original.fit = fitted, summary = NULL)
    ### inheriting from profile.glm is maybe not so smart but
    ### _currently_ works when calling MASS::confint.profile.glm
    class(val) <- c("profile.tram", "profile.glm", "profile")
    val
}

### score tests and confidence intervals
### hm, can't find such a generic
score_test <- function(object, ...)
    UseMethod("score_test")

score_test.default <- function(object, ...)
    stop("no score_test method implemented for class", class(object)[1])

score_test.tram <- function(object, parm = names(coef(object)), 
    alternative = c("two.sided", "less", "greater"), nullvalue = 0, 
    confint = TRUE, level = .95, maxsteps = 25, ...) {

    cf <- coef(object)
    stopifnot(all(parm %in% names(cf)))
    alternative <- match.arg(alternative)
    
    if (length(parm) > 1) {
        ret <- lapply(parm, score_test, object = object, 
                      alternative = alternative, 
                      nullvalue = nullvalue,
                      confint = confint, 
                      level = level, 
                      maxsteps = maxsteps, 
                      ...)
        names(ret) <- parm
        class(ret) <- "htests"
        return(ret)
    }

    m1 <- as.mlt(object)

    fx <- 0
    names(fx) <- parm
    off <- object$offset
    theta <- coef(as.mlt(object))
    theta <- theta[names(theta) != parm]
    m0 <- mlt(object$model, data = object$data, weights = object$weights,
              offset = off, scale = object$scale, fixed = fx,
              mltoptim = object$optim, theta = theta)

    cf <- coef(m1)
    X <- model.matrix(object)[, parm]

    sc <- function(b) {
        cf[] <- coef(update(m0, offset = off + b * X))
        cf[parm] <- b
        coef(m1) <- cf
        ### see Lehmann, Elements of Large-sample Theory,
        ### 1999, 539-540, for score tests in the presence

        ### of additional parameters
        sum(estfun(m1)[, parm]) * sqrt(vcov(m1)[parm, parm])
    }
    stat <- c("Z" = sc(nullvalue))
    pval <- switch(alternative, 
        "two.sided" = pnorm(-abs(stat)) * 2,
        "less" = pnorm(-abs(stat)),
        "greater" = pnorm(abs(stat)))

    if (confint) {
        alpha <- (1 - level)
        if (alternative == "two.sided") alpha <- alpha / 2

        Sci <- coef(object) + sqrt(vcov(object)[parm, parm]) * qp


if (FALSE) {
        Wci <- confint(object, level = 1 - alpha / 5)[parm,]
        grd <- seq(from = Wci[1], to = Wci[2], length.out = maxsteps)
        grd_sc <- sapply(grd, sc) 
        method <- "fmm"
        ### theoretically, sc is monotone (increasing or decreasing)
        ### but the optimiser might not always know this
        if (all(diff(grd_sc) < 0) || all(diff(grd_sc) > 0))
            method <- "hyman"
        s <- spline(x = grd, y = grd_sc, method = method)
        Sci <- approx(x = s$y, y = s$x, 
                      xout = qnorm(c(alpha, 1 - alpha)))$y
}
        est <- coef(object)[parm]
        attr(Sci, "conf.level") <- level
        if (alternative == "less")
            Sci[1] <- -Inf
        if (alternative == "greater")
            Sci[2] <- Inf
    }

    parameter <- paste(switch(class(object)[1], 
                                  "Colr" = "Log-odds ratio",
                                   "Coxph" = "Log-hazard ratio",
                                   "Lm" = "Standardised difference",
                                   "Lehmann" = "Lehmann parameter",
                                   "BoxCox" = "Standardised difference"))
    parameter <- paste(tolower(parameter), "for", parm)
    names(nullvalue) <- parameter

    ret <- list(statistic = stat,
                p.value = pval, 
                null.value = nullvalue, 
                alternative = alternative, 
                method = "Transformation Score Test",
                data.name = deparse(object$call))
    if (confint) {
        ret$conf.int <- Sci
        names(est) <- parameter
        ret$estimate <- est
    }
    class(ret) <- "htest"
    ret
}

confint.htest <- function(object, parm, level = .95, ...) {
    ci <- object$conf.int
    if (is.null(ci)) return(NULL)
    stopifnot(attr(ci, "conf.level") == level)
    return(ci)
}

confint.htests <- function(object, parm, level = .95, ...) {
    ret <- do.call("rbind", lapply(object, confint))
    rownames(ret) <- names(object)
    alt <- unique(sapply(object, function(x) x$alternative))
    stopifnot(length(alt) == 1)
    rn <- switch(alt, 
        "two.sided" = {
            a <- (1 - level)/2
            a <- c(a, 1 - a)
            paste(round(100 * a, 1), "%")
        }, 
        "less" = c("", paste(round(100 * level, 1), "%")),
        "greater" = c(paste(round(100 * level, 1), "%"), ""))
    colnames(ret) <- rn
    ret
}

### permutation score tests and confidence intervals
perm_test <- function(object, ...)
    UseMethod("perm_test")

perm_test.default <- function(object, ...)
    stop("no perm_test method implemented for class", class(object)[1])

perm_test.tram <- function(object, parm = names(coef(object)), 
    statistic = c("Score", "LinearScore", "Likelihood", "Wald"),
    alternative = c("two.sided", "less", "greater"), 
    nullvalue = 0, confint = FALSE, level = .95, 
    block_permutation = TRUE, maxsteps = 25, ...) {

    cf <- coef(object)
    stopifnot(all(parm %in% names(cf)))
    statistic <- match.arg(statistic, several.ok = TRUE)
    SCORE <- grep("Score", statistic)
    if (length(SCORE) > 0) statistic <- statistic[SCORE[1]]
    alternative <- match.arg(alternative)

    parameter <- paste(switch(class(object)[1], 
                                  "Colr" = "Log-odds ratio",
                                   "Coxph" = "Log-hazard ratio",
                                   "Lm" = "Standardised difference",
                                   "Lehmann" = "Lehmann parameter",
                                   "BoxCox" = "Standardised difference"))

    block <- NULL
    if (!is.null(object$model$bases$interacting) && block_permutation) {
        svar <- variable.names(object$model$bases$interacting)
        block <- object$data[, svar]
        if (is.data.frame(block))
            block <- do.call("interaction", block)
    }
    if (length(grep("Score", statistic)) > 0) {

        stopifnot(isTRUE(all.equal(nullvalue, 0)))

        if (length(parm) > 1) {
            ret <- lapply(parm, perm_test, object = object, 
                          alternative = alternative, 
                          nullvalue = nullvalue,
                          confint = confint, 
                          level = level, 
                          block_permutation = block_permutation,
                          ...)
            names(ret) <- parm
            class(ret) <- "htests"
            return(ret)
        }

        m1 <- as.mlt(object)

        fx <- 0
        names(fx) <- parm
        off <- object$offset
        theta <- coef(as.mlt(object))
        theta <- theta[names(theta) != parm]
        w <- object$weights
        m0 <- mlt(object$model, data = object$data, weights = object$weights,
                  offset = off, scale = object$scale, fixed = fx,
                  mltoptim = object$optim, theta = theta)

        cf <- coef(m1)
        X <- Xf <- model.matrix(object)[, parm]
        ### this is a hack and not really necessary but
        ### distribution = "exact" needs factors @ 2 levels
        ### use baseline level 1 such that p-values are increasing
        if (length(unique(X)) == 2)
            Xf <- relevel(factor(X, levels = sort(unique(X)), 
                          labels = 0:1), "1")

        cf[] <- c(coef(m0, fixed = FALSE), 0)
        coef(m1) <- cf
        ### resid is weighted, remove weights and feed them to coin
        r0 <- (resid(m1) / w) * sqrt(vcov(m1)[parm,parm])
        if (is.null(block)) {
            it0 <- coin::independence_test(r0 ~ Xf, teststat = "scalar", 
                alternative = alternative, weights = ~ w, ...)
        } else {
            it0 <- coin::independence_test(r0 ~ X | block, teststat = "scalar", 
                alternative = alternative, weights = ~ w, ...)
        }
        stat <- c("Z" = coin::statistic(it0, "standardized"))
        pval <- coin::pvalue(it0)

        if (confint) {

            alpha <- (1 - level)
            if (alternative == "two.sided") alpha <- alpha / 2

            ### we always have Prob(Q(alpha) <= S) >= alpha
            ### for alpha < .5, we need Prob(Q(alpha) <= S) <= alpha
            qa <- coin::qperm(it0, p = alpha)
            if (coin::pperm(it0, q = qa) > alpha - 1e3) {
                sprt <- coin::support(it0)
                if (!all(is.na(sprt))) {
                    qa <- max(sprt[sprt < qa])
                } else {
                    a <- alpha
                    while(coin::pperm(it0, qa) > alpha) {
                        a <- a - 1e-4
                        qa <- coin::qperm(it0, p = a)
                    }
                }
            }
            q1a <- coin::qperm(it0, p = 1 - alpha)
            qp <- c(qa, q1a)
            achieved <- coin::pperm(it0, q = qp)
            achieved <- 1 - (achieved[1] + (1 - achieved[2]))
            qp <- qp * sqrt(coin::variance(it0)) +
                coin::expectation(it0)

            est <- coef(object)[parm]

            if ("LinearScore" %in% statistic) {
                Sci <- coef(object) + sqrt(vcov(object)[parm, parm]) * qp
            } else {
                sc <- function(b) {
                    cf[] <- coef(update(m0, offset = off + b * X, 
                                        theta = theta))
                    cf[parm] <- b
                    coef(m1) <- cf
                    ### see Lehmann, Elements of Large-sample Theory,
                    ### 1999, 539-540, for score tests in the presence
                    ### of additional parameters
                    ### estfun is already weighted
                    sum(estfun(m1)[, parm]) * sqrt(vcov(m1)[parm, parm])
                }
                Wci <- confint(object, level = 1 - alpha / 5)[parm,]
                grd <- seq(from = Wci[1], to = Wci[2], length.out = maxsteps)
                grd_sc <- sapply(grd, sc) 
                method <- "fmm"
                ### theoretically, sc is monotone (increasing or decreasing)
                ### but the optimiser might not always know this
                if (all(diff(grd_sc) < 0) || all(diff(grd_sc) > 0))
                    method <- "hyman"
                s <- spline(x = grd, y = grd_sc, method = method)
                Sci <- approx(x = s$y, y = s$x, xout = qp)$y
            } 
            
            attr(Sci, "conf.level") <- level
            attr(Sci, "achieved.conf.level") <- achieved
            if (alternative == "less")
                Sci[1] <- -Inf
            if (alternative == "greater")
                Sci[2] <- Inf
        }

        distname <- switch(class(it0@distribution),
            "AsymptNullDistribution" = "Asymptotic",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution"  = "Exact"
        )

        parameter <- paste(tolower(parameter), "for", parm)
        names(nullvalue) <- parameter
  
        ret <- list(statistic = stat,
                    p.value = pval, 
                    null.value = nullvalue, 
                    alternative = alternative, 
                    method = paste(distname, "Permutation Transformation Score Test"),
                    data.name = deparse(object$call))
        if (confint) {
            ret$conf.int <- Sci
            names(est) <- parameter
            ret$estimate <- est
        }
        class(ret) <- "htest"
        return(ret)

    } else {

        if ("Wald" %in% statistic) {
            stopifnot(length(parm) == 1)
        } else {
            stopifnot(alternative == "two.sided")
        }

        if (length(nullvalue) != length(parm))
            nullvalue <- rep(nullvalue, length(parm))

        if (confint && !isTRUE(all.equal(nullvalue, coef(object)[parm],
                                         check.attributes = FALSE)))
            nullvalue <- coef(object)[parm]

        off <- object$offset + c(1, -1)[object$negative + 1L] *
            model.matrix(object)[, parm, drop = FALSE] %*% nullvalue

        theta <- coef(as.mlt(object), fixed = FALSE)
        thetafix <- coef(as.mlt(object), fixed = TRUE)
        fx <- NULL
        if (length(thetafix) > length(theta))
            fx <- thetafix[-names(theta)]
        m0 <- mlt(object$model, data = object$data, weights = object$weights,
                  offset = off, scale = object$scale, fixed = fx,
                  mltoptim = object$optim, theta = theta)
        
        if (length(list(...)) == 0) {
            nperm <- 1000
        } else {
            nperm <- sapply(list(...), function(x) {
                np <- get("nresample", environment(x))
                if (is.numeric(np)) return(np)
                return(NA)
            })
            if (!all(is.na(nperm))) nperm <- max(nperm, na.rm = TRUE)
        }

        ll <- numeric(nperm)
        cf <- matrix(0, nrow = nperm, ncol = length(parm))
        colnames(cf) <- parm

        idx <- perm <- 1:NROW(object$data)
        if (!is.null(block))
            sidx <- do.call("c", tapply(idx, block, function(i) i))
        for (b in 1:nperm) {
            if (!is.null(block)) {
                perm[sidx] <- do.call("c", tapply(idx, block, sample))
            } else {
                perm <- sample(idx)
            }

            um0 <- update(m0, perm = parm, permutation = perm)
            ll[b] <- -um0$value
            cf[b,] <- coef(um0)[parm]
        }

        if ("Likelihood" %in% statistic) {
            stat <- logLik(object)
            sull <- sort(unique(ll))
            s <- spline(x = sull, y = ecdf(ll)(sull), method = "hyman")
            pval <- NA
            if (!confint)
                pval <- approx(x = s$x, y = 1 - s$y, xout = stat)$y
            qlevel <- approx(x = s$y, y = s$x, xout = level)$y
            if (confint) {
                if(length(parm) > 1)
                    return(list(stat = stat, pval = pval, 
                                confcoef = t(t(cf[ll < qlevel, ]) + nullvalue)))
                conf.int <- range(cf[ll < qlevel, 1]) + nullvalue
            }
            retL <- list(statistic = stat,
                p.value = pval, 
                null.value = nullvalue, 
                alternative = alternative, 
                method = "Permutation Transformation Likelihood Test",
                data.name = deparse(object$call))
            if (confint) {
                retL$conf.int <- conf.int
                retL$estimate <- coef(object)[parm]
                names(retL$estimate) <- parameter
            }
            class(retL) <- "htest"
            if (length(statistic) == 1)
                return(retL)
        } 

        if ("Wald" %in% statistic) {
            stat <- coef(object)
            sucf <- sort(unique(cf))
            s <- spline(x = sucf, y = ecdf(cf)(sucf), method = "hyman")
            pval <- NA
            if (!confint) {
                pval <- approx(x = s$x, y = s$y, xout = stat)$y
                if (alternative == "greater") pval <- 1 - pval
                if (alternative == "two.sided") 
                    pval <- 2 * min(pval, 1 - pval)
            }
            alpha <- 1 - level
            if (alternative == "two.sided")
                alpha <- alpha / 2
            qlevel <- approx(x = s$y, y = s$x, xout = c(alpha, 1 - alpha))$y
            if (confint) {
                conf.int <- nullvalue + qlevel
                if (alternative == "less") confint[1] <- -Inf
                if (alternative == "greater") confint[1] <- -Inf
            }
            retW <- list(statistic = stat,
                p.value = pval, 
                null.value = nullvalue, 
                alternative = alternative, 
                method = "Permutation Transformation Wald Test",
                data.name = deparse(object$call))
            if (confint) {
                retW$conf.int <- conf.int
                retW$estimate <- coef(object)[parm]
                names(retW$estimate) <- parameter
            }
            class(retW) <- "htest"
            if (length(statistic) == 1)
                return(retW)
         }
         return(list(Likelihood = retL, Wald = retW))
    }
}
