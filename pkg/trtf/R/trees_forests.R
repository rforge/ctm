
.ctmfit <- function(object, parm, mltargs) {
    
    ctmobject <- object

    function(formula, data, weights, cluster, ctrl, ...) {
        weights <- model.weights(data)
        if (is.null(weights)) weights <- integer(0)
        f <- Formula(formula)    
        mf <- model.frame(formula = f, data = data, na.action = na.pass)
        if (ctrl$nmax < Inf) {
            bdr <- inum::inum(mf, total = TRUE, nmax = ctrl$nmax, complete.cases.only = TRUE,
                              as.interval = ctmobject$response)
            iy <- c(bdr)
            mf1 <- as.data.frame(bdr)
            attr(iy, "levels") <- 1:nrow(mf1)
            w <- libcoin::ctabs(iy, weights = weights)[-1L]
        } else {
            mf1 <- mf
            w <- weights
            if (length(w) == 0) w <- rep(1, nrow(mf1))
            iy <- NULL
        }
        mltargs$data <- mf1  
        mltargs$weights <- w     
        ctmobject <- do.call("mlt", mltargs)
        thetastart <- coef(ctmobject)

        function(subset = NULL, info = NULL, newweights = NULL, model = FALSE, 
                 estfun = TRUE, object = FALSE) {
            if (model) return(list(object = ctmobject, iy = iy))
            if (is.null(newweights)) {
                if (ctrl$nmax < Inf) {
                    w <- libcoin::ctabs(iy, weights = weights, subset = subset)[-1L]
                    subset <- NULL
                } else {
                    w[-subset] <- 0 ### mlt >= 1.0-3 allows subset but we
                    ### still need the weights to be zero for some operations below
                }
            } else {
                w <- newweights
            }
            if (!is.null(info$coef)) thetastart <- info$coef
            umod <- suppressWarnings(try(update(ctmobject, weights = w, subset = subset, theta = thetastart), silent = TRUE))
            if (inherits(umod, "try-error") || umod$convergence != 0) {
                umod <- suppressWarnings(try(update(ctmobject, weights = w, subset = subset), silent = TRUE))
                if (inherits(umod, "try-error") || umod$convergence != 0) {
                    mltargs$weights <- w
                    ### no subset allowed here, so used zero weights (see
                    ### above)!!!
                    umod <- try(do.call("mlt", mltargs))
                }
            }
            if (inherits(umod, "try-error")) {
                ### we badly need some estimate in each node, even if fitting
                ### fails
                if (!estfun) 
                    return(list(coef = thetastart, objfun = NA,
                                converged = FALSE))
                return(list(estfun = matrix(0, nrow = nrow(mf1), ncol = length(w)),
                            iy = iy, converged = FALSE))
            }
            ret <- NULL
            if (estfun) {
                ret <- estfun(umod)[, parm, drop = FALSE]
                if (!is.null(subset)) {
                    tmp <- matrix(0, nrow = length(w), 
                                  ncol = ncol(ret))
                    tmp[subset,] <- ret
                    ret <- tmp
                }
                ret <- ret / w ### we need UNWEIGHTED scores
                ret[w == 0,] <- 0
            }
            return(list(estfun = ret, index = iy, 
                        coefficients = coef(umod), objfun = logLik(umod), 
                        object = if (object) umod else NULL,
                        converged = isTRUE(all.equal(umod$convergence, 0))))
        }
    }
} 

trafotree <- function(object, parm = 1:length(coef(object)), mltargs = list(maxit = 10000), ...) {

    mltargs$model <- object
    args <- list(...)
    args$ytrafo <- .ctmfit(object, parm, mltargs)
    ret <- do.call("ctree", args)
    ret$model <- object
    ret$mltobj <- ret$trafo(model = TRUE, estfun = FALSE)
    ret$mltargs <- mltargs

    ### store coefs and logLik _outside_ tree
    ### <FIXME> this will cause problems with nodeprune </FIXME>
    nd <- predict(ret, type = "node")
    ret$models <- tapply(1:length(nd), factor(nd), function(i) 
        ret$trafo(i, estfun = FALSE)) ### note: trafo is (potentially) weighted
    ret$coef <- do.call("rbind", lapply(ret$models, function(x) x$coef))
    ret$logLik <- sapply(ret$models, function(x) x$objfun)

    class(ret) <- c("trafotree", class(ret))
    ret
}

traforest <- function(object, parm = 1:length(coef(object)), mltargs = list(maxit = 10000), ...) {

    mltargs$model <- object
    args <- list(...)
    args$ytrafo <- .ctmfit(object, parm, mltargs)
    ret <- do.call("cforest", args)
    ret$model <- object
    ret$mltargs <- mltargs
    ret$mltobj <- ret$trafo(model = TRUE, estfun = FALSE, object = TRUE)
    class(ret) <- c("traforest", class(ret))
    ret
}
