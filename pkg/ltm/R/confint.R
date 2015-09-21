
confint.ltm <- function(object, parm, level = 0.95, 
                        type = c("trafo", "distribution", "survivor", "cumhazard"), n = 20, ...) {

    cf <- coef(object)
    if (is.null(cf))
        parm <- data.frame(var = 1)
    if (missing(parm))
        parm <- names(cf)
    type <- match.arg(type)

    if (is.character(parm)) {
        K <- diag(length(cf))
        rownames(K) <- names(cf)
        K <- K[parm,,drop = FALSE]
        return(confint(glht(multcomp::parm(cf, vcov(object)), linfct = K), 
                       level = level, ...))
    } 

    ### confint.mlt?
    stopifnot(is.data.frame(parm))
    if (nrow(parm) > 1) {
        ret <- lapply(1:nrow(parm), function(i) 
            confint(object = object, parm = parm[i,,drop = FALSE], 
                    level = level, type = type, ...))
        return(ret)
    }

    stopifnot(nrow(parm) == 1)
    y <- object$response
    q <- mkgrid(object, n = n)[[y]]
    parm <- parm[rep(1, length(q)),,drop = FALSE]
    parm[[y]] <- q
    X <- model.matrix(object$model$model, data = parm)
    ci <- confint(glht(parm(coef(object, lp_only = FALSE), 
                            vcov(object, lp_only = FALSE)), 
                       linfct = X), ...)$confint
    if (type == "distribution") ci <- object$model$todistr$p(ci)
    if (type == "survivor") ci <- 1 - object$model$todistr$p(ci)
    if (type == "cumhazard") ci <- -log(1 - object$model$todistr$p(ci))
    ci <- cbind(q = q, ci)
    return(ci)
}
