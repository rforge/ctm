
confint.ltm <- function(object, parm, level = 0.95, method = c("univariate", "adjusted"), 
                        type = c("trafo", "distribution", "survivor"), n = 20, ...) {

    cf <- coef(object)
    method <- match.arg(method)
    if (method == "univariate") {
        calpha <- univariate_calpha()
    } else {
        calpha <- adjusted_calpha()
    }
    type <- match.arg(type)

    if (missing(parm)) {
        K <- diag(length(cf))
        rownames(K) <- names(cf)
        return(confint(glht(multcomp::parm(cf, vcov(object)), linfct = K), 
                       calpha = calpha, ...))
    } 

    stopifnot(is.data.frame(parm))
    stopifnot(nrow(parm) == 1)
    y <- object$response
    q <- mkgrid(object, n = n)[[y]]
    parm <- parm[rep(1, length(q)),,drop = FALSE]
    parm[[y]] <- q
    X <- model.matrix(object$model$model, data = parm)
    ci <- confint(glht(parm(coef(object, all = TRUE), 
                            vcov(object, all = TRUE)), 
                       linfct = X), calpha = calpha, ...)
    attr(ci, "q") <- q
    if (type == "trafo") return(ci)
    ### sapply(ci$confint...)
    if (type == "distribution") return(object$model$todistr$p(ci))
    if (type == "survivor") return(1 - object$model$todistr$p(ci))
}
