
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

    mlt::confband(as.mlt(object), newdata = parm, level = level, type = type, n = n, ...)
}
