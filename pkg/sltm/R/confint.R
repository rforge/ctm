
confint.sltm <- function(object, parm, level = 0.95, ...) {

    cf <- coef(object)
    if (is.null(cf))
        parm <- data.frame(var = 1)
    if (missing(parm))
        parm <- names(cf)
    stopifnot(all(parm %in% names(cf)))
    K <- diag(length(cf))
    rownames(K) <- names(cf)
    K <- K[parm,,drop = FALSE]
    return(confint(glht(multcomp::parm(cf, vcov(object)), linfct = K), 
                   level = level, ...))
}

confband.sltm <- function(object, newdata, level = 0.95,
                          type = c("trafo", "distribution", "survivor", "cumhazard"), 
                          n = 20, ...)
    mlt::confband(as.mlt(object), newdata = newdata, 
                  level = level, type = type, n = n, ...)
