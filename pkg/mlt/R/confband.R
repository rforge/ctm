
confband <- function(object, newdata, level = .95, ...)
    UseMethod("confband")

confband.mlt <- function(object, newdata, level = 0.95,
                        type = c("trafo", "distribution", "survivor", "cumhazard"), n = 20, ...) {

    stopifnot(!missing(newdata))
    stopifnot(is.data.frame(newdata))
    type <- match.arg(type)
    if (nrow(newdata) > 1) {
        ret <- lapply(1:nrow(newdata), function(i)
            confband(object = object, newdata = newdata[i,,drop = FALSE],
                     level = level, type = type, ...))
        return(ret)
    }

    stopifnot(nrow(newdata) == 1)
    y <- object$response
    q <- mkgrid(object, n = n)[[y]]
    newdata <- newdata[rep(1, length(q)),,drop = FALSE]
    newdata[[y]] <- q
    X <- model.matrix(object$model$model, data = newdata)
    ci <- confint(multcomp::glht(multcomp::parm(coef(object), vcov(object)),
                                 linfct = X), ...)$confint
    if (type == "distribution") ci <- object$model$todistr$p(ci)
    if (type == "survivor") ci <- 1 - object$model$todistr$p(ci)
    if (type == "cumhazard") ci <- -log(1 - object$model$todistr$p(ci))
    ci <- cbind(q = q, ci)
    return(ci)
}

