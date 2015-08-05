
predict.ltm <- function(object, newdata, type = c("lp", "trafo", "distribution", 
    "survivor", "density", "logdensity", "hazard", "loghazard", "cumhazard", "quantile"),
    q = NULL, n = 50, ...) {

    cf <- coef(object, lp_only = TRUE)
    object <- as.mlt(object)

    if (missing(newdata))
        newdata <- object$data
    stopifnot(is.data.frame(newdata))

    type <- match.arg(type)
    if (type == "lp") {
        if (!is.null(q)) warning("argument q ignored")
        if (is.null(cf)) return(NULL)
        return(drop(predict(object$model$model$bshifting, newdata = newdata, 
                            coef = cf)))
    }

    if (is.null(q))
        q <- mkgrid(object, n = n, bounds = object$bounds)[[object$response]]

    predict(object, newdata = newdata, type = type, q = q, ...)
}
