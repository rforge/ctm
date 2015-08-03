
predict.ltm <- function(object, newdata, type = c("lp", "trafo", "distribution", 
    "survivor", "density", "logdensity", "hazard", "loghazard", "cumhazard", "quantile"),
    q = NULL, n = 50, ...) {

    class(object) <- class(object)[-1]

    if (missing(newdata))
        newdata <- object$data
    stopifnot(is.data.frame(newdata))

    type <- match.arg(type)
    if (type == "lp") {
        if (!is.null(q)) warning("argument q ignored")
        q <- mkgrid(object, n = 1)[[object$response]] 
        return(drop(predict(object, newdata = newdata, 
                            type = "trafo", terms = "bshifting", q = q)))
    }

    if (is.null(q))
        q <- mkgrid(object, n = n)[[object$response]]

    predict(object, newdata = newdata, type = type, q = q, ...)
}
