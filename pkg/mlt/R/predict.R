
predict.mlt <- function(object, newdata = object$data, type = c("trafo", "distribution", "survivor", 
    "density", "logdensity", "hazard", "loghazard", "cumhazard", "quantile"), 
    terms = c("all", "bresponse", "binteracting", "bshifting"), q = NULL, p = NULL, n = 50,
    interpolate = TRUE, ...) {

    type <- match.arg(type)
    terms <- match.arg(terms)
    if (terms != "all") {
        stopifnot(type == "trafo")
    } else {
       terms <- NULL
    }

    if (type == "quantile")
        stopifnot(!is.null(p) && (min(p) > 0 & max(p) < 1))

    if (!is.data.frame(object))
        stopifnot(object$respone %in% names(object))

    ret <- switch(type, 
        "trafo" = tmlt(object = object, newdata = newdata, q = q, terms = terms, ...),
        "distribution" = pmlt(object = object, newdata = newdata, q = q),
        "survivor" = smlt(object = object, newdata = newdata, q = q),
        "density" = dmlt(object = object, newdata = newdata, q = q, log = FALSE),
        "logdensity" = dmlt(object = object, newdata = newdata, q = q, log = TRUE),
        "hazard" = hmlt(object = object, newdata = newdata, q = q, log = FALSE),
        "loghazard" = hmlt(object = object, newdata = newdata, q = q, log = TRUE),
        "cumhazard" = Hmlt(object = object, newdata = newdata, q = q),
        "quantile" = qmlt(object = object, newdata = newdata, n = n,
                          p = p, interpolate = interpolate))

    return(ret)
}
