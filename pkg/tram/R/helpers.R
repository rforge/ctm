
PI <- function(logOR, prob) {

    stopifnot(xor(missing(logOR), missing(prob)))

    if (missing(prob)) {
        OR <- exp(logOR)
        ret <- OR * (OR - 1 - logOR) / (OR - 1)^2
        ret[abs(logOR) < .Machine$double.eps] <- .5
        return(ret)
    }

    logOR <- 1:999 / 50
    s <- spline(x = logOR, y = PI(logOR = logOR), method = "hyman")
    wl5 <- (prob < .5 - .Machine$double.eps)
    wg5 <- (prob > .5 + .Machine$double.eps)
    ret <- numeric(length(prob))
    if (any(wl5))
        ret[wl5] <- - approx(x = s$y, y = s$x, xout = 1 - prob[wl5])$y
    if (any(wg5))
        ret[wg5] <- approx(x = s$y, y = s$x, xout = prob[wg5])$y
    ret
}


