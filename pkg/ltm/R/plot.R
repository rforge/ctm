
### plot.mlt?

plot.ltm <- function(x, newdata, type = c("distribution",
    "survivor", "density", "logdensity", "hazard", "loghazard", "cumhazard", "quantile", "trafo"),
    q = NULL, p = 1:9 / 10, n = 50, col = rgb(.1, .1, .1, .1), ...) {

    if (is.null(q))
        q <- mkgrid(x, n = n, bounds = x$bounds)[[x$response]]
    type <- match.arg(type)
    pr <- predict(x, newdata = newdata, type = type, q = q, p = p, ...)
    pr[!is.finite(pr)] <- NA
    rpr <- range(pr, na.rm = TRUE)
    ylim <- switch(type, "distribution" = c(0, 1),
                         "survivor" = c(0, 1),
                         "density" = c(0, rpr[2]),
                         "hazard" = c(0, rpr[2]),
                         "cumhazard" = c(0, rpr[2]),
                         rpr)
    if (type == "quantile")  q <- p
    
    plot(q, rep(rpr[1], length(q)), ylim = ylim, xlab = x$response,
         ylab = type, type = "n")
    if (is.double(q)) {
        out <- apply(pr, 2, function(y) lines(q, y, col = col))
    } else {
        out <- apply(pr, 2, function(y) points(q, y, col = col))
    }
}

