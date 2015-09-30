
plot.mlt <- function(x, newdata, type = c("distribution",
    "survivor", "density", "logdensity", "hazard", "loghazard", "cumhazard", "quantile", "trafo"),
    q = NULL, p = 1:(n-1) / n, n = 50, col = rgb(.1, .1, .1, .1), add = FALSE, ...) {

    if (is.null(q))
        q <- mkgrid(x, n = n)[[x$response]]
    type <- match.arg(type)
    pr <- predict(x, newdata = newdata, type = type, q = q, p = p, ...)
    pr[!is.finite(pr)] <- NA
    rpr <- range(pr, na.rm = TRUE)
    if (is.null(dim(pr))) pr <- matrix(pr, ncol = 1)
    ylim <- switch(type, "distribution" = c(0, 1),
                         "survivor" = c(0, 1),
                         "density" = c(0, rpr[2]),
                         "hazard" = c(0, rpr[2]),
                         "cumhazard" = c(0, rpr[2]),
                         rpr)
    if (type == "quantile")  q <- p
    if (length(col) == 1) col <- rep(col, ncol(pr))
    
    if (!add) {
        plot(unclass(q), rep(rpr[1], length(q)), ylim = ylim, xlab = x$response,
             ylab = type, type = "n", axes = FALSE)
        if (is.factor(q)) {
            axis(1, at = unclass(q), labels = levels(q))
        } else {
            axis(1)
        }
        axis(2)
        box()
    }
    y <- as.vars(x)[[x$response]]
    if (inherits(y, "continuous_var")) {
        for (i in 1:ncol(pr)) lines(q, pr[,i], col = col[i])
    } else {
        for (i in 1:ncol(pr)) lines(stepfun(q, c(ylim[1], pr[,i])), col = col[i])
    }
    invisible(pr)
}
