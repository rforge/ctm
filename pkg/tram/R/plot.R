
plot.tram <- function(x, newdata = x$data, 
    which = c("QQresiduals", "baseline", "distribution"), ...) {

    which <- match.arg(which)
    object <- as.mlt(x)

    if (which == "QQresiduals") {
        U <- predict(object, newdata = newdata, type = "distribution")
        qqplot(U, 1:length(U) / length(U), ...)
        qqline(U, distribution = qunif)
    }
    if (which == "baseline") {
        B <- predict(object, newdata = newdata, type = "trafo") - 
             predict(object, newdata = newdata, terms = "bshifting")
        y <- newdata[[variable.names(x, "response")]]
        plot(y, B, ...)
    }
    if (which == "distribution")
        plot(object, newdata = newdata, ...)
}
