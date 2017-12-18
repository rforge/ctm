
plot.tram <- function(x, newdata = model.frame(x), 
    which = c("QQ-PIT", "baseline only", "distribution"), ...) {

    which <- match.arg(which)
    object <- as.mlt(x)
    y <- newdata[[variable.names(x, "response")]]

    censored <- inherits(y, "Surv") || inherits(y, "response")

    if (which != "distribution" && censored)
        stop("Cannot compute in-sample ", which, " for censored responses")
        
    ret <- switch(which, 
        "QQ-PIT" = {
            U <- predict(object, newdata = newdata, type = "distribution")
            qqplot(U, 1:length(U) / length(U), ...)
            qqline(U, distribution = qunif)
        },
        "baseline only" = {
            scf <- object$shiftparm
            if (length(scf) > 0) {
                mobj <- as.mlt(object)
                cf <- coef(mobj)
                cf[scf] <- 0
                coef(mobj) <- cf
                plot(mobj, newdata = newdata, type = "trafo", ...)
            } else {
                plot(object, newdata = newdata, type = "trafo", ...)
            }
        },
        "distribution" = {
            plot(object, newdata = newdata, ...)
        })
}
