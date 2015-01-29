
.Normal <- function()
    list(p = pnorm, d = dnorm, dd = ddnorm)

.Logistic <- function()
    list(p = plogis, d = dlogis, dd = function(x)
         (2 * exp(-x)^2 / (1 + exp(-x))^3) - (exp(-x) / (1 + exp(-x))^2))

.MinExtrVal <- function()
    list(p = function(x) 1 - exp(-exp(x)),
         d = function(x, log = FALSE) {
             ret <- x - exp(x)
             if (!log) return(exp(ret))
             ret
         },
         dd = function(x)
             exp(x - exp(x)) * (1 - exp(x)))

Distr <- function(which = c("Normal", "Logistic", "MinExtrVal")) {
    which <- match.arg(which)
    do.call(paste(".", which, sep = ""), list())
}

