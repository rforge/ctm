
.Normal <- function()
    list(p = pnorm, d = dnorm, q = qnorm, 
         ### see also MiscTools::ddnorm
         dd = function(x) -dnorm(x = x) * x,
         ddd = function(x) dnorm(x = x) * (x^2 - 1), 
         dd2d = function(x) -x,
         name = "normal")

.Exponential <- function()
    list(p = pexp, d = dexp, q = qexp, 
         dd = function(x) -dexp(x = x),
         ddd = function(x) dexp(x = x),
         dd2d = function(x) -1,
         name = "exponential")

.Logistic <- function()
    list(p = plogis, d = dlogis, q = qlogis,
         dd = function(x) {
             ex <- exp(x)
             (ex - ex^2) / (1 + ex)^3
         },
         ddd = function(x) {
             ex <- exp(x)
             (ex - 4 * ex^2 + ex^3) / (1 + ex)^4
         },
         dd2d = function(x) {
             ex <- exp(x)
             (1 - ex) / (1 + ex)
         },
         name = "logistic")

### Gompertz distribution
.MinExtrVal <- function()
    list(p = function(x) 1 - exp(-exp(x)),
         q = function(p) log(-log1p(- p)),
         d = function(x, log = FALSE) {
             ret <- x - exp(x)
             if (!log) return(exp(ret))
             ret
         },
         dd = function(x) {
             ex <- exp(x)
             (ex - ex^2) / exp(ex)
         },
         ddd = function(x) {
             ex <- exp(x)
             (ex - 3*ex^2 + ex^3) / exp(ex)
         },
         dd2d = function(x)
             1 - exp(x),
         name = "minimum extreme value")

### Gumbel distribution
.MaxExtrVal <- function()
    list(p = function(x) exp(-exp(-x)),
         q = function(p) -log(-log(p)),
         d = function(x, log = FALSE) {
             ret <- - x - exp(-x)
             if (!log) return(exp(ret))
             ret
         },
         dd = function(x) {
             ex <- exp(-x)
             exp(-ex - x) * (ex - 1)
         },
         ddd = function(x) {
             ex <- exp(-x)
             exp(-x - ex) * (ex - 1)^2 - exp(-ex - 2 * x)
         },
         dd2d = function(x)
             exp(-x) - 1,
         name = "maximum extreme value")

### see 10.1080/15598608.2013.772835
.Logarithmic <- function(logrho = 0) {
    if (logrho < log(sqrt(.Machine$double.eps)))
        return(.MinExtrVal())
    list(p = .F <- function(x) 1 - (1 + exp(x + logrho))^(-exp(-logrho)),
         q = function(p)
             log((1 - p)^(-exp(logrho)) - 1) - logrho,
         d = .d <- function(x, log = FALSE) {
             ret <- x + (-exp(-logrho) - 1) * log(exp(x + logrho) + 1)
             if (!log) return(exp(ret))
             ret
         },
         dd = .dd <- function(x) {
             exlr <- exp(x + logrho)
             memlr <- -exp(-logrho)
             (memlr - 1) * (exlr + 1)^(memlr - 2) * exp(2 * x + logrho) + 
               exp(x) * (exlr + 1)^(memlr - 1)
         },
         ddd = function(x) {
             exlr <- exp(x + logrho)
             memlr <- -exp(-logrho)
             (memlr - 2) * (memlr - 1) * (exlr + 1)^(memlr - 3) * exp(3 * x + 2 * logrho) +
               3 * (memlr - 1) * (exlr + 1)^(memlr - 2) * exp(2 * x + logrho) +
               exp(x) * (exlr + 1)^(memlr - 1)
         },
         dd2d = function(x) {
             .dd(x) / .d(x)
         },
         name = paste0("logarithmic(rho = ", 
                       round(exp(logrho), 
                       options("digits")$digits), ")"))
}

.distr <- function(which = c("Normal", "Logistic", 
                             "MinExtrVal", "MaxExtrVal", "Exponential")) {
    which <- match.arg(which)
    do.call(paste(".", which, sep = ""), list())
}

