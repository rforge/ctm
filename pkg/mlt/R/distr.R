
.Normal <- function()
    list(parm = function(x) NULL,
         p = pnorm, d = dnorm, q = qnorm, 
         ### see also MiscTools::ddnorm
         dd = function(x) -dnorm(x = x) * x,
         ddd = function(x) dnorm(x = x) * (x^2 - 1), 
         dd2d = function(x) -x,
         call = ".Normal",
         name = "normal")

.Exponential <- function()
    list(parm = function(x) NULL,
         p = pexp, d = dexp, q = qexp, 
         dd = function(x) -dexp(x = x),
         ddd = function(x) dexp(x = x),
         dd2d = function(x) -1,
         call = ".Exponential",
         name = "exponential")

.Logistic <- function()
    list(parm = function(x) NULL,
         p = plogis, d = dlogis, q = qlogis,
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
         call = ".Logistic",
         name = "logistic")

### Gompertz distribution
.MinExtrVal <- function()
    list(parm = function(x) NULL,
         p = function(x) 1 - exp(-exp(x)),
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
         call = ".MinExtrVal",
         name = "minimum extreme value")

### Gumbel distribution
.MaxExtrVal <- function()
    list(parm = function(x) NULL,
         p = function(x) exp(-exp(-x)),
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
         call = ".MaxExtrVal",
         name = "maximum extreme value")

### see 10.1080/15598608.2013.772835
.GammaFrailty <- function(logrho = 0) {
    logrho <- pmax(logrho, log(sqrt(.Machine$double.eps)))
    list(parm = function() c("logrho" = logrho),
         p = function(x) 1 - (1 + exp(x + logrho))^(-exp(-logrho)),
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
         call = ".GammaFrailty",
         name = paste0("GammaFrailty(rho = ", 
                       round(exp(logrho), 
                       options("digits")$digits), ")"))
}

### see 10.1002/sim.687
.InvGaussFrailty <- function(logtheta = 0) {
    logtheta <- pmax(logtheta, log(sqrt(.Machine$double.eps)))
    list(parm = function() c("logtheta" = logtheta),
         p = function(x)
             1 - exp(- sqrt(4 * exp(logtheta) * (exp(logtheta) + exp(x))) + 2 * exp(logtheta)),
         q = function(p) {
             theta <- exp(logtheta)
             log((-log1p(- p) + 2 * theta)^2 / (4 * theta) - theta)
         },
         d = .d <- function(x, log = FALSE) {
             ret <- (logtheta - 2 * sqrt(exp(logtheta) * (exp(x) + exp(logtheta))) + x + 2 * exp(logtheta)) - 
                     .5 * (logtheta + log(exp(x) + exp(logtheta)))
             if (!log) return(exp(ret))
             ret
         },
         dd = .dd <- function(x) {
             theta <- exp(logtheta)
             stx <- sqrt(theta * (exp(x) + theta))
             ret <- theta * (1 - theta * exp(x) / stx) * exp(-2 * stx + x + 2 * theta)
             ret <- ret / stx
             ret <- ret - theta^2 * exp(-2 * stx + 2 * x + 2 * theta) / (2 * stx * stx^2)
             ret
         },
         ddd = function(x) {
             theta <- exp(logtheta)
             txt <- theta * (exp(x) + theta)
             stx <- sqrt(txt)
             ret <- 3 * theta^3 * exp(-2 * stx + 3 * x + 2 * theta) / (4 * txt^(5/2))
             ret <- ret - theta^2 * (2 - theta * exp(x) / stx) * exp(-2 * stx + 2 * x + 2 * theta) / (2 * txt^(3 / 2))
             ret <- ret - theta^2 * (1 - theta * exp(x) / stx) * exp(-2 * stx + 2 * x + 2 * theta) / (2 * txt^(3 / 2))
             ret <- ret + theta * (theta^2 * exp(2 * x) / (2 * txt^(3/2)) - theta * exp(x) / stx) * exp(-2 * stx + x + 2 * theta) / stx
             ret <- ret + theta * (1 - theta * exp(x) / stx)^2 * exp(-2 * stx + x + 2 * theta) / stx
             ret
         },
         dd2d = function(x) {
             .dd(x) / .d(x)
         },
         call = ".InvGaussFrailty",
         name = paste0("InvGaussFrailty(theta = ", 
                       round(exp(logtheta), 
                       options("digits")$digits), ")"))
}

### Cure rate models: 10.1002/sim.687
.Cure <- function(logitrho = 0, ..., frailty = .GammaFrailty) {
    f <- frailty(...)
    list(parm = function() c("logitrho" = logitrho, f$parm()),
         p = function(x) plogis(logitrho) + plogis(logitrho) * f$p(x),
         q = function(p)
             f$q(p / plogis(logisrho)),
         d = .d <- function(x, log = FALSE) {
             if (log)
                 return(plogis(logitrho, log.p = TRUE) + f$d(x, log = TRUE))
             plogis(logitrho) * f$d(x)
         },
         dd = .dd <- function(x)
             plogis(logitrho) * f$dd(x),
         ddd = function(x)
             plogis(logitrho) * f$ddd(x),
         dd2d = function(x)
             f$dd2d(x),
         call = ".CureRate",
         name = paste0("CureRate(rho = ", 
                       round(plogis(logitrho), 
                       options("digits")$digits), ")"))
}

.distr <- function(which = c("Normal", "Logistic", 
                             "MinExtrVal", "MaxExtrVal", "Exponential")) {
    which <- match.arg(which)
    do.call(paste(".", which, sep = ""), list())
}

