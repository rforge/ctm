
ltm <- function(formula, data, subset, weights, na.action = na.omit,
                method = c("logistic", "probit", "cloglog"),
                trafo = c("Bernstein", "Legendre", "fixed log", "log"),
                integerAsFactor = FALSE, order = 5, support = NULL, 
                negative_lp = TRUE, contrasts.arg = NULL, ...) {

    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(mf), 0)
    mf <- mf[c(1, m)]
 
    formula <- Formula::Formula(formula)

    stopifnot(length(formula)[1] == 1L)   
    stopifnot(length(formula)[2] %in% 1:2)

    strata <- FALSE
    if (length(formula)[2] == 2) strata <- TRUE

    mf$formula <- formula
    mf$drop.unused.levels <- FALSE
    mf$na.action <- na.action
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- terms(formula, data = data)
    mtX <- delete.response(terms(formula, data = data, rhs = 1L))
    if (strata) {
        stratum <- Formula::model.part(formula, mf, rhs = 2L)
        stopifnot(length(stratum) == 1)
        stopifnot(is.factor(stratum[[1]]) & nlevels(stratum[[1]]) > 1)
    }
    weights <- model.weights(mf)
    Y <- Formula::model.part(formula, mf, lhs = 1L)
    response <- names(Y)
    Y <- Y[[response]]
    if (inherits(Y, "Surv")) {
        Y <- mlt:::.Surv2R(Y)
        rclass <- "Surv"
    } else {
        rclass <- class(Y)[1]
    }
    if (inherits(Y, "response")) {
        if (is.null(support)) {
            tmp <- unlist(Y[, c("cleft", "exact", "cright")])
            support <- range(tmp[is.finite(tmp)], na.rm = TRUE)
            if (rclass == "Surv") support[1] <- 0
        rclass <- ifelse (all(is.na(Y$exact)), class(Y$cleft), class(Y$exact))
        }
    } else {
        if (is.null(support))
            support <- range(Y[is.finite(Y)])
    }
    rtype <- c("continuous", "discrete")[rclass %in% c("factor", "integer") + 1L]

    method <- match.arg(method)
    todistr <- switch(method, "logistic" = "Logistic",
                              "probit" = "Normal",
                              "cloglog" = "MinExtrVal")

    fixed <- NULL
    if (rtype == "discrete" & ifelse(rclass == "integer", integerAsFactor, TRUE)) {
        if (rclass == "integer") Y <- ordered(Y)
        nlev <- nlevels(Y)
        cntr <- list(function(n) contr.treatment(n, base = nlev))
        names(cntr) <- response
        rtrafo <- as.basis(as.formula(paste("~", response)), data = mf, remove_intercept = TRUE,
              contrasts.arg = cntr, ui = diff(diag(nlev - 1)), ci = 0)
    } else {
        trafo <- match.arg(trafo)
        rtrafo <- switch(trafo, 
            "fixed log" = log_basis(ui = "increasing", varname = response, support = support),
            "log" = log_basis(ui = "increasing", varname = response, support = support),
            "Bernstein" = Bernstein_basis(order = order, support = support,
                          ui = "increasing", varname = response),
            "Legendre" = Legendre_basis(order = order, support = support,
                          ui = "increasing", varname = response))
        if (trafo == "fixed log") {
            fixed <- 1
            names(fixed) <- names(model.matrix(rtrafo, data = mf))
        }
    }

    strafo <- xtrafo <- NULL
    if (strata)
        strafo <- as.basis(as.formula(paste("~ ", names(stratum), "- 1")), data = mf)
    xtrafo <- as.basis(mtX, data = mf, remove_intercept = TRUE, contrasts.arg = contrasts.arg,
                       negative = negative_lp)
    if (ncol(model.matrix(xtrafo, data = mf)) == 0) xtrafo <- NULL

    m <- model(rtrafo, interacting = strafo, shifting = xtrafo, todistr = todistr)

    ret <- mlt(m, data = mf, fixed = fixed, scale = TRUE, check = FALSE, ...)
    class(ret) <- c("ltm", class(ret))
    ret
}

model.frame.ltm <- function(object, ...)
    object$data

model.matrix.ltm <- function(object, data = model.frame(object), all = FALSE, ...) {
    if (all) {
        class(object) <- class(object)[-1]
        return(model.matrix(object, data = data, ...))
    }
    if (is.null(object$model$model$bshifting))
        return(NA)
    model.matrix(object$model$model$bshifting, data = data, ...)
}

coef.ltm <- function(object, all = FALSE, ...) {
    if (all) {
        class(object) <- class(object)[-1]
        return(coef(object, ...))
    }      
    xn <- colnames(model.matrix(object))
    if (length(xn) == 0) return(NA)
    class(object) <- class(object)[-1]
    coef(object)[xn]
}

vcov.ltm <- function(object, all = FALSE, ...) {
    if (all) {
        class(object) <- class(object)[-1]
        return(vcov(object, ...))
    }      
    xn <- colnames(model.matrix(object))
    if (length(xn) == 0) return(NA)
    class(object) <- class(object)[-1]
    vcov(object)[xn, xn, drop = FALSE]
}

AIC.ltm <- function(object, ..., k = 2) {
    class(object) <- class(object)[-1]
    AIC(object, ..., k = k)
}

logLik.ltm <- function(object, ...) {
    class(object) <- class(object)[-1]
    logLik(object, ...)
}

print.ltm <- function(x, ...) 
    print(cftest(x))

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

library("mlt")
library("Formula")
library("multcomp")
coef(a <- ltm(Sepal.Length ~ Sepal.Width | Species, data = iris, trafo = "log"))

### predict(a)

library("survival")
data("GBSG2", package = "TH.data")
storage.mode(GBSG2$time) <- "double"

b <- ltm(time ~ horTh , data = GBSG2)

predict(b)

predict(b, newdata = data.frame(horTh = unique(GBSG2$horTh)), q = 100:110, type = "trafo")

b <- ltm(time ~ 1, data = GBSG2)

predict(b, q = 100:110, type = "distribution")

b <- ltm(time ~ 1 | tgrade, data = GBSG2)

predict(b, newdata = data.frame(tgrade = unique(GBSG2$tgrade)), q = 100:110, type = "distribution")

cc <- ltm(Surv(time, cens) ~ horTh + menostat + pnodes | tgrade, data = GBSG2, method = "cloglog")

predict(cc, newdata = expand.grid(horTh = c("no", "yes"), menostat = "Pre",
                           pnodes = 100, tgrade = "II"), q = 100:101, type = "trafo")

d <- ltm(Surv(time, cens) ~ 1, data = GBSG2, method = "cloglog")
class(d) <- class(d)[-1]
cf <- coef(d)
v <- vcov(d)
prm <- parm(cf, v)
K <- diag(length(cf))
rownames(K) <- colnames(K) <- names(cf)
ci <- confint(glht(prm, linfct = K), calpha = qnorm(.975))

lwr <- d
class(lwr) <- class(lwr)[-1]
upr <- d
class(upr) <- class(upr)[-1]
coef(lwr) <- ci$confint[, "lwr"]
coef(upr) <- ci$confint[, "upr"]

tm <- 10:2700

s <- seq(from = min(GBSG2$time), to = max(GBSG2$time), length = length(cf))
plot(tm, predict(d, q = tm), type = "l", ylim = range(ci$confint))
points(s, cf, col = "red")

lines(tm, predict(lwr, q = tm))
points(s, coef(lwr), col = "blue")

lines(tm, predict(upr, q = tm))
points(s, coef(upr), col = "green")

plot(tm, 1 - predict(d, q = tm, type = "distr"), type = "l")

lines(tm, 1 - predict(lwr, q = tm, type = "distr"))

lines(tm, 1 - predict(upr, q = tm, type = "distr"))



plot(survfit(Surv(time, cens) ~ 1, data = GBSG2))
lines(tm, 1 - predict(d, q = tm, type = "distr"))
lines(tm, 1 - predict(lwr, q = tm, type = "distr"))
lines(tm, 1 - predict(upr, q = tm, type = "distr"))

