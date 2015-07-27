
ltm <- function(formula, data, subset, weights, na.action = na.omit,
                method = c("logistic", "probit", "cloglog"),
                trafo = c("fixed log", "log", "Bernstein", "Legendre"),
                order = 5, support = NULL, ...) {

    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(mf), 0)
    mf <- mf[c(1, m)]
 
    formula <- Formula::Formula(formula)

    if (length(formula)[2L] < 2L) {
        formula <- Formula::Formula(formula(Formula::as.Formula(formula(formula), ~ 0), 
                                    rhs = 2L:1L))
        xreg <- FALSE
    } else {
        if( length(formula)[2L] > 2L) {
            formula <- Formula::Formula(formula(formula, rhs = 1L:2L))
            warning("Formula must not have more than two RHS parts")
        }
        xreg <- TRUE
        strata <- (length(formula)[2] == 2)
    }
    mf$formula <- formula
    mf$drop.unused.levels <- FALSE
    mf$na.action <- na.action
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- terms(formula, data = data)
    if (xreg)
        mtX <- delete.response(terms(formula, data = data, rhs = 1L))
    if (strata) {
        stratum <- names(Formula::model.part(formula, mf, rhs = 2L))
        stopifnot(length(stratum) == 1)
    }

    Y <- Formula::model.part(formula, mf, lhs = 1L)
    response <- names(Y)
    Y <- Y[[response]]
    weights <- model.weights(mf)

    method <- match.arg(method)
    todistr <- switch(method, "logistic" = "Logistic",
                              "probit" = "Normal",
                              "cloglog" = "MinExtVal")
    if (is.factor(Y)) {
        nlev <- nlevels(Y)
        cntr <- list(function(n) contr.treatment(n, base = nlev))
        names(cntr) <- response
        rtrafo <- as.basis(as.formula(paste("~", response)), data = mf, remove_intercept = TRUE,
              contrasts.arg = cntr, ui = diff(diag(nlev - 1)), ci = 0)
    } else {
        support <- range(unlist(Y))
        trafo <- match.arg(trafo)
        rtrafo <- switch(trafo, 
            "fixed log" = log_basis(ui = "increasing", varname = response, support = support),
            "log" = log_basis(ui = "increasing", varname = response, support = support),
            "Bernstein" = Bernstein_basis(order = order, support = support,
                          ui = "increasing", varname = "response"),
            "Legendre" = Legendre_basis(order = order, support = support,
                          ui = "increasing", varname = "response"))
    }

    strafo <- xtrafo <- NULL
    if (strata)
        strafo <- as.basis(as.formula(paste("~ ", stratum, "- 1")), data = mf)
    if (xreg)
        xtrafo <- as.basis(mtX, data = mf, remove_intercept = TRUE)

    m <- model(rtrafo, interacting = strafo, shifting = xtrafo, todistr = todistr)
    mlt(m, data = mf, ...)
}

library("mlt")
library("Formula")
coef(ltm(Sepal.Length ~ Sepal.Width | Species, data = iris))
