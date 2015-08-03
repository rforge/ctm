
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

