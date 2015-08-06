
ltm <- function(formula, data, subset, weights, na.action = na.omit,
                method = c("logistic", "probit", "cloglog"),
                trafo = c("Bernstein", "Legendre", "fixed log", "log", "linear"),
                integerAsFactor = TRUE, order = 5, support = NULL, bounds = c(-Inf, Inf),
                negative_lp = TRUE, contrasts.arg = NULL, ...) {

    ### extract response and possibly a stratum and predictor variables
    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    formula <- Formula::Formula(formula)
    ## formula is `response | stratum ~ x'
    stopifnot(length(formula)[1] %in% 1:2)
    STRATA <- (length(formula)[1] == 2)
    stopifnot(length(formula)[2] == 1L)   
    mf$formula <- formula
    mf$drop.unused.levels <- FALSE
    mf$na.action <- na.action
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- terms(formula, data = data)
    mtX <- delete.response(terms(formula, data = data, lhs = 0, rhs = 1L))
    if (STRATA) {
        stratum <- Formula::model.part(formula, mf, lhs = 2L, rhs = 0)
        stopifnot(length(stratum) == 1)
        stopifnot(is.factor(stratum[[1]]) & nlevels(stratum[[1]]) > 1)
    }
    weights <- model.weights(mf)
    Y <- Formula::model.part(formula, mf, lhs = 1L, rhs = 0)
    response <- names(Y)
    Y <- Y[[response]]

    SURV <- inherits(Y, "Surv")
    if (!inherits(Y, "response"))
        RY <- R(Y, bounds = bounds)
    type <- attr(RY, "type")
    DISCRETE <- type != "double" 
    if (type == "integer") DISCRETE <- integerAsFactor
    
    fixed <- NULL
    if (DISCRETE) {
        if (type == "integer") Y <- mf[[response]] <- ordered(Y)
        nlev <- nlevels(Y)
        cntr <- list(function(n) contr.treatment(n, base = nlev))
        names(cntr) <- response
        rtrafo <- as.basis(as.formula(paste("~", response)), data = mf, remove_intercept = TRUE,
              contrasts.arg = cntr, ui = diff(diag(nlev - 1)), ci = rep(0, nlev - 2))
    } else {
        if (is.null(support)) {
            tmp <- unlist(RY[, c("cleft", "exact", "cright")])
            support <- range(tmp[is.finite(tmp)], na.rm = TRUE)
            interval <- quantile(tmp[is.finite(tmp)], na.rm = TRUE, prob = c(.025, .975))
            if (SURV) bounds[1] <- 0
        }
        if (is.null(match.call()$trafo)) trafo <- "Bernstein"
        trafo <- match.arg(trafo, several.ok = TRUE)
        rtrafo <- lapply(trafo, function(tr) switch(tr, 
            "fixed log" = log_basis(ui = "increasing", varname = response, support = support),
            "log" = log_basis(ui = "increasing", varname = response, support = support),
            "Bernstein" = Bernstein_basis(order = order, support = support, interval = interval,
                          ui = "increasing", varname = response),
            "Legendre" = Legendre_basis(order = order, support = support, interval = interval,
                          ui = "increasing", varname = response),
            "linear" = polynomial_basis(c(TRUE, TRUE), varname = response, 
                                        support = support, ci = c(-Inf, 0)))
        )
        if ("fixed log" %in% trafo) {
            fixed <- 1
            names(fixed) <- names(model.matrix(rtrafo[["fixed log"]], data = mf))
        }
        if (length(rtrafo) > 1) {
            rtrafo <- do.call("c", rtrafo)
        } else {
            rtrafo <- rtrafo[[1L]]
        }
    }

    strafo <- xtrafo <- NULL
    if (STRATA)
        strafo <- as.basis(as.formula(paste("~ ", names(stratum), "- 1")), data = mf)
    xtrafo <- as.basis(mtX, data = mf, remove_intercept = TRUE, contrasts.arg = contrasts.arg,
                       negative = negative_lp)
    if (ncol(model.matrix(xtrafo, data = mf)) == 0) xtrafo <- NULL

    method <- match.arg(method)
    todistr <- switch(method, "logistic" = "Logistic",
                              "probit" = "Normal",
                              "cloglog" = "MinExtrVal")

    m <- model(rtrafo, interacting = strafo, shifting = xtrafo, todistr = todistr)

    ret <- mlt(m, data = mf, fixed = fixed, bounds = bounds, scale = TRUE, check = FALSE, ...)
    ret$call <- match.call(expand.dots = TRUE)

    class(ret) <- c(ifelse(DISCRETE, "dltm", "cltm"), "ltm", class(ret))
    ret
}

