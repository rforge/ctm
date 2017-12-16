
tram_data <- function(formula, data, subset, weights, offset, cluster, na.action = na.omit) 
{

    ## set up model.frame() call
    mf <- match.call(expand.dots = FALSE)
    mf$na.action <- na.action ### evaluate na.action
    if(missing(data)) data <- environment(formula)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$dot <- "sequential"

    ## formula processing
    oformula <- as.formula(formula)
    formula <- Formula::as.Formula(formula)
    mf$formula <- formula
    npart <- length(formula)
    if(any(npart < 1L)) stop("'formula' must specify at least one left-hand and one right-hand side")
    if(any(npart > 2L)) stop("'formula' must specify at least one left-hand and one right-hand side")

    ## evaluate model.frame
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    ## extract terms in various combinations
    mt <- list(
      "all" = terms(formula, data = data,                        dot = "sequential"),
      "y"   = terms(formula, data = data, lhs = 1L, rhs = 0L,    dot = "sequential"),
      "s"   = if (npart[1] == 2) 
              terms(formula, data = data, lhs = 2L, rhs = 0L,    dot = "sequential"),
      "x"   = terms(formula, data = data, lhs = 0L, rhs = 1L,    dot = "sequential"),
      "z"   = if (npart[2] == 2) 
              terms(formula, data = data, lhs = 0L, rhs = 2L,    dot = "sequential")
    )

    stopifnot(is.null(mt$z))

    ### Surv(...) etc. will be altered anyway...
    # names(mf) <- make.names(names(mf))

    response <- mf[[1L]]
    weights <- model.weights(mf)
    offset <- model.offset(mf)
    cluster <- mf[["(cluster)"]]

    ret <- list(response = response, weights = weights, offset = offset, cluster =
                cluster, mf = mf, mt = mt)
    class(ret) <- "tram_data"
    ret
}

tram <- function(formula, data, subset, weights, offset, cluster, na.action = na.omit,
                 distribution = c("Normal", "Logistic", "MinExtrVal"),
                 transformation = c("discrete", "linear", "logarithmic", "smooth"),
                 LRtest = TRUE, 
                 prob = c(.1, .9), support = NULL, order = 6, negative = FALSE, scale =
                 TRUE, asFamily = FALSE, model_only = FALSE, ...) 
{

    if (!inherits(td <- formula, "tram_data")) {
        mf <- match.call(expand.dots = FALSE)
        m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
        mf <- mf[c(1L, m)]
        mf[[1L]] <- quote(tram_data)
        td <- eval(mf, parent.frame())
    } 

    rvar <- asvar(td$response, names(td$mf)[1L], prob = prob, support = support)
    rbasis <- mkbasis(rvar, transformation = transformation, order = order)

    iS <- NULL
    if (!is.null(td$mt$s)) 
        iS <- as.basis(formula(Formula(td$mt$s)[-3]), data = td$mf)
    iX <- NULL
    if (!is.null(td$mt$x)) 
        iX <- as.basis(td$mt$x, data = td$mf, remove_intercept = !asFamily, 
                       negative = negative)

    model <- ctm(response = rbasis, interacting = iS, shifting = iX, 
                 todistr = distribution, data = td$mf)

    if (model_only) return(model)

    if (asFamily) return(ctmFamily(model, td$mf))

    ret <- mlt(model, data = td$mf, weights = td$weights, offset = td$offset, 
               scale = scale, ...)
    ret$terms <- td$terms
    ret$cluster <- td$cluster
    ret$shiftcoef <- colnames(model.matrix(iX, data = td$mf))
    class(ret) <- c("tram", class(ret))

    if (LRtest & !is.null(iX)) {
        nullmodel <- ctm(response = rbasis, interacting = iS, 
                         todistr = distribution, data = td$mf)
        nullret <- mlt(nullmodel, data = td$mf, weights = td$weights, offset = td$offset, 
                       scale = scale, ...)
        nulllogLik <- logLik(nullret)
        fulllogLik <- logLik(ret)
        ret$LRtest <- c(LRstat = -2 * (nulllogLik - fulllogLik), 
                        df = attr(fulllogLik, "df") - attr(nulllogLik, "df"))
    }

    ret
}
