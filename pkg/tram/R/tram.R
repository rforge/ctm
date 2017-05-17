
### y | strata ~ x | z
### cluster: Hessian with cluster weights, sum-up
### variances

tram <- function(formula, data, subset, weights, offset,
                 order = 1, distribution = c("Normal", "Logistic", "MinExtrVal"),
                 algorithm = c("ML", "boosting", "tree", "forest"), ...) 
{

    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    names(mf)[names(mf) %in% c("subset", "weights", "offset")]  <- 
        paste("(", names(mf)[names(mf) %in% c("subset", "weights", "offset")], ")", sep = "")
    m <- match(c("formula", "data", "(subset)", "(weights)", "(offset)"),
               names(mf), 0)
    mf <- mf[c(1, m)]

    formula <- Formula::Formula(formula)
    mf$formula <- formula
    mf[[1]] <- quote(stats::get_all_vars)
    mf <- eval(mf, parent.frame())

    mfterms <- terms(formula, data = data, dot = "sequential") 
    ### there might be dots in formula, fdot
    ### is formula with dots replaced
    fdot <- attr(mfterms, "Formula_without_dot")
    if (!is.null(fdot)) formula <- fdot

    lfm <- length(formula)
    stopifnot(length(lfm) == 2)
    stopifnot(lfm[1] == 1)
    yfm <- formula(formula, lhs = 1, rhs = 0)

    xfm <- formula(formula, lhs = 0, rhs = 1)
    sfm <- NULL
    if (lfm[2] >= 2)
        sfm <- formula(formula, lhs = 0, rhs = 2)
    rfm <- NULL
    if (lfm[2] >= 3)
        rfm <- formula(formula, lhs = 0, rhs = 3)
    zfm <- NULL
    if (lfm[2] == 4)
        zfm <- formula(formula, lhs = 0, rhs = 4)

    return(list(xfm, sfm, rfm, zfm, mf))
}


tram(y ~ ., 
     data = data.frame(y = 1:3, x = 2:4, s = 3:5, r = 4:6, z = 5:7),
     weights = rep(1, 3),
     offset = 0,
     subset = z > 6)

