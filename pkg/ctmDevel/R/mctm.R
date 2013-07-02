	
upsilon <- function(y, ngrid = NULL) {

    if (is.factor(y)) {
        if (!is.ordered(y) && nlevels(y) > 2) 
            stop("multiclass ctm not implemented")
        return(list(y = (1:nlevels(y))[y], 
                    upsilon = 1:nlevels(y),
                    uorig = factor(levels(y), ordered = is.ordered(y))))
    } 
    upsilon <- sort(unique(y))
    upsilon <- c(upsilon[1] - (upsilon[2] - upsilon[1]), upsilon)
    if (!is.null(ngrid) && !is.ordered(y)) {
        if (length(ngrid) > 1) {
            upsilon <- ngrid
        } else {
            upsilon <- seq(from = min(upsilon), to = max(upsilon),
                           length = ngrid)
        }
    }
    return(list(y = y, upsilon = upsilon, uorig = upsilon))
}

### multivariate ctm
mctm <- function(formula, data, weights = NULL, constant = NULL,
                 varying = NULL, ngrid = NULL, fit = TRUE, offset = 0, 
                 ...) {

    yname <- all.vars(formula[[2]])
    response <- data[yname]
    ### cdf <- ecdf(response)
    for (y in yname) data[[y]] <- NULL

    yu <- lapply(response, upsilon, ngrid = ngrid)
    mresponse <- as.matrix(sapply(yu, function(x) x$y))
    uresponse <- do.call("expand.grid", 
                         lapply(yu, function(x) x$upsilon))
    oresponse <- do.call("expand.grid", 
                         lapply(yu, function(x) x$uorig))

    ### need inverse!
    ### offset <- family@cdf(uresponse)

    if (all(sapply(response, is.factor))) {
        uresponse <- uresponse[-nrow(uresponse),,drop = FALSE]
        oresponse <- oresponse[-nrow(oresponse),,drop = FALSE]
        for (y in yname) oresponse[[y]] <- oresponse[[y]][,drop = TRUE]
    }
    
    dresponse <- apply(uresponse, 1, function(u) {
	utmp <- matrix(u, nrow = nrow(mresponse), ncol = length(u),
                       byrow = TRUE)
        rowSums(mresponse <= utmp) == ncol(response)
    })
    dresponse <- factor(as.vector(dresponse))

    ### lhs may have multiple terms
    cfm <- paste(deparse(formula), collapse = "")
    cfm <- strsplit(cfm, "~")[[1]]
    xfm <- strsplit(cfm[2], "\\+")[[1]]
    yfm <- strsplit(cfm[1], "\\+")[[1]]
    tmp <- outer(xfm, yfm, function(x, y) paste(x, y, sep = "%O%"))
    xpart <- paste(as.vector(tmp), collapse = " + ")

    # yfm <- paste("%O%", as.character(formula)[2], sep = "")
    # fm <- paste("dresponse ~ ", paste(xfm, yfm, collapse = "+"))
    fm <- paste("dresponse ~ ", xpart)
    ### terms that depend on x only but not on y
    if (!is.null(constant)) {
        constant <- strsplit(constant, "\\+")[[1]]
        fm <- paste(fm, paste(constant, " %O% bols(ONEy, intercept = FALSE, df = 1)",
                        collapse = " + "),
                    sep = "+")
    }
    ### terms that depend on y only but not on x
    if (!is.null(varying)) {
        varying <- strsplit(varying, "\\+")[[1]]
        fm <- paste(fm, paste("bols(ONEx, intercept = FALSE, df = 1) %O%", varying,
                        collapse = " + "),
                    sep = "+")
    }
    fm <- formula(fm, env = new.env())
    for (y in yname) assign(y, oresponse[[y]], envir = environment(fm))
    ### ONEy is a constant on the lhs; same length as pseudo-response
    ### this is error prone
    assign("ONEy", rep(1.0, nrow(uresponse))) ###, environment(formula))
    assign("ONEx", rep(1.0, nrow(data))) ###, environment(formula))
    if (is.null(weights)) weights <- rep(1, nrow(data))
    w <- weights
    if (length(w) == nrow(data))
        w <- rep(w, nrow(uresponse))

    fitfct <- function(formula) {
        ret <- mboost(formula, data = data, weights = w, ...)
        class(ret) <- c("ctm", class(ret))
        ### reset weights for cvrisk etc., expanding works OK in bl_lin_matrix!
        ret$"(weights)" <- weights
        ret$ycdf <- yname
        ret$originalresponse <- response[,,drop = TRUE]
        ret$uresponse <- oresponse[,,drop = TRUE]
        ret$data <- data
        ret$call <- match.call()
        ret
    }
    if (fit) return(fitfct(fm))
    ret <- list(fitfct = fitfct, formula = fm)
    class(ret) <- "ctm_unfitted"
    ret
}

update.ctm_unfitted <- function(object, formula = object$formula) {
    environment(formula) <- environment(object$formula)
    return(object$fitfct(formula))
}

print.ctm_unfitted <- function(x, ...) {
    cat("\n")
    cat("Conditional Transformation Model:\n")
    print(x$formula)
    cat("\n")
    return(invisible(x))
}
