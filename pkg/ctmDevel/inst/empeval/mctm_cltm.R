mctm_cltm <- function(formula, data, weights = NULL, constant = NULL, varying = NULL, 
                      ngrid = NULL, no.variance = NULL, yname, 
                      fit = TRUE, offset = 0, ...){
    #bl <- eval(formula[[2]], envir = data)
    #yname <- all.vars(formula)[1]
    response <- data[yname]
    for (y in yname) data[[y]] <- NULL
    yu <- lapply(response, upsilon, ngrid = ngrid)
    mresponse <- as.matrix(sapply(yu, function(x) x$y))
    uresponse <- do.call("expand.grid", lapply(yu, function(x) x$upsilon))
    oresponse <- do.call("expand.grid", lapply(yu, function(x) x$uorig))
    if (all(sapply(response, is.factor))) {
        uresponse <- uresponse[-nrow(uresponse), , drop = FALSE]
        oresponse <- oresponse[-nrow(oresponse), , drop = FALSE]
        for (y in yname) oresponse[[y]] <- oresponse[[y]][, drop = TRUE]
    }
    dresponse <- apply(uresponse, 1, function(u) {
        utmp <- matrix(u, nrow = nrow(mresponse), ncol = length(u), 
            byrow = TRUE)
        rowSums(mresponse <= utmp) == ncol(response)
    })
    dresponse <- factor(as.vector(dresponse))

    if(no.variance) fm <- "dresponse ~ "
    else{
    cfm <- paste(deparse(formula), collapse = "")
    cfm <- strsplit(cfm, "~")[[1]]
    xfm <- strsplit(cfm[2], "\\+")[[1]]
    yfm <- strsplit(cfm[1], "\\+")[[1]]
    tmp <- outer(xfm, yfm, function(x, y) paste(x, y, sep = "%O%"))
    xpart <- paste(as.vector(tmp), collapse = " + ")
    fm <- paste("dresponse ~ ", xpart)}

    if (!is.null(constant)) {
        constant <- strsplit(constant, "\\+")[[1]]
        if(no.variance)
        fm <- paste(fm, paste(constant, " %O% bols(ONEy, intercept = FALSE, df = 1)",
                    collapse = " + "))
        else{
        fm <- paste(fm, paste(constant, " %O% bols(ONEy, intercept = FALSE, df = 1)",
                        collapse = " + "),
                    sep = "+")}
    }
    if (!is.null(varying)) {
        varying <- strsplit(varying, "\\+")[[1]]
        fm <- paste(fm, paste("bols(ONEx, intercept = FALSE, df = 1) %O%", 
            varying, collapse = " + "), sep = "+")
    }
    fm <- formula(fm, env = new.env())
    for (y in yname) assign(y, oresponse[[y]], envir = environment(fm))
    assign("ONEy", rep(1, nrow(uresponse)))
    assign("ONEx", rep(1, nrow(data)))
    if (is.null(weights)) 
        weights <- rep(1, nrow(data))
    w <- weights
    if (length(w) == nrow(data)) 
        w <- rep(w, nrow(uresponse))
    fitfct <- function(formula) {
        ret <- mboost(formula, data = data, weights = w, offset = offset, 
            ...)
        class(ret) <- c("ctm", class(ret))
        ret$"(weights)" <- weights
        ret$ycdf <- yname
        ret$originalresponse <- response[, , drop = TRUE]
        ret$uresponse <- oresponse[, , drop = TRUE]
        ret$data <- data
        ret$call <- match.call()
        ret
    }
    if (fit) 
        return(fitfct(fm))
    ret <- list(fitfct = fitfct, formula = fm)
    class(ret) <- "ctm_unfitted"
    ret
}

