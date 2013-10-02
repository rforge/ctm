
ctm <- function(formula, data, weights = NULL, constant = NULL, 
                varying = NULL, asis = NULL, monotone = FALSE, 
                ngrid = NULL, ...) {

    opt <- options(mboost_Xmonotone = monotone)
    yname <- all.vars(formula)[1]
    response <- data[[yname]]
    data[[yname]] <- NULL

    uresponse <- sort(unique(response))
    ### make sure that for each observation there is at least
    ### one pseudo response = 0 AND = 1
    if (is.numeric(uresponse)) {
        uresponse <- c(uresponse[1] - (uresponse[2] - uresponse[1]), 
                       uresponse)
        if (!is.null(ngrid))
            uresponse <- seq(from = min(uresponse), to = max(uresponse),
                             length = ngrid)
    } else {
        if (!is.ordered(response)) stop(sQuote("response"), " is neither ordered nor numeric")
        uresponse <- uresponse[-length(uresponse)]
    }
    dresponse <- factor(sapply(uresponse, function(r) response <= r))
    if (is.ordered(uresponse)) uresponse <- uresponse[,drop = TRUE]

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
    ### you know what you do, right?
    if (!is.null(asis)) 
        fm <- paste(fm, asis, sep = "+")

    fm <- as.formula(fm)
    assign(yname, uresponse)###, environment(formula))
    ### ONEy is a constant on the lhs; same length as pseudo-response
    ### this is error prone
    assign("ONEy", rep(1.0, length(uresponse))) ###, environment(formula))
    assign("ONEx", rep(1.0, nrow(data))) ###, environment(formula))
    if (is.null(weights)) weights <- rep(1, nrow(data))
    w <- weights
    if (length(w) == nrow(data))
        w <- rep(w, length(uresponse))
    ret <- mboost(fm, data = data, weights = w, offset = 0, ...)
    class(ret) <- c("ctm", class(ret))
    ### reset weights for cvrisk etc., expanding works OK in bl_lin_matrix!
    ret$"(weights)" <- weights
    ret$ycdf <- all.vars(formula)[1]
    ret$originalresponse <- response
    ret$uresponse <- uresponse
    ret$data <- data
    ret$call <- match.call()
    options(opt)
    ret
}

tune <- function(object, alpha = 0.05, mstopmax = 1000, 
                 mstep = 10, trace = TRUE) {

    mf <- object$data
    mf$y <- object$originalresponse
    names(mf)[names(mf) == "y"] <- object$ycdf
    fm <- object$call$family
    link <- link2dist(fm$link)
    while(mstop(object) < mstopmax) {
        z <- predict(object, newdata = mf, y = NULL)
        ### <FIXME> link$p is not available for "logit" </FIXME>
        p <- ks.test(z, link$p)$p.value
        if (trace) 
            cat("mstop: ", mstop(object), "; p = ", p, "\n")
        if (p > alpha) break()
        object[mstop(object) + mstep]
    }
    mstop(object)
}

predict.ctm <- function(object, newdata, y = object$uresponse, 
                           annotated = FALSE, ...) {

    class(object) <- class(object)[-1]
    if (is.null(y)) {
        n <- nrow(newdata)
        p <- predict(object, newdata = newdata, ...)
        indx <- (1:n * n) + 1:n - n
        if (!is.matrix(p))
            return(p[indx])
        return(p[indx, , drop = FALSE])
    }
    x <- newdata[, !(colnames(newdata) %in% object$ycdf), drop = FALSE]
    nd <- as.list(x)
    if (NCOL(y) == 1) {
        y <- data.frame(y)
        colnames(y) <- object$ycdf
    }
    for (yn in object$ycdf)
        nd[[yn]] <- y[[yn]]
    nd$ONEy <- rep(1, length(y))
    nd$ONEx <- rep(1, nrow(x))
    p <- predict(object, newdata = nd, ...)
    if (!annotated) return(p)
    indx <- rep(1:NROW(y), rep(nrow(newdata), NROW(y)))
    ytmp <- y[indx,,drop = FALSE]
    ret <- cbind(ytmp, data.frame(
                      ID = rep(1:nrow(newdata), NROW(y)),
                      p = as.vector(p)))
    ret
}

outrisk <- function(object, newdata) {

    nresponse <- newdata[[object$ycdf]]
    y <- object$uresponse
    ndresponse <- factor(sapply(y, function(r) nresponse <= r))
    f <- predict(object, newdata = newdata, y = y, aggregate = "cumsum")
    y <- object$family@check_y(ndresponse)
    apply(f, 2, function(x) object$family@risk(y, x))
}

qest <- function(x, prob = c(.1, .5, .9)) {
    y <- x[[1]]
    p <- x$p
    f <- approxfun(y, p)
    zf <- function(q, prob) (f(q) - prob)^2
    sapply(prob, function(p) optimize(zf, interval = range(y), prob = p)$minimum)
}

plot.ctm <- function(x, which = sort(unique(selected(object))), 
                        prob = TRUE, obs = TRUE, ...) {

    object <- x
    pfun <- function(which) {
        stopifnot(length(which) == 1)
        df <- model.frame(object, which = which)[[1]]
        df <- lapply(df, function(x) seq(from = min(x), 
                                         to = max(x), length = 50))
        df2 <- do.call("expand.grid", df)
        class(object) <- class(object)[-1]
        df2$pr <- predict(object, newdata = df, which = which)
        ret <- cbind(df2, which = names(df)[1])
        names(ret)[1:2] <- c("x", "y")
        ret
    }

    ret <- do.call("rbind", lapply(which, pfun))
    if (nlevels(ret$which) == 1 & prob) {
        ret$pr <- object$family@response(ret$pr)
        ret$pr <- ifelse(ret$pr < .5, ret$pr, 1 - ret$pr)
        at <- seq(from = .05, to = .5, by = .05)
    } else {
        ret$pr <- ifelse(ret$pr < 0, ret$pr, - ret$pr)
        at <- seq(from = 0, to = min(ret$pr), length = 10)
    }
    pf <- function(x, y, z, subscripts, ...) {
        xname <- as.character(unique(ret[subscripts, "which"]))
        xo <- model.frame(object, which = xname)[[1]][[xname]]
        yo <- object$originalresponse
        panel.levelplot.raster(x, y, z, subscripts, ...)
        panel.contourplot(x, y, z, subscripts, contour = TRUE, ...)
        if (obs)
            panel.xyplot(x = xo, y = yo, pch = 20)
    }
    levelplot(pr ~ x + y | which, data = ret, panel = pf,
              scales = list(x = list(relation = "free")), region = TRUE, 
              at = at, col.regions = grey.colors(length(at)), ...)
}

### not exported from mboost
link2dist <- function(link, choices = c("logit", "probit"), ...) {
    i <- pmatch(link, choices, nomatch = 0L, duplicates.ok = TRUE)
    if (i[1] == 1) return("logit")
    if (i[1] == 2) {
        ret <- list(p = pnorm, d = dnorm, q = qnorm)
        attr(ret, "link") <- link
        return(ret)
    }
    p <- get(paste("p", link, sep = ""))
    d <- get(paste("d", link, sep = ""))
    q <- get(paste("q", link, sep = ""))
    ret <- list(p = function(x) p(x, ...),
                d = function(x) d(x, ...),
                q = function(x) q(x, ...))
    attr(ret, "link") <- link
    ret
}
    
### rho(y, f) = (y - link(f))^2
ClassL2 <- function(link = c("norm", "logis", "t"), ...) {
    link <- link2dist(link, ...)
    biny <- function(y) {
        if (!is.factor(y))
            stop("response is not a factor but ",
                  sQuote("family = Class()"))
            if (nlevels(y) != 2)
                stop("response is not a factor at two levels but ",
                      sQuote("family = ClassL2()"))
        return(c(-1, 1)[as.integer(y)])
    }

    trf <- function(f) {
        thresh <- -link$q(.Machine$double.eps)
        pmin(pmax(f, -thresh), thresh)
    }

    return(Family(ngradient = function(y, f, w = 1) {
               y <- (y + 1) / 2
               p <- link$p(trf(f))
               (y - p) * link$d(trf(f))
           },
           loss = function(y, f) {
               p <- link$p(trf(f))
               y <- (y + 1) / 2
               (y - p)^2
           },
           offset = function(y, w) {
               p <- weighted.mean(y > 0, w)
               link$q(p)
           },
           response = function(f) {
               p <- link$p(trf(f))
               return(p)
           },
           rclass = function(f) (f > 0) + 1 ,
           check_y = biny,
           name = paste("Classification L2 loss --", 
                        attr(link, "link"), "Link")))
}

### rho(y, f) = |y - link(f)|
ClassL1 <- function(link = c("norm", "logis", "t"), ...) {
    link <- link2dist(link, ...)
    biny <- function(y) {
        if (!is.factor(y))
            stop("response is not a factor but ",
                  sQuote("family = Class()"))
            if (nlevels(y) != 2)
                stop("response is not a factor at two levels but ",
                      sQuote("family = ClassL1()"))
        return(c(-1, 1)[as.integer(y)])
    }

    trf <- function(f) {
        thresh <- -link$q(.Machine$double.eps)
        pmin(pmax(f, -thresh), thresh)
    }

    return(Family(ngradient = function(y, f, w = 1) {
               y <- (y + 1) / 2
               p <- link$p(trf(f))
               sign(y - p) * link$d(trf(f))
           },
           loss = function(y, f) {
               p <- link$p(trf(f))
               y <- (y + 1) / 2
               abs(y - p)
           },
           offset = function(y, w) {
               p <- weighted.mean(y > 0, w)
               link$q(p)
           },
           response = function(f) {
               p <- link$p(trf(f))
               return(p)
           },
           rclass = function(f) (f > 0) + 1 ,
           check_y = biny,
           name = paste("Classification L1 loss --", 
                        attr(link, "link"), "Link")))
}

