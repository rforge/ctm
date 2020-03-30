.add_confband <-function(object, smooth = FALSE, col = "black", lwd = 1, 
                         fill = "lightgrey", border = NA) {
    q <- object[, "q"]
    lwr <- object[, "lwr"]
    upr <-  object[, "upr"]
    est <- object[, "Estimate"]
    polyx <- c(q[1], rep(q[2:length(q)], each = 2), q[length(q)],
               rep(rev(q)[2:length(rev(q))], each = 2))
    polyy <-  c(rep(lwr[1:(length(lwr)-1)], each = 2), rev(upr)[1], 
                rep(rev(upr)[2:length(rev(upr))], each = 2), lwr[1])
    type <- "s"
    if (smooth) {
        polyx <- c(q, rev(q))
        polyy <- c(lwr, rev(upr))
        type <- "l"
    }
    polygon(x = polyx, y = polyy, col = fill, border = border)
    lines(x = q, y = est, type = type, col = col, lwd = lwd)
    if (!smooth) points(x = q, y = est, col = col, pch = 20, cex =.75)
}


plot.cotram <- function(x, newdata,
                        type = c("distribution", "survivor", "density", "logdensity",
                                 "cumhazard", "quantile", "trafo"),
                        confidence = c("none", "band"), level = 0.95, 
                        smooth = FALSE, q = NULL, K = 20, cheat = K, prob = 1:(10-1)/10,
                        col = "black", fill = "lightgrey", lty = 1, lwd = 1, add = FALSE, ...) {
    
    args <- list(...)
    y <- variable.names(x, "response")
    
    if (is.null(q)){
        q <- mkgrid(x, n = K)[[y]] - as.integer(x$plus_one)
        if (smooth)
            q <- seq(from = min(q), to = max(q), length.out = K)
    }
    
    type <- match.arg(type)
    
    pr <- predict(x, newdata = newdata, type = type, q = q , smooth = smooth,
                  prob = prob, K = K, ...) 
    pr[!is.finite(pr)] <- NA
    
    cb <- NULL
    confidence <- match.arg(confidence)
    # calpha <- switch(confidence, "none" = NULL,
    #                  "interval" = univariate_calpha(),
    #                  "band" = adjusted_calpha())
    confidence <- confidence != "none"
    
    if (confidence)
        cb <- confband(x, newdata = newdata, level = level,
                       type = type, K = K, cheat = cheat, 
                       smooth = smooth)
    
    rpr <- range(pr, na.rm = TRUE)
    if (is.null(dim(pr))) pr <- matrix(pr, ncol = 1)
    ylim <- switch(type, "distribution" = c(0, 1),
                   "survivor" = c(0, 1),
                   "density" = c(0, rpr[2]),
                   "hazard" = c(0, rpr[2]),
                   "cumhazard" = c(0, rpr[2]),
                   rpr)
    if (!is.null(args$ylim)) ylim <- args$ylim
    if (type == "quantile")  {
        q <- prob
        if (is.null(args$xlab)) args$xlab <- type
        if (is.null(args$ylab)) args$ylab <- y
    }
    if (length(col) == 1) col <- rep(col, ncol(pr))
    if (length(lty) == 1) lty <- rep(lty, ncol(pr))
    if (length(lwd) == 1) lwd <- rep(lwd, ncol(pr))
    
    if (!add) {
        args$x <- unclass(q)
        args$y <- rep(rpr[1], length(q))
        args$ylim <- ylim
        args$xlab <- ifelse(is.null(args$xlab), y, args$xlab)
        args$ylab <- ifelse(is.null(args$ylab), type, args$ylab)
        args$type <- "n"
        args$axes <- FALSE
        do.call("plot", args)
        if (is.factor(q)) {
            axis(1, at = unclass(q), labels = levels(q))
        } else {
            axis(1)
        }
        axis(2)
        box()
    }
    if (smooth) {
        for (i in 1:ncol(pr)) lines(q, pr[,i], col = col[i], lty = lty[i], lwd = lwd[i])
    } else {
        ty <- ifelse(!length(grep(type, "density")) > 0, "s", "h")
        for (i in 1:ncol(pr)){
            lines(q, pr[,i], col = col[i], lty = lty[i], type = ty)
            points(q, pr[,i], col = col[i], pch = 20, cex =.75)
        }
    }
    invisible(pr)
    
    if (confidence) {
        if (length(fill) != NROW(newdata)) 
            fill <- rep(fill, length.out = NROW(newdata))
        if (length(col) != NROW(newdata)) 
            col <- rep(col, length.out = NROW(newdata))
        
        if (is.matrix(cb)) {
            .add_confband(cb, smooth = smooth, col = col[1], lwd = lwd[1], fill = fill[1]) 
        } else {
            out <- lapply(1:length(cb), function(i) 
                .add_confband(cb[[i]], smooth = smooth, col = col[i], lwd = lwd[i], fill = fill[i]))
        }
    }
}

