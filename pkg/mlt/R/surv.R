
.Surv2R <- function(object) {

    type <- attr(object, "type")
    stopifnot(type %in% c("left", "right", "interval", 
                          "interval2", "counting"))
    status <- object[, "status"]

    ret <- switch(type,
        "right" = R(y = ifelse(status == 1, object[, "time"], NA),
                    cleft = ifelse(status != 1, object[, "time"], NA),
                    cright = ifelse(status != 1, Inf, NA)),
        "left" =  R(y = ifelse(status == 1, object[, "time"], NA),
                    cleft = ifelse(status != 1, -Inf, NA),
                    cright = ifelse(status != 1, object[, "time"], NA)),

        "interval2" = {
            ret <- cbind(left = object[, "time1"], 
                         right = object[, "time2"])
            ret$left[is.na(ret$left)] <- -Inf
            ret$right[is.na(ret$right)] <- Inf
            R(cleft = ret$left, cright = ret$right)
        },
        "interval" = {
            status <- factor(status, levels = 0:3, 
                             labels = c("right", "exact", "left", "interval"))
            tmp <- matrix(NA, nrow = nrow(object), ncol = 2)
            colnames(tmp) <- c("left", "right")
            for (s in levels(status)) {
                idx <- which(status == s)
                tmp[idx, ] <- switch(s, 
                    "right" = cbind(object[idx, "time1"], Inf),
                    "exact" = cbind(object[idx, "time1"], NA),
                    "left" = cbind(-Inf, object[idx, "time1"]),
                    "interval" = object[idx, c("time1", "time2")])
            }
            R(y = ifelse(is.na(tmp[, "right"]), tmp[, "left"], NA), 
              cleft = ifelse(is.na(tmp[, "right"]), NA, tmp[, "left"]), 
              cright = tmp[, "right"])
        },
        ### left truncation, right censoring
        "counting" = R(y = ifelse(status == 1, object[, "stop"], NA),
                       cleft = ifelse(status != 1, object[, "stop"], NA),
                       cright = ifelse(status != 1, Inf, NA),
                       tleft = object[, "start"])
    )
    ret
}

### handle exact integer / factor as interval censored
R <- function(y = NA, cleft = NA, cright = NA, 
              tleft = NA, tright = NA) {

    if (inherits(y, "Surv"))
        return(.Surv2R(y))

    exact <- y
    n <- c(length(exact), length(cleft), length(cright), 
           length(tleft), length(tright))    
    if (length(unique(n)) > 1) {
        stopifnot(any(n > 1))
        n1 <- n[n > 1]
        stopifnot(length(unique(n1)) == 1)
    }
    nm <- max(n)
    if (length(exact) == 1 && is.na(exact)) 
        exact <- rep(exact, nm)
    if (length(cleft) == 1 && is.na(cleft)) 
        cleft <- rep(cleft, nm)
    if (length(cright) == 1 && is.na(cright)) 
        cright <- rep(cright, nm)
    if (length(tleft) == 1 && is.na(tleft)) 
        tleft <- rep(tleft, nm)
    if (length(tright) == 1 && is.na(tright)) 
        tright <- rep(tright, nm)

    if (all(!is.na(exact))) {
        rank <- rank(exact, ties.method = "max")
        type <- "double"
        if (is.factor(exact)) {
            type <- class(exact)[1]
            if (!is.ordered(exact)) 
                warning("response is unordered factor; 
                         results may depend on order of levels")
            stopifnot(all(is.na(c(tleft, cleft, cright, tright))))
            lev <- levels(exact)
            cright <- as.ordered(exact)
            cright[cright == lev[nlevels(exact)]] <- NA
            cleft <- factor(unclass(exact) - 1, levels = 1:length(lev), 
                            labels = lev, exclude = 0, ordered = TRUE)
            exact <- NA
        }
        if (is.integer(exact)) {
            type <- "integer"
            stopifnot(all(is.na(c(tleft, cleft, cright, tright))))
            cright <- exact
            cleft <- exact - 1
            cleft[cleft < 0] <- NA
            exact <- NA
        }
    } else {
        if (any(!is.na(exact)) && !is.double(exact))
            stop("Mix of exact and censored observations currently only implemented for doubles")
        ### <FIXME> allow for censord factors and integers </FIXME>
        type <- "double"

        ### some meaningful ordering of observations
        tmpexact <- as.numeric(exact)
        tmpleft <- as.numeric(cleft)
        tmpright <- as.numeric(cright)
        tmpleft[is.na(tmpleft)] <- 
            pmin(tmpleft, tmpexact, tmpright, na.rm = TRUE)[is.na(tmpleft)]
        tmpright[is.na(tmpright)] <- 
            pmax(tmpleft, tmpexact, tmpright, na.rm = TRUE)[is.na(tmpright)]
        tmpexact[is.na(tmpexact)] <- 
            ((tmpright - tmpleft) / 2)[is.na(tmpexact)]
        rank <- rank(tmpexact, ties.method = "max")
    }

    ret <- data.frame(tleft = tleft, cleft = cleft, exact = exact,
                      cright = cright, tright = tright, rank = rank)
    attr(ret, "type") <- type
    class(ret) <- c("response", class(ret))
    ret
}

.exact <- function(object)
    !is.na(object$exact)

.cleft <- function(object)
    is.finite(object$cleft) 

.cright <- function(object)
    is.finite(object$cright)

.cinterval <- function(object)
    .cleft(object) | .cright(object)

.tleft <- function(object)
    is.finite(object$tleft) 

.tright <- function(object)
   is.finite(object$tright)

.tinterval <- function(object)
    .tleft(object) | .tright(object)

.mm_exact <- function(model, data, response, object) {

    e <- .exact(object)
    if (!any(e)) return(NULL)
    tmp <- data[e,,drop = FALSE]
    tmp[[response]] <- object$exact[e]
    Y <- model.matrix(model, data = tmp)
    deriv <- 1
    names(deriv) <- response
    Yprime <- model.matrix(model, data = tmp, deriv = deriv)

    trunc <- NULL
    if (any(.tinterval(object) & e)) {
        Ytleft <- matrix(-Inf, nrow = nrow(Y), ncol = ncol(Y))
        Ytright <- matrix(Inf, nrow = nrow(Y), ncol = ncol(Y))
        if (any(il <- (.tleft(object) & e))) {
            tmp <- data[il,]
            tmp[[response]] <- object$tleft[il]
            Ytleft[il,] <- model.matrix(model, data = tmp)
        }
        if (any(ir <- (.tright(object) & e))) {
            tmp <- data[ir,,drop = FALSE]
            tmp[[response]] <- object$tright[ir]
            Ytright[ir,] <- model.matrix(model, data = tmp)
        }
        trunc <- list(left = Ytleft, right = Ytright)
    }

    list(Y = Y, Yprime = Yprime, trunc = trunc, which = which(e))
}

.mm_interval <- function(model, data, response, object) {

    i <- .cinterval(object)
    if (!any(i)) return(NULL)
    tmpdata <- data[i,,drop = FALSE]
    object <- object[i,, drop = FALSE]

    Yleft <- NULL
    if (any(il <- .cleft(object))) {
        tmp <- tmpdata[il,,drop = FALSE]
        tmp[[response]] <- object$cleft[il]
        Ytmp <- model.matrix(model, data = tmp)
        Yleft <- matrix(-Inf, nrow = length(il), ncol = ncol(Ytmp))
        colnames(Yleft) <- colnames(Ytmp)
        attr(Yleft, "constraint") <- attr(Ytmp, "constraint")
        attr(Yleft, "Assign") <- attr(Ytmp, "Assign")
        Yleft[il,] <- Ytmp
    }

    Yright <- NULL
    if (any(ir <- .cright(object))) {
        tmp <- tmpdata[ir,, drop = FALSE]
        tmp[[response]] <- object$cright[ir]
        Ytmp <- model.matrix(model, data = tmp)
        Yright <- matrix(Inf, nrow = length(ir), ncol = ncol(Ytmp))
        colnames(Yright) <- colnames(Ytmp)
        attr(Yright, "constraint") <- attr(Ytmp, "constraint")
        attr(Yright, "Assign") <- attr(Ytmp, "Assign")
        Yright[ir,] <- Ytmp
    }

    if (is.null(Yright)) { 
        Yright <- matrix(Inf, nrow = nrow(Yleft), ncol = ncol(Yleft))
        colnames(Yright) <- colnames(Yleft)
        attr(Yright, "constraint") <- attr(Yleft, "constraint")
        attr(Yright, "Assign") <- attr(Yleft, "Assign")
    }
    if (is.null(Yleft)) {
        Yleft <- matrix(-Inf, nrow = nrow(Yright), ncol = ncol(Yright))
        colnames(Yleft) <- colnames(Yright)
        attr(Yleft, "constraint") <- attr(Yright, "constraint")
        attr(Yleft, "Assign") <- attr(Yright, "Assign")
    }

    trunc <- NULL
    if (any(.tinterval(object))) {
        Ytleft <- matrix(-Inf, nrow = nrow(Yleft), ncol = ncol(Yleft))
        Ytright <- matrix(Inf, nrow = nrow(Yleft), ncol = ncol(Yleft))
        colnames(Ytleft) <- colnames(Ytright) <- colnames(Yleft)
        if (any(il <- (.tleft(object)))) {
            tmp <- tmpdata[il,,drop = FALSE]
            tmp[[response]] <- object$cleft[il]
            Ytleft[il,] <- model.matrix(model, data = tmp)
        }
        if (any(ir <- (.tright(object)))) {
            tmp <- tmpdata[ir,,drop = FALSE]
            tmp[[response]] <- object$cleft[ir]
            Ytright[ir,] <- model.matrix(model, data = tmp)
        }
        trunc <- list(left = Ytleft, right = Ytright)
    }

    list(Yleft = Yleft, Yright = Yright, trunc = trunc, which = which(i))
}
