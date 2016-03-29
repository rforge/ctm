
.trparty <- function(object, part, data, parm, weights, modelsplit, 
                     control, FUN, ...) {

    if (missing(parm)) parm <- 1:length(coef(object))

    cfleft <- cfright <- coef(object)
    wleft <- wright <- 0
    trfo <- function(data, weights) {
        thetastart <- coef(object)
        if (max(abs(weights - wleft)) < .Machine$double.eps)
           thetastart <- cfleft
        if (max(abs(weights - wright)) < .Machine$double.eps)
           thetastart <- cfright
        estfun(update(object, weights = weights,
                      theta = thetastart))[, parm, drop = FALSE]
    }

    split <- function(x, response, weights, mb) {


        if (is.ordered(x)) x <- unclass(x)
        if (is.factor(x)) {
            psplits <- partykit:::mob_grow_getlevels(x)
            lmod <- rmod <- object
            ll <- numeric(NROW(psplits))
            lmax <- -Inf
            for (i in 1:NROW(psplits)) {
                ll[i] <- logLik(lmod <- update(object, weights = weights * (x %in% lev[psplits[i,]]),
                                               theta = coef(lmod))) +
                         logLik(rmod <- update(object, weights = weights * (x %in% lev[!psplits[i,]]),
                                               theta = coef(rmod)))
                if (lll[i] > lmax) {
                    lmax <- ll[i]  
                    cfleft <<- coef(lmod)
                    cfright <<- coef(rmod)
                    xsplit <- i
                }
            }
            return((!psplit[i,]) + 1)
        }
        ox <- order(x)
        w <- cumsum(weights[ox])
        ux <- x[ox][(w > mb) & (w < (sum(weights) - mb))]
        sux <- sort(unique(ux))
        suxl <- rev(sux[sux <= median(sux)])
        suxr <- sux[sux > median(sux)]

        lll <- numeric(length(suxl))
        lmax <- -Inf
        xsplit <- NA
        lmod <- rmod <- object
        for (i in 1:length(suxl)) {
            lll[i] <- logLik(lmod <- update(object, weights = weights * (x <= suxl[i]),
                                        theta = coef(lmod))) +
                      logLik(rmod <- update(object, weights = weights * (x > suxl[i]),
                                        theta = coef(rmod)))
            if (lll[i] > lmax) {
                lmax <- lll[i]
                cfleft <<- coef(lmod)
                cfright <<- coef(rmod)
                xsplit <-  suxl[i]
            }
        }
        llr <- numeric(length(suxr))
        lmod <- rmod <- object
        for (i in 1:length(suxr)) {
            llr[i] <- logLik(lmod <- update(object, weights = weights * (x <= suxr[i]),
                                            theta = coef(lmod))) +
                      logLik(rmod <- update(object, weights = weights * (x > suxr[i]),
                                            theta = coef(rmod)))
            if (llr[i] > lmax) {
                lmax <- llr[i]
                cfleft <<- coef(lmod)
                cfright <<- coef(rmod)
                xsplit <-  suxr[i]
            }
        }

        wleft <<- weights * (x <= xsplit)
        wright <<- weights * (x > xsplit)

        xsplit
    }


    mf <- model.frame(part, data = data)
    stopifnot(all(!names(mf) %in% variable.names(object)))
    if (missing(weights)) weights <- rep(1, NROW(mf))
    tmp <- trfo(data, weights)

    mf$response <- 0
    fm <- response ~ x
    fm[[3L]] <- part[[2L]]

    if (modelsplit) control$splitfun <- split
    ct <- FUN(fm, data = mf, ytrafo = trfo, weights = weights, 
              control = control, ...)

    ct$model <- object
    ct
}

trtree <- function(object, part, data, parm, weights, modelsplit = FALSE, 
                   control = ctree_control(), ...) {

    ret <- .trparty(object = object, part = part, data = data, parm = parm,
                    weights = weights, modelsplit = modelsplit,
                    control = control, FUN = ctree, ...)

    nf <- predict(ret, type = "node")
    if (missing(weights)) weights <- rep(1, length(nf))
    cf <- tapply(1:length(nf), nf, function(nd) {
        coef(update(object, 
                    weights = weights * (1:length(nf) %in% nd), theta = coef(object))) 
    })

    ret$coefficients <- cf
    class(ret) <- c("trtree", class(ret))
    ret
}

predict.trtree <- function(object, newdata, K = 20, q = NULL,
    type = c("node", "coef", "trafo", "distribution", "survivor", "density", 
             "logdensity", "hazard", "loghazard", "cumhazard", "quantile"), 
    ...) {

    type <- match.arg(type)
    if (missing(newdata)) {
        nf <- fitted(object)[["(fitted)"]]
    } else {
        tmp <- object
        class(tmp) <- class(tmp)[-1L]
        nf <- predict(tmp, newdata = newdata, type = "node")
    }
    if (type == "node") return(nf)
    nf <- factor(nf)
    if (type == "coef") return(do.call("cbind", coef(object)[nf]))
    
    mod <- object$model
    if (is.null(q))
        q <- mkgrid(mod, n = K)[[mod$response]]

    if (missing(newdata)) newdata <- mod$data

    ### <FIXME> need .R2vec??? </FIXME>
    pr <- .R2vec(predict(mod, newdata = newdata, q = q, type = type, ...))
    if (!is.matrix(pr))
        pr <- matrix(pr, nrow = NROW(pr), ncol = NROW(newdata))
    for (nd in levels(nf)) {
        i <- nf == nd
        coef(mod) <- object$coef[[nd]]
        pr[,i] <- .R2vec(predict(mod, newdata = newdata[i,], q = q, 
                                 type = type, ...))
    }
    pr
}

node_mlt <- function(obj, col = "black", bg = "white", ylines = 2,
                     id = TRUE, mainlab = NULL, gp = gpar(), K = 20, 
                     type = c("trafo", "distribution", "survivor", 
                              "density", "logdensity", "hazard", 
                              "loghazard", "cumhazard", "quantile"), ...)
{
    mod <- obj$model
    q <- mkgrid(mod, n = K)[[mod$response]]
    type <- match.arg(type)

    if (type %in% c("distribution", "survivor")) {
        yscale <- c(0, 1)
    } else {
        pr <- predict.trtree(obj, q = q, type = type, ...)
        yscale <- range(pr)
    }
    xscale <- range(q)

    ### panel function for ecdf in nodes
    rval <- function(node) {

        nid <- id_node(node)
        dat <- data_party(obj, nid)
        wn <- dat[["(weights)"]]

        cf <- obj$coefficients[[as.character(nid)]]
        coef(mod) <- cf
        y <- as.vector(predict(mod, q = q, type = type, ))

        ## set up plot
        q <- q - xscale[1]
        q <- q / xscale[2]
        y <- y - yscale[1]
        y <- y / yscale[2]

        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines, 1, 1),
                                         c("lines", "null", "lines")),
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"),
                           height = unit(1, "npc") - unit(2, "lines"),
                           name = paste("node_mlt", nid, sep = ""), gp = gp)

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = bg, col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
        if (is.null(mainlab)) {
          mainlab <- if(id) {
            function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
          } else {
            function(id, nobs) sprintf("n = %s", nobs)
          }
        }
        if (is.function(mainlab)) {
          mainlab <- mainlab(nid, sum(wn))
        }
        grid.text(mainlab)
        popViewport()

        plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                         xscale=xscale, yscale=yscale,
                         name = paste0("node_mlt", nid, "plot"),
                         clip = FALSE)

        pushViewport(plot)
        grid.xaxis()
        grid.yaxis()
        grid.rect(gp = gpar(fill = "transparent"))
        grid.clip()
        grid.lines(q, y, gp = gpar(col = col))
        upViewport(2)
    }

    return(rval)
}
class(node_mlt) <- "grapcon_generator"

plot.trtree <- function(x, ...) {
    class(x) <- class(x)[-1L]
    plot(x, terminal_panel = node_mlt, ...)
}

coef.trtree <- function(object, ...)
    object$coef

logLik.trtree <- function(object, newdata, ...) {

    if (missing(newdata)) newdata <- object$model$data

    nd <- factor(predict(object, newdata = newdata, type = "node", ...))
    mod <- mlt(object$model$model, data = newdata, dofit = FALSE)

    ll <- rep(0, nlevels(nd))
    for (i in 1:nlevels(nd)) {
        w <- rep(0, NROW(newdata))
        w[nd == levels(nd)[i]] <- 1
        ll[i] <- logLik(mod, parm = coef(object)[[levels(nd)[i]]], w = w)
    }
    ret <- sum(ll)
    attr(ret, "df") <- NA
    class(ret) <- "logLik"
    ret
}

simulate.trtree <- simulate.trforest
