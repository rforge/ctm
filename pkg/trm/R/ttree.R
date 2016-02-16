
ttree <- function(object, part, data, parm, weights, ...) {

    if (missing(parm)) parm <- 1:length(coef(object))

    trfo <- function(data, weights)
        estfun(update(object, weights = weights,
                      theta = coef(object)))[, parm, drop = FALSE]

    mf <- model.frame(part, data = data)
    stopifnot(all(!names(mf) %in% variable.names(object)))
    if (missing(weights)) weights <- rep(1, NROW(mf))
    tmp <- trfo(data, weights)

    mf$response <- 0
    fm <- response ~ x
    fm[[3L]] <- part[[2L]]

    ct <- ctree(fm, data = mf, ytrafo = trfo, weights = weights, 
                ...)

    nf <- predict(ct, type = "node")

    cf <- tapply(1:nrow(mf), nf, function(nd) {
        coef(update(object, 
                    weights = weights * (1:nrow(mf) %in% nd))) 
    })

    ct$coefficients <- cf
    ct$model <- object
    class(ct) <- c("ttree", class(ct))
    ct
}

predict.ttree <- function(object, newdata, K = 20, 
    type = c("node", "trafo", "distribution", 
             "survivor", "density", "logdensity", 
             "hazard", "loghazard", 
             "cumhazard", "quantile"), ...) {

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
    
    mod <- object$model
    q <- mkgrid(mod, n = K)[[mod$response]]

    if (missing(newdata)) newdata <- mod$data

    pr <- predict(mod, newdata = newdata, q = q, type = type, ...)
    if (!is.matrix(pr))
        pr <- matrix(pr, nrow = length(q), ncol = NROW(newdata))
    for (nd in levels(nf)) {
        i <- nf == nd
        coef(mod) <- object$coef[[nd]]
        pr[,i] <- predict(mod, newdata = newdata[i,], q = q, 
                          type = type, ...)
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
        pr <- predict(obj, q = q, type = type, ...)
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
        q <- q - min(q)
        q <- q / max(q)
        y <- y - min(y)
        y <- y / max(y)

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

plot.ttree <- function(x, ...)
    partykit:::plot.constparty(x, terminal_panel = node_mlt, ...)
