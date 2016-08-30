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
        pr <- predict.trafotree(obj, q = q, type = type, ...)
        yscale <- range(pr)
    }
    xscale <- range(q)

    ### panel function for ecdf in nodes
    rval <- function(node) {

        nid <- id_node(node)
        dat <- data_party(obj, nid)
        wn <- dat[["(weights)"]]   

        cf <- obj$coef[as.character(nid),]
        coef(mod) <- cf
        y <- as.vector(predict(mod, newdata = data.frame(1), q = q, type = type))

        ## set up plot
        q <- q - xscale[1]
        q <- q / diff(xscale)
        y <- y - yscale[1]
        y <- y / diff(yscale)

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

plot.trafotree <- function(x, type = "trafo", ...) {
    class(x) <- class(x)[-1L]
    tp <- function(...) node_mlt(..., type = type)
    class(tp) <- class(node_mlt)
    plot(x, terminal_panel = tp, ...)
}

