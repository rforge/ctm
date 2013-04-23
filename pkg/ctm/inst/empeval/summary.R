
load("ret.Rda")
library("lattice")

### temporary results
ret <- ret[!sapply(ret, is.null)]

model <- names(ret[[1]][[1]])
id <- 1:length(ret[[1]])
chr <- names(ret[[1]][[1]][[1]])

p <- unlist(ret)
x <- expand.grid(                        
                 chr = chr,
                 model = model, 
                 id = id,
                 pnon = names(ret))
x <- x[1:length(p),]
x$p <- p
levels(x$chr)[3] <- "Med."
levels(x$chr) <- paste(levels(x$chr), "MAD")
pnon <- as.integer(levels(x$pnon))
levels(x$pnon) <- paste("p = ", pnon)

x <- subset(x, model %in% c("boosting", "kernel", "gamlss"))
x <- subset(x, chr %in% c("Min. MAD", "Med. MAD", "Max. MAD"))
ylim <- tapply(1:nrow(x), x$chr[,drop = TRUE], function(i) range(x[i,"p"]))
ylim <- ylim[rep(1:length(ylim), rep(nlevels(x$pnon), length(ylim)))]

pfun <- function(x, y, subscripts, ...) {
    panel.abline(h = median(y[x == "boosting"]), col = "lightgrey")
    panel.bwplot(x, y, ...)
}

bwplot(p ~ model | pnon + chr, data = x, panel = pfun,
    ylab = "MAD",
    scales = list(y = list(relation = "free", limits = ylim, tck = 0.2),
                  x = list(rot = 45)), fill="grey",
    par.settings = list(plot.symbol = list(col=1,pch=20, cex=0.7),
                        box.rectangle = list(col=1),
                        box.umbrella = list(lty=1, col=1)),
    strip= strip.custom(bg="white"))

load("qret.Rda")

pfun <- function(x, y, subscripts, ...) {
    panel.xyplot(x = x, y = y, ...)
    panel.abline(a = 0, b = 1, col = "lightgrey")
}

qctm <- cbind(q[, c("mb", "mt", "tau")], model = "CTM")
qq <- cbind(q[, c("mq", "mt", "tau")], model = "AQR")
names(qctm) <- names(qq) <- c("estimate", "true", "tau", "model")
q <- rbind(qctm, qq)

xyplot(estimate ~true | tau * model, data = q, panel = pfun, layout = c(3, 2),
    xlab = "True quantile", ylab = "Estimated quantile",
    par.settings = list(plot.symbol = list(col=1,pch=20, cex=0.7),
                        box.rectangle = list(col=1),
                        box.umbrella = list(lty=1, col=1)),
    strip= strip.custom(bg="white"))
