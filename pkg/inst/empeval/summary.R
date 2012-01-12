### R code from vignette source 'empeval.Rnw'

###################################################
### code chunk number 1: sim-data
###################################################
source("setup.R")
load("ret.Rda")

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
x$model <- x$model[, drop = TRUE]
levels(x$model) <- c("CTM", "Kernel", "GAMLSS")
x <- subset(x, chr %in% c("Min. MAD", "Med. MAD", "Max. MAD"))
ylim <- tapply(1:nrow(x), x$chr[,drop = TRUE], function(i) range(x[i,"p"]))
ylim <- ylim[rep(1:length(ylim), rep(nlevels(x$pnon), length(ylim)))]


###################################################
### code chunk number 2: sim-plot
###################################################
pfun <- function(x, y, subscripts, ...) {
    panel.abline(h = median(y[x == "CTM"]), col = "lightgrey")
    panel.bwplot(x, y, ...)
}

print(bwplot(p ~ model | pnon + chr, data = x, panel = pfun,
    ylab = "MAD",
    scales = list(y = list(relation = "free", limits = ylim, tck = 0.2),
                  x = list(rot = 45)), fill="grey"))


