## ----setup, echo = FALSE, results = "hide", message = FALSE--------------
set.seed(290875)

dir <- system.file("applications", package = "tbm")

sapply(c("tram", "survival", "tbm", "mboost", "lattice", "latticeExtra"), library, char = TRUE)
library("animation")

trellis.par.set(list(plot.symbol = list(col=1,pch=20, cex=0.7),
                     box.rectangle = list(col=1),
                     box.umbrella = list(lty=1, col=1),
                     strip.background = list(col = "white")))
ltheme <- canonical.theme(color = FALSE)     ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)

knitr::opts_chunk$set(echo = TRUE, results = 'markup', error = FALSE,
                      warning = FALSE, message = FALSE,
                      tidy = FALSE, cache = FALSE, size = "small",
                      fig.width = 6, fig.height = 4, fig.align = "center",
                      out.width = NULL, ###'.6\\linewidth', 
                      out.height = NULL,
                      fig.scap = NA)
knitr::render_sweave()  # use Sweave environments
knitr::set_header(highlight = '')  # do not \usepackage{Sweave}
## R settings
options(prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE)  # JSS style
options(width = 75)
library("colorspace")
col <- diverge_hcl(2, h = c(246, 40), c = 96, l = c(65, 90))
fill <- diverge_hcl(2, h = c(246, 40), c = 96, l = c(65, 90), alpha = .3)

## ----citation, echo = FALSE----------------------------------------------
year <- substr(packageDescription("tbm")$Date, 1, 4)
version <- packageDescription("tbm")$Version

## ----results, echo = FALSE-----------------------------------------------
ds <- c("Beetles Exctinction Risk",
        "Birth Weight Prediction",
        "Body Fat Mass",
        "CAO/ARO/AIO-04 DFS",
        "Childhood Malnutrition",
        "Head Circumference",
        "PRO-ACT ALSFRS",
        "PRO-ACT OS")
names(ds) <- c("ex_beetles.rda", "ex_fetus.rda", "ex_bodyfat.rda", "ex_CAO.rda", 
        "ex_india.rda", "ex_heads.rda", "ex_ALSFRS.rda", "ex_ALSsurv.rda")
m <- c(expression(paste("L ", bold(vartheta)(bold(x)))),
       expression(paste("L ", beta(bold(x)))),
       expression(paste("N ", bold(vartheta)(bold(x)))),
       expression(paste("N ", beta(bold(x)))),
       expression(paste("T ", bold(vartheta)(bold(x)))),
       expression(paste("T ", beta(bold(x)))),
       expression(paste("Tree ", bold(vartheta)(bold(x)))),
       expression(paste("F ", bold(vartheta)(bold(x)))))
names(m) <- c("glm_ctm", "glm_tram", "gam_ctm", "gam_tram", 
              "tree_ctm", "tree_tram", "trtf_tree", "trtf_forest")
mor <- c("gam_ctm", "glm_ctm", "tree_ctm", "gam_tram", "glm_tram", 
        "tree_tram", "trtf_tree", "trtf_forest")

pit <- function(object, data) {
    cf <- predict(object, newdata = data, coef = TRUE)
    tmp <- object$model
    ret <- numeric(NROW(data))
    for (i in 1:NROW(data)) {
        coef(tmp) <- cf[i,]
        ret[i] <- predict(tmp, newdata = data[i,,drop = FALSE], type = "distribution")
    }
    ret
}

pitplot <- function(object, data, main = "") {
    u <- pit(object, data)
    plot(1:length(u) / length(u), sort(u), xlab = "Theoretical Uniform Quantiles", 
         ylab = "PIT Quantiles", main = main, ylim = c(0, 1), xlim = c(0, 1))
    abline(a = 0, b = 1, col = "red")
} 

plot.ctm <- function(x, which = sort(unique(selected(object))), 
                     obs = TRUE, ...) {

    object <- x
    q <- mkgrid(x$model, n = 100)
    X <- model.matrix(x$model$model, data = q)
    pfun <- function(which) {
        stopifnot(length(which) == 1)
        df <- model.frame(object, which = which)[[1]]
        df <- lapply(df, function(x) seq(from = min(x), 
                                         to = max(x), length = 50))
        df2 <- do.call("expand.grid", c(df, q))
        class(object) <- class(object)[-1]
        tmp <- predict(object, newdata = df, which = which)
        CS <- diag(ncol(tmp[[1]]))
        CS[upper.tri(CS)] <- 1
        cf <- tmp[[1]] %*% CS
        df2$pr <- as.vector(t(X %*% t(cf)))
        ret <- cbind(df2, which = names(df)[1])
        names(ret)[1:2] <- c("x", "y")
        ret
    }

    ret <- do.call("rbind", lapply(which, pfun))
    at <- quantile(ret$pr, prob = 0:10/10)

    pf <- function(x, y, z, subscripts, ...) {
        xname <- as.character(unique(ret[subscripts, "which"]))
        xo <- model.frame(object, which = xname)[[1]][[xname]]
        yo <- object$response
        panel.levelplot.raster(x, y, z, subscripts, ...)
        panel.contourplot(x, y, z, subscripts, contour = TRUE, ...)
        if (obs)
            panel.xyplot(x = xo, y = yo, pch = 20)
    }
    levelplot(pr ~ x + y | which, data = ret, panel = pf,
              scales = list(x = list(relation = "free")), region = TRUE, 
              at = at, col.regions = grey.colors(length(at)), ...)
}




## ----bodyfat-results, echo = FALSE, fig.width = 7, fig.height = 6--------
load(file.path(dir, file <- "ex_bodyfat.rda"))
ylim <- c(0, 1.7)
x <- risk[,mor] - ll0
med <- apply(x, 2, median, na.rm = TRUE)
boxplot(x, axes = FALSE, main = ds[file], outline = FALSE, var.width = TRUE,
        ylab = "Out-of-sample log-likelihood (centered)", col = col[(med == max(med)) + 1],
        ylim = ylim)
abline(v = c(3.5, 6.5), col = "lightgrey")
abline(h = 0, col = "lightgrey")
abline(h = max(med), col = col[2])
axis(1, at = 1:ncol(x), labels = m[colnames(x)])
axis(2)
box()

## ----bodyfat-setup, echo = FALSE-----------------------------------------

data("bodyfat", package = "TH.data")
ldata <- bodyfat

xvars <- colnames(ldata)
xvars <- xvars[xvars != "DEXfat"]

ldata[xvars] <- lapply(xvars, function(x) as.numeric(scale(ldata[[x]])))


fm_gam <- c("ctm" = as.formula(paste("DEXfat ~ ", 
                paste("bbs(", xvars, ")", collapse = "+"))),
            "stm" = as.formula(paste("DEXfat ~ ", 
                paste("bols(", xvars, ", intercept = FALSE)", collapse = "+"), "+",
                paste("bbs(", xvars, ", center = TRUE, df = 1)", collapse = "+"))))
fm_glm <- c("ctm" = as.formula(paste("DEXfat ~ ",
                paste("bols(", xvars, ")", collapse = "+"))),
            "stm" = as.formula(paste("DEXfat ~ ", 
                paste("bols(", xvars, ", intercept = FALSE)", collapse = "+"))))

fm_tree <- DEXfat ~ .

mf <- model.frame(fm_tree, data = ldata)
xf <- sum(sapply(mf[,-1], is.factor))
N <- dim(mf)[1]
p <- dim(mf)[2] - 1
vars <- names(mf)[-1]

## ----bodyfat-model0, cache = TRUE----------------------------------------
m_mlt <- BoxCox(DEXfat ~ 1, data = ldata, prob = c(.1, .99))
logLik(m_mlt)

## ----bodyfat-null, echo = FALSE------------------------------------------
plot(as.mlt(m_mlt), newdata = data.frame(1), type = "density", K = 200,
     xlab = "Body Fat Mass (in kilogram)", ylab = "Density", col = "black")

## ----bodyfat-model, cache = TRUE-----------------------------------------
fm_gam[["ctm"]]
bm <- ctmboost(m_mlt, formula = fm_gam[["ctm"]], data = ldata,
               method = quote(mboost::mboost))[101]

## ----bodyfat-dens, echo = FALSE------------------------------------------
q <- mkgrid(m_mlt, n = 200)[[1]]


bm[101]
d0 <- predict(bm, newdata = ldata, type = "density", q = q)
r <- risk(bm)
saveLatex(
for (m in 0:100)  {
    dev.hold()
    layout(matrix(1:2, ncol = 2))

    if (m == 0) {
        plot(as.mlt(m_mlt), newdata = data.frame(1), q = q, type = "density",
             col = rgb(.1, .1, .1, .6), lty = 1, xlab = "Body Fat Mass (in kilogram)", 
             ylab = "Density", ylim = c(0, max(d0)))
        ll <- exp(m_mlt$logliki(coef(as.mlt(m_mlt)), rep(1, nrow(ldata))))
        points(ldata$DEXfat, ll, pch = 19, col = col[1])
        plot(-r, type = "l", xlab = "Boosting Iteration", ylab = "Log-Likelihood")
        ani.pause(1)
    } else {
        matplot(x = q, y = d <- predict(bm[m], newdata = ldata, type = "density", q = q), 
                type = "l", col = rgb(.1, .1, .1, .6), lty = 1, xlab = "Body Fat Mass (in kilogram)", 
                ylab = "Density", ylim = c(0, max(d0)))
        j <- sapply(ldata[, "DEXfat"], function(y) which.min((q - y)^2))
        points(ldata[,"DEXfat"], sapply(1:length(j), function(i) d[j[i], i]), pch = 19, col = col[1])
        plot(-r, type = "l", xlab = "Boosting Iteration", ylab = "Log-Likelihood")
        points(m + 1, -rev(risk(bm[m]))[1], pch = 19, col = col[1])
        ani.pause(.1)
    }
}, ani.height = 380, ani.width = 580, img.name = "fit_movie"
)
