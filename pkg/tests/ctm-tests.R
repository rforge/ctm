
library("mboost")
library("lattice")
library("ctm")

trf <- function(formula, data, grid = 50) {

    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$formula <- formula
    mf$drop.unused.levels <- FALSE
    mf$na.action <- na.pass
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    response <- names(mf)[1]
    y <- mf[[1]]
    if (inherits(y, "Surv")) {
        mf <- mf[-1]	
        mf$time <- y[,1]
        mf$event <- as.logical(y[,2])
        response <- "time"
        y <- y[ ,1]
    } else {
        mf$event <- TRUE
    }
    if (!(is.ordered(y) || is.numeric(y)))
        stop("response is not ordered")
    if (any(is.na(y)))
        stop("response must not contain missing values")
    ngrid <- 50
    if (!is.null(grid) & length(grid) == 1) {
        ngrid <- grid
        grid <- NULL
    }
    if (is.null(grid)) {
        if (is.ordered(y))
            grid <- levels(y)[-nlevels(y)]
        if (is.numeric(y))
            grid <- seq(from = min(y), to = max(y), length = ngrid)
    }

    n <- nrow(mf)
    p <- length(grid)
    mf$obs <- factor(paste("ID", 1:nrow(mf), sep = ""))  
    mf <- mf[rep(1:n, rep(p, n)),]
    mf$grid <- rep(grid, n)
    mf$ytmp <- factor(mf[[response]] <= grid)
    ### remove entries after censoring time
    ### mf <- subset(mf, !(ytmp == "TRUE" & !mf[["event"]]))
    if (is.ordered(y))   
        mf$grid <- ordered(mf$grid, levels = grid)
    mf
}


x <- runif(200)
y <- rnorm(length(x), mean = sin(x * 3) + x, sd = .25 + .25 * (1 - x^2))
df <- data.frame(x = x, y = y)
yy <- sort(unique(y))
yy <- yy[-length(yy)]

xg <- seq(from = min(x), to = max(x), length = 50)
yg <- seq(from = min(yy), to = max(yy), length = 60)
nd <- expand.grid(x = xg, y = yg)
nd$ptrue <- pnorm(nd$y, mean = sin(nd$x * 3) + nd$x, 
                  sd = .25 + .25 * (1 - nd$x^2))


wireframe(ptrue ~ x + y, data = nd)

## df <- df[order(df$y),]
plot(y ~ x, data = df)

system.time(m1 <- ctm(bbs(y, df = sqrt(6)) ~ bbs(x, df = sqrt(6)), 
             data = df, weights = rep(1, nrow(df)),
             family = Binomial()))

nd$p <- p1 <- predict(m1, newdata = data.frame(x = xg), y = yg,
                      type = "response")
wireframe(p ~ x + y, data = nd)
sapply(xg, function(a) min(diff(subset(nd, x == a)$p)))


df2 <- trf(y ~ x, data = df, grid = yy)
df2$y <- NULL
names(df2)[names(df2) == "grid"] <- "y"
system.time(m2 <- mboost(ytmp ~ bbs(y, df = sqrt(6)) %X% bbs(x, df = sqrt(6)), 
            data = df2, family = Binomial(), offset = 0))
nd$p <- p2 <- predict(m2, newdata = nd, type = "response")
wireframe(p ~ x + y, data = nd)
sapply(xg, function(a) min(diff(subset(nd, x == a)$p)))


system.time(m3 <- mboost(ytmp ~ bspatial(y, x),
            data = df2, family = Binomial(), offset = 0))
nd$p <- p3 <- predict(m3, newdata = nd, type = "response")
wireframe(p ~ x + y, data = nd)
sapply(xg, function(a) min(diff(subset(nd, x == a)$p)))

max(abs(p1 - p2))
max(abs(p1 - p3))

cf1 <- coef(m1)[[1]]
cf2 <- coef(m2)[[1]]
cf3 <- coef(m3)[[1]]

max(abs(cf1 - cf2))
max(abs(cf1 - cf3))

nd$pp <- ifelse(nd$ptrue < .5, nd$ptrue, 1 - nd$ptrue)
levelplot(pp ~ x + y, data = nd)

plot(m1)

data("bodyfat", package = "mboost")

nm <- names(bodyfat)
nm <- nm[nm != "DEXfat"]
bl <- paste("bbs(", nm, ", df = 2.5)", collapse = "+")
fm <- as.formula(paste("bbs(DEXfat, df = 2.5) ~ ", bl))
mod <- ctm(fm, data = bodyfat, family = Binomial())

plot(mod, which = sort(unique(selected(mod))))

mod1 <- mboost(DEXfat ~ ., data = bodyfat)
layout(matrix(1:9, nc = 3))
plot(mod1, which = sort(unique(selected(mod))))
