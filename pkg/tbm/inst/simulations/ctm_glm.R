
library("tram")
library("tbm")

source("setup.R")

myFUN <- function(ldata, lweights, model = c("normal", "logistic", "minextrval"), order) {

    bctrl <- boost_control(mstop = 100, risk = "oob", nu = 0.1)

    model <- match.arg(model)
    m0 <- switch(model, 
        "normal" = BoxCox(y ~ 1, data = ldata, order = order, support = sup, bounds = bds),
        "logistic" = Colr(y ~ 1, data = ldata, order = order, support = sup, bounds = bds),
        "minextrval" = Coxph(y ~ 1, data = ldata, order = order, support = sup, bounds = bds))

    fm <- "y ~ bols(x1) + bols(x2) + bols(x2, by = x1)"
    nm <- colnames(ldata)
    nx <- nm[grep("^nx", nm)]
    if (length(nx) > 0)
        fm <- paste(fm, "+", paste("bols(", nx, ")", collapse = "+"))
    fm <- as.formula(fm)
    l1 <- ctmboost(m0, formula = fm,
                   data = ldata, method = quote(mboost:::mboost), 
                   weights = lweights, control = bctrl)
    while(!(which.min(risk(l1)) < mstop(l1)) && mstop(l1) < 2000)
        l1 <- l1[2 * mstop(l1)]
    mstop <- which.min(risk(l1))
    l1 <- ctmboost(m0, formula = fm,
                   data = ldata, method = quote(mboost:::mboost), 
                   control = bctrl)[mstop]
    return(l1)
}

ret <- c()

for (NOBS in tNOBS) {
    for (PNON in tPNON) {
        for (TD in tTD) {
            for (OR in tOR * order) {

                FUN <- function(...) myFUN(..., model = TD, order = OR)

                source("run.R", echo = TRUE)

                res$model <- "ctm_glm"
                res$PNON <- PNON
                res$NOBS <- NOBS
                res$order <- OR
                res$todistr <- TD
                ret <- rbind(ret, res)
            }
        } 
    }
}

save(ret, file = "ctm_glm.rda")
sessionInfo()
