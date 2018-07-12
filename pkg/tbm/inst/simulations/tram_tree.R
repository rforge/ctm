
library("tram")
library("tbm")
library("partykit")

source("setup.R")

myFUN <- function(ldata, lweights, model = c("normal", "logistic", "minextrval"), order) {

    ### we use large trees, so train with nu = 0.01 instead of nu 0.1
    bctrl <- boost_control(mstop = 100, risk = "oob", nu = 0.01)

    model <- match.arg(model)
    m0 <- switch(model, 
        "normal" = BoxCox(y ~ 1, data = ldata, order = order, support = sup, bounds = bds),
        "logistic" = Colr(y ~ 1, data = ldata, order = order, support = sup, bounds = bds),
        "minextrval" = Coxph(y ~ 1, data = ldata, order = order, support = sup, bounds = bds))

    fm <- "y ~ x1 + x2"
    nm <- colnames(ldata)
    nx <- nm[grep("^nx", nm)]
    if (length(nx) > 0)
        fm <- paste(fm, "+", paste(nx, collapse = "+"))
    fm <- as.formula(fm)
    tctrl <- ctree_control(minsplit = 4, minbucket = 2, mincriterion = 0,
                           maxdepth = 6, splittest = FALSE, 
                           testtype = "Teststatistic")
    l1 <- stmboost(m0, formula = fm,
                   data = ldata, method = quote(mboost::blackboost), 
                   weights = lweights, control = bctrl, tree_control = tctrl)
    while(!(which.min(risk(l1)) < mstop(l1))  && mstop(l1) < 2000)
        l1 <- l1[2 * mstop(l1)]
    mstop <- which.min(risk(l1))
    l1 <- stmboost(m0, formula = fm,
                   data = ldata, method = quote(mboost::blackboost), 
                   control = bctrl, tree_control = tctrl)[mstop]
    return(l1)
}

ret <- c()

for (NOBS in tNOBS) {
    for (PNON in tPNON) {
        for (TD in tTD) {
            for (OR in tOR * order) {

                FUN <- function(...) myFUN(..., model = TD, order = OR)

                source("run.R", echo = TRUE)

                res$model <- "tram_tree"
                res$PNON <- PNON
                res$NOBS <- NOBS
                res$order <- OR
                res$todistr <- TD
                ret <- rbind(ret, res)
            }
        } 
    }
}

save(ret, file = "tram_tree.rda")
sessionInfo()
