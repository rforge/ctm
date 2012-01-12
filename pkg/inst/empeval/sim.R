
source("simfuns.R")

set.seed(290875)

sfun <- function(seed, pnon) {

    set.seed(seed)
    xdf <- dgp(c(100, 100), pnon = pnon)
    test <- dgp(pnon = pnon)

    boost <- Fboost(xdf)
    boostl <- Fboost(xdf, ylin = TRUE)
    kern <- Fnp(xdf)
    glss <- Fgamlss(xdf)

    pt <- pdf(truth, test)
    pb <- pdf(boost, test)
    pbl <- pdf(boostl, test)
    pk <- pdf(kern, test)
    pg <- pdf(glss, test)

    boost[20]
    tune(boost, alpha = 0.05, mstopmax = 1200)
    pbt <- pdf(boost, test)

    boostl[20]
    tune(boostl, alpha = 0.05, mstopmax = 1200)
    pblt <- pdf(boostl, test)

    list(boosting = summary(colMeans(abs(pt - pb))),
         boostingt = summary(colMeans(abs(pt - pbt))),
         boostinglin = summary(colMeans(abs(pt - pbl))),
         boostinglint = summary(colMeans(abs(pt - pblt))),
         kernel = summary(colMeans(abs(pt - pk))),
         gamlss = summary(colMeans(abs(pt - pg)))
    )
}

pnon <- c(0, 1, 2, 3, 4, 5)

ret <- vector(mode = "list", length = length(pnon))
names(ret) <- pnon

for (p in pnon) {

    ### 5 x 20 cores 
    ### make sure to set up the seeds before 
    ### distributing the computations
    seeds <- matrix(round(runif(100) * 10000), nrow = 5)
    for (j in 1:nrow(seeds)) {
        cat("j: ", j, " pnon: ", p, "\n")
        ret[[as.character(p)]] <- c(ret[[as.character(p)]],
                                    mclapply(seeds[j,], sfun, pnon = p))
        save(ret, file = "ret.Rda")
    }
    save(ret, file = "ret.Rda")
}

save(ret, file = "ret.Rda")
