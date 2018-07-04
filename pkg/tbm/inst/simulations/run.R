
res <- expand.grid(i = 1:NSIM, m = 1:length(model))

simfun <- function(m, i) {

    cat("m: ", m, " i: ", i, "\n")

    ret <- c("boost_ll" = NA, "true_ll" = NA, "mlt_ll" = NA, "mstop" = NA)

    if (PNON > 0) {
        ### mean centered
        Xnon <- matrix(runif(nrow(d) * PNON), ncol = PNON) - .5
        colnames(Xnon) <- paste("nx", 1:PNON, sep = "")
        d <- cbind(d, Xnon)
    }

    idx <- sample(1:nrow(d), 2 * NOBS, replace = FALSE)
    ind <- d[idx,,drop = FALSE]
    ind$y <- y[[m]][[i]][idx]

    test <- d
    test$y <- y[[m]][[NSIM + i]]

    t0 <- mlt(model[[m]], data = test, theta = coef(model[[m]]), 
              dofit = FALSE)
    coef(t0) <- coef(model[[m]])
    print(ret["true_ll"] <- logLik(t0))

    t1 <- mlt(model[[m]], data = ind, theta = coef(model[[m]]))
    print(ret["mlt_ll"] <- logLik(t1, newdata = test))

    w <- rep(c(1, 0), each = NOBS)

    ### mean centering for boosting
    ind$x2 <- ind$x2 - .5
    test$x2 <- test$x2 - .5

    l1 <- try(FUN(ind, w))
    if (!inherits(l1, "try-error")) {
        print(ret["boost_ll"] <- logLik(l1$model, newdata = test, 
                                        parm = coef(l1, newdata = test)))
        ret["mstop"] <- mstop(l1)
    }

    cat("----\n")

    return(ret)

}

mc.reset.stream()
x <- mclapply(1:NROW(res), 
    function(r) simfun(res[r, "m"], res[r, "i"]), mc.cores = MC_CORES)

res <- cbind(res, do.call("rbind", x))
