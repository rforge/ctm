
set.seed(290875)

ret <- TRUE

if (require("rms") && require("tram")) {

    x <- runif(2000)
    y <- round(rnorm(length(x), mean = 2 * x, sd = .25), 1)
    d <- data.frame(y = y, x = x)
    d$Ry <- with(d, R(y, as.R.ordered = TRUE))

    mP <- Polr(Ry ~ x, data = d)
    mM <- Polr(Ry ~ x, data = d, sparse_nlevels = 2)
    mO <- orm(y ~ x, data = d)

    tol <- 1e-4
    ret <- 
      isTRUE(all.equal(coef(mP), coef(mM), tol = tol)) && 
      isTRUE(all.equal(coef(mP), coef(mO)["x"], tol = tol)) && 
      isTRUE(all.equal(coef(as.mlt(mP)), coef(as.mlt(mM)), tol = tol)) &&
      isTRUE(all.equal(rev(coef(as.mlt(mP)))[-1L], 
                       -rev(coef(mO))[-1L], tol = tol, 
                       check.attributes = FALSE)) &&
      isTRUE(all.equal(logLik(mP), logLik(mM), tol = tol)) &&
      isTRUE(all.equal(c(logLik(mP)), c(logLik(mO)), tol = tol)) &&

      isTRUE(all.equal(c(vcov(mP)), as.numeric(vcov(mM)), tol = tol)) &&
      isTRUE(all.equal(c(vcov(mP)), vcov(mO)["x", "x"], tol = tol))

}

ret
