library("cotram")

set.seed(29)

## 

## dgp
dgp <- function(n = 200){
  x <- runif(n)
  yd <- as.integer(rnbinom(n, mu = exp(.5 + .8 * x), size = 10))
  yn <- as.numeric(yd)
  yd.p1 <- yd + 1L
  data.frame(x = x, yd = yd, yn = yn, yd.p1 = yd.p1)
}

df <- dgp()
trainID <- sample(1:nrow(df), size = 0.25 * nrow(df))
df.train <- df[trainID,]
df.test <- df[-trainID,]

## test model
m1d <- cotram(yd ~ x, data = df.train, method = "cloglog")

m1n <- cotram(yn ~ x, data = df.train, method = "cloglog")

m3d <- cotram(yd ~ x , data = df.train,
              log_first = FALSE, method = "cloglog")

m2 <- Coxph(yd.p1 ~ x, data = df.train,
             support = m1d$support,
             bounds = m1d$bounds,
             log_first = TRUE)

# compare coefficients
cf1n <- coef(as.mlt(m1n))
cf1d <- coef(as.mlt(m1d))
stopifnot(all.equal(cf1n, cf1d, check.attributes = FALSE))

c2 <- as.mlt(m2)
cf2 <- coef(c2)
cf2[c2$shiftcoef] <- -cf2[c2$shiftcoef]
stopifnot(all.equal(cf1n, cf2, check.attributes = FALSE))

# compare likelihoods
l1d <- m1d$logliki(coef(as.mlt(m1d)), rep(1, length(df.train$yd)))
l1n <- m1d$logliki(coef(as.mlt(m1n)), rep(1, length(df.train$yn)))
stopifnot(all.equal(l1d, l1n))

l2 <- m2$logliki(coef(as.mlt(m2)), rep(1, length(df.train$yd.p1)))
stopifnot(all.equal(l1d, l2))

# logLik for newdata
stopifnot(isTRUE(logLik(m1d) == logLik(m1d, newdata = df.train)))
stopifnot(isTRUE(logLik(m1d, newdata = df.test) == logLik(m2, newdata = df.test)))
stopifnot(isTRUE(logLik(m1d, newdata = df.test) == logLik(m1n, newdata = df.test)))
