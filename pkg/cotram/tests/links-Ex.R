library("cotram")

set.seed(29)

## dgp
dgp <- function(n = 200){
  x <- runif(n)
  y <- as.integer(rnbinom(n, mu = exp(.5 + .8 * x), size = 10))
  y.p1 <- y + 1L
  data.frame(x = x, y = y, y.p1 = y.p1)
}

df <- dgp()

# cloglog link
m1a <- cotram(y ~ x, data = df, method = "cloglog")

m1b <- Coxph(y.p1 ~ x, data = df,
             support = m1a$support,
             bounds = m1a$bounds,
             log_first = TRUE)

cf1a <- coef(as.mlt(m1a))

c1b <- as.mlt(m1b)
cf1b <- coef(c1b)
cf1b[c1b$shiftcoef] <- -cf1b[c1b$shiftcoef]
stopifnot(all.equal(cf1a, cf1b, check.attributes = FALSE))


# logit link
m2a <- cotram(y ~ x, data = df, method = "logit")

m2b <- Colr(y.p1 ~ x, data = df,
             support = m2a$support,
             bounds = m2a$bounds,
             log_first = TRUE)

cf2a <- coef(as.mlt(m2a))

c2b <- as.mlt(m2b)
cf2b <- coef(c2b)
cf2b[c2b$shiftcoef] <- -cf2b[c2b$shiftcoef]
stopifnot(all.equal(cf2a, cf2b, check.attributes = FALSE))

# loglog link
m3a <- cotram(y ~ x, data = df, method = "loglog")

m3b <- Lehmann(y.p1 ~ x, data = df,
            support = m3a$support,
            bounds = m3a$bounds,
            log_first = TRUE)

cf3a <- coef(as.mlt(m3a))

c3b <- as.mlt(m3b)
cf3b <- coef(c3b)
stopifnot(all.equal(cf3a, cf3b, check.attributes = FALSE))

# probit link
m4a <- cotram(y ~ x, data = df, method = "probit")

m4b <- BoxCox(y.p1 ~ x, data = df,
               support = m4a$support,
               bounds = m4a$bounds,
               log_first = TRUE)

cf4a <- coef(as.mlt(m4a))

c4b <- as.mlt(m4b)
cf4b <- coef(c4b)
stopifnot(all.equal(cf4a, cf4b, check.attributes = FALSE))