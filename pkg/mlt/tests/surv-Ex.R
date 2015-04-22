
library("mlt")
library("survival")
set.seed(29)

y <- rexp(10)
d <- rep(c(TRUE, FALSE), length = length(y))
### right
(s <- Surv(y, d, type = "right"))
attr(s, "type")
cbind(unclass(s), mlt:::.Surv2matrix(s))

### right censored, left truncated
(s <- Surv(y, y + 1, d))
attr(s, "type")
cbind(unclass(s), mlt:::.Surv2matrix(s))

### left censored
(s <- Surv(y, d, type = "left"))
attr(s, "type")
cbind(unclass(s), mlt:::.Surv2matrix(s))

### interval
dd <- rep(0:3, length = length(y))
(s <- Surv(y, y + 1, dd, type = "interval"))
attr(s, "type")
cbind(unclass(s), mlt:::.Surv2matrix(s))

### interval2
(s <- Surv(y, y + 1, type = "interval2"))
attr(s, "type")
cbind(unclass(s), mlt:::.Surv2matrix(s))

(s <- Surv(y, ifelse(d, Inf, y + 1), type = "interval2"))
attr(s, "type")
cbind(unclass(s), mlt:::.Surv2matrix(s))

(s <- Surv(y, ifelse(d, NA, y + 1), type = "interval2"))
attr(s, "type")
cbind(unclass(s), mlt:::.Surv2matrix(s))

### new interface
tm <- runif(10)
e <- sample(c(TRUE, FALSE), length(tm), replace = TRUE)
s <- Surv(tm, e)
cbind(s, as.data.frame(mlt:::.Surv2matrix(s)))
cbind(s, mlt:::.Surv2R(s))

s <- Surv(tm, e, type = "left")
cbind(s, as.data.frame(mlt:::.Surv2matrix(s)))
cbind(s, mlt:::.Surv2R(s))

s <- Surv(tm, tm + 1, type = "interval2")
cbind(s, as.data.frame(mlt:::.Surv2matrix(s)))
cbind(s, mlt:::.Surv2R(s))

s <- Surv(tm, tm + 1, sample(0:3, length(tm), replace = TRUE), type = "interval")
cbind(s, as.data.frame(mlt:::.Surv2matrix(s)))
cbind(s, mlt:::.Surv2R(s))

s <- Surv(tm, tm + 1, e, type = "counting")
cbind(s, as.data.frame(mlt:::.Surv2matrix(s)))
cbind(s, mlt:::.Surv2R(s))

R(gl(3, 2, ordered = TRUE))
R(0:10)

