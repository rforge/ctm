
library("mlt")

x <- rexp(100)
sr <- Surv(x, sample(c(TRUE, FALSE), length(x), replace = TRUE))
sl <- Surv(x, sample(c(TRUE, FALSE), length(x), replace = TRUE), type = "left")
si <- Surv(x, x + 1, sample(0:3, length(x), replace = TRUE), type = "interval")

mlt:::.Surv2matrix(sr)
mlt:::.Surv2matrix(sl)
mlt:::.Surv2matrix(si)


