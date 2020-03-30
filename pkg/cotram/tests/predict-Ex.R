
library("cotram")

set.seed(25)

yr <- rpois(1000, lambda = 2)
y <- 0:10
yy <- seq(from = min(y), to = max(y), length.out = 50)

pr <- 1:9/10

## cotram: log_first = FALSE
m0 <- cotram(yr ~ 1, log_first = FALSE, extrapolate = TRUE, prob = .99)


stopifnot(all.equal(predict(m0, newdata = data.frame(1), q = 0, type = "density"), 
                    predict(m0, newdata = data.frame(1), q = 0, type = "distribution")))

stopifnot(!isTRUE(all.equal(predict(m0, newdata = data.frame(1), q = 1, type = "density"), 
                            predict(m0, newdata = data.frame(1), q = 1, type = "distribution"))))

## cotram: log_first = TRUE
m0_lf <- cotram(yr ~ 1, log_first = TRUE, extrapolate = TRUE, prob = .99)

stopifnot(all.equal(predict(m0_lf, newdata = data.frame(1), q = 0, type = "density"), 
                    predict(m0_lf, newdata = data.frame(1), q = 0, type = "distribution")))

stopifnot(!isTRUE(all.equal(predict(m0_lf, newdata = data.frame(1), q = 1, type = "density"), 
                            predict(m0_lf, newdata = data.frame(1), q = 1, type = "distribution"))))

stopifnot(all.equal(predict(m0, newdata = data.frame(y), type = "density"), 
                    predict(m0, newdata = data.frame(1), q = y, type = "density")))

stopifnot(all.equal(predict(m0, newdata = data.frame(yr = y), type = "density"), 
                    predict(m0, newdata = data.frame(1), q = y, type = "density")))

# mt <- Colr(yr ~ 1, log_first = FALSE, extrapolate = TRUE, bounds = m0$bounds, support = m0$support)
# predict(mt, newdata = data.frame(1), type = "quantile", prob = pr)
# predict(mt, newdata = data.frame(1), q = y, type = "quantile", prob = pr)

if (FALSE) {
        layout(matrix(c(1:4), nrow = 2, byrow = TRUE))
        
        main <- "predict - density: log_first = FALSE"
        plot(y, predict(m0, newdata = data.frame(1), q = y , type = "density"), 
             type = "h", main = main)
        points(y, predict(m0, newdata = data.frame(1), q = y, type = "density"), pch = 20)
        points(y, dpois(y, lambda = 2), col = "blue")
        lines(yy, predict(m0, newdata = data.frame(1), q = yy, type = "density", smooth = TRUE), col = "grey")
        lines(yy, predict(m0, newdata = data.frame(1), type = "density", smooth = TRUE), 
             col = "red")
        plot(m0, newdata = data.frame(1), type = "density", smooth = TRUE, 
              col = "red", add = TRUE)
        plot(m0, newdata = data.frame(1), type = "density", 
             col = "red", add = TRUE)
        
        main <- "predict - distribution: log_first = FALSE"
        plot(y, predict(m0, newdata = data.frame(1), q = y, type = "distribution"),
             type = "s", ylim = c(0, 1), main = main)
        lines(yy, predict(m0, newdata = data.frame(1), q = yy, type = "distribution", smooth = TRUE), col = "grey")
        points(y, ppois(y, lambda = 2), col = "blue")
        points(predict(m0, newdata = data.frame(1), q = y, prob = pr, type = "quantile", 
                       smooth = TRUE), pr, col = "darkgreen")
        legend("bottomright", legend = c("ppois", "quantiles"), pch = 1, 
               col = c("blue", "green"))
        plot(m0, newdata = data.frame(1), type = "distribution", 
             col = "red", smooth = TRUE, add = TRUE)
        plot(m0, newdata = data.frame(1), q = y, type = "distribution", 
             col = "red", add = TRUE)
        
        main <- "predict - density: log_first = TRUE"
        plot(y, predict(m0_lf, newdata = data.frame(1), q = y, type = "density"), 
             type = "h", main = main)
        points(y, predict(m0_lf, newdata = data.frame(1), q = y, type = "density"), pch = 20)
        points(y, dpois(y, lambda = 2), col = "blue")
        lines(yy, predict(m0_lf, newdata = data.frame(1), type = "density",  q = yy, smooth = TRUE))
        plot(m0_lf, newdata = data.frame(1), type = "density", 
             col = "red", smooth = TRUE, add = TRUE, ylim = c(0, .3))
        plot(m0_lf, newdata = data.frame(1), q = y, type = "density", 
             col = "red", add = TRUE)
        
        main <- "predict - distribution: log_first = TRUE"
        plot(y, predict(m0_lf, newdata = data.frame(1), q = y, type = "distribution"),
             type = "s", ylim = c(0, 1), main = main)
        points(y, ppois(y, lambda = 2), col = "blue")
        lines(yy, predict(m0_lf, newdata = data.frame(1), q = yy, type = "distribution", smooth = TRUE))
        points(predict(m0_lf, newdata = data.frame(1), 
                       smooth = TRUE, prob = pr, type = "quantile"), pr, col = "green")
        legend("bottomright", legend = c("ppois", "quantiles"), pch = 1, 
               col = c("blue", "green"))
        plot(m0_lf, newdata = data.frame(1), type = "distribution", q = 0:3,
             col = "red", smooth = TRUE, add = TRUE)
        plot(m0_lf, newdata = data.frame(1), type = "distribution", 
             col = "red", add = TRUE)
}


## Check settings predict-function with tram
dgp <- function(n = 200){
        x <- runif(n)
        y <- as.integer(rnbinom(n, mu = exp(.5 + .8 * x), size = 10))
        y.p1 <- y + 1L
        data.frame(x = x, y = y, y.p1 = y.p1)
}
df <- dgp()
m1 <- cotram(y ~ x, data = df, method = "cloglog")

m2 <- Coxph(y.p1 ~ x, data = df,
            support = m1$support,
            bounds = m1$bounds,
            log_first = TRUE)

q <- m1$count_response$support - 1 
q2 <- m1$count_response$support

nd <- data.frame(y = q,
                 y.p1 = q + 1, 
                 x = seq(0, 1, length.out = length(q)))

stopifnot(all.equal(
        predict(m1, newdata = nd), 
        (-1) * m1$negative * predict(m2, newdata = nd)
))

stopifnot(all.equal(
        predict(m1, newdata = nd, type = "distribution"), 
        predict(m2, newdata = nd, type = "distribution"), 
        check.attributes = FALSE
))


stopifnot(all.equal(
        predict(m1, newdata = nd, q = q, type = "density"), 
        predict(m2, newdata = nd, q = q2, type = "distribution") - 
        predict(m2, newdata = nd, q = q2 - 1, type = "distribution"),
        check.attributes = FALSE
))

nd.m1 <- data.frame(y.p1 = nd$y.p1 - 1 , x = nd$x)
stopifnot(all.equal(
        predict(m1, newdata = nd, type = "density"), 
        predict(m2, newdata = nd, type = "distribution") - 
                predict(m2, newdata = nd.m1, type = "distribution"),
        check.attributes = FALSE
))

stopifnot(all.equal(
        predict(m1, newdata = nd[, "x", drop = FALSE], q = q, type = "density"), 
        predict(m2, newdata = nd[, "x", drop = FALSE], type = "distribution", q = q2) - 
                predict(m2, newdata = nd[, "x", drop = FALSE], type = "distribution", q = q2 - 1),
        check.attributes = FALSE
))

