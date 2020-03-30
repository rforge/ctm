
library("cotram")

set.seed(25)


## dgp
dgp <- function(n = 200, seed = 25){
        set.seed(seed)
        x <- runif(n, 0, 20)
        y <- as.integer(rnbinom(n, mu = exp(.2 + .1 * x), size = 3))
        data.frame(x = x, y)
}

y <- 1:10
m1 <- cotram(y ~ x, data = dgp())


## Confband for grid of counts
confband(m1, type = "distribution", newdata = model.frame(m1)[3,])

## Confband for K grid points
confband(m1, type = "distribution", newdata = model.frame(m1)[3, ],
                 smooth = TRUE, K = 40)


if (FALSE){
        layout(matrix(1:2, nrow = 1))
        type = "trafo"
        nd <- model.frame(m1)[3,]
        cb <- confband(m1, type =  type, newdata = nd)
        plot(m1, type = type, newdata = nd, 
             confidence = "band", col = "red", ylim = c(-2, 15))
        lines(x = cb[, "q"], y = cb[, "lwr"], type = "s")
        lines(x = cb[, "q"], y = cb[, "upr"], type = "s")
        
        cb.s <- confband(m1, type = type, newdata = nd, 
                       smooth = TRUE)
        plot(m1, type = type, newdata = nd, 
             confidence = "band", col = "red", smooth = TRUE, ylim = c(-2, 15))
        lines(x = cb.s[, "q"], y = cb.s[, "lwr"], type = "l")
        lines(x = cb.s[, "q"], y = cb.s[, "upr"], type = "l")
y}


