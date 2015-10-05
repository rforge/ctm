
library("mlt")
set.seed(29)

### Nadja Klein
dat <- data.frame(matrix(rnorm(300),ncol=3))
names(dat) <- c("y","x1","x2")
### set-up conditional transformation model for conditional
y <- numeric_var("y", support = c(min(dat$y), max(dat$y)), bounds = c(-Inf, Inf))
x1 <- numeric_var("x1", support = c(min(dat$x1), max(dat$x1)), bounds = c(-Inf, Inf)) 
x2 <- numeric_var("x2", support = c(min(dat$x2), max(dat$x2)), bounds = c(-Inf, Inf)) 
ctmm2 <- ctm(response = Bernstein_basis(y, order = 4, ui = "increasing"),
              interacting = c(x1=Bernstein_basis(x1, order = 3),
                              x2=Bernstein_basis(x2, order= 3)))
### fit model
mltm2 <- mlt(ctmm2, data = dat, check = FALSE)
(p <- predict(mltm2, newdata = data.frame(x1=0, x2 = 0), q = mkgrid(mltm2, n = 10)[["y"]]))
### plot data
plot(mltm2,newdata=expand.grid(x1=0:1, x2 = 0:1))
