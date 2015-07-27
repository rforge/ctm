
library("mlt")
set.seed(29)

n <- 20
### design has rank deficit; just for interface checking
### we need something better!
d <- data.frame(x1 = 1:n, x2 = 1:n + 1, y = rnorm(n))
m <- model(polynomial_basis(c(TRUE, TRUE), varname = "y",
           ci = c(-Inf, 0), support = range(d$y)),
           shift = ~ x1 + x2)
mod <- mlt(m, data = d)

p <- tmlt(mod, newdata = d)

(p0 <- predict(mod$model$model, 
    newdata = expand.grid(d), coef = coef(mod)))
(p1 <- tmlt(mod, newdata = as.list(d)))
(p2 <- tmlt(mod, newdata = d, q = d$y[1]))

max(abs(p0 - as.vector(p1)))

all.equal(p1[cbind(1:n, 1:n, 1), drop = TRUE],
          drop(p2))

all.equal(p1[cbind(1:n, 1:n, 1:n), drop = TRUE],
          drop(p), check.attributes = FALSE)

qmlt(mod, newdata = list(x1 = 1:3, x2 = 2:3), p = c(.25, .5))

simulate(mod, nsim = 1, seed = 291, interpolate = FALSE)

d$y <- gl(3, 1, ordered = TRUE)[rep(1:3, length = n)]

r <- as.basis(~ y, data = d, remove_intercept = TRUE,
              contrasts.arg = list(y = function(n)
                  contr.treatment(n, base = 3)),
              ui = diff(diag(2)), ci = 0)

mod2 <- mlt(model(r, shift = ~ x1 + x2), data = d)

tmlt(mod2, q = unique(d$y))

qmlt(mod2, p = 1:9 / 10)

simulate(mod2, nsim = 3, seed = 29)

dmlt(mod2, q = unique(d$y))

dmlt(mod2, list(y = unique(d$y), x1 = 1:3, x2 = 2:3))

