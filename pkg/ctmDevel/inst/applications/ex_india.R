
source("setup.R")
set.seed(29)

load("india.Rda")

### only stunting and district
kids <- kids[, c("stunting", "distH")]

### extract lambda from univariate baselearners
### not really what is needed but a good guess
a <- with(kids, bbs(stunting, df = 4))$dpp(rep(1, nrow(kids)))
lambda1 <- extract(a, what = "lambda")
b <- with(kids, bmrf(distH, bnd = nb, df = 4))$dpp(rep(1, nrow(kids)))
lambda2 <- extract(b, what = "lambda")

mod <- ctm(bbs(stunting, lambda = lambda1) ~ bmrf(distH, bnd = nb, lambda = lambda2), 
              data = kids,
              family = Binomial(link = "probit"), control = boost_control(nu = 0.4, 
              mstop = 100, trace = TRUE), ngrid = 50)

nd <- data.frame(distH = unique(kids$distH))
nd$ID <- 1:nrow(nd)
x <- predict(mod, newdata = nd, type = "response", y = mod$uresponse, annotated = TRUE)
x <- merge(nd, x, by = "ID")
x$ID <- NULL

save(x, file = "ex_india.Rda")

