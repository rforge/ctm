
library("ctm")   

set.seed(29)

load("india.Rda")

### only stunting and district
kids <- kids[, c("stunting", "distH")]

mod <- ctm(bbs(stunting, lambda = 100) ~ bmrf(distH, bnd = nb, lambda = 100), 
              data = kids,
              family = Binomial(link = "probit"), control = boost_control(nu = 0.4, 
              mstop = 100, trace = TRUE), ngrid = 50)

nd <- data.frame(distH = unique(kids$distH))
nd$ID <- 1:nrow(nd)
x <- predict(mod, newdata = nd, type = "response", y = mod$uresponse, annotated = TRUE)
x <- merge(nd, x, by = "ID")
x$ID <- NULL

save(x, file = "ex_india.Rda")

