
library("ctm")   

data("db", package = "gamlss.data")
db$lage <- with(db, age^(1/3))

mod <- ctm(bbs(head, df = 2.5) ~ bbs(lage, df = 3, knots = 50), data = db,
              family = Binomial(link = "probit"), monotone = FALSE, ngrid = 50,
              control = boost_control(nu = .5, trace = TRUE))

l <- with(db, seq(from = min(lage), to = max(lage), length = 100))
pr <- predict(mod, newdata = data.frame(lage = l), type = "response",
              annotate = TRUE)
pr$lage <- l^3
pr$cut <- factor(pr$lage > 2.5)
levels(pr$cut) <- c("Age < 2.5 yrs", "Age > 2.5 yrs")

save(db, pr, file = "ex_heads.Rda")
