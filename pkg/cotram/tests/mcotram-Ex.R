### test for mcotram
library("cotram")
library("mvtnorm")
########################################################
################ constant Lambda - data ################
########################################################
set.seed(123)
N <- 1000
J <- 2 ## define dimension

Lambda <- diag(J)
Lambda[lower.tri(Lambda)] <- (value <- 0.9)  ## choose 0 for uncorrelated responses
Sigma <- tcrossprod(solve(Lambda))

p <- 2
x <- matrix(runif(N * p), ncol = p)

beta <- c(1, -1, -.5, .5) ## depends on dimension J and number of predictors p

ly <- x %*% matrix(beta, nrow = p) + rmvnorm(N, sigma = Sigma)
ly_marg <- pnorm(ly)
y <- qbinom(ly_marg, size = 10, prob = 0.3) 

## or take different distributions for the two margins
# y <- ly
# y[,1] <- qbinom(ly_marg[,1], size = 15, prob = 0.2) 
# y[,2] <- qbinom(ly_marg[,2], size = 10, prob = 0.3) 

d <- data.frame(y = y, x = x)

##################################################
###### constant lambda, marginal models ~ 1 ######
##################################################
## here we don't expect a great performance as the marginal models are
## very simple

## marginal cotram models with log_first = TRUE (default)
u1 <- cotram(y.1 ~ 1, data = d, method = "probit")
u2 <- cotram(y.2 ~ 1, data = d, method = "probit")

## joint models with different orders of the marginals, constant lambdas
uc1 <- mcotram(u1, u2, data = d)
uc2 <- mcotram(u2, u1, data = d)
## these log-likelihoods are expected to be very close, but not equal
logLik(uc1)
logLik(uc2)

## check gradient
all.equal(uc1$sc(uc2$par), numDeriv::grad(uc1$ll, uc2$par),
          tol = 1e-6, check.attributes = FALSE)

## marginal cotram models with log_first = FALSE
u1l <- cotram(y.1 ~ 1, data = d, method = "probit", log_first = FALSE)
u2l <- cotram(y.2 ~ 1, data = d, method = "probit", log_first = FALSE)
uc1l <- mcotram(u1l, u2l, data = d)
uc2l <- mcotram(u2l, u1l, data = d)
logLik(uc1l)
logLik(uc2l)
all.equal(uc1l$sc(uc2l$par), numDeriv::grad(uc1l$ll, uc2l$par),
          tol = 1e-6, check.attributes = FALSE)

## compare coefficients
## log_first = TRUE
# coef(uc1)
# coef(uc2)
value ## recall original value for Lambda
coef(uc1, type = "Lambda")[1,] 
coef(uc2, type = "Lambda")[1,]
cov2cor(Sigma)[1, 2] ## original value for correlation
coef(uc1, type = "Corr")[1,]
coef(uc2, type = "Corr")[1,]

## log_first = FALSE
# coef(uc1l)
# coef(uc2l)
coef(uc1l, type = "Lambda")[1,]
coef(uc2l, type = "Lambda")[1,]

coef(uc1l, type = "Corr")[1,]
coef(uc2l, type = "Corr")[1,]

## transformation theory check
## for constant lambdas, the transformed coefficients should agree with the
## coefficients of the model with the different order

## here, the results are not great because the marginal models are too simple

## log_first = TRUE
coef(uc1)[-c(1:7, 15)] / sqrt(coef(uc1, type = "Sigma")$diagonal[1,2])
coef(uc2)[1:7]

coef(uc2)[-c(1:7, 15)] / sqrt(coef(uc2, type = "Sigma")$diagonal[1,2])
coef(uc1)[1:7]

## log_first = FALSE?
coef(uc1l)[-c(1:7, 15)] / sqrt(coef(uc1l, type = "Sigma")$diagonal[1,2])
coef(uc2l)[1:7]

coef(uc2l)[-c(1:7, 15)] / sqrt(coef(uc2l, type = "Sigma")$diagonal[1,2])
coef(uc1l)[1:7]

### predict accuracy
cbind(u1 = predict(u1, newdata = d[1:3, ], type = "distribution"),
      uc1 = predict(uc1, marginal = 1, newdata = d[1:3, ], type = "distribution"),
      uc2 = predict(uc2, marginal = 2, newdata = d[1:3, ], type = "distribution"))
cbind(u2 = predict(u2, newdata = d[1:3, ], type = "distribution"),
      uc1 = predict(uc1, marginal = 2, newdata = d[1:3, ], type = "distribution"),
      uc2 = predict(uc2, marginal = 1, newdata = d[1:3, ], type = "distribution"))
## log_first = TRUE seems to perform slightly better

cbind(u1l = predict(u1l, newdata = d[1:3, ], type = "distribution"),
      uc1l = predict(uc1l, marginal = 1, newdata = d[1:3, ], type = "distribution"),
      uc2l = predict(uc2l, marginal = 2, newdata = d[1:3, ], type = "distribution"))
cbind(u2l = predict(u2l, newdata = d[1:3, ], type = "distribution"),
      uc1l = predict(uc1l, marginal = 2, newdata = d[1:3, ], type = "distribution"),
      uc2l = predict(uc2l, marginal = 1, newdata = d[1:3, ], type = "distribution"))

## QQ-plot to check that marginal transformations are standard normal:
## plots don't look great since the marginal models are too simple
nd <- d
par(mfrow = c(2, 2))
d1 <- predict(uc1, marginal = 1, newdata = d, type = "trafo")
qqnorm(d1, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 12, marginal = 1")
qqline(d1, col = "red")

d2 <- predict(uc1, marginal = 2, newdata = d, type = "trafo")
qqnorm(d2, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 12, marginal = 2")
qqline(d2, col = "red")

d3 <- predict(uc2, marginal = 2, newdata = d, type = "trafo")
qqnorm(d3, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 21, marginal = 2")
qqline(d3, col = "red")

d4 <- predict(uc2, marginal = 1, newdata = d, type = "trafo")
qqnorm(d4, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 21, marginal = 2")
qqline(d4, col = "red")

## log_first = FALSE
d1 <- predict(uc1l, marginal = 1, newdata = d, type = "trafo")
qqnorm(d1, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 12, marginal = 1")
qqline(d1, col = "red")

d2 <- predict(uc1l, marginal = 2, newdata = d, type = "trafo")
qqnorm(d2, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 12, marginal = 2")
qqline(d2, col = "red")

d3 <- predict(uc2l, marginal = 2, newdata = d, type = "trafo")
qqnorm(d3, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 21, marginal = 2")
qqline(d3, col = "red")

d4 <- predict(uc2l, marginal = 1, newdata = d, type = "trafo")
qqnorm(d4, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 21, marginal = 2")
qqline(d4, col = "red")

#### continuous approximation
## data - 0.5 (and + 1 depending on log_first)
yt <- y - 0.5*(y > 0)
ytt <- yt + 1
d <- data.frame(y = y, yt = yt, ytt = ytt, x = x)

## univariate BoxCox models, log_first = TRUE
b1 <- BoxCox(ytt.1 ~ 1, data = d, log_first = TRUE)
b2 <- BoxCox(ytt.2 ~ 1, data = d, log_first = TRUE)

b1l <- BoxCox(yt.1 ~ 1, data = d, log_first = FALSE)
b2l <- BoxCox(yt.2 ~ 1, data = d, log_first = FALSE)

## multivariate mmlt models
mb1 <- mmlt(b1, b2, data = d)
mb2 <- mmlt(b2, b1, data = d)
logLik(mb1)
logLik(mb2)

mb1l <- mmlt(b1l, b2l, data = d)
mb2l <- mmlt(b2l, b1l, data = d)
logLik(mb1l)
logLik(mb2l)

## compare coef
# coef(mb1)
# coef(mb2)
coef(mb1, type = "Lambda")[1,]
coef(mb2, type = "Lambda")[1,]
coef(mb1, type = "Corr")[1,]
coef(mb2, type = "Corr")[1,]

## log_first = FALSE
# coef(mb1l)
# coef(mb2l)
coef(mb1l, type = "Lambda")[1,]
coef(mb2l, type = "Lambda")[1,]
coef(mb1l, type = "Corr")[1,]
coef(mb2l, type = "Corr")[1,]

## transformation theory check
## for constant lambdas, the transformed coefficients should agree with the
## coefficients of the model with the different order

## log_first = TRUE
coef(mb1)[-c(1:7, 15)] / sqrt(coef(mb1, type = "Sigma")$diagonal[1,2])
coef(mb2)[1:7]

coef(mb2)[-c(1:7, 15)] / sqrt(coef(mb2, type = "Sigma")$diagonal[1,2])
coef(mb1)[1:7]

## log_first = FALSE
## results not great but probably because the marginal models are too simple
coef(mb1l)[-c(1:7, 15)] / sqrt(coef(mb1l, type = "Sigma")$diagonal[1,2])
coef(mb2l)[1:7]

coef(mb2l)[-c(1:7, 15)] / sqrt(coef(mb2l, type = "Sigma")$diagonal[1,2])
coef(mb1l)[1:7]

### predict accuracy
cbind(b1 = predict(b1, newdata = d[1:3, ], type = "distribution"),
      mb1 = predict(mb1, marginal = 1, newdata = d[1:3, ], type = "distribution"),
      mb2 = predict(mb2, marginal = 2, newdata = d[1:3, ], type = "distribution"))
cbind(b2 = predict(b2, newdata = d[1:3, ], type = "distribution"),
      mb1 = predict(mb1, marginal = 2, newdata = d[1:3, ], type = "distribution"),
      mb2 = predict(mb2, marginal = 1, newdata = d[1:3, ], type = "distribution"))

cbind(b1l = predict(b1l, newdata = d[1:3, ], type = "distribution"),
      mb1l = predict(mb1l, marginal = 1, newdata = d[1:3, ], type = "distribution"),
      mb2l = predict(mb2l, marginal = 2, newdata = d[1:3, ], type = "distribution"))
cbind(b2l = predict(b2l, newdata = d[1:3, ], type = "distribution"),
      mb1l = predict(mb1l, marginal = 2, newdata = d[1:3, ], type = "distribution"),
      mb2l = predict(mb2l, marginal = 1, newdata = d[1:3, ], type = "distribution"))

## QQ-plot to check that marginal transformations are standard normal:
nd <- d
par(mfrow = c(2, 2))
d1 <- predict(mb1, marginal = 1, newdata = d, type = "trafo")
qqnorm(d1, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 12, marginal = 1")
qqline(d1, col = "red")

d2 <- predict(mb1, marginal = 2, newdata = d, type = "trafo")
qqnorm(d2, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 12, marginal = 2")
qqline(d2, col = "red")

d3 <- predict(mb2, marginal = 2, newdata = d, type = "trafo")
qqnorm(d3, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 21, marginal = 2")
qqline(d3, col = "red")

d4 <- predict(mb2, marginal = 1, newdata = d, type = "trafo")
qqnorm(d4, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 21, marginal = 2")
qqline(d4, col = "red")

## log_first = FALSE
d1 <- predict(mb1l, marginal = 1, newdata = d, type = "trafo")
qqnorm(d1, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 12, marginal = 1")
qqline(d1, col = "red")

d2 <- predict(mb1l, marginal = 2, newdata = d, type = "trafo")
qqnorm(d2, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 12, marginal = 2")
qqline(d2, col = "red")

d3 <- predict(mb2l, marginal = 2, newdata = d, type = "trafo")
qqnorm(d3, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 21, marginal = 2")
qqline(d3, col = "red")

d4 <- predict(mb2l, marginal = 1, newdata = d, type = "trafo")
qqnorm(d4, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 21, marginal = 2")
qqline(d4, col = "red")

##########################################################
###### constant lambda, marginal models ~ x.1 + x.2 ######
##########################################################
## cotram models for the marginals
u1 <- cotram(y.1 ~ x.1 + x.2, data = d, method = "probit")
u2 <- cotram(y.2 ~ x.1 + x.2, data = d, method = "probit")

## joint models with different orders of the marginals, constant lambdas
uc1 <- mcotram(u1, u2, data = d)
uc2 <- mcotram(u2, u1, data = d)
logLik(uc1)
logLik(uc2) ## these log-likelihoods are expected to be very close, but not equal

## check gradient
all.equal(uc1$sc(uc2$par), numDeriv::grad(uc1$ll, uc2$par),
          tol = 1e-6, check.attributes = FALSE)

## marginal cotram models with log_first = FALSE
# uc1 <- uc2 <- 0
u1l <- cotram(y.1 ~ x.1 + x.2, data = d, method = "probit", log_first = FALSE)
u2l <- cotram(y.2 ~ x.1 + x.2, data = d, method = "probit", log_first = FALSE)
uc1l <- mcotram(u1l, u2l, data = d)
uc2l <- mcotram(u2l, u1l, data = d)
logLik(uc1l)
logLik(uc2l)
all.equal(uc1l$sc(uc2l$par), numDeriv::grad(uc1l$ll, uc2l$par),
          tol = 1e-6, check.attributes = FALSE)

## compare coefficients
## recall beta:
beta
## log_first = TRUE
coef(uc1)
coef(uc2)
value
coef(uc1, type = "Lambda")[1,] 
coef(uc2, type = "Lambda")[1,]
cov2cor(Sigma)[1, 2]
coef(uc1, type = "Corr")[1,]
coef(uc2, type = "Corr")[1,]

## log_first = FALSE
beta
coef(uc1l)
coef(uc2l)
coef(uc1l, type = "Lambda")[1,]
coef(uc2l, type = "Lambda")[1,]
coef(uc1l, type = "Corr")[1,]
coef(uc2l, type = "Corr")[1,]

## transformation theory check
## for constant lambdas, the transformed coefficients should agree with the
## coefficients of the model with the different order

## log_first = TRUE
## seems to agree better than log_first = FALSE
coef(uc1)[-c(1:9, 19)] / sqrt(coef(uc1, type = "Sigma")$diagonal[1,2])
coef(uc2)[1:9]

coef(uc2)[-c(1:9, 19)] / sqrt(coef(uc2, type = "Sigma")$diagonal[1,2])
coef(uc1)[1:9]

## log_first = FALSE
coef(uc1l)[-c(1:9, 19)] / sqrt(coef(uc1l, type = "Sigma")$diagonal[1,2])
coef(uc2l)[1:9]

coef(uc2l)[-c(1:9, 19)] / sqrt(coef(uc2l, type = "Sigma")$diagonal[1,2])
coef(uc1l)[1:9]

### predict accuracy
cbind(u1 = predict(u1, newdata = d[1:3, ], type = "distribution"),
      uc1 = predict(uc1, marginal = 1, newdata = d[1:3, ], type = "distribution"),
      uc2 = predict(uc2, marginal = 2, newdata = d[1:3, ], type = "distribution"))
cbind(u2 = predict(u2, newdata = d[1:3, ], type = "distribution"),
      uc1 = predict(uc1, marginal = 2, newdata = d[1:3, ], type = "distribution"),
      uc2 = predict(uc2, marginal = 1, newdata = d[1:3, ], type = "distribution"))

cbind(u1l = predict(u1l, newdata = d[1:3, ], type = "distribution"),
      uc1l = predict(uc1l, marginal = 1, newdata = d[1:3, ], type = "distribution"),
      uc2l = predict(uc2l, marginal = 2, newdata = d[1:3, ], type = "distribution"))
cbind(u2l = predict(u2l, newdata = d[1:3, ], type = "distribution"),
      uc1l = predict(uc1l, marginal = 2, newdata = d[1:3, ], type = "distribution"),
      uc2l = predict(uc2l, marginal = 1, newdata = d[1:3, ], type = "distribution"))

## QQ-plot to check that marginal transformations are standard normal:
nd <- d
par(mfrow = c(2, 2))
d1 <- predict(uc1, marginal = 1, newdata = d, type = "trafo")
qqnorm(d1, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 12, marginal = 1")
qqline(d1, col = "red")

d2 <- predict(uc1, marginal = 2, newdata = d, type = "trafo")
qqnorm(d2, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 12, marginal = 2")
qqline(d2, col = "red")

d3 <- predict(uc2, marginal = 2, newdata = d, type = "trafo")
qqnorm(d3, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 21, marginal = 2")
qqline(d3, col = "red")

d4 <- predict(uc2, marginal = 1, newdata = d, type = "trafo")
qqnorm(d4, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 21, marginal = 2")
qqline(d4, col = "red")

## log_first = FALSE
d1 <- predict(uc1l, marginal = 1, newdata = d, type = "trafo")
qqnorm(d1, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 12, marginal = 1")
qqline(d1, col = "red")

d2 <- predict(uc1l, marginal = 2, newdata = d, type = "trafo")
qqnorm(d2, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 12, marginal = 2")
qqline(d2, col = "red")

d3 <- predict(uc2l, marginal = 2, newdata = d, type = "trafo")
qqnorm(d3, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 21, marginal = 2")
qqline(d3, col = "red")

d4 <- predict(uc2l, marginal = 1, newdata = d, type = "trafo")
qqnorm(d4, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 21, marginal = 2")
qqline(d4, col = "red")

#### continuous approximation
## univariate BoxCox models, log_first = TRUE
b1 <- BoxCox(ytt.1 ~ x.1 + x.2, data = d, log_first = TRUE)
b2 <- BoxCox(ytt.2 ~ x.1 + x.2, data = d, log_first = TRUE)

b1l <- BoxCox(yt.1 ~ x.1 + x.2, data = d, log_first = FALSE)
b2l <- BoxCox(yt.2 ~ x.1 + x.2, data = d, log_first = FALSE)

## multivariate mmlt models
mb1 <- mmlt(b1, b2, data = d)
mb2 <- mmlt(b2, b1, data = d)
logLik(mb1)
logLik(mb2)

mb1l <- mmlt(b1l, b2l, data = d)
mb2l <- mmlt(b2l, b1l, data = d)
logLik(mb1l)
logLik(mb2l)

## compare coef
## log_first = TRUE
beta
coef(mb1)
coef(mb2)
value
coef(mb1, type = "Lambda")[1,]
coef(mb2, type = "Lambda")[1,]
cov2cor(Sigma)[1, 2]
coef(mb1, type = "Corr")[1,]
coef(mb2, type = "Corr")[1,]

## log_first = FALSE
beta
coef(mb1l)
coef(mb2l)
coef(mb1l, type = "Lambda")[1,]
coef(mb2l, type = "Lambda")[1,]
coef(mb1l, type = "Corr")[1,]
coef(mb2l, type = "Corr")[1,]

## transformation theory check
## for constant lambdas, the transformed coefficients should agree with the
## coefficients of the model with the different order
coef(mb1)[-c(1:9, 19)] / sqrt(coef(mb1, type = "Sigma")$diagonal[1,2])
coef(mb2)[1:9]

coef(mb2)[-c(1:9, 19)] / sqrt(coef(mb2, type = "Sigma")$diagonal[1,2])
coef(mb1)[1:9]

## doesn't seem to be the case for log_first = FALSE?
coef(mb1l)[-c(1:9, 19)] / sqrt(coef(mb1l, type = "Sigma")$diagonal[1,2])
coef(mb2l)[1:9]

coef(mb2l)[-c(1:9, 19)] / sqrt(coef(mb2l, type = "Sigma")$diagonal[1,2])
coef(mb1l)[1:9]

### predict accuracy
cbind(b1 = predict(b1, newdata = d[1:3, ], type = "distribution"),
      mb1 = predict(mb1, marginal = 1, newdata = d[1:3, ], type = "distribution"),
      mb2 = predict(mb2, marginal = 2, newdata = d[1:3, ], type = "distribution"))
cbind(b2 = predict(b2, newdata = d[1:3, ], type = "distribution"),
      mb1 = predict(mb1, marginal = 2, newdata = d[1:3, ], type = "distribution"),
      mb2 = predict(mb2, marginal = 1, newdata = d[1:3, ], type = "distribution"))

cbind(b1l = predict(b1l, newdata = d[1:3, ], type = "distribution"),
      mb1l = predict(mb1l, marginal = 1, newdata = d[1:3, ], type = "distribution"),
      mb2l = predict(mb2l, marginal = 2, newdata = d[1:3, ], type = "distribution"))
cbind(b2l = predict(b2l, newdata = d[1:3, ], type = "distribution"),
      mb1l = predict(mb1l, marginal = 2, newdata = d[1:3, ], type = "distribution"),
      mb2l = predict(mb2l, marginal = 1, newdata = d[1:3, ], type = "distribution"))

## QQ-plot to check that marginal transformations are standard normal:
nd <- d
par(mfrow = c(2, 2))
d1 <- predict(mb1, marginal = 1, newdata = d, type = "trafo")
qqnorm(d1, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 12, marginal = 1")
qqline(d1, col = "red")

d2 <- predict(mb1, marginal = 2, newdata = d, type = "trafo")
qqnorm(d2, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 12, marginal = 2")
qqline(d2, col = "red")

d3 <- predict(mb2, marginal = 2, newdata = d, type = "trafo")
qqnorm(d3, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 21, marginal = 2")
qqline(d3, col = "red")

d4 <- predict(mb2, marginal = 1, newdata = d, type = "trafo")
qqnorm(d4, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 21, marginal = 2")
qqline(d4, col = "red")

## log_first = FALSE
d1 <- predict(mb1l, marginal = 1, newdata = d, type = "trafo")
qqnorm(d1, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 12, marginal = 1")
qqline(d1, col = "red")

d2 <- predict(mb1l, marginal = 2, newdata = d, type = "trafo")
qqnorm(d2, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 12, marginal = 2")
qqline(d2, col = "red")

d3 <- predict(mb2l, marginal = 2, newdata = d, type = "trafo")
qqnorm(d3, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 21, marginal = 2")
qqline(d3, col = "red")

d4 <- predict(mb2l, marginal = 1, newdata = d, type = "trafo")
qqnorm(d4, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 21, marginal = 2")
qqline(d4, col = "red")

#############################################################
###### lambda = lambda(x), marginal models ~ x.1 + x.2 ######
#############################################################

## data
## need to simulate data with lambdas that depend on x.1, x.2
set.seed(2701)
N <- 1000
J <- 2
x <- matrix(runif(N * J), ncol = J)

beta <- c(1, -1, -.5, .5)
# beta <- rep(0, 4)

xb <- x %*% matrix(beta, nrow = 2)
ly <- xb

lb <- (off_l <- 0.4) + x %*% (cf_l <- c(1.3, 0.6))
# plot(lb)
cr <- numeric(N)

for (i in 1:N) {
   L_i <- diag(J)
   L_i[2, 1] <- lb[i]
   S_i <- tcrossprod(solve(L_i))
   ly[i, ] <- ly[i, ] + rmvnorm(1, sigma = S_i)
   cr[i] <- cov2cor(S_i)[1, 2]
}

ly_marg <- pnorm(ly)
y <- qbinom(ly_marg, size = 10, prob = 0.3) 

## or take different distributions for the two margins
# y <- ly
# y[,1] <- qbinom(ly_marg[,1], size = 15, prob = 0.2) 
# y[,2] <- qbinom(ly_marg[,2], size = 10, prob = 0.3) 

# d <- data.frame(y = y, x = x)
yt <- y - 0.5*(y > 0)
ytt <- yt + 1
d <- data.frame(y = y, yt = yt, ytt = ytt, x = x)
# hist(d$y.1, breaks = 50)
# hist(d$y.2, breaks = 50)
# sort(unique(y)) ## we only have counts
# plot(y)

##################################
## cotram models for the marginals
u1 <- cotram(y.1 ~ x.1 + x.2, data = d, method = "probit")
u2 <- cotram(y.2 ~ x.1 + x.2, data = d, method = "probit")

## joint models with different orders of the marginals, constant lambdas
uc1 <- mcotram(u1, u2, data = d, formula = ~ x.1 + x.2)
uc2 <- mcotram(u2, u1, data = d, formula = ~ x.1 + x.2)
logLik(uc1)
logLik(uc2) ## these log-likelihoods will not be the same

## check gradient
all.equal(uc1$sc(uc2$par), numDeriv::grad(uc1$ll, uc2$par),
          tol = 1e-6, check.attributes = FALSE)

## marginal cotram models with log_first = FALSE
# uc1 <- uc2 <- 0
u1l <- cotram(y.1 ~ x.1 + x.2, data = d, method = "probit", log_first = FALSE)
u2l <- cotram(y.2 ~ x.1 + x.2, data = d, method = "probit", log_first = FALSE)
uc1l <- mcotram(u1l, u2l, data = d, formula = ~ x.1 + x.2)
uc2l <- mcotram(u2l, u1l, data = d, formula = ~ x.1 + x.2)
logLik(uc1l)
logLik(uc2l)
all.equal(uc1l$sc(uc2l$par), numDeriv::grad(uc1l$ll, uc2l$par),
          tol = 1e-6, check.attributes = FALSE)

## compare coefficients
## recall true coefs:
beta
off_l
cf_l

R <- numeric(3)
for (i in 1:3) {
   L <- diag(2)
   L[2,1] <- lb[i]
   S <- tcrossprod(solve(L))
   R[i] <- cov2cor(S)[2, 1]
}

## log_first = TRUE
coef(uc1)
coef(uc2)
cbind(lb[1:3], coef(uc1, type = "Lambda")[1:3,], coef(uc2, type = "Lambda")[1:3,])
cbind(orig = R, uc1 = coef(uc1, type = "Corr")[1:3,], 
      uc2 = coef(uc2, type = "Corr")[1:3,])

## log_first = FALSE
coef(uc1l)
coef(uc2l)
cbind(orig = lb[1:3], uc1l = coef(uc1l, type = "Lambda")[1:3,], 
      uc2l = coef(uc2l, type = "Lambda")[1:3,])
cbind(orig = R, uc1 = coef(uc1l, type = "Corr")[1:3,], 
      uc2 = coef(uc2l, type = "Corr")[1:3,])

### predict accuracy
cbind(u1 = predict(u1, newdata = d[1:3, ], type = "distribution"),
      uc1 = predict(uc1, marginal = 1, newdata = d[1:3, ], type = "distribution"),
      uc2 = predict(uc2, marginal = 2, newdata = d[1:3, ], type = "distribution"))
cbind(u2 = predict(u2, newdata = d[1:3, ], type = "distribution"),
      uc1 = predict(uc1, marginal = 2, newdata = d[1:3, ], type = "distribution"),
      uc2 = predict(uc2, marginal = 1, newdata = d[1:3, ], type = "distribution"))

cbind(u1l = predict(u1l, newdata = d[1:3, ], type = "distribution"),
      uc1l = predict(uc1l, marginal = 1, newdata = d[1:3, ], type = "distribution"),
      uc2l = predict(uc2l, marginal = 2, newdata = d[1:3, ], type = "distribution"))
cbind(u2l = predict(u2l, newdata = d[1:3, ], type = "distribution"),
      uc1l = predict(uc1l, marginal = 2, newdata = d[1:3, ], type = "distribution"),
      uc2l = predict(uc2l, marginal = 1, newdata = d[1:3, ], type = "distribution"))

## QQ-plot to check that marginal transformations are standard normal:
nd <- d
par(mfrow = c(2, 2))
d1 <- predict(uc1, marginal = 1, newdata = d, type = "trafo")
qqnorm(d1, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 12, marginal = 1")
qqline(d1, col = "red")

d2 <- predict(uc1, marginal = 2, newdata = d, type = "trafo")
qqnorm(d2, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 12, marginal = 2")
qqline(d2, col = "red")

d3 <- predict(uc2, marginal = 2, newdata = d, type = "trafo")
qqnorm(d3, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 21, marginal = 2")
qqline(d3, col = "red")

d4 <- predict(uc2, marginal = 1, newdata = d, type = "trafo")
qqnorm(d4, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 21, marginal = 2")
qqline(d4, col = "red")

## log_first = FALSE
d1 <- predict(uc1l, marginal = 1, newdata = d, type = "trafo")
qqnorm(d1, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 12, marginal = 1")
qqline(d1, col = "red")

d2 <- predict(uc1l, marginal = 2, newdata = d, type = "trafo")
qqnorm(d2, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 12, marginal = 2")
qqline(d2, col = "red")

d3 <- predict(uc2l, marginal = 2, newdata = d, type = "trafo")
qqnorm(d3, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 21, marginal = 2")
qqline(d3, col = "red")

d4 <- predict(uc2l, marginal = 1, newdata = d, type = "trafo")
qqnorm(d4, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 21, marginal = 2")
qqline(d4, col = "red")

#### continuous approximation
## univariate BoxCox models, log_first = TRUE
b1 <- BoxCox(ytt.1 ~ x.1 + x.2, data = d, log_first = TRUE)
b2 <- BoxCox(ytt.2 ~ x.1 + x.2, data = d, log_first = TRUE)

b1l <- BoxCox(yt.1 ~ x.1 + x.2, data = d, log_first = FALSE)
b2l <- BoxCox(yt.2 ~ x.1 + x.2, data = d, log_first = FALSE)

## multivariate mmlt models
mb1 <- mmlt(b1, b2, data = d, formula = ~ x.1 + x.2)
mb2 <- mmlt(b2, b1, data = d, formula = ~ x.1 + x.2)
logLik(mb1)
logLik(mb2) ## these log-likelihoods are *not* expected to be the same

mb1l <- mmlt(b1l, b2l, data = d, formula = ~ x.1 + x.2)
mb2l <- mmlt(b2l, b1l, data = d, formula = ~ x.1 + x.2)
logLik(mb1l)
logLik(mb2l)

## compare coefficients
## recall true coefs:
beta
off_l
cf_l

## log_first = TRUE
coef(mb1)
coef(mb2)
cbind(lb[1:3], coef(mb1, type = "Lambda")[1:3,], coef(mb2, type = "Lambda")[1:3,])
cbind(orig = R, mb1 = coef(mb1, type = "Corr")[1:3,], 
      mb2 = coef(mb2, type = "Corr")[1:3,])

## log_first = FALSE
coef(mb1l)
coef(mb2l)
cbind(lb[1:3], coef(mb1l, type = "Lambda")[1:3,], coef(mb2l, type = "Lambda")[1:3,])
cbind(orig = R, mb1l = coef(mb1l, type = "Corr")[1:3,], 
      mb2l = coef(mb2l, type = "Corr")[1:3,])

### predict accuracy
cbind(b1 = predict(b1, newdata = d[1:3, ], type = "distribution"),
      mb1 = predict(mb1, marginal = 1, newdata = d[1:3, ], type = "distribution"),
      mb2 = predict(mb2, marginal = 2, newdata = d[1:3, ], type = "distribution"))
cbind(b2 = predict(b2, newdata = d[1:3, ], type = "distribution"),
      mb1 = predict(mb1, marginal = 2, newdata = d[1:3, ], type = "distribution"),
      mb2 = predict(mb2, marginal = 1, newdata = d[1:3, ], type = "distribution"))

cbind(b1l = predict(b1l, newdata = d[1:3, ], type = "distribution"),
      mb1l = predict(mb1l, marginal = 1, newdata = d[1:3, ], type = "distribution"),
      mb2l = predict(mb2l, marginal = 2, newdata = d[1:3, ], type = "distribution"))
cbind(b2l = predict(b2l, newdata = d[1:3, ], type = "distribution"),
      mb1l = predict(mb1l, marginal = 2, newdata = d[1:3, ], type = "distribution"),
      mb2l = predict(mb2l, marginal = 1, newdata = d[1:3, ], type = "distribution"))

## QQ-plot to check that marginal transformations are standard normal:
nd <- d
par(mfrow = c(2, 2))
d1 <- predict(mb1, marginal = 1, newdata = d, type = "trafo")
qqnorm(d1, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 12, marginal = 1")
qqline(d1, col = "red")

d2 <- predict(mb1, marginal = 2, newdata = d, type = "trafo")
qqnorm(d2, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 12, marginal = 2")
qqline(d2, col = "red")

d3 <- predict(mb2, marginal = 2, newdata = d, type = "trafo")
qqnorm(d3, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 21, marginal = 2")
qqline(d3, col = "red")

d4 <- predict(mb2, marginal = 1, newdata = d, type = "trafo")
qqnorm(d4, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 21, marginal = 2")
qqline(d4, col = "red")

## log_first = FALSE
d1 <- predict(mb1l, marginal = 1, newdata = d, type = "trafo")
qqnorm(d1, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 12, marginal = 1")
qqline(d1, col = "red")

d2 <- predict(mb1l, marginal = 2, newdata = d, type = "trafo")
qqnorm(d2, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 12, marginal = 2")
qqline(d2, col = "red")

d3 <- predict(mb2l, marginal = 2, newdata = d, type = "trafo")
qqnorm(d3, pch = 19, col = rgb(.1, .1, .1, .1), main = "order = 21, marginal = 2")
qqline(d3, col = "red")

d4 <- predict(mb2l, marginal = 1, newdata = d, type = "trafo")
qqnorm(d4, pch = 19, col = rgb(.1, .1, .1, .1), main = "oder = 21, marginal = 2")
qqline(d4, col = "red")




## add tests for diag = TRUE in future?
## add J = 3?
