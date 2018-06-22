
source("setup.R")

### Paul Eilers: Using 1/2 is better
data("db", package = "gamlss.data")
db$lage <- with(db, age^(1/3))

ldata <- db

m_mlt <- BoxCox(head ~ 1, data = ldata)
ll0 <- logLik(m_mlt) / nrow(ldata)

fm_gam <- c("ctm" = head ~ bbs(lage),
            "tram" = head ~ bols(lage, intercept = FALSE) + bbs(lage, center = TRUE, df = 1))
fm_glm <- c("ctm" = head ~ bols(lage),
            "tram" = head ~ bols(lage, intercept = FALSE))

fm_tree <- head ~ .

### no need to adapt here
fd <- cv(weights(m_mlt), type = "subsampling", B = B, prob = .75)
bctrl <- boost_control(mstop = M, trace = TRUE)

(m_glm <- FUN(m_mlt, fm_glm, ldata, control = bctrl, folds = fd))
(m_gam <- FUN(m_mlt, fm_gam, ldata, control = bctrl, folds = fd))
(m_tree <- FUN(m_mlt, fm_tree, ldata, control = bctrl, method =
              quote(mboost::blackboost), folds = fd))

tctrl <- ctree_control(saveinfo = FALSE, alpha = .01,
                       minbucket = length(coef(as.mlt(m_mlt))) * 2)
fctrl <- ctree_control(saveinfo = FALSE, alpha = 1,
                       minsplit = 50, minbucket = 25, nmax = c("yx" = Inf, "z" = 100))
r_trtf <- FUN2(m_mlt, fm_tree, ldata, tcontrol = tctrl, fcontrol = fctrl, fd)

r_glm <- m_glm$risk
r_gam <- m_gam$risk
r_tree <- m_tree$risk

colnames(r_glm) <- paste("glm", colnames(r_glm), sep = "_")
colnames(r_gam) <- paste("gam", colnames(r_gam), sep = "_")
colnames(r_tree) <- paste("tree", colnames(r_tree), sep = "_")
colnames(r_trtf) <- paste("trtf", colnames(r_trtf), sep = "_")

risk <- cbind(r_glm, r_gam, r_tree, r_trtf)

ll0 <- numeric(ncol(fd))
for (i in 1:ncol(fd)) {
    w <- fd[,i]
    ll0[i] <- logLik(update(m_mlt, theta = coef(as.mlt(m_mlt)), weights = w), w = 1 - w) / sum(1 - w)
}

save(risk, ll0, file = "ex_heads.rda")

warnings()
sessionInfo()
