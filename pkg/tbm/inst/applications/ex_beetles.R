
source("setup.R")

dir <- system.file("rda", package = "TH.data")
load(file.path(dir, "beetles.rda"))

ldata <- bmodel
ldata$TS <- NULL

m_mlt <- Polr(RL ~ 1, data = ldata)
ll0 <- logLik(m_mlt) / nrow(ldata)

X <- model.matrix(~ species - 1, data = ldata)
ldata$species <- NULL

nvars <- c("mean_body_size", "mean_elev", "flowers",
          "niche_decay", "niche_diam", "niche_canopy", "distribution")
fvars <- colnames(ldata)[!colnames(ldata) %in% c("RL", nvars)]
xvars <- c(nvars, fvars)

ldata[nvars] <- lapply(nvars, function(x) scale(ldata[[x]]))

fm_gam <- c(
"ctm" = as.formula(paste("RL ~ buser(X, K = Ainv, df = 2) + ",
    paste("bols(", xvars, ", df = 2)", collapse = "+"), "+",
    paste("bbs(", nvars, ", center = TRUE, df = 2)", collapse = "+"))),
"tram" = as.formula(paste("RL ~ buser(X, K = Ainv, df = 1) + ",
    paste("bols(", xvars, ", intercept = FALSE, df = 1)", collapse = "+"), "+",
    paste("bbs(", nvars, ", center = TRUE, df = 1)", collapse = "+"))))

fm_glm <- c(
"ctm" = as.formula(paste("RL ~ buser(X, K = Ainv, df = 2) +",
    paste("bols(", xvars, ", df = 2)", collapse = "+"))),
"tram" = as.formula(paste("RL ~ buser(X, K = Ainv, df = 1) +",
    paste("bols(", xvars, ", intercept = FALSE, df = 1)", collapse = "+"))))

fm_tree <- as.formula(paste("RL ~ ", paste(xvars, collapse = "+")))

fm_tree <- RL ~ . 

### no need to adapt here
fd <- cv(weights(m_mlt), type = "subsampling", B = B, prob = .75, strata = ldata$RL)
bctrl <- boost_control(mstop = M, trace = TRUE)

(m_glm <- FUN(m_mlt, fm_glm, ldata, control = bctrl, folds = fd))
(m_gam <- FUN(m_mlt, fm_gam, ldata, control = bctrl, folds = fd))
(m_tree <- FUN(m_mlt, fm_tree, ldata, control = bctrl, method =
               quote(mboost::blackboost), folds = fd))

tctrl <- ctree_control(saveinfo = FALSE, alpha = .05,
                       minbucket = length(coef(as.mlt(m_mlt))) * 2)
fctrl <- ctree_control(saveinfo = FALSE, alpha = 1,
                       minsplit = 5, minbucket = 2, nmax = c("yx" = Inf, "z" = 100))
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

save(risk, ll0, file = "ex_beetles.rda")

warnings()
sessionInfo()
