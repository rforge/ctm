
files <- list.files(pattern = "ex.*rda")

ret <- c()

pdf("summary.pdf", width = 10)
for (f in files) {

   load(f)
   boxplot(risk - ll0, main = f, outline = FALSE)
   df <- as.data.frame(risk - ll0)
   df$problem <- f
   ret <- rbind(ret, df)
}

for (f in files) {

   load(f)
   boxplot(risk[,-7] - ll0, main = f, outline = FALSE)
}

dev.off()

ret$problem <- as.factor(ret$problem)

save(ret, file = "summary.rda")

x <- ret[, c("gam_ctm", "glm_ctm", "tree_ctm", "gam_tram", "glm_tram", "tree_tram", "trtf_tree", "trtf_forest")]

tab <- (tapply(1:NROW(ret), ret$problem, function(i) {
    apply(x[i,], 2, quantile, prob = c(.5, .25, .75), na.rm = TRUE)
}))

tab <- lapply(tab, function(x) x)

tab <- tab[c("ex_beetles.rda", "ex_fetus.rda", "ex_bodyfat.rda", "ex_CAO.rda", 
             "ex_india.rda", "ex_heads.rda", "ex_ALSFRS.rda", "ex_ALSsurv.rda")]

names(tab) <- c("Beetles Exctinction Risk",
                "Birth Weight Prediction",
                "Body Fat Mass",
                "CAO/ARO/AIO-04 DFS",
                "Childhood Malnutrition",
                "Head Circumference",
                "PRO-ACT ALSFRS",
                "PRO-ACT OS")

tab <- lapply(tab, function(x) {
    i <- which(x[1,] == max(x[1,]))
    med <- formatC(round(x[1,], 2), format = "f", digits = 2)
    med[i] <- paste("\\textbf{", med[i], "}", sep = "")
    q1 <- formatC(round(x[2,], 2), format = "f", digits = 2)
    q3 <- formatC(round(x[3,], 2), format = "f", digits = 2)
    ret <- rbind(paste(paste(med, collapse = " & "), "\\\\"),
                 paste(paste(paste("(", q1, ",", q3, ")", sep = ""), 
                 collapse = " & "), "\\\\"))
    ret
})

tab <- lapply(names(tab), function(n) {
    rbind(paste(n, "& &", tab[[n]][1,]),
          paste("& &", tab[[n]][2,]))
})

dump("tab", file = "tab.R")
