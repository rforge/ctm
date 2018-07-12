	
library("lattice")

ctm <- list.files(pattern = "ctm.*rda")
tram <- list.files(pattern = "tram.*rda")

out <- c()

for (f in c(ctm, tram)) {

    load(f)
  
    out <- rbind(out, ret)

}

out$m <- factor(out$m)
out$PNON <- factor(out$PNON)
out$NOBS <- factor(out$NOBS)
out$model <- factor(out$model)
out$todistr <- factor(out$todistr)
out$order <- factor(out$order)

out$boost_ll[grep("rror", out$boost_ll)] <- NA
out$boost_ll <- as.double(out$boost_ll)
out$true_ll <- as.double(out$true_ll)
out$mlt_ll <- as.double(out$mlt_ll)
out$mstop <- as.double(out$mstop)

out$todistr <- factor(as.character(out$todistr), levels = c("normal", "logistic", "minextrval"), 
                      labels = c("Normal", "Logistic", "Gompertz"))

i <- with(out, interaction(model, m, todistr, order, PNON, NOBS))
cf <- tapply(out$boost_ll - out$true_ll, i, median, na.rm = TRUE)

x <- xtabs(~ model + m + todistr + order + PNON + NOBS, data = out)

tmp <- matrix(cf, nrow = nlevels(out$model))
tmp <- apply(tmp, 2, function(x) {
    i <- which.max(x)
    x <- round(x)
    x[i] <- paste("\\textbf{", x[i], "}", sep = "")
    x
})

x[] <- tmp

x <- aperm(x, c(2, 1, 3:6))

write.csv(stats:::format.ftable(ftable(x, col.vars = 1:2), quote=FALSE),
file = "stab.tex")

