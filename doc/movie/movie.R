
library("mlt")
library("colorspace")

yvar <- numeric_var("y", support = c(-2, 2), add = c(-1, 1))
By <- Bernstein_basis(yvar, order = 3, ui = "increasing")
yctm <- ctm(By, todistr = "Min")

ygrid <- as.data.frame(mkgrid(yctm, n = 100))
cf <- seq(from = -3, to = 3, length.out = length(coef(yctm))) 

dcf <- diff(cf[-length(cf)])
n <- 5	

ss <- lapply(dcf, function(m) {
    ret <- seq(from = m, to = 0, length.out = n)
    c(ret, rev(ret)[-c(1, length(rev))])
})

ss <- expand.grid(ss)

mm <- model.matrix(yctm, data = ygrid)

col <- rainbow_hcl(n = length(coef(yctm)))


for (i in 1:nrow(ss)) {

png(paste(sprintf("pic_%02d",i), ".png", sep = ""))
layout(matrix(1:3, nrow = 1))

    coef(yctm) <- c(as.vector(cumsum(c(cf[1], ss[i,]))), cf[length(cf)])

    yl <- range(cf)
    plot(ygrid$y, mm[,1] * coef(yctm)[1], type = "l", ylim = yl, col = col[1])
    legend("topright", col = col, lty = 1, legend = 1:length(coef(yctm)))
    for (j in 2:length(coef(yctm)))
        lines(ygrid$y, mm[,j] * coef(yctm)[j], type = "l", col = col[j])

    plot(yctm, newdata = data.frame(1), q = ygrid[[1]], type = "trafo", ylim = c(-3, 3), col = "black")
    plot(yctm, newdata = data.frame(1), q = ygrid[[1]], type = "density", col = "black", ylim = c(0, 1))
dev.off()
}


    

