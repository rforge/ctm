\name{trafotree}
\alias{trafotree}
\title{Transformation Trees}
\description{
  Partitioned transformation models
}
\usage{
trafotree(object, parm = 1:length(coef(object)), mltargs = list(maxit = 10000), ...)
}
\arguments{
  \item{object}{an object of class \code{\link[mlt]{ctm}} specifying the
                abstract model to be partitioned.}
  \item{parm}{parameters of \code{object} those corresponding score is
              used for finding partitions.}
  \item{mltargs}{arguments to \code{\link[mlt]{mlt}} for fitting the
                 transformation models.}
  \item{\dots}{arguments to \code{\link[partykit]{ctree}}, at least
               \code{formula} and \code{data}.}
}
\details{
  Conditional inference trees are used for partitioning likelihood-based transformation
  models as described in Hothorn and Zeileis (2017). The method can be seen
  in action in Hothorn (2018) and the corresponding code is available as
  \code{demo("BMI")}. \code{demo("applications")} performs transformation
  tree analyses for some standard benchmarking problems.
}
\value{
  An object of class \code{trafotree} with corresponding \code{plot}, \code{logLik} and
  \code{predict} methods.}
\references{
  Torsten Hothorn and Achim Zeileis (2017). Transformation Forests.
  \url{https://arxiv.org/abs/1701.02110}.

  Torsten Hothorn (2018). Top-Down Transformation Choice. \emph{Statistical
  Modelling}, \url{https://arxiv.org/abs/1706.08269}.
}
\examples{

### Example: Stratified Medicine Using Partitioned Cox-Models
### A combination of <DOI:10.1515/ijb-2015-0032> and <arXiv:1701.02110>
### based on infrastructure in the mlt R add-on package described in
### https://cran.r-project.org/web/packages/mlt.docreg/vignettes/mlt.pdf

library("trtf")
library("survival")
### German Breast Cancer Study Group 2 data set
data("GBSG2", package = "TH.data")

### set-up Cox model with overall treatment effect in hormonal therapy
yvar <- numeric_var("y", support = c(100, 2000), bounds = c(0, Inf))
By <- Bernstein_basis(yvar, order = 5, ui = "incre")
m <- ctm(response = By, shifting = ~ horTh, todistr = "MinExt", data = GBSG2)
GBSG2$y <- with(GBSG2, Surv(time, cens))

### overall log-hazard ratio
coef(cmod <- mlt(m, data = GBSG2))["horThyes"]
### roughly the same as 
coef(coxph(y ~ horTh, data = GBSG2))

### partition the model, ie both the baseline hazard function AND the
### treatment effect
(part_cmod <- trafotree(m, formula = y ~ horTh | age + menostat + tsize + 
    tgrade + pnodes + progrec + estrec, data = GBSG2))

### compare the log-likelihoods
logLik(cmod)
logLik(part_cmod)

### stronger effects in nodes 2 and 4 and no effect in node 5
coef(part_cmod)[, "horThyes"]

### plot the conditional survivor functions
nd <- data.frame(horTh = sort(unique(GBSG2$horTh)))
col <- c("cadetblue3", "chartreuse4")
layout(matrix(1:3, ncol = 3))
cf <- coef(part_cmod)
tmp <- sapply(rownames(cf), function(node) {
    coef(cmod) <- cf[node,]
    plot(cmod, newdata = nd, type = "survivor", col = col, 
         xlab = "Time", ylab = "Probability", lwd = 2,
         main = paste("Node", node))
    legend("topright", col = col, legend = levels(nd$horTh), 
           lty = 1, title = "Hormonal Therapy", bty = "n", lwd = 2)
})
}
\keyword{trees}
