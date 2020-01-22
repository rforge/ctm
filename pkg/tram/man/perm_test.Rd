\name{perm_test}
\alias{perm_test}
\alias{perm_test.tram}
\title{
  Permutation Transformation Score Tests 
}
\description{
  P-values for a parameter in a linear transformation
  model obtained from by the permutation score test principle
}
\usage{
perm_test(object, ...)
\method{perm_test}{tram}(object, parm = names(coef(object)), 
    alternative = c("two.sided", "less", "greater"), 
    block_permutation = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{\link{tram}}}
  \item{parm}{a vector of names of parameters to be tested.
              These parameters must be present in \code{object}.}
  \item{alternative}{a character string specifying the alternative hypothesis,
          must be one of \code{"two.sided"} (default), \code{"greater"} 
          or \code{"less"}.}
  \item{block_permutation}{a logical indicating wheather stratifying
      variables shall be interpreted as blocks defining admissible
      permutations.}
  \item{\dots}{additional arguments to \code{\link[coin]{independence_test}}.}
}
\details{
  A permuation score test for one single parameters in the linear
  predictor of \code{object} is computed. This parameters must be present
  in \code{object}.
}
\value{
  An object of class \code{htest} or a list thereof. See \code{\link{Coxph}}
  for an example.
}
\examples{

  ## Tritiated Water Diffusion Across Human Chorioamnion
  ## Hollander and Wolfe (1999, p. 110, Tab. 4.1)
  diffusion <- data.frame(
      pd = c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46,
             1.15, 0.88, 0.90, 0.74, 1.21),
      age = factor(rep(c("At term", "12-26 Weeks"), c(10, 5)))
  )

  ### plot the two quantile functions
  boxplot(pd ~ age, data = diffusion)

  ### the Wilcoxon rank sum test, with a confidence interval
  ### for a median shift
  wilcox.test(pd ~ age, data = diffusion, conf.int = TRUE, exact = TRUE)

  ### a corresponding transformation model with a log-odds ratio
  ### difference parameter
  md <- Colr(pd ~ age, data = diffusion)

  ### assess model fit by plotting estimated distribution fcts
  agef <- sort(unique(diffusion$age))
  col <- c("black", "darkred")
  plot(as.mlt(md), newdata = data.frame(age = agef),
       type = "distribution", col = col)
  legend("bottomright", col = col, lty = 1, legend = levels(agef), 
         bty = "n", pch = 19)
  ## compare with ECDFs: not too bad (but not good, either)
  npfit <- with(diffusion, tapply(pd, age, ecdf))
  lines(npfit[[1]], col = col[1])
  lines(npfit[[2]], col = col[2])

  ### Wald confidence interval
  confint(md)

  ### Likelihood confidence interval
  confint(profile(md))

  ### Score confidence interval
  confint(score_test(md))

  ### exact permutation score test
  (pt <- perm_test(md, distribution = "exact"))

  ### compare with probabilistic indices obtained from asht::wmwTest
  if (require("asht")) {
      (wt2 <- wmwTest(pd ~ I(relevel(age, "At term")), 
                      data = diffusion, method = "exact.ce"))
      ### as log-odds ratios
      PI(prob = wt2$conf.int)
      PI(prob = wt2$estimate)
       
  }
}
