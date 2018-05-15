\name{stmboost}
\alias{stmboost}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stmboost(model, formula, data = list(), weights = NULL, method = quote(mboost::mboost), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (model, formula, data = list(), weights = NULL, method = quote(mboost::mboost), 
    ...) 
{
    basedata <- model$data
    if (inherits(model, "tram")) 
        model <- as.mlt(model)
    mf <- match.call(expand.dots = TRUE)
    mf$model <- mf$method <- NULL
    if (missing(data)) 
        data <- environment(formula)
    tmp <- model$model$bases
    if (!is.null(tmp$shifting)) {
        tmp$shifting <- c(int = intercept_basis(), shift = tmp$shifting)
    }
    else {
        tmp$shifting <- intercept_basis()
    }
    td <- model$model$todistr$name
    td <- switch(td, normal = "Normal", logistic = "Logistic", 
        `minimum extreme value` = "MinExtrVal")
    myctm <- ctm(tmp$response, tmp$interacting, tmp$shifting, 
        todistr = td)
    mf$family <- .stmFamily(myctm, basedata, weights)
    mf[[1L]] <- method
    ret <- eval(mf, parent.frame())
    ret$model <- mlt(myctm, data = basedata, dofit = FALSE)
    class(ret) <- c("stmboost", "tbm", class(ret))
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
