
cotram <- function(formula, data, method = c("logit", "cloglog", "loglog", "probit"),
                   log_first = TRUE, plus_one = log_first, prob = 0.9,
                   subset, weights, offset, cluster, na.action = na.omit, ...)
{
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(tram_data)
    td <- eval(mf, parent.frame())
    
    method <- match.arg(method)
    distribution <- c("logit" = "Logistic", "probit" = "Normal", 
                      "loglog" = "MaxExtrVal", "cloglog" = "MinExtrVal")
    distribution <- distribution[method]
    name <- c("logit" = "Odds", "loglog" = "Reverse Time Hazards",
              "cloglog" = "Hazards Cox")
    
    stopifnot(inherits(td$response, "response") || is.numeric(td$response))
    
    ## check whether response is positive integer
    if (any(td$response < 0))
        stop("response is not a positive number")
    if(!all(td$response %% 1 == 0))
        stop("response is not an integer number")
    
    ## as.integer for correct likelihood
    td$response <- as.integer(td$response)
    td$mf[,td$rname] <- as.integer(td$mf[,td$rname])
    
    ## y + 1 for log_first
    td$response <- td$response + as.integer(plus_one)
    td$mf[,td$rname] <- td$mf[,td$rname] + as.integer(plus_one)
    
    # support & bounds
    support <- c(-.5 + as.integer(log_first), quantile(td$response, probs = prob))
    bounds <- c(-.9 + as.integer(log_first), Inf)
    
    ret <- tram(td, transformation = "smooth", distribution = distribution, 
                log_first = log_first, support = support, bounds = bounds, 
                negative = TRUE, ...)
    if (!inherits(ret, "mlt")) return(ret)
    
    ret$call <- match.call(expand.dots = TRUE)
    ret$plus_one <- plus_one
    ret$support <- support
    ret$bounds <- bounds
    if (method != "probit") {
        ret$tram <- paste(ifelse(is.null(td$terms$s), "", "(Stratified)"), 
                          "Discrete", name[method], "Count Transformation Model")
    } else {
        ret$tram <- paste(ifelse(is.null(td$terms$s), "", "(Stratified)"),
                          "Transformed Counts Probit Transformation Model")
    }
    ret$count_response <- numeric_var(td$rname, support = min(td$response):max(td$response))
    class(ret) <- c("cotram", class(ret))
    ret
}
