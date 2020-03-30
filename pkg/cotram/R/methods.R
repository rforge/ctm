
mkgrid.cotram <- function(object, n, ...)
  mkgrid(object$count_response, n = n, ...)

as.mlt.cotram <- function(object) {
  class(object) <- class(object)[-(1:2)]
  return(object)
}

logLik.cotram <- function(object, parm = coef(as.mlt(object), fixed = FALSE), newdata,...){
  response <- variable.names(object, "response")
  
  if(!missing(newdata)){
    ## check whether response is non-negative integer
    if (any(newdata[response] < 0))
      stop("response is non-positive")
    if(isTRUE(newdata[response] %% 1 == 0))
      stop("response is non-integer")
    newdata[,response] <- as.integer(newdata[,response]) + as.integer(object$plus_one)
    return(logLik(as.mlt(object), parm = parm, newdata = newdata, ...))
  }
  logLik(as.mlt(object), parm = parm, ...)
}
