
pextreme <- function(x){
  answer <- 1 - exp(-exp(x))
  answer
}

dextreme <- function(x){
  answer <- exp(x - exp(x))
  answer
}

qextreme <- function(p){
  q <- log(-log(1 - p))
  q
}

library("VGAM")
Binomial_extreme <- function(){
  biny <- function(y){
          if(!is.factor(y)) 
            stop("response is not a factor but ", sQuote("family = Binomial()"))
          if(nlevels(y) != 2) 
            stop("response is not a factor at two levels but ", 
                sQuote("family = Binomial()"))
    return(c(-1, 1)[as.integer(y)])
  }
  trf <- function(f) {
    pmax(pmin(f, 2), -7)
  }
  return(Family(
      ngradient = function(y, f, w = 1){
          trf <- function(f){
                 pmax(pmin(f, 2), -7)
          }
          y <- (y + 1)/2
          p <- pextreme(trf(f))
          d <- dextreme(trf(f))
          d * (y/p - (1 - y)/(1 - p))
      }, 
      loss = function(y, f){
          trf <- function(f) {
            pmax(pmin(f, 2), -7)
          }
          p <- pextreme(trf(f))
          y <- (y + 1)/2
          -y * log(p) - (1 - y) * log(1 - p)
      }, 
     offset = function(y, w){
         p <- weighted.mean(y > 0, w)
         qextreme(p)
     }, 
     response = function(f){
         p <- pextreme(trf(f))
         return(p)
     }, 
  rclass = function(f) (f > 0) + 1, check_y = biny, 
  name = paste("Negative Binomial Likelihood -- Minimum-Extremevalue-Link")))
}
