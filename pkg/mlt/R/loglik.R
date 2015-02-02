
.mlt_loglik_exact <- function(d, mm, mmprime, offset = 0, mmtrunc = NULL) {
    trunc <- function(beta) {
        if (is.null(mmtrunc)) return(0)
        if (is.null(mmtrunc$right)) return(.log(1 - d$p(mmtrunc$left %*% beta)))
        if (is.null(mmtrunc$left)) return(.log(d$p(mmtrunc$right %*% beta)))
        return(.log(d$p(mmtrunc$right %*% beta) - d$p(mmtrunc$left %*% beta)))
    }
    function(beta)
        d$d(offset + mm %*% beta, log = TRUE) + .log(mmprime %*% beta) - trunc(beta)
}

.mlt_score_exact <- function(d, mm, mmprime, offset = 0, mmtrunc = NULL) {                          
    trunc <- function(beta) {
        if (is.null(mmtrunc)) return(0)
        if (is.null(mmtrunc$right)) {
            mmtb <- mmtrunc$left %*% beta + offset
            return(- 1 / (1 - d$p(mmtb)) * d$d(mmtb) * mmtrunc)
        }
        if (is.null(mmtrunc$left)) {
            mmtb <- mmtrunc$right %*% beta + offset
            return(1 / d$p(mmtb) * d$d(mmtb) * mmtrunc)
        }
        mmtbl <- drop(mmtrunc$left %*% beta) + offset
        mmtbr <- drop(mmtrunc$right %*% beta) + offset
        return(1 / (d$p(mmtbr) - d$p(mmtbl)) * 
               (d$d(mmtbr) * mmtrunc$right - d$d(mmtbl) * mmtrunc$left))
    }
    function(beta) {
        mmb <- drop(mm %*% beta) + offset
        d$dd(mmb) / d$d(mmb) * mm + (1 / drop(mmprime %*% beta)) * mmprime - 
            trunc(beta)
    }
}

### if (finite) fun(X %*% beta + offset) else value
.dealinf <- function(X, beta, offset, fun, value, Xmult = FALSE) {
    OK <- is.finite(X[,1])
    tmp <- X[OK,] %*% beta
    ret <- numeric(nrow(X))
    ret[OK] <- fun(offset[OK] + tmp)
    ret[!OK] <- value
    if (Xmult) {
        X[!OK,] <- 0
        ret <- ret * X
    }
    ret
}

.log <- function(x)
    log(pmax(.Machine$double.eps, x))

.mlt_loglik_interval <- function(d, mml, mmr, offset = 0, mmtrunc = NULL) {
    trunc <- function(beta) {
        if (is.null(mmtrunc)) return(0)
        if (is.null(mmtrunc$right)) return(.log(1 - d$p(mmtrunc$left %*% beta)))
        if (is.null(mmtrunc$left))return(.log(d$p(mmtrunc$right %*% beta)))
        return(.log(d$p(mmtrunc$right %*% beta) - d$p(mmtrunc$left %*% beta)))
    }
    function(beta)
        .log(.dealinf(mmr, beta, offset, d$p, 1) - 
             .dealinf(mml, beta, offset, d$p, 0)) - trunc(beta)
}

.mlt_score_interval <- function(d, mml, mmr, offset = 0, mmtrunc = NULL) {
    trunc <- function(beta) {
        if (is.null(mmtrunc)) return(0)
        if (is.null(mmtrunc$right)) {
            mmtb <- mmtrunc$left %*% beta
            return(- 1 / (1 - d$p(mmtb)) * d$d(mmtb) * mmtrunc)
        }
        if (is.null(mmtrunc$left)) {
            mmtb <- mmtrunc$right %*% beta
            return(1 / d$p(mmtb) * d$d(mmtb) * mmtrunc)
        }
        mmtbl <- drop(mmtrunc$left %*% beta) + offset
        mmtbr <- drop(mmtrunc$right %*% beta) + offset
        return(1 / (d$p(mmtbr) - d$p(mmtbl)) * 
               (d$d(mmtbr) * mmtrunc$right - d$d(mmtbl) * mmtrunc$left))
    }
    function(beta) {
        1 / (.dealinf(mmr, beta, 0, d$p, 1) - 
             .dealinf(mml, beta, 0, d$p, 0)) *
        (.dealinf(mmr, beta, 0, d$d, 0, Xmult = TRUE) - 
         .dealinf(mml, beta, 0, d$d, 0, Xmult = TRUE))
    }
}
