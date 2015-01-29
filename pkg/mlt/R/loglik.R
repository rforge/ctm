
.mlt_loglik_exact <- function(d, mm, mmprime, offset = 0, mmtrunc = NULL) {
    trunc <- function(beta) {
        if (is.null(mmtrunc)) return(0)
        if (is.null(mmtrunc$right)) return(log(1 - d$p(mmtrunc$left %*% beta)))
        if (is.null(mmtrunc$left))return(log(d$p(mmtrunc$right %*% beta)))
        return(log(d$p(mmtrunc$right %*% beta) - d$p(mmtrunc$left %*% beta)))
    }
    function(beta)
        d$d(offset + mm %*% beta, log = TRUE) + log(mmprime %*% beta) - trunc(beta)
}

.mlt_score_exact <- function(d, mm, mmprime, offset = 0, mmtrunc = NULL) {                          
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
        mmtbr <- mmtrunc$right %*% beta

    }
    function(beta) {
        mmb <- drop(mm %*% beta) + offset
        d$dd(mmb) / d$d(mmb) * mm + (1 / drop(mmprime %*% beta)) * mmprime - 
            trunc(beta)
    }
}

.mlt_loglik_interval <- function(d, mml, mmr, offset = 0, mmtrunc = NULL) {
    trunc <- function(beta) {
        if (is.null(mmtrunc)) return(0)
        if (is.null(mmtrunc$right)) return(log(1 - d$p(mmtrunc$left %*% beta)))
        if (is.null(mmtrunc$left))return(log(d$p(mmtrunc$right %*% beta)))
        return(log(d$p(mmtrunc$right %*% beta) - d$p(mmtrunc$left %*% beta)))
    }
    function(beta)
        log(ifelse(!is.finite(mmr[,1]), 1, d$p(offset + mmr %*% beta)) - 
            ifelse(!is.finite(mml[,1]), 0, d$p(offset + mml %*% beta))) - trunc(beta)
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
        mmtbr <- mmtrunc$right %*% beta

    }
    function(beta) {
        mmrb <- drop(mmr %*% beta) + offset
        mmlb <- drop(mml %*% beta) + offset
        1 / (ifelse(!is.finite(mmrb), 1, d$p(mmrb)) - 
             ifelse(!is.finite(mmlb), 1, d$p(mmlb))) *
        (ifelse(!is.finite(mmrb), 0, d$d(mmrb) * mmr) -  
         ifelse(!is.finite(mmlb), 0, d$d(mmlb) * mml))
    }
}
