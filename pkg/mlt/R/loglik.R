
.mlt_loglik <- function(d, mm, mmprime, offset = 0, mmtrunc = NULL) {
    trunc <- function(beta) {
        if (is.null(mmtrunc)) return(0)
        if (is.null(mmtrunc$right)) return(log(1 - d$p(mmtrunc$left %*% beta)))
        if (is.null(mmtrunc$left))return(log(d$p(mmtrunc$right %*% beta)))
        return(log(d$p(mmtrunc$right %*% beta) - d$p(mmtrunc$left %*% beta)))
    }
    function(beta)
        d$d(offset + mm %*% beta, log = TRUE) + log(mmprime %*% beta) - trunc(beta)
}

.mlt_score <- function(d, mm, mmprime, offset = 0, mmtrunc = NULL) {                          
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

