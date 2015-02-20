
### if (finite) fun(X %*% beta + offset) else value
.dealinf <- function(X, beta, offset, fun, value, Xmult = FALSE) {
    if (is.null(X)) return(value)
    OK <- is.finite(X[,1])
    if (all(!OK)) return(rep(value, nrow(X)))
    tmp <- X[OK,] %*% beta
    ret <- numeric(nrow(X))
    ret[OK] <- fun(offset[OK] + tmp)
    ret[!OK] <- value
    if (Xmult) {
        X[!OK,] <- 0
        ret <- ret * X
    }
    return(ret)
}

.log <- function(x)
    log(pmax(.Machine$double.eps, x))

.trunc_loglik <- function(beta, d, offset = 0, mmtrunc) {
    if (is.null(mmtrunc)) return(0)
    return(.log(.dealinf(mmtrunc$right, beta, offset, d$p, 1) - 
                .dealinf(mmtrunc$left, beta, offset, d$p, 0)))
}

.trunc_score <- function(beta, d, offset = 0, mmtrunc) {
    if (is.null(mmtrunc)) return(0)
    return((.dealinf(mmtrunc$right, beta, offset, d$d, 0, TRUE) - 
            .dealinf(mmtrunc$left, beta, offset, d$d, 0, TRUE)) / 
           (.dealinf(mmtrunc$right, beta, offset, d$p, 1) - 
            .dealinf(mmtrunc$left, beta, offset, d$p, 0)))
}

.trunc_hessian <- function(beta, d, offset = 0, mmtrunc, w = 1) {
    if (is.null(mmtrunc)) return(0)
    Fr <- .dealinf(mmtrunc$right, beta, offset, d$p, 1)
    Fl <- .dealinf(mmtrunc$left, beta, offset, d$p, 0)
    fr <- .dealinf(mmtrunc$right, beta, offset, d$d, 0)
    fl <- .dealinf(mmtrunc$left, beta, offset, d$d, 0)
    dfr <- .dealinf(mmtrunc$right, beta, offset, d$dd, 0)
    dfl <- .dealinf(mmtrunc$left, beta, offset, d$dd, 0)
    Frl <- Fr - Fl
    w1 <- fr / sqrt(Frl^2) * sqrt(w)
    w2 <- fl / sqrt(Frl^2) * sqrt(w)
    w3 <- dfr / Frl * w
    w4 <- dfl / Frl * w
    if (is.null(mmtrunc$right)) 
        mmtrunc$right <- matrix(0, nrow = nrow(mmtrunc$left), 
                                   ncol = ncol(mmtrunc$left))
    if (is.null(mmtrunc$left)) 
        mmtrunc$left <- matrix(0, nrow = nrow(mmtrunc$right), 
                                  ncol = ncol(mmtrunc$right))
    return((crossprod(mmtrunc$right * w1) - 
            crossprod(mmtrunc$right * w1, mmtrunc$left * w2) - 
            crossprod(mmtrunc$left * w2, mmtrunc$right * w1) + 
            crossprod(mmtrunc$left * w2))
           - (crossprod(mmtrunc$right * w3, mmtrunc$right) - 
              crossprod(mmtrunc$left * w4, mmtrunc$left)))
}

.mlt_loglik_exact <- function(d, mm, mmprime, offset = 0, mmtrunc = NULL) {
    function(beta)
        return(d$d(offset + mm %*% beta, log = TRUE) + 
               .log(mmprime %*% beta) - .trunc_loglik(beta, d, offset, mmtrunc))
}

.mlt_score_exact <- function(d, mm, mmprime, offset = 0, mmtrunc = NULL) {                          
    function(beta) {
        mmb <- drop(mm %*% beta) + offset
        return(d$dd(mmb) / d$d(mmb) * mm + (1 / drop(mmprime %*% beta)) * mmprime - 
               .trunc_score(beta, d, offset, mmtrunc))
    }
}

.mlt_hessian_exact <- function(d, mm, mmprime, offset = 0, mmtrunc = NULL, w = 1) {
    function(beta) {
        mmb <- drop(mm %*% beta) + offset
        if (length(w) != length(mmb)) w <- rep(w, length(mmb))
        w1 <- -(d$ddd(mmb) / d$d(mmb) - (d$dd(mmb) / d$d(mmb))^2) * w
        w2 <- w / (drop(mmprime %*% beta)^2)
        return(crossprod(mm * w1, mm) + crossprod(mmprime * w2, mmprime) 
               - .trunc_hessian(beta, d, offset, mmtrunc, w)) ### - or + ???
    }
}

.mlt_loglik_interval <- function(d, mml, mmr, offset = 0, mmtrunc = NULL) {
    function(beta)
        .log(.dealinf(mmr, beta, offset, d$p, 1) - 
             .dealinf(mml, beta, offset, d$p, 0)) - 
        .trunc_loglik(beta, d, offset, mmtrunc)
}

.mlt_score_interval <- function(d, mml, mmr, offset = 0, mmtrunc = NULL) {
    function(beta) {
        (.dealinf(mmr, beta, offset, d$d, 0, Xmult = TRUE) -
         .dealinf(mml, beta, offset, d$d, 0, Xmult = TRUE)) / 
        (.dealinf(mmr, beta, offset, d$p, 1) - 
         .dealinf(mml, beta, offset, d$p, 0)) 
    }
}

.mlt_hessian_interval <- function(d, mml, mmr, offset = 0, mmtrunc = NULL, w = 1) {
    function(beta) {
        Fr <- .dealinf(mmr, beta, offset, d$p, 1)
        Fl <- .dealinf(mml, beta, offset, d$p, 0)
        fr <- .dealinf(mmr, beta, offset, d$d, 0)
        fl <- .dealinf(mml, beta, offset, d$d, 0)
        dfr <- .dealinf(mmr, beta, offset, d$dd, 0)
        dfl <- .dealinf(mml, beta, offset, d$dd, 0)
        if (length(w) != length(Fr)) w <- rep(w, length(Fr))
        Frl <- Fr - Fl
        w1 <- dfr / Frl * w
        w2 <- dfl / Frl * w
        w3 <- fr / Frl * sqrt(w)
        w4 <- fl / Frl * sqrt(w)
        mmr[!is.finite(mmr)] <- 0
        mml[!is.finite(mml)] <- 0
        W3 <- mmr * w3
        W4 <- mml * w4
        return(-(crossprod(mmr * w1, mmr) - crossprod(mml * w2, mml) - 
                 (crossprod(W3) - crossprod(W3, W4) - crossprod(W4, W3) + crossprod(W4)) - 
                 .trunc_hessian(beta, d, offset, mmtrunc, w)))
    }
}
