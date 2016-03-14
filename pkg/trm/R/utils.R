
.R2vec <- function(object) {
    if (!inherits(object, "response")) return(object)
    ex <- object$exact
    le <- object$cleft
    ri <- object$cright
    ex[is.na(ex)] <- 0
    le[is.na(le) | !is.finite(le)] <- 0
    ri[is.na(ri) | !is.finite(ri)] <- 0
    ex + (le + (ri - le) / 2)
}
