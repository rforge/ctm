
### cbind of multiple _named_ basi(e)s objects
c.basis <- function(..., recursive = FALSE) {

    stopifnot(!recursive)

    bases <- list(...)
    stopifnot(all(sapply(bases, inherits, what = "basis")))
    bnames <- names(bases)
    stopifnot(length(unique(bnames)) == length(bnames))

    varnames <- sapply(bases, varnames)
    stopifnot(all(!is.null(varnames)))
    support <- lapply(bases, support)

    ret <- function(data, model.matrix = TRUE, deriv = NULL, integrate = NULL) {
        if (!is.null(deriv)) {
            stopifnot(length(deriv) == 1)
            if (!names(deriv) %in% varnames) deriv <- NULL
        }
        if (!is.null(integrate)) {
            stopifnot(length(integrate) == 1)
            if (!names(integrate) %in% varnames) integrate <- NULL
        }
        ret <- lapply(bnames, function(b) {
            thisargs <- list()
            if (!is.null(deriv)) {
                if (names(deriv) %in% varnames(bases[[b]])) {
                    thisargs$deriv <- deriv
                } else {
                    X <- model.matrix(bases[[b]], data)
                    X[] <- 0
                    return(X)
                }
            }
            if (!is.null(integrate)) {
                if (names(integrate) %in% varnames(bases[[b]]))
                    thisargs$integrate <- integrate
            }
            thisargs$object <- bases[[b]]
            thisargs$data <- data
            X <- do.call("model.matrix", thisargs)
            X <- do.call("model.matrix", thisargs)
            attr(X, "Assign") <- rbind(attr(X, "Assign"), b)
            X
        })
        if (!model.matrix) return(ret)
        ui <- do.call("bdiag", lapply(ret, function(r)
                      attr(r, "constraint")$ui))
        ci <- do.call("c", lapply(ret, function(r)
                      attr(r, "constraint")$ci))
        if (length(bases) > 1) {
            a <- lapply(ret, function(r) matrix(attr(r, "Assign"), 
                                                ncol = ncol(r)))
            mr <- max(sapply(a, NROW))
            for (i in 1:length(a)) 
                a[[i]] <- a[[i]][rep(1:nrow(a[[i]]), length = mr),,
                                 drop = FALSE]
            ret <- do.call("cbind", ret)
            attr(ret, "Assign") <- do.call("cbind", a)
        } else {
            ret <- ret[[1]]
        }
        attr(ret, "constraint") <- list(ui = ui, ci = ci)
        return(ret )
    }
    attr(ret, "varnames") <- varnames
    attr(ret, "support") <- support
    class(ret) <- c("cbind_bases", "bases", "basis", class(ret))
    ret
}
