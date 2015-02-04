
.is.formula <- function(x) {
    if (is.null(x)) return(FALSE)
    inherits(x, "formula")
}

.is.Surv <- function(x)
    inherits(x, "Surv")
