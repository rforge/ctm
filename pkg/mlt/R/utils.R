
.is.formula <- function(x) {
    if (is.null(x)) return(FALSE)
    inherits(x, "formula")
}

.is.Surv <- function(x)
    inherits(x, "Surv")

.type_of_response <- function(y) {
    if (!is.null(dim(y))) {
        if (.is.Surv(y)) return("survival")
        return(NA)
    }
    if (storage.mode(y) == "double") return("double")
    if (is.integer(y)) return("integer")
    if (is.ordered(y)) return("ordered")
    if (is.factor(y)) return("unordered")
    return(NA)
}
