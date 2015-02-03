
.is.formula <- function(x) {
    if (is.null(x)) return(FALSE)
    inherits(x, "formula")
}
