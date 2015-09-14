
.var <- function(name, desc = NULL) {
    ret <- list(name = name, desc = desc)
    class(ret) <- "var"
    ret
}

factor_var <- function(name, desc = NULL, levels) {
    ret <- .var(name = name, desc = desc)
    ret$support <- factor(levels)
    class(ret) <- c("factor_var", class(ret))
    ret
}

ordered_var <- function(name, desc = NULL, levels) {
    ret <- factor_var(name = name, desc = desc, levels = levels)
    ret$support <- as.ordered(ret$support)
    class(ret) <- c("ordered_var", class(ret))
    ret
}

numeric_var <- function(name, desc = NULL, unit = NULL,
                        support = c(0.0, 1.0), bounds = NULL) {
    ret <- .var(name = name, desc = desc)
    ret$unit <- unit
    stopifnot(length(support) >= 2L)
    stopifnot(all(is.finite(support)))
    stopifnot(is.integer(support) || is.double(support))
    ret$support <- support
    discrete <- is.integer(support) || (length(support) > 2L)
    if (discrete) {
        stopifnot(is.null(bounds))
        class(ret) <- c("discrete_var", "numeric_var", class(ret))
        return(ret)
    }
    if (is.null(bounds))
        bounds <- c(-Inf, Inf)
    stopifnot(bounds[1] <= min(support))
    stopifnot(max(support) <= bounds[2])
    ret$bounds <- bounds
    class(ret) <- c("continuous_var", "numeric_var", class(ret))
    ret
}

c.var <- function(...) {
    ret <- list(...)
    names(ret) <- sapply(ret, variable.names)
    stopifnot(all(sapply(ret, function(x) inherits(x, "var"))))
    class(ret) <- "vars"
    ret
}
    
variable.names.var <- function(object, ...)
    object$name

variable.names.vars <- function(object, ...)
    sapply(object, variable.names)

desc <- function(x)
    UseMethod("desc")

desc.var <- function(x)
    x$desc

desc.vars <- function(x)
    sapply(x, desc)

units <- function(x)
    UseMethod("units")

units.numeric_var <- function(x)
    x$unit

units.var <- function(x)
    return(NA)

units.vars <- function(x)
    sapply(x, units)

support <- function(object)
    UseMethod("support")

support.var <- function(object)
    object$support

support.vars <- function(object)
   lapply(object, support)

levels.factor_var <- function(x)
    levels(support(x))

levels.discrete_var <- function(x)
    support(x)

levels.var <- function(x)
    return(NA)

levels.vars <- function(x)
    lapply(x, levels)

bounds <- function(object)
    UseMethod("bounds")

bounds.continuous_var <- function(object) 
    object$bounds

bounds.discrete_var <- function(object)
    return(range(support(object)))

bounds.ordered_var <- function(object) {
    f <- support(object)
    return(f[c(1, nlevels(f))])
}

bounds.vars <- function(object)
    lapply(object, bounds)
    
bounds.default <- function(object)
    return(NA)

is.bounded <- function(object)
    UseMethod("is.bounded")

is.bounded.continuous_var <- function(object)
    any(is.finite(bounds(object)))

is.bounded.var <- function(object)
    return(TRUE)

is.bounded.vars <- function(object)
    sapply(object, is.bounded)

mkgrid <- function(object, ...)
    UseMethod("mkgrid")

mkgrid.var <- function(object, ...)
    return(support(object))

mkgrid.continuous_var <- function(object, n = 2, ...) {
    s <- support(object)
    stopifnot(n > 0)
    if (n == 1L) return(diff(s))    
    return(seq(from = s[1], to = s[2], length.out = n))
}
    
mkgrid.vars <- function(object, ...)
    lapply(object, mkgrid, ...)

as.data.frame.vars <- function(x, row.names = NULL, optional = FALSE, n = 1L, ...) {
    g <- mkgrid(x, n = n)
    len <- max(sapply(g, length))
    as.data.frame(lapply(g, function(x) rep_len(x, length.out = len)))
}

as.vars <- function(data) {
    stopifnot(is.data.frame(data))
    v <- lapply(colnames(data), function(x) {
        if (is.ordered(data[[x]])) return(ordered_var(x, levels = levels(data[[x]])))
        if (is.factor(data[[x]])) return(factor_var(x, levels = levels(data[[x]])))
        if (is.integer(data[[x]])) {
            s <- sort(unique(data[[x]]))
        } else {
            s <- range(data[[x]])
        }
        return(numeric_var(x, support = s))
    })
    return(do.call("c", v))
}

check <- function(object, data)
    UseMethod("check")

check.ordered_var <- function(object, data) {
    v <- variable.names(object)
    stopifnot(v %in% colnames(data))
    is.ordered(data[[v]]) && all.equal(levels(data[[v]]), levels(object))
}

check.factor_var <- function(object, data) {
    v <- variable.names(object)
    stopifnot(v %in% colnames(data))
    is.factor(data[[v]]) && all.equal(levels(data[[v]]), levels(object))
}

check.discrete_var <- function(object, data) {
    v <- variable.names(object)
    stopifnot(v %in% colnames(data))
    all(data[[v]] %in% support(object))
}

check.continuous_var <- function(object, data) {
    v <- variable.names(object)
    stopifnot(v %in% colnames(data))
    b <- bounds(object)
    min(data[[v]], na.rm = TRUE) >= b[1] && 
    max(data[[v]], na.rm = TRUE) <= b[2]
}

check.vars <- function(object, data)
    all(sapply(object, check, data = data))
    
