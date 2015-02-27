
model <- function(response, interacting = NULL, shifting = NULL, 
                  remove_intercept = FALSE,
                  todistr = c("Normal", "Logistic", "MinExtrVal")) {

    ### generate() will not work due to missing data
    if (.is.formula(response)) 
        response <- as.basis(response)
    if (.is.formula(interacting)) 
        interacting <- as.basis(interacting, 
                                remove_intercept = !remove_intercept)
    if (.is.formula(shifting)) 
        shifting <- as.basis(shifting, remove_intercept = TRUE)

    if (is.character(todistr))
        todistr <- .distr(todistr)

    if (!is.null(interacting)) {
        if (!(varnames(response) %in% varnames(interacting)))
            interacting <- b(iresponse = response, iinteracting = interacting)
    }

    if (is.null(interacting) && is.null(shifting)) {
        stopifnot(!remove_intercept)
        mod <- c(bresponse = response)
    } else if (!is.null(interacting) && is.null(shifting)) {
        if (remove_intercept)
            mod <- c(binteracting = interacting)
        else
            mod <- c(bresponse = response, binteracting = interacting)
    } else if (is.null(interacting) && !is.null(shifting)) {
        stopifnot(!remove_intercept)
        mod <- c(bresponse = response, bshifting = shifting)
    } else {
        if (remove_intercept)
            mod <- c(binteracting = interacting, bshifting = shifting)
        else
            mod <- c(bresponse = response, binteracting = interacting, 
                     bshifting = shifting)
    }
    ret <- list(model = mod, response = varnames(response), 
                todistr = todistr, remove_intercept = remove_intercept)
    class(ret) <- "model"
    return(ret)
}

model.matrix.model <- function(object, data, ...) {
    if (is.null(object$model$binteracting) ||
        is.null(object$model$bresponse))
        return(model.matrix(object$model, data = data, ...))
    X <- model.matrix(object$model, data = data, ...)
    ui <- attr(X, "constraint")$ui
    ci <- attr(X, "constraint")$ci
    Xy <- model.matrix(object$model$bresponse, data = data, ...)
    uiy <- attr(Xy, "constraint")$ui
    ciy <- attr(Xy, "constraint")$ci
    if (!any(is.finite(ciy))) return(X)
    uiy <- uiy[is.finite(ciy),,drop = FALSE]
    a <- attr(X, "Assign")
    i <- grep("bresponse", apply(a, 2, paste, collapse = "-"))
    ui[is.finite(ci),i] <- uiy[rep(1:nrow(uiy), sum(is.finite(ci)) / nrow(uiy)),]
    attr(X, "constraint") <- list(ui = ui, 
        ci = attr(X, "constraint")$ci)
    X
}



varnames.model <- function(x)
    varnames(x$model)

