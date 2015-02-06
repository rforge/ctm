
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
    ret <- list(model = mod, response = varnames(response), todistr = todistr)
    class(ret) <- "model"
    return(ret)
}

model.matrix.model <- function(object, data, ...)
    model.matrix(object$model, data = data, ...)

varnames.model <- function(x)
    varnames(x$model)
