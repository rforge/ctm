
model <- function(response, interacting = NULL, shifting = NULL, 
                  todistr = c("Normal", "Logistic", "MinExtrVal"), 
                  sumconstr = inherits(interacting, c("formula", "formula_basis"))) {

    ### mkgrid() will not work due to missing data
    if (.is.formula(response)) 
        response <- as.basis(response)
    if (.is.formula(interacting)) 
        interacting <- as.basis(interacting)
    if (.is.formula(shifting)) 
        shifting <- as.basis(shifting, remove_intercept = TRUE)

    if (is.character(todistr))
        todistr <- .distr(todistr)

    if (!is.null(interacting))
        interacting <- b(iresponse = response, iinteracting = interacting, 
                         sumconstr = sumconstr)

    if (is.null(interacting) && is.null(shifting)) {
        mod <- c(bresponse = response)
    } else if (!is.null(interacting) && is.null(shifting)) {
        mod <- c(binteracting = interacting)
    } else if (is.null(interacting) && !is.null(shifting)) {
        mod <- c(bresponse = response, bshifting = shifting)
    } else {
        mod <- c(binteracting = interacting, bshifting = shifting)
    }
    ret <- list(model = mod, response = variable.names(response), 
                todistr = todistr)
    class(ret) <- "model"
    return(ret)
}

model.matrix.model <- function(object, data, ...)
    return(model.matrix(object$model, data = data, ...))

variable.names.model <- function(object, ...)
    variable.names(object$model)

