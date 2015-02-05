
model <- function(response, interacting = NULL, shifting = NULL, 
                  todistr = c("Normal", "Logistic", "MinExtrVal")) {

    if (.is.formula(response)) 
        response <- as.basis(response)
    if (.is.formula(interacting)) 
        interacting <- as.basis(interacting)
    if (.is.formula(shifting)) 
        shifting <- as.basis(shifting, remove_intercept = TRUE)

    if (is.character(todistr))
        todistr <- .distr(todistr)

    if (is.null(shifting)) {
        if (!is.null(interacting)) {
            mod <- b(bresponse = response, binteracting = interacting)
        } else {
            mod <- c(bresponse = response)
        }   
    } else {
        if (!is.null(interacting)) {
            mod <- c(b(bresponse = response, binteracting = interacting), 
                    bshifting = shifting)
        } else {
            mod <- c(bresponse = response, bshifting = shifting)
        }
    }
    ret <- list(model = mod, response = varnames(response), todistr = todistr)
    class(ret) <- "model"
    return(ret)
}

model.matrix.model <- function(object, data, ...)
    model.matrix(object$model, data = data, ...)

varnames.model <- function(x)
    varnames(x$model)
