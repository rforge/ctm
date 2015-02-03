
model <- function(response = NULL, interacting = NULL, shifting = NULL) {

    if (.is.formula(response)) 
        response <- as.basis(response)
    if (.is.formula(interacting)) 
        interacting <- as.basis(interacting)
    if (.is.formula(shifting)) 
        shifting <- as.basis(shifting, remove_intercept = TRUE)

    if (is.null(shifting)) {
        if (!is.null(interacting)) {
            ret <- b(bresponse = response, binteracting = interacting)
        } else {
            ret <- c(bresponse = response)
        }   
    } else {
        if (!is.null(interacting)) {
            ret <- c(b(bresponse = response, binteracting = interacting), 
                    bshifting = shifting)
        } else {
            ret <- c(bresponse = response, bshifting = shifting)
        }
    }
    attr(ret, "response") <- varnames(response)
    return(ret)
}
