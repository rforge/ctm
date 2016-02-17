
checkGH <- function(object) {

    ### check gradient  and hessian
    gr <- numDeriv::grad(object$loglik, coef(object), weights = weights(object))
    s <- Gradient(object)
    cat("Gradient")
    print(all.equal(gr, s, check.attributes = FALSE))

    H1 <- numDeriv::hessian(object$loglik, coef(object), weights = weights(object))
    H2 <- Hessian(object)
    cat("Hessian:")
    print(all.equal(H1, H2, check.attributes = FALSE))
}
