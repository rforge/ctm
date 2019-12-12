# mlrMBO tramnet interface functions for model based optimization of
# the tuning parameters based on cross-validated log-likelihoods

mbo_tramnet <- function(object, fold = 2, n_design = 5, n_iter = 5, minlambda = 0,
                        maxlambda = 16, minalpha = 0, maxalpha = 1,
                        folds = NULL, learner = "regr.km", pred.type = "se",
                        opt_crit = makeMBOInfillCritEI(), noisy = FALSE,
                        obj_type = c("lasso", "ridge", "elnet"), verbose = TRUE,
                        ...) {
  if (length(list(...)) > 0) {
    warning("Additional arguments ignored.")
  }
  stopifnot(inherits(object, "tramnet"))
  df <- .get_tramnet_data(object)
  rsp <- variable.names(object$model, "response")
  n <- nrow(df)
  if (!noisy) {
    if (is.null(folds) & !is.null(fold)) {
      folds <- sample(rep(1:fold, ceiling(n/fold)), n)
    } else {
      folds <- round(folds)
      fold <- max(folds)
    }
  } else {
    folds <- NULL
  }
  obj_type <- match.arg(obj_type)
  obj_fun <- switch(obj_type,
                    "elnet" = elnet_obj(object = object, minlambda = minlambda,
                                        maxlambda = maxlambda, minalpha = minalpha,
                                        maxalpha = maxalpha, folds = folds,
                                        noisy = noisy, fold = fold),
                    "lasso" = lasso_obj(object = object, minlambda = minlambda,
                                        maxlambda = maxlambda, folds = folds,
                                        noisy = noisy, fold = fold),
                    "ridge" = ridge_obj(object = object,minlambda = minlambda,
                                        maxlambda = maxlambda, folds = folds,
                                        noisy = noisy, fold = fold))
  des <- generateDesign(n = n_design, par.set = getParamSet(obj_fun),
                        fun = randomLHS)
  my_learner <- makeLearner(learner, predict.type = pred.type)
  control <- makeMBOControl()
  control <- setMBOControlTermination(control, iters = n_iter)
  control <- setMBOControlInfill(control, crit = opt_crit)
  ret <- mbo(obj_fun, design = des, learner = my_learner,
             control = control, show.info = verbose)
  return(ret)
}

# Elastic net objective function for model based optimization

elnet_obj <- function(object, minlambda = 0, maxlambda = 16, minalpha = 0,
                      maxalpha = 1, folds, noisy = FALSE, fold) {
  ret <- makeSingleObjectiveFunction(
    name = "elnet_obj",
    fn = function(lmb, alp = 1) {
      -2*cvl_tramnet(object, folds = folds, lambda = lmb, fold = fold,
                     alpha = alp)[["logLik_tab"]][["sum_logLik"]][1]
    },
    par.set = makeParamSet(
      makeNumericVectorParam("lmb", len = 1, lower = minlambda, upper = maxlambda),
      makeNumericVectorParam("alp", len = 1, lower = minalpha, upper = maxalpha)
    ),
    minimize = TRUE, noisy = noisy
  )
  return(ret)
}

# Lasso objective function for model based optimization

lasso_obj <- function(object, minlambda = 0, maxlambda = 16, folds, noisy = FALSE,
                      fold) {
  ret <- makeSingleObjectiveFunction(
    name = "lasso_obj",
    fn = function(lmb) {
      -2*cvl_tramnet(object, folds = folds, lambda = lmb, fold = fold,
                     alpha = 1)[["logLik_tab"]][["sum_logLik"]][1]
    },
    par.set = makeParamSet(
      makeNumericVectorParam("lmb", len = 1, lower = minlambda, upper = maxlambda)
    ),
    minimize = TRUE, noisy = noisy
  )
  return(ret)
}

# Ridge objective function for model based optimization

ridge_obj <- function(object, minlambda = 0, maxlambda = 16, folds, noisy = FALSE,
                      fold) {
  ret <- makeSingleObjectiveFunction(
    name = "ridge_obj",
    fn = function(lmb) {
      -2*cvl_tramnet(object, folds = folds, lambda = lmb, fold = fold,
                     alpha = 0)[["logLik_tab"]][["sum_logLik"]][1]
    },
    par.set = makeParamSet(
      makeNumericVectorParam("lmb", len = 1, lower = minlambda, upper = maxlambda)
    ),
    minimize = TRUE, noisy = noisy
  )
  return(ret)
}

# Fit recommended regularized tram based on model based optimization output

mbo_recommended <- function(mbo_obj, m0, x, ...) {
  rec_lambda <- ifelse(is.null(mbo_obj[["x"]][["lmb"]]),
                       0, mbo_obj[["x"]][["lmb"]])
  rec_alpha <- ifelse(is.null(mbo_obj[["x"]][["alp"]]),
                      .get_alpha_from_obj(mbo_obj),
                      mbo_obj[["x"]][["alp"]])
  ret <- tramnet(m0, x = x, lambda = rec_lambda, alpha = rec_alpha, ...)
  return(ret)
}

# Helper Functions

.get_alpha_from_obj <- function(mbo_obj) {
  nm <- attr(mbo_obj$final.opt.state$opt.problem$fun, "name")
  ifelse(nm == "lasso_obj", 1, 0)
}
