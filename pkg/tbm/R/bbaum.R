
### the classical tree-based baselearner; stumps by default
### (also fits an additive model)
bbaum <- function(..., mtry = Inf,
    tree_controls = partykit::ctree_control(stump = TRUE,
                                            mincriterion = 0,
                                            saveinfo = FALSE)) {

    cll <- match.call()
    cll[[1]] <- as.name("bbaum")

    ctrl <- tree_controls
    mf <- list(...)
    if (length(mf) == 1 && is.data.frame(mf[[1]])) {
        mf <- mf[[1]]
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) as.character(x))
    }

    ret <- list(model.frame = function() return(mf),
                get_call = function(){
                    cll <- deparse(cll, width.cutoff=500L)
                    if (length(cll) > 1)
                        cll <- paste(cll, collapse="")
                    cll
                },
                get_names = function() colnames(mf),
                set_names = function(value) {
                    if(length(value) != length(colnames(mf)))
                        stop(sQuote("value"), " must have same length as ",
                             sQuote("colnames(mf)"))
                    for (i in 1:length(value)){
                        cll[[i+1]] <<- as.name(value[i])
                    }
                    attr(mf, "names") <<- value
                })
    class(ret) <- "blg"


    ret$dpp <- function(weights) {

        ### construct design matrix etc.
        y <- vector(length = nrow(mf), mode = "numeric")
        ### name for working response (different from any x)
        rname <- paste(sample(LETTERS, 25, replace = TRUE), collapse = "")
        fm <- as.formula(paste(rname, " ~ ", paste(colnames(mf), collapse = "+")))
        df <- mf
        df[[rname]] <- y
        d <- extree_data(fm, data = df, yx = "none")
        ytrafo <- function(subset, weights, info, estfun, object, ...) 
            list(estfun = Y, unweighted = TRUE) 
        mymf <- model.frame(d)
        tree_controls$update <- FALSE
        subset <- which(weights > 0)

        fitfun <- function(y) {
            if (!is.matrix(y)) y <- matrix(y, ncol = 1)
 
            assign("Y", y, envir = environment(ytrafo))
            tree <- extree_fit(data = d, trafo = ytrafo, converged = TRUE, 
                               partyvars = d$variables$z, subset = subset, 
                               weights = weights, ctrl = tree_controls, 
                               doFit = TRUE)$node
            where <- factor(fitted_node(tree, mymf))
            coef <- do.call("rbind", tapply(1:NROW(y), where, function(i)
                colSums(y[i,,drop = FALSE] * weights[i]) / sum(weights[i]), 
                simplify = FALSE))

            fitted <- function()
                return(coef[unclass(where),,drop = FALSE])

            predict <- function(newdata = NULL) {
                if (is.null(newdata)) newdata <- mymf
                vmatch <- match(names(mymf), names(newdata))
                wh <- factor(fitted_node(tree, newdata, vmatch = vmatch), 
                             levels = levels(where), labels = levels(where))
                return(coef[unclass(wh),,drop = FALSE])
            }

            ret <- list(model = tree, fitted = fitted, predict = predict)
            class(ret) <- c("bm_tree", "bm")
            ret
        }

        predict <- function(bm, newdata = NULL, aggregate = c("sum", "cumsum", "none")) {
            aggregate <- match.arg(aggregate)

            if (is.null(newdata)) 
                newdata <- mymf 

            pr <- 0
            for (i in 1:length(bm)) {
                pri <- bm[[i]]$predict(newdata)
                if (aggregate == "sum") {
                    pr <- pr + pri
                } else {
                    if (i > 1) {
                        pr <- cbind(pr, pri)
                    } else {
                        pr <- pri
                    }
                    if (aggregate == "cumsum")
                        if (i > 1) pr[,i] <- pr[,i] + pr[,i-1]
                }
            }
            return(pr)
        }

        ret <- list(fit = fitfun, predict = predict)
        class(ret) <- c("bl_tree", "bl")
        ret
    }
    return(ret)
}

schwarzboost <- function (formula, data = list(), weights = NULL, na.action = na.pass, 
    offset = NULL, family = Gaussian(), control = boost_control(), 
    oobweights = NULL, tree_controls = partykit::ctree_control(teststat = "quad", 
        testtype = "Teststatistic", mincriterion = 0, maxdepth = 2, 
        saveinfo = FALSE), ...) 
{
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights", "na.action"), 
        names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    response <- model.response(mf)
    weights <- model.weights(mf)
    mf <- mf[, -1, drop = FALSE]
    mf$"(weights)" <- NULL
    bl <- list(bbaum(mf, tree_controls = tree_controls))
    ret <- mboost_fit(bl, response = response, weights = weights, 
        offset = offset, family = family, control = control, 
        oobweights = oobweights, ...)
    ret$call <- cl
    ret$rownames <- rownames(mf)
    class(ret) <- c("blackboost", class(ret))
    ret
}

