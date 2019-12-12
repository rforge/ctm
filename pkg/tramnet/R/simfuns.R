# Functions for comparing tramnet vs tram

.compare <- function(tram, tramnet, fm, tramtime, tramnettime) {
  stopifnot(inherits(tram, "tram") && inherits(tramnet, "tramnet"))
  cfx_tram <- coef(tram, with_baseline = TRUE)
  cfx_tramnet <- coef(tramnet, with_baseline = TRUE, tol = 0)
  logLik_tram <- logLik(tram)
  logLik_tramnet <- logLik(tramnet)
  tramsize <- object.size(tram)
  tramnetsize <- object.size(tramnet)
  ret <- list(
    model = tram$tram,
    formula = fm,
    coefficients = .cmp(cfx_tram, cfx_tramnet),
    logLiks = .cmp(logLik_tram, logLik_tramnet),
    time = .cmp(tramtime, tramnettime),
    size = .cmp(tramsize, tramnetsize)
  )
  return(ret)
}

.cmp <- function(x, y) {
  ret <- cbind(x, y, x - y, (x - y) / x)
  colnames(ret) <- c("tram", "tramnet", "diff", "rel_diff")
  return(ret)
}

.runFUN <- function(modFUN, stratum, shift, response, data, ...) {
  if (is.null(stratum)) {
    if (is.null(shift)) {
      fm1 <- as.formula(
        paste(response, "~ 1")
      )
    } else {
      fm1 <- as.formula(
        paste(response, "~", paste(shift, collapse = " + "))
      )
    }
    fm3 <- as.formula(
      paste(response, "~ 1")
    )
  } else {
    if (is.null(shift)) {
      fm1 <- as.formula(
        paste(response, "| 0 + ",
              paste(stratum, collapse = "+"),
              "~ 1")
      )
    } else {
      fm1 <- as.formula(
        paste(response, "| 0 + ",
              paste(stratum, collapse = "+"),
              "~",
              paste(shift, collapse = "+"))
      )
    }
    fm3 <- as.formula(
      paste(response, "| 0 + ",
            paste(stratum, collapse = "+"),
            "~ 1"
      )
    )
  }
  if (is.null(shift)) {
    fm2 <- ~ 1
    x <- matrix(0, nrow = nrow(data))
  } else {
    fm2 <- as.formula(
      paste("~", paste(shift, collapse = "+"))
    )
    x <- model.matrix(fm2, data = data)[, -1, drop = FALSE]
  }

  trmt <- system.time(trm <- modFUN(fm1, data = data, ...))
  tmp <- modFUN(fm3, data = data, ...)
  trntt <- system.time(trnt <- tramnet(tmp, x = x, lambda = 0, alpha = 0,
                                       check_dcp = FALSE, solver = "ECOS"))
  return(.compare(tram = trm, tramnet = trnt, fm = fm1,
                  tramtime = trmt, tramnettime = trntt))
}
