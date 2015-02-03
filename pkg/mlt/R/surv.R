
.Surv2matrix <- function(object) {

    type <- attr(object, "type")
    stopifnot(type %in% c("left", "right", "interval", 
                          "interval2", "counting"))
    status <- object[, "status"]

    ret <- switch(type,
        "right" = cbind(left = object[, "time"],
                        right = ifelse(status == 1, NA, Inf)),
        "left" = cbind(left = ifelse(status == 1, NA, -Inf),
                       right = object[, "time"]),
        "interval2" = {
            ret <- cbind(left = object[, "time1"], 
                         right = object[, "time2"])
            ret$left[is.na(ret$left)] <- -Inf
            ret$right[is.na(ret$right)] <- Inf
        },
        "interval" = {
            status <- factor(status, levels = 0:3, 
                             labels = c("right", "exact", "left", "interval"))
            tmp <- matrix(NA, nrow = nrow(object), ncol = 2)
            colnames(tmp) <- c("left", "right")
            for (s in levels(status)) {
                idx <- which(status == s)
                tmp[idx, ] <- switch(s, 
                    "right" = cbind(object[idx, "time1"], Inf),
                    "exact" = cbind(object[idx, "time1"], NA),
                    "left" = cbind(-Inf, object[idx, "time1"]),
                    "interval" = object[idx, c("time1", "time2")])
            }
            return(tmp)
        },
        ### left truncation, right censoring
        "counting" = cbind(left = object[, "stop"],
                           right = ifelse(status == 1, NA, Inf),
                           lefttrunc = object[, "start"])
    )
    ### right truncation, left censoring? truncation and interval?
    ret
}
