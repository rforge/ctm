
.Surv2matrix <- function(object) {

    type <- attr(object, "type")
    stopifnot(type %in% c("left", "right", "interval"))
    status <- object[, "status"]

    ret <- switch(type,
        "right" = cbind(left = object[, "time"],
                        right = ifelse(status == 1, NA, Inf)),
        "left" = cbind(left = ifelse(status == 1, NA, -Inf),
                       right = object[, "time"]),
        "interval" = {
            status <- factor(object[, "status"], levels = 0:3, 
                             labels = c("right", "exact", "left", "interval"))
            tmp <- matrix(NA, nrow = nrow(object), ncol = 2)
            colnames(tmp) <- c("left", "right")
            for (s in levels(status)) {
                idx <- which(status == s)
                tmp[idx, ] <- switch(s, 
                    "left" = cbind(-Inf, object[idx, "time1"]),
                    "exact" = cbind(object[idx, "time1"], NA),
                    "right" = cbind(object[idx, "time1"], Inf),
                    "interval" = object[idx, c("time1", "time2")])
            }
            return(tmp)
        }
    )
    ret
}
