#' Derive WoE, IV, and KS based on the pooled binning outcome
#'
#' The function \code{gen_woe2} calculates WoE, IV, and KS statistics based on 
#' the pooled binning outcome from the function \code{pool_bin}. 
#'
#' @param tbl A input dataframe
#' @param cut A list of numeric values as cut points
#'
#' @return A dataframe
#' @noRd

gen_woe2 <- function(tbl, cut) {
  d1 <- tbl[order(tbl$bin), ]

  gdist <- (d1$freq - d1$bads) / (sum(d1$freq) - sum(d1$bads))
  bdist <- d1$bads / sum(d1$bads)
  adist <- d1$freq / sum(d1$freq)

  d1$rate <- round(d1$bads / d1$freq, 4)
  d1$woe  <- round(log(bdist / adist), 4)
  d1$woe2 <- round(log(bdist / gdist), 4)
  d1$iv   <- round((bdist - gdist) * d1$woe2, 4)
  d1$ks   <- round(abs(cumsum(bdist) - cumsum(gdist)) * 100, 2)
  d1$rule <- NA

  for (i in seq(nrow(d1))) {
    if (d1[i, ]$bin == 0) {
      d1[i, ]$rule <- paste("is.na($X$)", sep = '')
    } else if (d1[i, ]$bin == 1) {
      d1[i, ]$rule <- ifelse(d1[i, ]$miss == 0, 
                             paste("$X$ <= ", cut[1], sep = ''), 
                             paste("$X$ <= ", cut[1], " | is.na($X$)", sep = ''))
    } else if (d1[i, ]$bin == max(d1$bin)) {
      d1[i, ]$rule <- ifelse(d1[i, ]$miss == 0,
                             paste("$X$ > ", cut[length(cut)], sep = ''),
                             paste("$X$ > ", cut[length(cut)], " | is.na($X$)", sep = ''))
    } else {
      d1[i, ]$rule <- paste("$X$ > ", cut[d1[i, ]$bin - 1], " & $X$ <= ", cut[d1[i, ]$bin], sep = '')
    }
  }
 
  keep <- c("bin", "freq", "miss", "bads", "rate", "woe", "iv", "ks", "rule")
  rownames(d1) <- NULL

  return(d1[, keep])
}
