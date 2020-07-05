#' Perform WoE transformation of a numeric variable
#'
#' The function \code{cal_woe} applies the WoE transformation to a numeric 
#' vector based on the binning outcome from a binning function, e.g. qtl_bin()
#' or iso_bin().
#'
#' @param x   A numeric vector that will be transformed to WoE values. 
#' @param bin A list with the binning outcome from the binning function, 
#'            e.g. qtl_bin() or iso_bin()
#'
#' @return A numeric vector with WoE transformed values.
#'
#' @examples
#' data(hmeq)
#' bin_out <- qtl_bin(hmeq$DEROG, hmeq$BAD)
#' cal_woe(hmeq$DEROG[1:10], bin_out)

cal_woe <- function(x, bin) {
  cut <- sort(c(bin$cut, -Inf, Inf))
  cat <- Reduce(c, lapply(x, function(x_) ifelse(is.na(x_), 0, findInterval(x_, cut, left.open = T))))
  tbl <- bin$tbl[, c("bin", "woe")]

  d1 <- data.frame(i = seq(length(x)), x = x, bin = cat)

  d2 <- merge(x = d1, y = tbl, by = "bin", all.x = TRUE)
  d2$woe <- ifelse(is.na(d2$woe), 0, d2$woe)

  return(d2[order(d2$i), ]$woe)
}

