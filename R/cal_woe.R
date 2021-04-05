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
  tbl <- bin$tbl[, c("bin", "woe")]
  cut <- sort(c(bin$cut, -Inf, Inf))
  cat <- Reduce(c, lapply(unique(x), function(x_) ifelse(is.na(x_), 0, findInterval(x_, cut, left.open = TRUE))))

  d1 <- data.frame(i = seq(length(x)), x = x)
  d2 <- data.frame(x = unique(x), bin = cat)
  d3 <- merge(x = merge(x = d1, y = d2, by = "x", all.x = TRUE), y = tbl, by = "bin", all.x = TRUE)
  d3$woe <- ifelse(is.na(d3$woe), 0, d3$woe)

  return(d3[order(d3$i), ]$woe)
}

