#' Monotonic binning based on isotonic regression
#'
#' The function \code{iso_bin} implements the monotonic binning based on 
#' the isotonic regression.
#'
#' @param x A numeric vector 
#' @param y A numeric vector with 0/1 binary values
#'
#' @return A list of binning outcomes, including a numeric vector with cut
#'         points and a dataframe with binning summary
#'
#' @examples
#' data(hmeq)
#' iso_bin(hmeq$DEROG, hmeq$BAD)

iso_bin <- function(x, y) {
  x_ <- x[!is.na(x)]
  y_ <- y[!is.na(x)]

  odx <- x_[order(x_)]
  ody <- y_[order(x_)]
  spc <- cor(odx, ody, method = "spearman")

  d1 <- with(isoreg(odx, spc / abs(spc) * ody), data.frame(x = x, y = y, cat = yf))

  l1 <- lapply(split(d1, d1$cat), 
               function(d) list(rate = abs(round(mean(d$y), 8)), maxx = max(d$x)))

  l2 <- l1[Reduce(c, lapply(l1, function(l) l$rate > 0 & l$rate < 1))]

  l3 <- sort(Reduce(c, lapply(l2, function(l) l$maxx)))[-length(l2)]

  l4 <- manual_bin(x_, y_, l3)

  return(list(cut = l3, tbl = gen_woe(add_miss(l4, x, y), l3)))
}
