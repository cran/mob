#' Monotonic binning based on generalized boosted model
#'
#' The function \code{gbm_bin} implements the monotonic binning based on 
#' the generalized boosted model (GBM).
#'
#' @param x A numeric vector 
#' @param y A numeric vector with 0/1 binary values
#'
#' @return A list of binning outcomes, including a numeric vector with cut
#'         points and a dataframe with binning summary
#'
#' @examples
#' data(hmeq)
#' gbm_bin(hmeq$DEROG, hmeq$BAD)

gbm_bin <- function(x, y) {
  x_ <- x[!is.na(x)]
  y_ <- y[!is.na(x)]

  spc <- cor(x_, y_, method = "spearman")

  set.seed(1)
  m_ <- gbm::gbm(y ~ x1 + x2, distribution = "bernoulli", data = data.frame(y = y_, x1 = x_, x2 = x_), 
                 var.monotone = c(spc / abs(spc), spc / abs(spc)),
                 bag.fraction = 1, n.minobsinnode = round(length(x_) / 100), n.trees = 500)

  d1 <- data.frame(y = y_, x = x_, cat = gbm::predict.gbm(m_, n.trees = m_$n.trees, type = "response"))

  l1 <- lapply(split(d1, d1$cat), 
               function(d) list(rate = abs(round(mean(d$y), 8)), maxx = max(d$x)))

  l2 <- l1[Reduce(c, lapply(l1, function(l) l$rate > 0 & l$rate < 1))]

  l3 <- sort(Reduce(c, lapply(l2, function(l) l$maxx)))[-length(l2)]

  l4 <- manual_bin(x_, y_, l3)

  return(list(cut = l3, tbl = gen_woe(add_miss(l4, x, y), l3)))
}

