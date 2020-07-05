#' Monotonic binning based on decision tree model
#'
#' The function \code{arb_bin} implements the monotonic binning based on 
#' the decision tree.
#'
#' @param x A numeric vector 
#' @param y A numeric vector with 0/1 binary values
#'
#' @return A list of binning outcomes, including a numeric vector with cut
#'         points and a dataframe with binning summary
#'
#' @examples
#' data(hmeq)
#' arb_bin(hmeq$DEROG, hmeq$BAD)

arb_bin <- function(x, y) {
  x_ <- x[!is.na(x)]
  y_ <- y[!is.na(x)]

  spc <- cor(x_, y_, method = "spearman")

  set.seed(1)
  m_ <- Rborist::Rborist(cbind(x_, x_), y_, noValidate = T, nTree = 1, regMono = c(spc / abs(spc), spc / abs(spc)),
                         ctgCensus = "prob", nSamp = length(x_) , withRepl = F, minNode = round(length(x_) / 100))

  d1 <- data.frame(y = y_, x = x_, cat = predict(m_, newdata = cbind(x_, x_), ctgCensus = "prob")$yPred)

  l1 <- lapply(split(d1, d1$cat), 
               function(d) list(rate = abs(round(mean(d$y), 8)), maxx = max(d$x)))

  l2 <- l1[Reduce(c, lapply(l1, function(l) l$rate > 0 & l$rate < 1))]

  l3 <- sort(Reduce(c, lapply(l2, function(l) l$maxx)))[-length(l2)]

  l4 <- manual_bin(x_, y_, l3)

  return(list(cut = l3, tbl = gen_woe(add_miss(l4, x, y), l3)))
}
