#' Monotonic binning based on k-means clustering
#'
#' The function \code{kmn_bin} implements the monotonic binning based on 
#' the k-means clustering
#'
#' @param x A numeric vector 
#' @param y A numeric vector with 0/1 binary values
#'
#' @return A list of binning outcomes, including a numeric vector with cut
#'         points and a dataframe with binning summary
#'
#' @examples
#' data(hmeq)
#' kmn_bin(hmeq$DEROG, hmeq$BAD)

kmn_bin <- function(x, y) {
  x_ <- x[!is.na(x)]
  y_ <- y[!is.na(x)]
  n_ <- 2:max(2, min(50, length(unique(x_)) - 1))

  set.seed(1)
  c1 <- lapply(n_, function(n) kmeans(x_, centers = n, nstart = 10, iter.max = 500, algorithm = "MacQueen"))

  c2 <- lapply(c1, function(c) sort(Reduce(c, lapply(split(x_, c$cluster), max))))
  
  p_ <- unique(append(lapply(c2, function(c) c[-length(c)]), list(median(x_[y_ == 1]))))

  l1 <- lapply(p_, function(p) list(cut = p, out = manual_bin(x_, y_, p)))
 
  l2 <- lapply(l1[order(Reduce(c, lapply(l1, function(l) -length(l$cut))))], 
               function(l) list(cut  = l$cut, 
                                minr = min(l$out$bads / l$out$freq), 
                                maxr = max(l$out$bads / l$out$freq), 
                                scor = round(cor(l$out$bin, l$out$bads / l$out$freq, method = "spearman"), 8)))

  l3 <- l2[Reduce(c, lapply(l2, function(l) abs(l$scor) == 1 & l$minr > 0 & l$maxr < 1))][[1]]

  l4 <- l1[Reduce(c, lapply(l1, function(l) identical(l$cut, l3$cut)))][[1]]$out

  return(list(cut = l3$cut, tbl = gen_woe(add_miss(l4, x, y), l3$cut)))
}
