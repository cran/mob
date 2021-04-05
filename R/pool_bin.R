#' Monotonic binning for the pool data
#'
#' The function \code{pool_bin} implements the monotonic binning for the pool data
#' based on the generalized boosted model (GBM).
#'
#' @param x   A numeric vector
#' @param num A numeric vector with integer values for numerators to calculate bad rates
#' @param den A numeric vector with integer values for denominators to calculate bad rates
#'
#' @return A list of binning outcomes, including a numeric vector with cut
#'         points and a dataframe with binning summary
#'
#' @examples
#' data(hmeq)
#' df <- rbind(Reduce(rbind, 
#'                    lapply(split(hmeq, floor(hmeq$CLAGE)),
#'                           function(d) data.frame(AGE = unique(floor(d$CLAGE)),
#'                                                  NUM = sum(d$BAD),
#'                                                  DEN = nrow(d)))),
#'             data.frame(AGE = NA, 
#'                        NUM = sum(hmeq[is.na(hmeq$CLAGE), ]$BAD),
#'                        DEN = nrow(hmeq[is.na(hmeq$CLAGE), ])))
#' pool_bin(df$AGE, df$NUM, df$DEN) 

pool_bin <- function(x, num, den) {
  x_ <- x[!is.na(x)]
  n_ <- num[!is.na(x)]
  d_ <- den[!is.na(x)]
  y_ <- n_ / d_

  spc <- cor(x_, y_, method = "spearman")

  set.seed(1)
  m_ <- gbm::gbm(y ~ x1 + x2, distribution = "gaussian", data = data.frame(y = y_, x1 = x_, x2 = x_),
                 weight = d_, var.monotone = c(spc / abs(spc), spc / abs(spc)),
                 bag.fraction = 1, n.minobsinnode = round(length(x_) / 100), n.trees = 100)

  d1 <- data.frame(y = y_, x = x_, n = n_, d = d_, 
                   cat = gbm::predict.gbm(m_, n.trees = m_$n.trees, type = "response"))

  l1 <- lapply(split(d1, d1$cat),
               function(d) list(rate = abs(round(sum(d$n) / sum(d$d), 8)), maxx = max(d$x)))

  l2 <- l1[Reduce(c, lapply(l1, function(l) l$rate > 0 & l$rate < 1))]

  l3 <- sort(Reduce(c, lapply(l2, function(l) l$maxx)))[-length(l2)]

  d2 <- data.frame(x = x_, n = n_, d = d_, cut = findInterval(x_, sort(c(l3, -Inf, Inf)), left.open = T))

  d3 <- Reduce(rbind, lapply(split(d2, d2$cut),
                             function(d) data.frame(bin  = d$cut[1],
                                                    freq = sum(d$d),
                                                    miss = 0,
                                                    bads = sum(d$n),
                                                    minx = min(d$x),
                                                    maxx = max(d$x))))
  d4 <- d3[order(d3$bads / d3$freq), ]

  if (length(x[is.na(x)]) > 0) {
    m_ <- list(bin = 0, freq = sum(den[is.na(x)]), miss = sum(den[is.na(x)]), 
               bads = sum(num[is.na(x)]), minx = NA, maxx = NA)
    if (m_$bads == 0 | m_$bads == m_$freq) {
      r_ <- ifelse(m_$bads == 0, 1, nrow(d4))
      d4[r_, ]$freq <- d4[r_, ]$freq + m_$freq
      d4[r_, ]$miss <- m_$freq
      d4[r_, ]$bads <- d4[r_, ]$bads + m_$bads
    } else {
      d4 <- rbind(d4, data.frame(m_))
    }
  }

  return(list(cut = l3, tbl = gen_woe2(d4, l3)))
}
