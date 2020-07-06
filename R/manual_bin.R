#' Binning by the cut points
#'
#' The function \code{manual_bin} discretizes the x vector and summarizes over 
#' the y vector based on the discretization result. It is an utility function
#' and is not supposed to be called directly by the user.
#'
#' @param x   A numeric vector.
#' @param y   A numeric vector with 0/1 binary values.
#' @param cut A numeric vector of cut points for discretize x.
#'
#' @return A data frame to summarize the binning outcome
#' @noRd

manual_bin <- function(x, y, cut) {
  d1 <- data.frame(x = x, y = y, cut = findInterval(x, sort(c(cut, -Inf, Inf)), left.open = T))

  return(
    Reduce(rbind, 
           lapply(split(d1, d1$cut), 
                  function(x) 
                    data.frame(bin  = x$cut[1],
                               freq = nrow(x),
                               miss = 0,
                               bads = sum(x$y),
                               minx = min(x$x),
                               maxx = max(x$x)))))
}
