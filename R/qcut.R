#' Discretizing a numeric vector
#'
#' The function \code{qcut} discretizes a numeric vector into N pieces based 
#' on quantiles.
#'
#' @param x A numeric vector.
#' @param n An integer indicating the number of categories to discretize.
#'
#' @return A numeric vector to divide the vector x into n categories.
#'
#' @examples
#' x <- 1:10
#' # [1]  1  2  3  4  5  6  7  8  9 10
#' v <- qcut(1:10, 4)
#' # [1] 3 5 8
#' findInterval(x, sort(c(v, -Inf, Inf)), left.open = TRUE)
#' # [1] 1 1 1 2 2 3 3 3 4 4

qcut <- function(x, n) {
  return(unique(stats::quantile(x, probs = seq(0, 1, 1 / n)[2:n], na.rm = TRUE, type = 1)))
}

