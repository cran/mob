#' Append missing value category
#'
#' The function \code{add_miss} appends missing value category, if any,
#' to the binning outcome. It is an utility function and is not supposed to
#' be called directly by the user.
#'
#' @param d A data frame
#' @param x The x vector used in binning functions
#' @param y The y vector used in binning functions
#'
#' @return A data frame
#' @noRd

add_miss <- function(d, x, y) {
  d1 <- d[order(d$bads / d$freq), ] 

  if (length(x[is.na(x)]) > 0) {
    m_ <- list(bin = 0, freq = length(y[is.na(x)]), miss = length(y[is.na(x)]), bads = sum(y[is.na(x)]),
               minx = NA, maxx = NA)
    if (m_$bads == 0 | m_$bads == m_$freq) {
      r_ <- ifelse(m_$bads == 0, 1, nrow(d1))
      d1[r_, ]$freq <- d1[r_, ]$freq + m_$freq
      d1[r_, ]$miss <- m_$freq
      d1[r_, ]$bads <- d1[r_, ]$bads + m_$bads
    } else {
      d1 <- rbind(d1, data.frame(m_))
    }
  }

  return(d1)
}
