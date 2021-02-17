#' Apply WoE transformations to vectors in dataframe 
#'
#' The function \code{batch_woe} applies WoE transformations to vectors 
#' in the dataframe.
#'
#' @param xs      A dataframe with numeric vectors to discretize.
#' @param bin_out A binning output from the function batch_bin(). 
#'
#' @return A dataframe with identical headers as the input xs. However, values
#'         of each variable have been transformed to WoE values.
#'
#' @examples
#' data(hmeq)
#' bin_out <- batch_bin(hmeq$BAD, hmeq[, c('DEROG', 'DELINQ')])$bin_out
#' head(batch_woe(hmeq[, c('DEROG', 'DELINQ')], bin_out))

batch_woe <- function(xs, bin_out) {

  xlst <- names(bin_out)

  woe_out <- Reduce(data.frame, 
                    lapply(xlst, function(x) cal_woe(xs[, x], bin_out[[x]])))
  
  colnames(woe_out) <- xlst
  return(woe_out)
}
