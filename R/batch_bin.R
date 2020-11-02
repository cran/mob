#' Apply monotonic binning to all vectors in dataframe
#'
#' The function \code{batch_bin} applies multiple binning algorithms in
#' batch to each vector in the dataframe.
#'
#' @param y      A numeric vector with 0/1 binary values.
#' @param xs     A dataframe with numeric vectors to discretize.
#' @param method A integer from 1 to 7 referring to implementations below:
#'               1. Implementation of iso_bin()  2. Implementation of qtl_bin()
#'               3. Implementation of bad_bin()  4. Implementation of rng_bin()
#'               5. Implementation of gbm_bin()  6. Implementation of kmn_bin()
#'               7. Implementation of arb_bin()
#'
#' @return A list of binning outcomes with 2 dataframes:
#'           bin_sum: A dataframe of binning summary.
#'           bin_out: A list of binning output from binning functions,
#'                    e.g. qtl_bin().
#'
#' @examples
#' data(hmeq)
#' batch_bin(hmeq$BAD, hmeq[, c('DEROG', 'DELINQ')])

batch_bin <- function(y, xs, method = 1) {

  bin_fn <- switch(method, 
                   "1" = iso_bin, "2" = qtl_bin, "3" = bad_bin, "4" = rng_bin,
                   "5" = gbm_bin, "6" = kmn_bin, "7" = arb_bin)

  bin_out <- lapply(xs, function(x) bin_fn(x, y))

  bin_sum <- Reduce(rbind, 
                    lapply(names(bin_out), 
                           function(n) data.frame(var  = n, 
                                                  nbin = nrow(bin_out[[n]]$tbl),
                                                  freq = sum(bin_out[[n]]$tbl$freq),
                                                  bads = sum(bin_out[[n]]$tbl$bads),  
                                                  miss = sum(bin_out[[n]]$tbl$miss),  
                                                  iv   = sum(bin_out[[n]]$tbl$iv),  
                                                  ks   = max(bin_out[[n]]$tbl$ks))))

  return(list(bin_sum = bin_sum, bin_out = bin_out))
}

