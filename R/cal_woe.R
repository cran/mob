#' Perform WoE transformation of a numeric variable 
#'
#' The function \code{cal_woe} performs the WoE transformation of a numeric variable based on the output 
#' specification from a binning function, e.g. qtl_bin() or iso_bin().
#'
#' @param data  A input dataframe
#' @param xname The name string of X with numeric values to which the WoE is applied
#' @param spec  The output table from the binning function, e.g. qtl_bin() or iso_bin()
#'
#' @return A list of WoE transformation outputs, including a dataframe with the transformed variable and a PSI summary
#'
#' @examples
#' data(hmeq)
#' bin_out <- qtl_bin(hmeq, BAD, DEROG)
#' cal_woe(hmeq, "DEROG", bin_out$df)

cal_woe <- function(data, xname, spec) {
  wname <- paste("woe", xname, sep = ".")
  calc <- function(i) {
    s <- spec[i, ]
    if (length(with(data, which(eval(parse(text = gsub("$X", xname, s$rule, fixed = T)))))) == 0) {
      return()
    } else {
      d <- data[with(data, which(eval(parse(text = gsub("$X", xname, s$rule, fixed = T))))), ]
      d$woe.bin <- s$bin
      return(within(d, assign(wname, s$woe)))
    }
  }
  df1 <- Reduce(rbind, Map(calc, seq(nrow(spec))))
  sm1 <- Reduce(rbind,
           Map(function(x) data.frame(bin      = unique(x$woe.bin), 
                                      cal_freq = nrow(x),
                                      cal_dist = round(nrow(x) / nrow(df1), 4), 
                                      cal_woe  = mean(x[[wname]])),
             split(df1, df1$woe.bin, drop = T)))
  sm2 <- merge(spec[, c("bin", "rule", "dist", "woe")], sm1, by = c("bin"), all = T)
  sm2$psi <- round((sm2$cal_dist - sm2$dist) * log(sm2$cal_dist / sm2$dist), 4)
  woe_out <- list(df = df1, psi = sm2)
  class(woe_out) <- "psi"
  return(woe_out)
}

print.psi <- function(x) {
  print(x$psi)
}
