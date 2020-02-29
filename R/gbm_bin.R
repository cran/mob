#' Monotonic binning based on the generalized boosted regression model (GBM) 
#'
#' The function \code{gbm_bin} implements the monotonic binning based on the generalized boosted 
#' regression model (GBM) by calling the gbm library.
#'
#' @param data A input dataframe
#' @param y    The name of Y with 0/1 binary values
#' @param x    The name of X with numeric values
#'
#' @return A list of binning outcomes, including a list of cut points and a summary dataframe
#'
#' @examples
#' data(hmeq)
#' gbm_bin(hmeq, BAD, DEROG)

gbm_bin <- function(data, y, x) {
  yname <- deparse(substitute(y))
  xname <- deparse(substitute(x))
  df1 <- subset(data, !is.na(data[[xname]]) & data[[yname]] %in% c(0, 1), select = c(xname, yname))

  if (length(unique(df1[[xname]])) == 1) {
    stop(paste("there is only a single value in", xname), call. = F)
  } else if (length(unique(df1[[xname]])) == 2) {
    return(list(df   = manual_bin(data, yname, xname, cuts = min(unique(df1[[xname]]))),
                cuts = min(unique(df1[[xname]]))))
  } else {
    df2 <- data.frame(y = df1[[yname]], x = df1[[xname]], x2 = df1[[xname]])
    spc <- cor(df2[, 2], df2[, 1], method = "spearman", use = "complete.obs")
    mdl <- gbm::gbm(y ~ x + x2, distribution = "bernoulli", data = df2, var.monotone = c(spc / abs(spc), spc / abs(spc)), 
                  bag.fraction = 1, n.minobsinnode = round(nrow(df2) / 100))
    df3 <- data.frame(y = df2$y, x = df2$x, yhat = gbm::predict.gbm(mdl, n.trees = mdl$n.trees, type = "response"))
    df4 <- Reduce(rbind, 
             lapply(split(df3, df3$yhat), 
               function(x) data.frame(maxx = max(x$x), yavg = mean(x$y), yhat = round(mean(x$yhat), 8))))
    df5 <- df4[order(df4$maxx), ]
    h <- ifelse(df5[["yavg"]][1] %in% c(0, 1), 2, 1)
    t <- ifelse(df5[["yavg"]][nrow(df5)] %in% c(0, 1), 2, 1)
    cuts <- df5$maxx[h:max(h, (nrow(df5) - t))]
    return(list(df   = manual_bin(data, yname, xname, cuts = cuts), 
                cuts = cuts))  
  }
}
