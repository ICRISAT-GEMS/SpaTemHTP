########################
# outliers_det_boxplot #
########################

#' Sequential outlier detection using Boxplot
#'
#' Detect outliers using boxplot method (Sun and Genton, 2011). For each day,
#' the 25% quantile (QR1), 75% quantile (QR3), and 50% inter-quantile range (IQR)
#' are calculated. The observations below QR1 - 1.5 x IQR or above QR3 + 1.5 x IQR
#' are considered as outliers. Outliers are replaced by missing value.
#'
#' @param data \code{data.frame} of dimension (N_genotype * N_replicate) x N_days
#' containing the measured phenotypic values.
#' 
#' @param plot \code{Logical} value indicating if the boxplot of each day should
#' be plotted. Default = \code{FALSE}.
#' 
#' @return Return:
#' 
#' \code{data.frame} with outlying values put as NA.
#'
#' @author Soumyashree Kar, Vincent Garin
#'
#' @references
#' 
#' Sun, Y. and Genton, M.G. (2011). Functional boxplots. Journal of
#' Computational and Graphical Statistics, 20(2), pp.316-334
#'
#' @examples
#'
#'
#' data(SG_PH_data)
#' 
#' data <- outliers_det_boxplot(data = SG_PH_data[, 6:28])
#'
#'
#' @import ggplot2
#' @import mice
#' @import SpATS
#' @import VIM
#'
#' @export
#'


outliers_det_boxplot <- function(data, plot = TRUE) {
  
  op <- data
  
  for(i in 1:ncol(data))
  {
    df <- data[,i]
    ol=boxplot(data[i], plot = plot, show.names = T)$out
    which(df %in% ol)
    df[df %in% ol] <- NA
    op[,i] <- df
    
  }
  return(op)
}
