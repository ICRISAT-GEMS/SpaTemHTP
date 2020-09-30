###########
# TS_plot #
###########

#' Plot G-BLUEs time series data
#' 
#' Plot G-BLUEs time series data.
#' 
#' The function plots (G-BLUEs) times series data organized by genotype using the
#' ggplot2 package. Each line represent a specifice genotype or individual.
#' Using argument 'trend', it is possible to provided an averaged trend over the genotypes
#' or any other trend that will be plot on the top of the time series.
#' 
#' @param data_TS \code{numeric} \code{matrix} of dimension N_genotype 
#' x N_days containing the G-BLUEs time series. Such an object can be obtained
#' with the function \code{\link{SpaTemHTP_proc}}.
#' 
#' @param trend \code{numeric} \code{matrix} of dimension n_days containing the
#' values of the trend to be ploted on the top of the time series Default = NULL.
#' 
#' @param main \code{Character} string representing the title of the graph. Default =  'G-BLUEs TS'.
#'
#' @return
#' 
#' Plot of the G-BLUEs time series
#'
#' @author Soumyashree Kar, Vincent Garin
#' 
#' @seealso
#' 
#' \code{\link{SpaTemHTP_proc}}
#'
#' @examples
#'
#' data(SG_PH_data)
#' data <- SG_PH_data
#' geno_id <- data$genotype
#' 
#' g_means <- matrix(NA, 384, 23)
#' 
#' for(i in 1:23){
#'   
#'   g_av <- tapply(X = data[, i + 5], INDEX = geno_id,
#'                  FUN = function(x) mean(x, na.rm = TRUE))
#'   g_means[, i] <- g_av
#'   
#' }
#' 
#' rownames(g_means) <- names(g_av)
#' 
#' plot_TS(data_TS = g_means, main = 'Raw data genotype means')
#' 
#'
#' @export
#'


plot_TS <- function(data_TS, trend = NULL, main = 'G-BLUEs TS'){
  
  n_days <- dim(data_TS)[2]
  
  dt <- data.frame(geno = rep(rownames(data_TS), n_days),
                   day = rep(1:n_days, each = dim(data_TS)[1]),
                   trait = c(data_TS))
  
  if(is.null(trend)){
    
    plot <- ggplot(data = dt, aes(x = day, y = trait, group = geno)) +
      geom_line(aes(group = geno)) + ggtitle(main)
    
    print(plot)
    
  } else {
    
    
    dt_trend <- data.frame(trait = trend, day = 1:n_days, geno = 'G_av')
    
    plot <- ggplot(data = dt, aes(x = day, y = trait, group = geno)) +
      geom_line(aes(group = geno)) + geom_line(data = dt_trend, color = 'red', size = 1) +
      ggtitle(main)
    
    print(plot)
    
  }
  
  
}