###############
# plot_TS_raw #
###############

#' Plot raw phenotypic time series data
#' 
#' Plot raw phenotypic time series data
#' 
#' @param data_TS \code{numeric} \code{matrix} of dimension (N_genotype x N_rep) 
#' x N_days containing the phenotype information. the different columns
#' represent the different time points (e.g. days, hours) of measurement.
#' 
#' @param geno_id \code{character} vector representing the genotype indicator
#' eventually replicated a number of time.
#' 
#' @param rep_av \code{Logical} value indicating if the value should be averaged
#' over replicates of no. Default = TRUE.
#' 
#' @param main \code{Character} string representing the title of the graph.
#' Default =  'TS'.
#'
#' @return
#' 
#' Plot of the G-BLUEs time series. The optional trend will be ploted in red.
#' The heritability will be ploted in dashed blue and time window limits in black.
#'
#' @author Vincent Garin
#' 
#' @examples
#'
#' data(SG_PH_data)
#' data_TS <- SG_PH_data[, 6:28]
#' geno_id <- SG_PH_data$genotype
#' 
#' plot_TS_raw(data_TS, geno_id, rep_av = FALSE)
#' 
#' plot_TS_raw(data_TS, geno_id, rep_av = TRUE,
#' main = 'TS averaged over replicates')
#' 
#'
#' @export
#'

plot_TS_raw <- function(data_TS, geno_id, rep_av = TRUE, main = 'TS'){
  
  
  n_days <- dim(data_TS)[2]
  geno_names <- unique(geno_id)
  n_geno <- length(geno_names)
  
  if(rep_av){
    
    g_means <- matrix(NA, n_geno, n_days)
    
    for(i in 1:n_days){
      
      g_av <- tapply(X = data_TS[, i], INDEX = geno_id,
                     FUN = function(x) mean(x, na.rm = TRUE))
      g_means[, i] <- g_av
      
    }
    
    rownames(g_means) <- names(g_av)
    
    plot_TS(data_TS = g_means, main = main)
    
  } else {

    rownames(data_TS) <- paste0('g', 1:dim(data_TS)[1])
    colnames(data_TS) <- paste0('d', 1:dim(data_TS)[2])
    
    plot_TS(data_TS = as.matrix(data_TS), main = main)
    
  }
  
  
}