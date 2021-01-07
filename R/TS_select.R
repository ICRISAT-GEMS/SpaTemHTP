#############
# TS_select #
#############

#' Select a part of the (G-BLUEs) time series
#'
#' Return a selection of a part of the G-BLUEs time series showing some optimal
#' properties in terms of internal consistency and heritability.
#'
#' If \code{method = 'CPA'}, the functions determine time windows using a change
#' point analysis and return the time window with the largest heritability. If
#' \code{method = 'h2_opt'}, the function return the TS data at the date where
#' the heritability is the maximum +/- \code{interval} days.
#'
#' @param data \code{Numeric} \code{Matrix} of dimension (N_genotype x N_days)
#' containing the G-BLUEs of a particular trait. The data colnmanes must be
#' in Date compatible format 'dd-mm-yyyy'.
#'
#' @param h2 \code{Numeric} \code{vector} of dimension (N_days)
#' containing the heritability estimates on each day for the trait
#' 
#' @param method \code{Character string} indicating which method to choose
#' to select the elements of the time series. Must be either 'CPA', or 'h2_opt'.
#' Default = 'CPA'.
#' 
#' @param interval \code{Numeric} value defining the number of day before and
#' after the day with the largest heritability that are returned. Default = 0.
#'
#' @return
#' 
#' List with the following items:
#' 
#' if \code{method = 'CPA'}
#' 
#' \item{TS_sel}{\code{Numeric matrix} with the  genotype trait values
#' within the OTW with the highest heritability.}
#' 
#' if \code{method = 'h2_opt'}
#' 
#' \item{TS_sel}{\code{Vector} of the genotype trait values at the date with
#' the highest heritability or \code{numeric matrix} with the  genotype trait
#' values at the day with the highest heritability with optionally +/-
#' \code{interval} days.}
#' 
#' @author Vincent Garin, Soumyashree Kar
#'
#' @examples
#'
#' data(SG_PH_data)
#' SG_PH_data$col_f <- factor(SG_PH_data$col)
#' SG_PH_data$row_f <- factor(SG_PH_data$row)
#' 
#' SG_PH_data$rep <- factor(SG_PH_data$rep)
#' SG_PH_data$block <- factor(SG_PH_data$block)
#' 
#' exp_des_data = SG_PH_data[, c("row", "col", "row_f", "col_f","genotype",
#'                               "rep", "block")]
#'
#' \dontrun{
#' 
#' op <- SpaTemHTP_proc(exp_des_data, pheno_data = SG_PH_data[, 6:28],
#'                      out_det = TRUE, miss_imp = TRUE, sp_adj = TRUE,
#'                      random = ~ rep +  rep:block + row_f + col_f,
#'                      h2_comp = TRUE, plot = TRUE)
#'                      
#' data <- op$G_BLUES
#' 
#' # make sure data colnmanes are dd-mm-yyyy Date format compatible
#' dates <- substr(colnames(data), 2, nchar(colnames(data)))
#' dates <- str_replace_all(string = dates, pattern = '\\.', replacement = '-')
#' colnames(data) <- dates
#' 
#' h2 <- op$h2
#' 
#' TS_sel <- TS_select(data = data, h2 = h2)
#'
#'}
#'
#' @export
#'


TS_select <- function(data, h2, method = 'CPA', interval = 0){
  
  
  # check the formats
  ###################
  
  # data format
  
  if(!is.matrix(data)){
    
    stop('data is not a matrix.')
    
  }
  
  if(!is.numeric(data)){
    
    stop('data is not numeric.')
    
  }
  
  # data colnames format
  
  if(method == 'CPA'){
    
    if(any(!IsDate(colnames(data)))){
      
      stop('colnames data are not all in dd-mm-yyyy Date format')
      
    }
    
  }
  
  
  
  # h2 format
  
  if(!is.vector(h2)){
    
    stop('h2 is not a vector.')
    
  }
  
  if(!is.numeric(h2)){
    
    stop('h2 is not numeric.')
    
  }
  
  if(!(method %in% c('CPA', 'h2_opt'))){
    
    stop('Method should be equal to "CPA" or "h2_opt".')
    
  }
  
  
  ######
  
  if(method == 'CPA'){
    
    TS_sel <- cpa_getOTW_2(data = data, h2 = h2)
    TS_sel <- TS_sel$OTW_data_opt
    
  } else if (method == 'h2_opt'){
    
    max_h2 <- which.max(h2)
    sel_vect <- (max_h2-interval):(max_h2+interval)
    
    TS_sel <- data[, sel_vect]
    
  }
  
  return(TS_sel)
  
}