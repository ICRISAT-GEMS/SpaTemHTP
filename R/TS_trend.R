############
# TS_trend #
############

#' Time series trend modelling
#' 
#' Fit linear and logistic trends to G-BLUEs time series.
#' 
#' The function fit a linear and/ or a logistic function to the genotype
#' best linear unbiased estimates time series obtained with the function
#' \code{\link{SpaTemHTP_proc}}. The logistic function is fitted using
#' the function drm from package drc.
#' 
#' @param data \code{numeric} \code{matrix} of dimension N_genotype 
#' x N_days containing the G-BLUEs time series. Such an object can be obtained
#' with the function \code{\link{SpaTemHTP_proc}}.
#' 
#' @param linear \code{Logical} value specifying if the linear trend should be
#' fitted. Default = TRUE.
#' 
#' @param logistic \code{Logical} value specifying if the logistic trend should be
#' fitted. Default = TRUE.
#'
#' @return Return:
#' 
#' \item{lin_res}{\code{list} containing the parameters of the G-BLUEs TS
#' linear trends.}
#' 
#' \item{lin_R2}{\code{Vector} of r squared goodness of fit statistic
#' for each G-BLUEs TS modeled with the linear trend.}
#' 
#' \item{log_res}{\code{list} containing the parameters of the G-BLUEs TS
#' linear trends.}
#' 
#' \item{log_R2}{\code{Vector} of r squared goodness of fit statistic
#' for each G-BLUEs TS modeled with the logistic trend.}
#' 
#' \item{log_val}{\code{matrix} of predicted values according to the
#' fitted logistic trend for each G-BLUEs TS.}
#' 
#'
#' @author Soumyashree Kar, Vincent Garin
#' 
#' @seealso
#' 
#' \code{\link{SpaTemHTP_proc}}
#' 
#' @references
#' 
#' Ritz, C., Baty, F., Streibig, J. C., Gerhard, D. (2015) Dose-Response
#' Analysis Using R PLOS ONE, 10(12), e0146021
#'
#' @examples
#'
#' data(SG_PH_data)
#' 
#' SG_PH_data$col_f <- factor(SG_PH_data$col)
#' SG_PH_data$row_f <- factor(SG_PH_data$row)
#' 
#' SG_PH_data$rep <- factor(SG_PH_data$rep)
#' SG_PH_data$block <- factor(SG_PH_data$block)
#' 
#' exp_des_data = SG_PH_data[, c("row", "col", "row_f", "col_f","genotype",
#' "rep", "block")]
#' 
#' \dontrun{
#' 
#' G_BLUEs <- SpaTemHTP_proc(exp_des_data, pheno_data = SG_PH_data[, 6:8],
#'                           out_det = TRUE, miss_imp = TRUE, sp_adj = TRUE,
#'                           random = ~ rep +  rep:block + row_f + col_f,
#'                           plot = TRUE)
#'                           
#' G_TS_trend <- TS_trend(G_BLUEs)
#' 
#' }
#'
#' @export
#'




TS_trend <- function(data, linear = TRUE, logistic = TRUE){
  
  # check data format
  ###################
  
  if(!is.matrix(data)){stop('The data object is not a matrix')}
  
  if(!is.numeric(data)){stop('The data matrix is not numeric')}
  
  n_geno <- dim(data)[1]
  geno_id <- rownames(data)
  days <- 1:dim(data)[2]
  n_days <- length(days)
  
  # linear trend
  
  if(linear){
    
    lin_res <- vector(mode = 'list', length = n_geno)
    lin_R2 <- rep(NA, n_geno)
    
    for(i in 1:n_geno){
      
      d_i <- data.frame(tr = data[i, ], day = days)
      m <- lm(tr ~ day, data = d_i)
      
      lin_res[[i]] <- m$coefficients 
      lin_R2[i] <- summary(m)$r.squared
      
    }
    
    names(lin_res) <- geno_id
    names(lin_R2) <- geno_id
    
  } else{
    
    lin_res <- NULL
    lin_R2 <- NULL
    
  }
  
  # logistic
  
  if(logistic){
    
    log_res <- vector(mode = 'list', length = n_geno)
    log_R2 <- rep(NA, n_geno)
    log_val <- matrix(NA, n_geno, n_days)
    
    for(i in 1:n_geno){
      
      d_i <- data.frame(tr = data[i, ], day = days)
      ml <- drm(tr ~ day, data = d_i, fct = L.4(), type = "continuous")
      
      log_res[[i]] <- ml$fit$par
      names(log_res[[i]]) <- c('b', 'c', 'd', 'e')
      
      log_R2[i] <- tryCatch(cor(ml$predres[, 1], d_i$tr)^2,
                            error = function(e) NA)
      
      if(length(ml$predres[, 1]) == n_days){log_val[i, ] <- ml$predres[, 1]}
      
    }
    
    names(lin_res) <- geno_id
    names(lin_R2) <- geno_id
    rownames(log_val) <- geno_id
    
  } else {
    
    log_res <- NULL
    log_R2 <- NULL
    log_val <- NULL
    
  }
  
  
  return(list(lin_res = lin_res, lin_R2 = lin_R2, log_res = log_res,
              log_R2 = log_R2, log_val = log_val))
  
}