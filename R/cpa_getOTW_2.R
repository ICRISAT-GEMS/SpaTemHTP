################
# cpa_getOTW_2 #
################

#' Obtain the Optimal Time Window (OTW)
#'
#' Perform Change Point Analysis (CPA) on genotype adjusted means or G-BLUEs
#' time series of a trait
#'
#' Entire temporal data set is partitioned into windows based on the differences
#' in the distribution of trait heritability and phenotypic separability using a
#' change point analysis procedure defined by Matteson et al. (2014).
#'
#'
#' @param data \code{Numeric} \code{Matrix} of dimension (N_genotype x N_days)
#' containing the G-BLUEs of a particular trait. The data colnmanes must be
#' in Date compatible format 'dd-mm-yyyy'.
#'
#' @param h2 \code{Numeric} \code{vector} of dimension (N_days)
#' containing the heritability estimates on each day for the trait
#'
#' @return
#' 
#' List with the following items:
#' 
#' \item{change_points}{\code{Numeric vector} indicating the location of the
#' change points on the time series.}
#' 
#' \item{change_points_dates}{\code{Character vector} dates of the time series
#' change points.}
#' 
#' \item{OTW_data}{\code{List} of \code{Numeric matrices} representing the
#' genotype trait values within the OTW.}
#' 
#' \item{OTW_data_opt}{\code{Numeric matrix} with the  genotype trait values
#' within the OTW with the highest heritability.}
#' 
#' @author Soumyashree Kar, Vincent Garin
#'
#' @references
#'
#' Matteson, D.S. and James, N.A. (2014). A nonparametric approach for multiple
#' change point analysis of multivariate data. Journal of the American
#' Statistical Association, 109(505), pp.334-345.
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
#' OTW <- cpa_getOTW_2(data = data, h2 = h2)
#'
#'}
#'
#' @export
#'


cpa_getOTW_2 <- function(data, h2){
  
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
  
  if(any(!IsDate(colnames(data)))){

    stop('colnames data are not all in dd-mm-yyyy Date format')

  }
  
  # h2 format
  
  if(!is.vector(h2)){
    
    stop('h2 is not a vector.')
    
  }
  
  if(!is.numeric(h2)){
    
    stop('h2 is not numeric.')
    
  }
  
  ######
  
  blueKKmeans <- kkmeans(x = as.matrix(data), centers = 3, 
                         kernel = "polydot", alg= "kkmeans", p=1,
                         na.action=na.omit)
  blue.res <- as.data.frame(centers(blueKKmeans))
  
  colnames(blue.res) <- colnames(data)
  
  clust_diff <- matrix(nrow = 1, ncol = ncol(blue.res))
  
  for(i in 1:ncol(blue.res))
  {
    # clust_diff[1,i]<-(abs(blue.res[1,i]-blue.res[2,i]) + abs(blue.res[1,i]-blue.res[3,i]) + abs(blue.res[3,i]-blue.res[2,i]))
    clust_diff[1,i]<-(sum(dist(blue.res[,i], method = "euclidean")))
  }
  
  colnames(clust_diff) <- colnames(blue.res)
  
  
  ### Start CPA for OTW identification
  
  ip.cpa<-as.data.frame(cbind(t(clust_diff), h2))
  colnames(ip.cpa)[1]<-"BLUE_CD"
  colnames(ip.cpa)[2]<-"h2"
  
  ip.cpa.ts <- xts(ip.cpa, order.by=as.Date(rownames(ip.cpa), "%d-%m-%Y"))
  
  ecp.ph <- e.cp3o(Z=ip.cpa.ts, K=4, minsize=3, alpha=1, verbose=FALSE)
  
  change_points <- ecp.ph$estimates
  change_points_dates <- rownames(ip.cpa)[change_points]
  
  # define the OTW section
  v1 <- c(1, change_points+1)
  v2 <- c(change_points, length(h2))
  
  seg <- vector(mode = 'list', length = length(v1))
  OTW <- vector(mode = 'list', length = length(seg))
  
  for(i in 1:length(v1)){
    
    seg[[i]] <- v1[i]:v2[i]
    OTW[[i]] <- data[, seg[[i]]]
    if(which.max(h2) %in% seg[[i]]){ seg_opt <- i; OTW_opt <- data[, seg[[i]]] }
    
  }
  
  res <- list(change_points = change_points,
              change_points_dates = change_points_dates, OTW_data = OTW,
              OTW_data_opt = OTW_opt)
  
  return(res)
  
}

