##############
# cpa_getOTW #
##############

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
#' @return Return:
#'
#' \code{Matrix} Genotype adjusted means within the Optimal Time Window or
#' the duration within an experiment with maximum genotypic and phenotypic diversity.
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
#' OTW <- cpa_getOTW(data = data, h2 = h2)
#'
#'}
#'
#' @export
#'

cpa_getOTW <- function(data, h2){
  
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
  
  cpaRES <- OTW <- list()
  
  cpaip <- data
  
  deflt <- 3
  
  blueKKmeans<-kkmeans(as.matrix(cpaip), deflt, 
                       kernel="polydot",alg="kkmeans",p=1, na.action=na.omit)
  blue.res<-as.data.frame(centers(blueKKmeans))
  
  colnames(blue.res)<-colnames(cpaip)
  
  clust_diff<-matrix(nrow = 1, ncol = ncol(blue.res))
  
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
  head(ip.cpa)
  
  ip.cpa.ts <- xts(ip.cpa, order.by=as.Date(rownames(ip.cpa), "%d-%m-%Y"))
  
  ecp.ph<-e.cp3o(Z=ip.cpa.ts, K=4, minsize=3, alpha=1, verbose=FALSE)
  
  E<-ecp.ph$estimates
  
  dates<-rownames(ip.cpa)
  bp.date<-xts(dates[E], order.by = as.Date(dates[E], "%d-%m-%Y"))
  
  cpaRES[[i]] <- list(ip.cpa.ts=ip.cpa.ts, bp.date=bp.date)
  
  range01 <- function(x) {(x-min(x))/(max(x) - min(x))}
  
  # Find OTW
  TWs <- length(bp.date)+1
  TWmetr <- as.data.frame(matrix(NA, nrow = 3, ncol = TWs))
  colnames(TWmetr) <- paste0("TW-", 1:TWs)
  rownames(TWmetr) <- c("medCD", "slpCD", "medH2")
  
  for(i in 1:TWs){
    
    if(i == 1){
      r.ind <- which(dates %in% bp.date[i])
      tmp.tw <- ip.cpa.ts[1:r.ind, ]
      
    } else if (i == TWs){
      r.ind <- which(dates %in% bp.date[i-1])
      tmp.tw <- ip.cpa.ts[r.ind:nrow(ip.cpa.ts), ]
      
    } else {
      r.ind1 <- which(dates %in% bp.date[i-1])
      r.ind2 <- which(dates %in% bp.date[i])
      tmp.tw <- ip.cpa.ts[(r.ind1) : (r.ind2-1), ]
      
    } # end if-else
    
    # get the median of cluster-distance
    TWmetr[1 ,i] <- round(mean(tmp.tw$BLUE_CD), 2)
    # get the slope of cluster-distance
    l.mod <- lm(tmp.tw$BLUE_CD ~ c(1:dim(tmp.tw)[1]))
    l.mod.st <- summary(l.mod)
    TWmetr[2 ,i] <- l.mod.st$coefficients[2, 1]
    # get the median of heritability
    TWmetr[3 ,i] <- round(mean(tmp.tw$h2), 2)
    
  } # end for loop
  
  TWmetr.sc <- as.data.frame(t(apply(TWmetr, 1, range01)))
  OTWid <- which.max(apply(TWmetr.sc[c(1,3), ], 2, sum))
  
  BLUEs <- data
  
  if (OTWid == 1) {
    c2 <- bp.date[(OTWid)]
    col2 <- which(colnames(BLUEs) %in% as.character(c2))
    OTW[[i]] <- BLUEs[ ,1:(col2)]
    
  } else if (OTWid == length(bp.date)) {
    c1 <- bp.date[(OTWid)]
    col1 <- which(colnames(BLUEs) %in% as.character(c1))
    OTW[[i]] <- BLUEs[ ,(col1:ncol(BLUEs))]
    
  } else {
    c1 <- bp.date[(OTWid-1)]
    c2 <- bp.date[(OTWid)] 
    
    col1 <- which(colnames(BLUEs) %in% as.character(c1))
    col2 <- which(colnames(BLUEs) %in% as.character(c2))
    
    OTW[[i]] <- BLUEs[ ,col1:(col2-1)]
    
  }
  
  return(OTW)
  
}

