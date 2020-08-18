################
# miss_imp_PMM #
################

#' Missing values imputation using PMM
#' 
#' Impute the missing values using predictive mean matching
#' 
#' Missing values are imputed sequentially from the first to the
#' last day.
#'
#' @param data \code{data.frame} of dimension (N_genotype * N_replicate) x N_days
#' containing the measured phenotypic values.
#' 
#' @return Return:
#' 
#' \code{data.frame} with missing values imputed.
#'
#' @author Vincent Garin, Soumyashree Kar
#'
#' @references
#' 
#' Rubin, D. B. (1986). Statistical matching using file concatenation with
#' adjusted weights and multiple imputations. Journal of Business & Economic
#' Statistics, 4(1), 87-94.
#'
#' @examples
#'
#'
#' data(SG_PH_data)
#' 
#' data <- outliers_det_boxplot(data = SG_PH_data[, 6:28])
#' 
#' \dontrun{
#' 
#' data <- miss_imp_PMM(data = data)
#' 
#' }
#' 
#' 
#'
#' @export
#'


miss_imp_PMM <- function(data){
  #### Perform MICE and then repeat the below steps #####
  mice.ip<-data
  mice_plot <- aggr(mice.ip, col=c('navyblue','yellow'),
                    numbers=TRUE, sortVars=FALSE,
                    labels=names(data), cex.axis=.7,
                    gap=3, ylab=c("Missing data","Pattern"))
  # Now, impute the missing values.
  imputed_Data <- mice(mice.ip, m=5, maxit = 10, method = 'pmm', seed = 500,
                       printFlag = FALSE)
  # summary(imputed_Data)
  #check imputed values
  OP<-mice::complete(imputed_Data, 5, include = FALSE)
  #mice_plot <- aggr(OP, col=c('navyblue','yellow'),
  #                  numbers=TRUE, sortVars=FALSE,
  #                  labels=names(OP), cex.axis=.7,
  #                  gap=3, ylab=c("Missing data","Pattern"))
  # DF <- OP # This is the imputed data to be used further 
  
  return(OP)
}
