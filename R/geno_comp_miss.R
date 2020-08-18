##################
# geno_comp_miss #
##################

#' Identify the genotype with missing values
#'
#' @param d \code{data.frame} containing the measured phenotypic values.
#'
#' @return Return:
#' 
#' \code{character} vector with genotype having only missing values.
#'
#' @export
#'

geno_comp_miss <- function(d){
  
  geno_id <- unique(as.character(d$genotype))
  
  prob_geno <- c()
  
  for(k in 1:length(geno_id)){
    
    if(all(is.na(d[d$genotype == geno_id[k], 'pheno']))){
      
      prob_geno <- c(prob_geno, geno_id[k])
      
    }
    
  }
  
  return(prob_geno)
  
}