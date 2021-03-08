####################
# SpATS_model_comp #
####################

SpATS_model_comp <- function(spatial, fixed, random, data, genotype.as.random){
  
  if (spatial == 'SAP'){
    
    m <- tryCatch(SpATS(response = "pheno", genotype = "genotype",
                        geno.decomp = NULL,
                        genotype.as.random = genotype.as.random,
                        spatial = ~SAP(col, row, nseg = c(20,20)),
                        fixed = fixed,
                        random = as.formula(random),
                        data = data,
                        control = list(maxit = 50, tolerance = 1e-06, monitoring = 1)),
                  error = function(e) NULL)
    
  } else if (spatial == 'PSANOVA') {
    
    m <- tryCatch(SpATS(response = "pheno", genotype = "genotype",
                        geno.decomp = NULL,
                        genotype.as.random = genotype.as.random,
                        spatial = ~PSANOVA(col, row, nseg = c(20,20)),
                        fixed = fixed,
                        random = as.formula(random),
                        data = data,
                        control = list(maxit = 50, tolerance = 1e-06, monitoring = 1)),
                  error = function(e) NULL)
    
  } else if (spatial == 'SAP_if_PSANOVA_fail') {
    
    m <- tryCatch(SpATS(response = "pheno", genotype = "genotype",
                        geno.decomp = NULL,
                        genotype.as.random = genotype.as.random,
                        spatial = ~PSANOVA(col, row, nseg = c(20,20)),
                        fixed = fixed,
                        random = as.formula(random),
                        data = data,
                        control = list(maxit = 50, tolerance = 1e-06, monitoring = 1)),
                  error = function(e) NULL)
    
    if(is.null(m)){
      
      m <- tryCatch(SpATS(response = "pheno", genotype = "genotype",
                          geno.decomp = NULL,
                          genotype.as.random = genotype.as.random,
                          spatial = ~SAP(col, row, nseg = c(20,20)),
                          fixed = fixed,
                          random = as.formula(random),
                          data = data,
                          control = list(maxit = 50, tolerance = 1e-06, monitoring = 1)),
                    error = function(e) NULL)
      
    }
    
    
  } else if (spatial == 'PSANOVA_if_SAP_fail'){
    
    m <- tryCatch(SpATS(response = "pheno", genotype = "genotype",
                        geno.decomp = NULL,
                        genotype.as.random = genotype.as.random,
                        spatial = ~SAP(col, row, nseg = c(20,20)),
                        fixed = fixed,
                        random = as.formula(random),
                        data = data,
                        control = list(maxit = 50, tolerance = 1e-06, monitoring = 1)),
                  error = function(e) NULL)
    
    if(is.null(m)){
      
      m <- tryCatch(SpATS(response = "pheno", genotype = "genotype",
                          geno.decomp = NULL,
                          genotype.as.random = genotype.as.random,
                          spatial = ~PSANOVA(col, row, nseg = c(20,20)),
                          fixed = fixed,
                          random = as.formula(random),
                          data = data,
                          control = list(maxit = 50, tolerance = 1e-06, monitoring = 1)),
                    error = function(e) NULL)
      
    }
    
  }
  
  return(m)
  
}
  
  