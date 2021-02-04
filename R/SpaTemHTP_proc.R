##################
# SpaTemHTP_proc #
##################

#' HTP data processing
#' 
#' Calculate sequentially the genotype adjusted means (genotype best unbiased
#' linear estimates, G-BLUEs).
#' 
#' The function iterates over the different time points (e.g days) of the
#' experiment and calculate the G-BLUEs according to the following model:
#' pheno = Int + row(random) + col(random) + genotype(fixed) + f(row, col) + e,
#' where f(row, col) represent the spatial surface modeled using a
#' 2-D P-spline approach as proposed by Rodgriguez-Alvarez et al. (2018).
#' 
#' The user can run a different model using the fixed and random arguments
#' that specifies the fixed and random part of the mixed model used to
#' calculate the genotypes BLUEs. The spatially adjusted model are fitted
#' using function from the SpATS package (Rodgriguez-Alvarez et al., 2018).
#' 
#' If \code{single_mixed_model = TRUE}, the function calculates a single-step mixed
#' model where outliers are iteratively removed based on the model residuals
#' and the missing values imputed during the estimation procedure. The outliers
#' are detected using a Grubb Test (p<0.05, default).
#' 
#' If \code{h2_comp = TRUE}, the function calculate the heritability from the
#' computed SpATS model. If \code{h2_comp_alt = TRUE}, if the heritability
#' computation fail with the SpATS model, the function calculates the
#' heritability using a classicle linear mixed model (pheno = rep (F) + 
#' row (R) + col (R) + genotype (R) + e). The formula used in the following
#' h2 = Sgeno / (Sgeno + (Se/n_rep)). 
#' 
#' @param exp_des_data \code{data.frame} of dimension (N_genotype * N_replicate)
#' x N_variable containing the experimental design information. It must include:
#' a) a 'genotype' column representing the line phenotyped;
#' b) numeric 'row' and 'col' column representing the row and column informaiton,
#' c) the same row and column information into factor columns named 'row_f' and
#' 'col_f'. Other variables like replicate or block can be introduced to be used
#' in the spatially adjusted mixed model computation. The user must set those
#' extra variable in the correct format (generally factor).
#' 
#' @param pheno_data \code{data.frame} of dimension (N_genotype * N_replicate)
#' x N_days containing the measured phenotypic values.
#' 
#' @param out_det \code{Logical} value specifying if outlier detection should
#' be performed on the phenotypic data. Default = TRUE.
#' 
#' @param miss_imp \code{Logical} value specifying if missing value imputation
#' should be performed on the phenotypic data. Default = TRUE.
#' 
#' @param sp_adj \code{Logical} value specifying if a mixed model with spatial
#' adjustment (SpATS model) should be used to calculate the genotype BLUEs.
#' Default = TRUE.
#' 
#' @param single_mixed_model \code{Logical} value indicating if a 'single-step'
#' mixed model should be calculated. See Details for more explanations.
#' Default = FALSE.
#' 
#' @param out_p_val \code{Numeric} value indicating the signficance threshold
#' for outliers detection. Default = 0.05.
#' 
#' @param print_day \code{Logical} value indicating if the day progression should
#' be printed. Default = TRUE.
#' 
#' @param fixed Optional right hand formula object specifying the fixed effects
#' of the SpATS model. Default = NULL.
#' 
#' @param random Optional right hand formula object specifying the random effects
#' of the SpATS model. Default = ~ row_f + col_f.
#' 
#' @param h2_comp \code{Logical} value indicating if the should calculate the daily
#' genotypic heritability. Default = TRUE.
#' 
#' @param h2_comp_alt \code{Logical} value indicating if the should calculate the daily
#' genotypic heritability. Default = TRUE.
#' 
#' @param plot \code{Logical} value specifying if a time series plot of the
#' G-BLUEs should be produced. Default = TRUE.
#' 
#' @param plot_out \code{Logical} value specifying if the boxplots of the outlier
#' detection should be saved at the specified location (\code{res_folder}).
#' Default = TRUE.
#' 
#' @param plot_miss \code{Logical} value specifying if the plot of the missing
#' values imputation should be saved at the specified location (\code{res_folder}).
#' Default = TRUE.
#' 
#' @param plot_SpATS \code{Logical} value specifying if the plot of the
#' spatial adjustment should be saved at the specified location (\code{res_folder}).
#' Default = TRUE.
#' 
#' @param res_folder \code{Character} string specifying the path where the plot
#' results will be saved. Default = NULL.
#' 
#' @return Return:
#' 
#' If sp_adj = FALSE, the code to calculate the genotype BLUEs using a mixed
#' model without spatial adjustment is not available, so we return
#' the matrix of data after eventual processing operation (outlier detection,
#' missing value imputation). Those data can be used with another software
#' (e.g. Genstat) to calculate the genotype BLUEs without spatial adjustment.
#' 
#' If sp_adj = TRUE, list containing the following objects
#' 
#' \item{G_BLUES}{\code{Matrix} of the genotype BLUEs with the genotype in row
#' and the day (or measurement time) in column.}
#' 
#' \item{G_BLUES_stdev}{\code{Matrix} of the genotype BLUEs standard deviation
#' with the genotype in row and the day (or measurement time) in column.}
#' 
#' \item{h2}{\code{Vector} with the daily genotypic heritability.}
#' 
#'
#' @author Soumyashree Kar, Vincent Garin
#' 
#' @references
#' 
#' Maria Xose Rodriguez-Alvarez, Martin P. Boer, Fred A. van Eeuwijk, Paul
#' H.C. Eilers (2018). Correcting for spatial heterogeneity in plant breeding
#' experiments with P-splines. Spatial Statistics 23 52 - 71
#' URL https://doi.org/10.1016/j.spasta.2017.10.003
#'
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
#' G_BLUEs <- SpaTemHTP_proc(exp_des_data, pheno_data = SG_PH_data[, 6:28],
#'                           out_det = TRUE, miss_imp = TRUE, sp_adj = TRUE,
#'                           random = ~ rep +  rep:block + row_f + col_f,
#'                           plot = TRUE, plot_out = FALSE, plot_miss = FALSE,
#'                           plot_SpATS = FALSE)
#' 
#' }
#'
#' @export
#'


SpaTemHTP_proc <- function(exp_des_data, pheno_data, out_det = TRUE,
                           miss_imp = TRUE, sp_adj = TRUE,
                           single_mixed_model = FALSE,
                           out_p_val = 0.05, print_day = TRUE,
                           fixed = NULL, random = ~ row_f + col_f,
                           h2_comp = TRUE, h2_comp_alt = TRUE, plot = TRUE,
                           plot_out = TRUE, plot_miss = TRUE,
                           plot_SpATS = TRUE, res_folder = NULL) {
  
  # Check if specified column were included in the experimental design
  # with the right format.
  ################
  
  if(!('genotype' %in% colnames(exp_des_data))){
    
    stop('There is no column labelled genotype in exp_des_data')
  }
  
  if(!('col' %in% colnames(exp_des_data))){
    
    stop('There is no column labelled col in exp_des_data')
  }
  
  if(!is.numeric(exp_des_data$col)){
    
    stop('The col information in exp_des_data must be numeric')
  }
  
  if(!('row' %in% colnames(exp_des_data))){
    
    stop('There is no column labelled row in exp_des_data')
  }
  
  if(!is.numeric(exp_des_data$row)){
    
    stop('The row information in exp_des_data must be numeric')
  }
  
  if(!('col_f' %in% colnames(exp_des_data))){
    
    stop('There is no column labelled col_f in exp_des_data')
  }
  
  if(!is.factor(exp_des_data$col_f)){
    
    stop('The col_f information in exp_des_data must be factor')
  }
  
  if(!('row_f' %in% colnames(exp_des_data))){
    
    stop('There is no column labelled row_f in exp_des_data')
  }
  
  if(!is.factor(exp_des_data$row_f)){
    
    stop('The row_f information in exp_des_data must be factor')
  }
  
  if(plot_out | plot_miss | plot_SpATS){
    
    if(is.null(res_folder)){
      
      stop('No folder path was specified for argument res_folder.')
      
    }
    
    if(file.exists(res_folder)){
      
      res_dir <- file.path(res_folder, 'plots')
      dir.create(res_dir)
      
      
    } else { stop('The path specified in res_fold argument is not valid.') }
    
  }
  
  ##############
  
  # Keep day identifiers and genotype names
  day_id <- colnames(pheno_data)
  n_days <- length(day_id)
  
  geno_id <- unique(exp_des_data$genotype)
  geno_id <- as.character(geno_id)
  n_geno <- length(geno_id)
  n_rep <- length(unique(exp_des_data$rep))
  
  # Space to store the results
  
  G_BLUES_mat <- matrix(NA, nrow = n_geno, ncol = n_days)
  G_BLUES_stdev_mat <- matrix(NA, nrow = n_geno, ncol = n_days)
  h2 <- rep(NA, n_days)
  
  if(!single_mixed_model){ # Regular pipeline procedure
    
    # Outliers detection
    if(out_det) {
      
      if(plot_out){
        
        pdf(file.path(res_dir, 'plot_outlier.pdf'))
        pheno_data <- outliers_det_boxplot(data = pheno_data, plot = TRUE)
        dev.off()
        
      } else{
        pheno_data <- outliers_det_boxplot(data = pheno_data, plot = FALSE)
      }
      
    }
    
    # Missing values imputation
    if(miss_imp){
      
      if(plot_miss){
        
        pdf(file.path(res_dir, 'plot_missing.pdf'))
        pheno_data <- miss_imp_PMM(data = pheno_data, plot = TRUE)
        dev.off()
        
      } else{
        pheno_data <- miss_imp_PMM(data = pheno_data, plot = FALSE)
      }
      
    }
    
    # Genotype BLUEs computation
    
    if(!sp_adj){
      
      # Mixed model without spatial adjustment not available now
      
      data <- cbind.data.frame(exp_des_data, pheno_data)
      
      return(data)
      
    } else {
      
      # SpATS model with spatial adjustment
      
      # transform variable into factor
      
      exp_des_data$genotype <- factor(exp_des_data$genotype)
      
      if(!is.null(fixed)){fixed <- as.formula(fixed)}
      
      day_ind <- 1
      h2_lm <- rep(FALSE, n_days)
      
      if(plot_SpATS){pdf(file.path(res_dir, 'plot_SpATS.pdf'))}
      
      for(i in 1:dim(pheno_data)[2]){
        
        data <- cbind.data.frame(exp_des_data, pheno_data[, i])
        colnames(data)[dim(data)[2]] <- 'pheno'
        
        m <- tryCatch(SpATS(response = "pheno", genotype = "genotype",
                            geno.decomp = NULL, genotype.as.random = FALSE,
                            spatial = ~PSANOVA(col, row, nseg = c(20,20)),
                            fixed = fixed,
                            random = as.formula(random),
                            data = data,
                            control = list(maxit = 50, tolerance = 1e-06, monitoring = 1)),
                      error = function(e) NULL)
        
        if(!is.null(m)){
          
          pred <- predict(m, which = 'genotype')
          BLUE_i <- pred$predicted.values
          BLUE_stdev_i <- pred$standard.errors
          names(BLUE_i) <- as.character(pred$genotype)
          names(BLUE_stdev_i) <- as.character(pred$genotype)
          
          BLUE_i <- BLUE_i[geno_id]
          G_BLUES_mat[, i] <- BLUE_i
          
          BLUE_stdev_i <- BLUE_stdev_i[geno_id]
          G_BLUES_stdev_mat[, i] <- BLUE_stdev_i
          
          if(plot_SpATS){plot(m)}
          
        }
        
        if(h2_comp){ # option to calculate heritability
          
          m_h <- tryCatch(SpATS(response = "pheno", genotype = "genotype",
                                geno.decomp = NULL, genotype.as.random = TRUE,
                                spatial = ~PSANOVA(col, row, nseg = c(20,20)),
                                fixed = fixed, random = as.formula(random),
                                data = data,
                                control = list(maxit = 50, tolerance = 1e-06,
                                               monitoring = 1)),
                          error = function(e) NULL)
          
          if(!is.null(m_h)){
            
            h2[i] <- getHeritability(m_h)
            
          } else {
            
            if(h2_comp_alt){
              
              m_h <- tryCatch(lmer(formula = pheno ~  rep + (1|row_f) +
                                     (1|col_f) + (1|genotype),
                                   REML = TRUE, data = data), error = function(e) NULL)
              
              if(!is.null(m_h)){
                
                remat <- summary(m_h)
                Se <- remat$sigma^2
                V_var <- remat$varcor
                Sg <- V_var$genotype[1, 1]
                h2[i] <- (Sg)/(Sg + (Se/(n_rep)))
                h2_lm[i] <- TRUE
                
              }
              
            }
            
          }
          
        }
        
        if(print_day){print(paste('day', day_ind))}
        
        day_ind <- day_ind + 1
        
        
        
      }
      
      if(plot_SpATS){dev.off()}
      
      rownames(G_BLUES_mat) <- geno_id
      colnames(G_BLUES_mat) <- day_id
      
      rownames(G_BLUES_stdev_mat) <- geno_id
      colnames(G_BLUES_stdev_mat) <- day_id
      
      if(plot){ # plot the G-BLUEs time series
        
        dt <- data.frame(geno = rep(geno_id, n_days),
                         day = rep(1:n_days, each = n_geno),
                         trait = c(G_BLUES_mat))
        
        plot <- ggplot(data = dt, aes(x = day, y = trait, group = geno)) +
          geom_point() + geom_line(aes(group = geno))
        
        print(plot)
        
      }
      
      if(any(h2_lm)){
        n_h2_lm <- sum(h2_lm)
        warn_mess <- paste0(n_h2_lm, '/', n_days, ' h2 computation were performed using a classical linear mixed model and not the SpATS model.')
        warning(warn_mess)
      }
      
      return(list(G_BLUES = G_BLUES_mat, G_BLUES_stdev = G_BLUES_stdev_mat, h2 = h2))
      
    }
    
  } else { # Single-step mixed model
    
    
    # transform variable into factor
    
    exp_des_data$genotype <- factor(exp_des_data$genotype)
    
    if(!is.null(fixed)){fixed <- as.formula(fixed)}
    
    day_ind <- 1
    h2_lm <- rep(FALSE, n_days)
    
    if(plot_SpATS){pdf(file.path(res_dir, 'plot_SpATS.pdf'))}
    
    for(i in 1:dim(pheno_data)[2]){
      
      data <- cbind.data.frame(exp_des_data, pheno_data[, i])
      colnames(data)[dim(data)[2]] <- 'pheno'
      
      p_val <- 0
      
      while(p_val < out_p_val){
        
        # compute the mixed model
        
        m <- tryCatch(SpATS(response = "pheno", genotype = "genotype",
                            geno.decomp = NULL, genotype.as.random = FALSE,
                            spatial = ~PSANOVA(col, row, nseg = c(20,20)),
                            fixed = fixed,
                            random = as.formula(random),
                            data = data,
                            control = list(maxit = 50, tolerance = 1e-06, monitoring = 1)),
                      error = function(e) NULL)
        
        # test null
        
        if(!is.null(m)){
          
          m_resid <- m$residuals
          
          test <- grubbs.test(m_resid, type = 10, two.sided = TRUE) # Grubb test
          p_val <- test$p.value
          
          if(p_val < out_p_val){ # test if p_val is lower than ... (presence of an outlier)
            
            # remove the outlying value from the data
            
            pos_max_res <- which(abs(m_resid) == max(abs(m_resid), na.rm = TRUE))
            
            data$pheno[pos_max_res] <- NA
            
            
          }
          
        } else{
          
          break()
          
        }
        
      }
      
      
      if(!is.null(m)){
        
        pred <- predict(m, which = 'genotype')
        BLUE_i <- pred$predicted.values
        BLUE_stdev_i <- pred$standard.errors
        names(BLUE_i) <- as.character(pred$genotype)
        names(BLUE_stdev_i) <- as.character(pred$genotype)
        
        BLUE_i <- BLUE_i[geno_id]
        G_BLUES_mat[, i] <- BLUE_i
        
        BLUE_stdev_i <- BLUE_stdev_i[geno_id]
        G_BLUES_stdev_mat[, i] <- BLUE_stdev_i
        
        if(plot_SpATS){plot(m)}
        
      }
      
      if(h2_comp){ # option to calculate heritability
        
        m_h <- tryCatch(SpATS(response = "pheno", genotype = "genotype",
                              geno.decomp = NULL, genotype.as.random = TRUE,
                              spatial = ~PSANOVA(col, row, nseg = c(20,20)),
                              fixed = fixed, random = as.formula(random),
                              data = data,
                              control = list(maxit = 50, tolerance = 1e-06,
                                             monitoring = 1)),
                        error = function(e) NULL)
        
        if(!is.null(m_h)){
          
          h2[i] <- getHeritability(m_h)
          
        } else {
          
          if(h2_comp_alt){
            
            m_h <- tryCatch(lmer(formula = pheno ~  rep + (1|row_f) +
                                   (1|col_f) + (1|genotype),
                                 REML = TRUE, data = data), error = function(e) NULL)
            
            if(!is.null(m_h)){
              
              remat <- summary(m_h)
              Se <- remat$sigma^2
              V_var <- remat$varcor
              Sg <- V_var$genotype[1, 1]
              h2[i] <- (Sg)/(Sg + (Se/(n_rep)))
              h2_lm[i] <- TRUE
              
            }
            
          }
          
        }
        
      }
      
      if(print_day){print(paste('day', day_ind))}
      
      day_ind <- day_ind + 1
      
    }
    
    if(plot_SpATS){dev.off()}
    
    rownames(G_BLUES_mat) <- geno_id
    colnames(G_BLUES_mat) <- day_id
    
    rownames(G_BLUES_stdev_mat) <- geno_id
    colnames(G_BLUES_stdev_mat) <- day_id
    
    if(plot){ # plot the G-BLUEs time series
      
      dt <- data.frame(geno = rep(geno_id, n_days),
                       day = rep(1:n_days, each = n_geno),
                       trait = c(G_BLUES_mat))
      
      plot <- ggplot(data = dt, aes(x = day, y = trait, group = geno)) +
        geom_point() + geom_line(aes(group = geno))
      
      print(plot)
      
    }
    
    if(any(h2_lm)){
      n_h2_lm <- sum(h2_lm)
      warn_mess <- paste0(n_h2_lm, '/', n_days, ' h2 computation were performed using a classical linear mixed model and not the SpATS model.')
      warning(warn_mess)
    }
    
    return(list(G_BLUES = G_BLUES_mat, G_BLUES_stdev = G_BLUES_stdev_mat, h2 = h2))
    
  }
  
}