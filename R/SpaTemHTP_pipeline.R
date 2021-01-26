######################
# SpaTemHTP_pipeline #
######################

#' Pipeline function for HTP data
#' 
#' Pipeline function performing different level of data treatment on high
#' throughput phenotyping (HTP) time series data.
#' 
#' The function perform different operations to progressively enrich the data
#' in information content. The user can select the amount of treatment he wants
#' to apply on the data by selecting among the following options:
#' 
#' \enumerate{
#' 
#' \item{}{Raw data with experimetal design information.}
#' 
#' \item{}{Raw data with outliers detected (using \code{\link{outliers_det_boxplot}})
#' and experimental design information.}
#' 
#' \item{}{Raw data with missing values imputed after outliers detection
#' (using \code{\link{miss_imp_PMM}}) and experimental design information.}
#' 
#' \item{}{Genotype adjusted means (BLUEs) time series using the SpATS model for
#' spatial correction after outliers detection and imputation.}
#' 
#' \item{}{Selection of an optimal section or time point in the whole genotype
#' BLUEs time series according to an heritability criteria or change point
#' analysis (\code{\link{TS_select}}).}
#' 
#' \item{}{Further analysis of the time series fitting a logistic curve to the
#' genotype BLUEs TS.}
#' 
#' The two last options (selection on the time series and further modelling of
#' the time series) are conditional on the calculation of the genotype adjusted
#' means time series (option 4).
#' 
#' }
#' 
#' @param exp_id \code{Character} string indicating the name of the experiment.
#' Default = 'exp_x'.
#' 
#' @param trait_id \code{Character} string indicating the name of the trait
#' analyzed. Default = 'trait_i'
#' 
#' @param out_loc \code{Character} string indicating the path location where
#' results will be saved. Defaut is the working directory
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
#' @param raw_data \code{Logical} value specifying if the user wants the raw data
#' to be returned. Default = TRUE.
#' 
#' @param raw_data_out_det \code{Logical} value specifying if the user wants
#' the raw data after outlier detection to be returned. Default = FALSE.
#' 
#' @param raw_data_imput \code{Logical} value specifying if the user wants
#' the raw data after imputation to be returned. Default = FALSE.
#' 
#' @param raw_data_out_det_imput \code{Logical} value specifying if the user wants
#' the raw data after outliers detection and imputation to be returned.
#' Default = TRUE.
#' 
#' @param G_BLUES_TS \code{Logical} value specifying if the user wants
#' to calculate the genotypes adjusted means (BLUEs) time series after
#' spatial adjustment. Default = TRUE.
#' 
#' @param G_BLUES_TS_sel \code{Logical} value specifying if the user wants
#' to select the day with the largest h2 on the genotypes adjusted means (BLUEs)
#' time series after. Default = TRUE.
#' 
#' @param G_BLUES_TS_log_curve \code{Logical} value specifying if the user wants
#' to perform a logistic curve fitting on the G-BLUES time series Default = TRUE.
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
#' @param fixed Optional right hand formula object specifying the fixed effects
#' of the SpATS model. Default = NULL.
#' 
#' @param random Optional right hand formula object specifying the random effects
#' of the SpATS model. Default = ~ row_f + col_f.
#' 
#' @return
#' 
#' For each chosen options, the function will save the produced data in
#' a folder created at the specified location. Will also be added.
#' 
#' ... (develop further)
#'
#' @author Vincent Garin, Subhash Degala
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
#' pheno_data <- SG_PH_data[, 6:28]
#' 
#' \dontrun{
#' 
#' out_loc <- getwd() # specify a directory where the results will be saved
#' 
#' results <- SpaTemHTP_pipeline(exp_id = 'Exp_XX', trait_id = 'trait_1',
#' out_loc = out_loc, exp_des_data = exp_des_data, pheno_data = pheno_data,
#' random = ~ rep +  rep:block + row_f + col_f)
#' 
#' }
#'
#' @export
#'

# exp_id = 'exp_x'
# trait_id = 'trait_i'
# out_loc <- 'G:/PhD/Test'
# 
# data(SG_PH_data)
# 
# SG_PH_data$col_f <- factor(SG_PH_data$col)
# SG_PH_data$row_f <- factor(SG_PH_data$row)
# 
# SG_PH_data$rep <- factor(SG_PH_data$rep)
# SG_PH_data$block <- factor(SG_PH_data$block)
# 
# exp_des_data = SG_PH_data[, c("row", "col", "row_f", "col_f","genotype",
# "rep", "block")]
# 
# pheno_data <- SG_PH_data[, 6:28]
# 
# raw_data <- TRUE
# raw_data_out_det <- FALSE
# 
# out_det = TRUE
# miss_imp = TRUE
# sp_adj = TRUE
# single_mixed_model = FALSE
# out_p_val = 0.05
# fixed = NULL
# random = ~ row_f + col_f

SpaTemHTP_pipeline <- function(exp_id = 'exp_x', trait_id = 'trait_i',
                               out_loc = NULL,
                               exp_des_data, pheno_data,
                               raw_data = TRUE,
                               raw_data_out_det = FALSE,
                               raw_data_imput = FALSE,
                               raw_data_out_det_imput = TRUE,
                               G_BLUES_TS = TRUE,
                               G_BLUES_TS_sel = TRUE,
                               G_BLUES_TS_log_curve = TRUE,
                               out_det = TRUE, miss_imp = TRUE, sp_adj = TRUE,
                               single_mixed_model = FALSE, out_p_val = 0.05,
                               fixed = NULL, random = ~ row_f + col_f) {
  
  
  ### Check data format
  
  # is there a genotype column
  
  if (!('genotype' %in% colnames(exp_des_data))){
    
    stop('There is not column called genotype in exp_des_data.')
    
  }
  
  ### Check consistency in the command
  
  if(!G_BLUES_TS & G_BLUES_TS_sel){
    
    stop('To perform the selection of the G-BLUEs time series G_BLUES_TS must be TRUE')
    
  }
  
  if(!G_BLUES_TS & G_BLUES_TS_log_curve){
    
    stop('To perform the logistic curve analysis the option G_BLUES_TS must be TRUE')
    
  }
  
  
  ### create folder to store the results
  
  if (is.null(out_loc)){
    out_loc <- getwd() } else {
      
      if(!file.exists(out_loc)){stop('The ouput file location path is misspecified')}
      
    }
  
  fold_loc <- file.path(out_loc, paste0(exp_id, '_', trait_id))
  dir.create(fold_loc)
  
  # create an empty list to store the results
  
  res_list <- list()
  i <- 1
  
  ### Raw data
  
  if(raw_data){
    
    raw_data <- data.frame(exp_des_data, pheno_data, stringsAsFactors = FALSE)
    
    # save the data in .csv file
    write.csv(x = raw_data, file = file.path(fold_loc, 'raw_data.csv'),
              row.names = FALSE)
    
    # store the data in the list
    res_list[[i]] <- raw_data
    names(res_list)[i] <- 'raw_data'
    i <- i + 1
    
    # plot data
    
    geno_id <- exp_des_data$genotype
    
    jpeg(filename = file.path(fold_loc, 'raw_data_plot.jpeg'), width = 960,
         height = 960)
    plot_TS_raw(data_TS = pheno_data, geno_id = geno_id, rep_av = FALSE,
                main = 'TS plot of the replicated genotype data')
    dev.off()
    
    
  }
  
  ### raw data out detect
  
  if(raw_data_out_det){
    
    pheno_data_out <- outliers_det_boxplot(data = pheno_data, plot = FALSE)
    
    raw_data_out <- data.frame(exp_des_data, pheno_data_out,
                               stringsAsFactors = FALSE)
    
    # save the data in .csv file
    write.csv(x = raw_data_out, file = file.path(fold_loc, 'raw_data_out_det.csv'),
              row.names = FALSE)
    
    # store the data in the list
    res_list[[i]] <- raw_data_out
    names(res_list)[i] <- 'raw_data_out'
    i <- i + 1
    
    # plot data
    
    geno_id <- exp_des_data$genotype
    
    jpeg(filename = file.path(fold_loc, 'raw_data_out_det_plot.jpeg'), width = 960,
         height = 960)
    plot_TS_raw(data_TS = pheno_data_out, geno_id = geno_id, rep_av = FALSE,
                main = 'TS plot of the replicated genotype data after outlier detection')
    dev.off()
    
    
  }
  
  ### raw data imputation
  
  if(raw_data_imput){
    
    pheno_data_imp <- miss_imp_PMM(data = pheno_data, plot = FALSE)
    
    raw_data_imp <- data.frame(exp_des_data, pheno_data_imp,
                               stringsAsFactors = FALSE)
    
    # save the data in .csv file
    write.csv(x = raw_data_imp, file = file.path(fold_loc, 'raw_data_imp.csv'),
              row.names = FALSE)
    
    # store the data in the list
    res_list[[i]] <- raw_data_imp
    names(res_list)[i] <- 'raw_data_imp'
    i <- i + 1
    
    # plot data
    
    geno_id <- exp_des_data$genotype
    
    jpeg(filename = file.path(fold_loc, 'raw_data_imp_plot.jpeg'), width = 960,
         height = 960)
    plot_TS_raw(data_TS = pheno_data_imp, geno_id = geno_id, rep_av = FALSE,
                main = 'TS plot of the replicated genotype data after imputation')
    dev.off()
    
    
  }
  
  ### raw data outliers detection and imputation
  
  if(raw_data_out_det_imput){
    
    pheno_data_out <- outliers_det_boxplot(data = pheno_data, plot = FALSE)
    pheno_data_out_imp <- miss_imp_PMM(data = pheno_data_out, plot = FALSE)
    
    raw_data_out_imp <- data.frame(exp_des_data, pheno_data_out_imp,
                               stringsAsFactors = FALSE)
    
    # save the data in .csv file
    write.csv(x = raw_data_out_imp, file = file.path(fold_loc, 'raw_data_out_imp.csv'),
              row.names = FALSE)
    
    # store the data in the list
    res_list[[i]] <- raw_data_out_imp
    names(res_list)[i] <- 'raw_data_out_imp'
    i <- i + 1
    
    # plot data
    
    geno_id <- exp_des_data$genotype
    
    jpeg(filename = file.path(fold_loc, 'raw_data_out_imp_plot.jpeg'), width = 960,
         height = 960)
    plot_TS_raw(data_TS = pheno_data_out_imp, geno_id = geno_id, rep_av = FALSE,
                main = 'TS plot of the replicated genotype data after outliers detection and imputation')
    dev.off()
    
    
  }
  
  ### G-BLUES time series computation
  
  if(G_BLUES_TS){
    
    G_BLUEs <- SpaTemHTP_proc(exp_des_data, pheno_data,
                              out_det = out_det, miss_imp = miss_imp,
                              sp_adj = sp_adj,
                              single_mixed_model = single_mixed_model,
                              out_p_val = out_p_val, fixed = fixed,
                              random = random, h2_comp = TRUE,
                              print_day = FALSE, plot = FALSE)
    
    G_BLUES_TS_data <- G_BLUEs$G_BLUES
    G_BLUES_stdev <- G_BLUEs_std
    h2 <- G_BLUEs$h2
    
    # save the data in .csv file
    write.csv(x = G_BLUES_TS_data, file = file.path(fold_loc, 'G_BLUES_TS_data.csv'))
    
    write.csv(x = G_BLUES_stdev, file = file.path(fold_loc, 'G_BLUES_stdev.csv'))
    
    write.csv(x = h2, file = file.path(fold_loc, 'heritability_values.csv'),
              row.names = FALSE)
    
    # store the data in the list
    res_list[[i]] <- G_BLUES_TS_data
    names(res_list)[i] <- 'G_BLUES_TS_data'
    i <- i + 1
    
    res_list[[i]] <- h2
    names(res_list)[i] <- 'h2'
    i <- i + 1
    
    # plot data
    
    jpeg(filename = file.path(fold_loc, 'G_BLUES_TS_data_plot.jpeg'), width = 960,
         height = 960)
    
    plot_TS(data_TS = G_BLUES_TS_data, h2 = h2,
            main = 'Genotype BLUEs TS plot with time point having the highest h2')
    
    dev.off()
    
  }
  
  ### Selection on the G-BLUEs time series
  
  if(G_BLUES_TS_sel){
    
    G_BLUES_val_h2 <- TS_select(data = G_BLUES_TS_data, h2 = h2, method =  'h2_opt')
    
    # save the data in .csv file
    write.csv(x = G_BLUES_val_h2, file = file.path(fold_loc, 'G_BLUES_values_highest_h2.csv'),
              row.names = FALSE)
    
    res_list[[i]] <- G_BLUES_val_h2
    names(res_list)[i] <- 'G_BLUES_val_h2'
    i <- i + 1
    
    
  }
  
  ### logistic curve fitting on the G-BLUES time series
  
  if(G_BLUES_TS_log_curve){
    
    trend_an <- TS_trend(data = G_BLUES_TS_data)
    G_BLUES_TS_log_val <- trend_an$log_val
    
    # save the data in .csv file
    write.csv(x = G_BLUES_TS_log_val, file = file.path(fold_loc, 'G_BLUES_TS_log_val.csv'),
              row.names = FALSE)
    
    # store the data in the list
    res_list[[i]] <- G_BLUES_TS_log_val
    names(res_list)[i] <- 'G_BLUES_TS_log_val'
    i <- i + 1
    
    # plot data
    
    jpeg(filename = file.path(fold_loc, 'G_BLUES_TS_log_val_curve_plot.jpeg'), width = 960,
         height = 960)
    
    plot_TS(data_TS = G_BLUES_TS_log_val, main = 'G-BLUEs logistic curve fit plot')
    
    dev.off()
    
    ####
    
  }
  
  save(res_list, file = file.path(fold_loc, 'res_list.RData'))
  
  return(res_list)
  
}