% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpaTemHTP_pipeline.R
\name{SpaTemHTP_pipeline}
\alias{SpaTemHTP_pipeline}
\title{Pipeline function for HTP data}
\usage{
SpaTemHTP_pipeline(
  exp_id = "exp_x",
  trait_id = "trait_i",
  out_loc = NULL,
  exp_des_data,
  pheno_data,
  raw_data = TRUE,
  raw_data_out_det = FALSE,
  raw_data_imput = FALSE,
  raw_data_out_det_imput = TRUE,
  G_BLUES_TS = TRUE,
  G_BLUES_TS_sel = TRUE,
  G_BLUES_TS_log_curve = TRUE,
  out_det = TRUE,
  miss_imp = TRUE,
  sp_adj = TRUE,
  single_mixed_model = FALSE,
  out_p_val = 0.05,
  fixed = NULL,
  random = ~row_f + col_f
)
}
\arguments{
\item{exp_id}{\code{Character} string indicating the name of the experiment.
Default = 'exp_x'.}

\item{trait_id}{\code{Character} string indicating the name of the trait
analyzed. Default = 'trait_i'}

\item{out_loc}{\code{Character} string indicating the path location where
results will be saved. Defaut is the working directory}

\item{exp_des_data}{\code{data.frame} of dimension (N_genotype * N_replicate)
x N_variable containing the experimental design information. It must include:
a) a 'genotype' column representing the line phenotyped;
b) numeric 'row' and 'col' column representing the row and column informaiton,
c) the same row and column information into factor columns named 'row_f' and
'col_f'. Other variables like replicate or block can be introduced to be used
in the spatially adjusted mixed model computation. The user must set those
extra variable in the correct format (generally factor).}

\item{pheno_data}{\code{data.frame} of dimension (N_genotype * N_replicate)
x N_days containing the measured phenotypic values.}

\item{raw_data}{\code{Logical} value specifying if the user wants the raw data
to be returned. Default = TRUE.}

\item{raw_data_out_det}{\code{Logical} value specifying if the user wants
the raw data after outlier detection to be returned. Default = FALSE.}

\item{raw_data_imput}{\code{Logical} value specifying if the user wants
the raw data after imputation to be returned. Default = FALSE.}

\item{raw_data_out_det_imput}{\code{Logical} value specifying if the user wants
the raw data after outliers detection and imputation to be returned.
Default = TRUE.}

\item{G_BLUES_TS}{\code{Logical} value specifying if the user wants
to calculate the genotypes adjusted means (BLUEs) time series after
spatial adjustment. Default = TRUE.}

\item{G_BLUES_TS_sel}{\code{Logical} value specifying if the user wants
to select the day with the largest h2 on the genotypes adjusted means (BLUEs)
time series after. Default = TRUE.}

\item{G_BLUES_TS_log_curve}{\code{Logical} value specifying if the user wants
to perform a logistic curve fitting on the G-BLUES time series Default = TRUE.}

\item{out_det}{\code{Logical} value specifying if outlier detection should
be performed on the phenotypic data. Default = TRUE.}

\item{miss_imp}{\code{Logical} value specifying if missing value imputation
should be performed on the phenotypic data. Default = TRUE.}

\item{sp_adj}{\code{Logical} value specifying if a mixed model with spatial
adjustment (SpATS model) should be used to calculate the genotype BLUEs.
Default = TRUE.}

\item{single_mixed_model}{\code{Logical} value indicating if a 'single-step'
mixed model should be calculated. See Details for more explanations.
Default = FALSE.}

\item{out_p_val}{\code{Numeric} value indicating the signficance threshold
for outliers detection. Default = 0.05.}

\item{fixed}{Optional right hand formula object specifying the fixed effects
of the SpATS model. Default = NULL.}

\item{random}{Optional right hand formula object specifying the random effects
of the SpATS model. Default = ~ row_f + col_f.}
}
\value{
For each chosen options, the function will save the produced data in
a folder created at the specified location. Will also be added.

... (develop further)
}
\description{
Pipeline function performing different level of data treatment on high
throughput phenotyping (HTP) time series data.
}
\details{
The function perform different operations to progressively enrich the data
in information content. The user can select the amount of treatment he wants
to apply on the data by selecting among the following options:

\enumerate{

\item{}{Raw data with experimetal design information.}

\item{}{Raw data with outliers detected (using \code{\link{outliers_det_boxplot}})
and experimental design information.}

\item{}{Raw data with missing values imputed after outliers detection
(using \code{\link{miss_imp_PMM}}) and experimental design information.}

\item{}{Genotype adjusted means (BLUEs) time series using the SpATS model for
spatial correction after outliers detection and imputation.}

\item{}{Selection of an optimal section or time point in the whole genotype
BLUEs time series according to an heritability criteria or change point
analysis (\code{\link{TS_select}}).}

\item{}{Further analysis of the time series fitting a logistic curve to the
genotype BLUEs TS.}

The two last options (selection on the time series and further modelling of
the time series) are conditional on the calculation of the genotype adjusted
means time series (option 4).

}
}
\examples{

data(SG_PH_data)

SG_PH_data$col_f <- factor(SG_PH_data$col)
SG_PH_data$row_f <- factor(SG_PH_data$row)

SG_PH_data$rep <- factor(SG_PH_data$rep)
SG_PH_data$block <- factor(SG_PH_data$block)

exp_des_data = SG_PH_data[, c("row", "col", "row_f", "col_f","genotype",
"rep", "block")]

pheno_data <- SG_PH_data[, 6:28]

\dontrun{

out_loc <- getwd() # specify a directory where the results will be saved

results <- SpaTemHTP_pipeline(exp_id = 'Exp_XX', trait_id = 'trait_1',
out_loc = out_loc, exp_des_data = exp_des_data, pheno_data = pheno_data,
random = ~ rep +  rep:block + row_f + col_f)

}

}
\author{
ICRISAT GEMS team
}
