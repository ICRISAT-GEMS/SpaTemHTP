% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outliers_det_boxplot.R.R
\name{outliers_det_boxplot}
\alias{outliers_det_boxplot}
\title{Sequential outlier detection using Boxplot}
\usage{
outliers_det_boxplot(data, plot = TRUE)
}
\arguments{
\item{data}{\code{data.frame} of dimension (N_genotype * N_replicate) x N_days
containing the measured phenotypic values.}

\item{plot}{\code{Logical} value indicating if the boxplot of each day should
be plotted. Default = \code{FALSE}.}
}
\value{
Return:

\code{data.frame} with outlying values put as NA.
}
\description{
Detect outliers using boxplot method (Sun and Genton, 2011). For each day,
the 25\% quantile (QR1), 75\% quantile (QR3), and 50\% inter-quantile range (IQR)
are calculated. The observations below QR1 - 1.5 x IQR or above QR3 + 1.5 x IQR
are considered as outliers. Outliers are replaced by missing value.
}
\examples{


data(SG_PH_data)

data <- outliers_det_boxplot(data = SG_PH_data[, 6:28])

}
\references{
Sun, Y. and Genton, M.G. (2011). Functional boxplots. Journal of
Computational and Graphical Statistics, 20(2), pp.316-334
}
\author{
Soumyashree Kar, Vincent Garin
}
