SpaTempHTP: Processing and data analysis of high throughput phenotyping data
====
  
  
## Overview
  
  SpaTempHTP is an R package with functionalities to process and perform data analyses using time series data from (outdoors) high throughput phenotyping platforms. The main task of the package is to calculate genotype best linear unbiased estimates time series (G-BLUEs). Those G-BLUEs time series can be further analysed by fitting linear or logistic curves. It is also possible perform a change point analysis of the G-BLUEs time series to determine the sections of the time series with similar behaviour and find the one where genetic variance and heritability are maximized.

## Installation


devtools::install_github("ICRISAT-GEMS/SpaTemHTP")


## Usage

Data and scripts to illustrate the functioning of SpaTempHTP can be found on the following repository https://github.com/ICRISAT-GEMS/SpaTemHTP_Validation
