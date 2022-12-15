# =====================================================
# Overview - Packages that need to be installed
# =====================================================
# simIReff
# VineCopula
# rio (also have to rio::install_formats() to completely install rio)
# rmarkdown
# ggplot2
# ggplotify
# ggpubr
# gridExtra
# glue
# emmeans
# forcats
# moments
# devtools
# Rcpp
# stringr
# dplyr
# tidyr
# doParallel
# Hmisc
# conquer
# devEMF



# =====================================================
# Useful Information
# =====================================================
# If you have problems running my code due to the packages used,
# I have produced two .txt files in "/other" with all the versions and package names on:
# i) my Windows machine and ii) the TU Delft Cluster (CentOS-7)
# This was done via "print_all_package_versions()" below:
print_all_package_versions <- function(dir = 'output', output_path = 'all_package_versions_WINDOWS.txt') {
  dir.create(dir, recursive = TRUE)
  
  ip = as.data.frame(installed.packages()[,c(1,3:4)])
  ip = ip[is.na(ip$Priority),1:2,drop=FALSE]
  ip = data.frame(Package=sprintf('%25-s', ip$Package),
                  Version=sprintf('%25-s', ip$Version))
  write.table(ip, paste0(dir, '/', output_path), row.names=FALSE, quote=FALSE)
}
print_all_package_versions()


# =====================================================
# Windows - How to install packages
# =====================================================
install.packages('simIReff', type = 'binary')
install.packages('VineCopula', type = 'binary')
install.packages('rio', type = 'binary')
install.packages('pillar', type = 'binary')
library(rio)
install_formats(type = 'binary') # need to do this for rio to be completely installed
install.packages('rmarkdown', type = 'binary')
install.packages('ggplot2', type = 'binary')
install.packages('ggplotify', type = 'binary')
install.packages('ggpubr', type = 'binary')
install.packages('gridExtra', type = 'binary')
install.packages('glue', type = 'binary')
install.packages('emmeans', type = 'binary')
install.packages('forcats', type = 'binary')
install.packages('moments', type = 'binary')
install.packages('devtools', type = 'binary')
install.packages('Rcpp', type = 'binary')
install.packages('stringr', type = 'binary')
install.packages('dplyr', type = 'binary')
install.packages('plyr', type = 'binary')
install.packages('tidyr', type = 'binary')
install.packages('doParallel', type = 'binary')
install.packages('Hmisc', type = 'binary')
install.packages('conquer', type = 'binary')
install.packages('devEMF', type = 'binary')


# =====================================================
# TU Delft Cluster (CentOS-7) - How to install packages
# =====================================================
# Need to run this first, so that GCC is a higher version for later installs
# module use /opt/insy/modulefiles
# module avail
# module load devtoolset/7

install.packages('simIReff')
install.packages('VineCopula')
install.packages('devtools')
install.packages('rio')
library(rio)
install_formats()
# If the above line fails to install package "arrow", do the following fix:
# START FIX
Sys.setenv("NOT_CRAN" = TRUE) 
Sys.setenv("LIBARROW_BINARY" = TRUE)
Sys.setenv("LIBARROW_BUILD" = FALSE)
library(devtools)
install_version("arrow", version = "4.0.0") # for me it was downloaded from: https://mirror.lyrahosting.com/CRAN/src/contrib/Archive/arrow/arrow_4.0.0.tar.gz
# END fix
install.packages('rmarkdown')
install.packages('ggplot2')
install.packages('ggplotify')
# Do the following before installing "ggpubr", if you are getting errors with installing "ggpubr"
install_version("nloptr", version = "1.2.2.2") # for me it was downloaded from: https://mirror.lyrahosting.com/CRAN/src/contrib/Archive/nloptr/nloptr_1.2.2.2.tar.gz
install.packages('ggpubr') 
install.packages('gridExtra')
install.packages('glue')
install.packages('emmeans')
install.packages('forcats')
install.packages('moments')
install.packages('Rcpp')
install.packages('stringr')
install.packages('dplyr')
install.packages('plyr')
install.packages('tidyr')
install.packages('doParallel')
install.packages('Hmisc')
install.packages('conquer')
install.packages('devEMF')


# =====================================================
# Remember to always run all R scripts from the base directory
# =====================================================


# =====================================================
# Need to run this first to convert the Terabyte collection to "our" format
# =====================================================
source('R/utils.R')
read_terabyte2006('data/terabyte2006/', 'data/')




