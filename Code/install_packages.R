# Source of this script template: https://github.com/Lagkouvardos/Rhea/blob/master/install_packages.R
#' This script installs all required libraries automatically
#' Please install libraries manually if it was not possible to install an library automatically
#' Missing libraries are listed in missing_packages.txt
#' To install an library please use following two command:
#' install.packages("name of the missing library")
#' library("name of the missing library)

##################################################################################
######             Set parameters in this section manually                  ######
##################################################################################

#' Please set the directory of the present script as the working folder
#' Note: the path is denoted by forward slash "/"
source("https://bioconductor.org/biocLite.R")


######                  NO CHANGES ARE NEEDED BELOW THIS LINE               ######
##################################################################################
######                             Main Script                              ###### 
##################################################################################
###################       Load all required libraries     ########################

# Check if required packages are already installed, and install if missing
packages <- c('readxl','rhdf5','plotly','doParallel','foreach','plyr','RColorBrewer','ggplot2','pheatmap','bindrcpp','usethis','devtools',
                   'data.table','GSVA','RUVSeq','edgeR','limma','EDASeq','ShortRead','GenomicAlignments','SummarizedExperiment',
                   'DelayedArray','matrixStats','Rsamtools','GenomicRanges','GenomeInfoDb','Biostrings','XVector','IRanges',
                   'S4Vectors','BiocParallel','Biobase','BiocGenerics','dplyr','gplots','Rfast','httr')

# Function to check whether the package is installed
InsPack <- function(pack)
{
  if ((pack %in% installed.packages()) == FALSE) {
    biocLite(pack,suppressUpdates = T)
  } 
}



# Applying the installation on the list of packages
lapply(packages, InsPack)

# Make the libraries
lib <- lapply(packages, require, character.only = TRUE)

# check these. Manually install any packages present in these vectors.
not_installed <- which(lib==FALSE)
missing_packages <- lapply(not_installed, function(x) print(packages[x])) 



#################################################################################
######                           End of Script                             ######
#################################################################################                           