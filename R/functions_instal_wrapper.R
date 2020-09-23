#' @export
octadPackagesRequirements=function(){
bioconductor_packages=c('edgeR','RUVSeq','DESeq2','limma','rhdf5','artMS')
if (length(setdiff(bioconductor_packages, rownames(installed.packages()))) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(setdiff(bioconductor_packages, rownames(installed.packages())))
}


packages=c('magrittr','dplyr','ggplot2','doParallel','foreach','lme4','Rfast','plotly','reshape2')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}}
