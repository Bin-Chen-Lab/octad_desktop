# Web version: 
http://octad.org/

# How to Install
Before library installation install required Bioconductor and CRAN packages through this code:
```r
bioconductor_packages=c('edgeR','RUVSeq','DESeq2','limma','rhdf5','artMS')

#For R version 3.5> use BiocManager to install required bioconductor packages: 
if (length(setdiff(bioconductor_packages, rownames(installed.packages()))) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(setdiff(bioconductor_packages, rownames(installed.packages())))
}

#For R version <3.5 use the BiocInstaller to install required bioconductor packages: 
source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite(bioconductor_packages)

packages=c('magrittr','dplyr','ggplot2','doParallel','foreach','lme4','Rfast','httr')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
```

Next, install the octad.db, package with all required files for computation acalilable via link [octad.db](https://chenlab-data-public.s3.amazonaws.com/octad/octad.db_0.99.0.tar.gz%3Fdl%3D0)
```
install.packages("octad.db_0.99.0.tar", repos = NULL, type="source")
```
Finally, install the package:
```
devtools::install_github('Bin-Chen-Lab/octad_desktop')
```
It takes a few minutes to install the package and verify files. Afterward, the pipeline will be ready to run. 

# Additional data and tutorials
<li>By default, octad package uses expression data for 978 genes from the LINCS dataset. However, it can influence the result and we advice using whole octad database. To obtatin whole results for DE, downloading of the additional OCTAD database [octad.counts.and.tpm.h5](https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tpm.h5) from the AWS link is required.</li> 
<li>The tutorial available via following link:
[Tutorial](https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad_tutorial.pdf) </li> 

# Usage and examples
The several examples listed in the file [octad_example.R](https://github.com/Bin-Chen-Lab/octad_desktop/blob/master/octad_example.R) :
<li>Example 1. liver hepatocellular carcinoma vs adjacent reference tissues;</li> 
<li>Example 2. breast cancer invasive carcinoma with PIK3 mutation vs reference tissues;</li> 
<li>Example 3. lung adenocarcinoma with amplified MYC gene vs reference tissues;</li> 
<li>Example 4. Primary breast cancer invasive carcinoma vs metastatic breast cancer invasive carcinoma using only 978 genes expression data for DE;</li> 
<li>Example 5. Compute sRGES score using GEO obtained dataset</li> 



# Contacts and citation
If you use our work, please cite us: https://www.biorxiv.org/content/10.1101/821546v1
OCTAD pipeline was developed by Bin Chen laboratory.
Examples and questions can be addressed to Evgenii Chekalin, PhD, chekali1@msu.edu or Bin Chen, PhD, PI, bin.chen@hc.msu.edu


Pipeline for Open Cancer TherApeutic Discovery. Based on paper by Bin Chen, PhD using public data to repurpose drugs for Liver cancers.
https://www.nature.com/articles/ncomms16022
https://www.gastrojournal.org/article/S0016-5085(17)30264-0/abstract
http://www.dahshu.org/events/JournalClub/BinChen_Sep26_BigDataAnalytics.pdf

