# Web version: http://octad.org/

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
By default, octad package uses expression data for 978 genes from the LINCS datased. However, it can influence the result and we advice using whole octad database. Whole octad database contains expression information for over 60 000 transcrips for all samples from the octad dataset. The dataset is available via link:
[Whole octad dataset](https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tpm.h5)
The tutorial available via following link:
[Octad tutorial](https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad_tutorial.pdf)

# Example
The several examples listed in the file octad_example.R:
<li>Example 1. liver hepatocellular carcinoma vs adjacent reference tissues;</li> 
<li>Example 2. breast cancer invasive carcinoma with PIK3 mutation vs reference tissues;</li> 
<li>Example 3. lung adenocarcinoma with amplified MYC gene vs reference tissues;</li> 
<li>Example 4. Primary breast cancer invasive carcinoma vs metastatic breast cancer invasive carcinoma using only 978 genes expression data for DE;</li> 
<li>Example 5. Compute sRGES score using GEO obtained dataset</li> 











# Overview
Several typical usages of the pipeline are sourced in the octad_examples.R, where users can run drug prediction for hepatocellular carcinoma, PI3KCA mutated breast invasive carcinoma, lung adenocarcinoma with MYC amplification and comparison of the metastatic vs primary breast invasive carcinoma.

To obtatin whole results for DE, downloading of the additional OCTAD database octad.counts.and.tpm.h5 from the AWS link is required:

The dataset and code can be also accessed by the link: https://www.synapse.org/#!Synapse:syn22101254/files/
https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tpm.h5

octad.counts.and.tpm.h5 contains expression of all transcripts, while octad.LINCS.h5, used as a default in this package, only contains expression of 978 genes profiled in the LINCS database.
















Web version: http://octad.org/

# Requisite data files can be downloaded here:

In the parent folder for this pipeline, create a folder named "data". Place all files from this download into this data folder.
to download the files, using the url prefix (https://s3-us-west-2.amazonaws.com/chenlab-data-public/octad) followed by file name. For example, https://s3-us-west-2.amazonaws.com/chenlab-data-public/octad/metadata.RData
Files: 
<li>metadata.RData</li>
<li>CCLE_OCTAD.RData</li>
<li>encoderDF_AEmodel1.RData</li>
<li>random_gsea_score.RData</li>
<li>cmpd_sets_chembl_targets.RData</li>
<li>cmpd_sets_ChemCluster_targets.RData</li>
<li>cmpd_sets_chembl_targets.RData</li>
<li>cmpd_sets_mesh.RData</li>
<li>cmpd_sets_sea_targets.RData</li>
<li>repurposing_drugs_20170327.csv </li>
<li>octad.h5</li>
<li>lincs_sig_info.csv</li>
<li>lincs_signatures_cmpd_landmark.RData</li>
<li>lincs_signatures_cmpd_landmark_debug.RData</li>
<li>lincs_signatures_cmpd_landmark_symbol.RData</li>
<li>all_TCGA_CNAdata.RData</li>
</ul>

# Navigating this repositorty:
<ul>
  <li>Code: In this folder you will find all necessary R code to run the desktop version of the drug discovery pipeline</li>
  <li>example_results: Example results from breast_cancer_lumA.Rmd (found in /Code/).</li>
  <li>deprecated: Old code and files kept for internal reference.</li>
</ul>

# Contacts and citation
If you use our work, please cite us: https://www.biorxiv.org/content/10.1101/821546v1
OCTAD pipeline was developed by Bin Chen laboratory.


Pipeline for Open Cancer TherApeutic Discovery. Based on paper by Bin Chen, PhD using public data to repurpose drugs for Liver cancers.
https://www.nature.com/articles/ncomms16022
https://www.gastrojournal.org/article/S0016-5085(17)30264-0/abstract
http://www.dahshu.org/events/JournalClub/BinChen_Sep26_BigDataAnalytics.pdf




