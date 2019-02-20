Codes for OCTAD desktop 
Updated 02/20/19

All files can be found on the Chen Lab RNA-Seq server: 
/home/ubuntu/chenlab_v2/proj/PatrickN/desktop_pipeline_v0.1/

This version includes some new data files, so be sure to download those from the above directory. 
The pipeline is configured to run within the folder structure as found in the above directory.

To run:
- Start a new RStudio project under the installation folder (like the above directory)
- source('install.packages.R')
- Open breast_cancer_lumA.Rmd
- Edit line 15 "base.folder"
- Knit (takes about 10-15 minutes)
- If this is able to complete, then try new combinations of case_id and normal_id (lines 154-191)
- Also try experimenting with compound enrichment 
line 477: ```{r EnrichmentAnalysis, eval=TRUE, message=FALSE, warning=FALSE, echo=FALSE, results='hide'}

Please contact Patrick: patrick.newbury@hc.msu.edu with any issues


Changelog:
- 02/20/19: Updated core_functions.R to properly output drug_sensitivity_insilico_results.txt. 
It had been putting out a blank .txt.
