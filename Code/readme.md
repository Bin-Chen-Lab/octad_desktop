## Codes for OCTAD desktop 


All files can be found on the Chen Lab RNA-Seq server: 
/home/ubuntu/chenlab_v2/proj/PatrickN/desktop_pipeline_v0.1/

This version includes some new data files, so be sure to download those from the above directory. 
The pipeline is configured to run within the folder structure as found in the above directory.

Please contact Patrick: patrick.newbury@hc.msu.edu with any issues

### To run:
- Start a new RStudio project under the installation folder (like the above directory)
- source('install_packages.R')
- Open breast_cancer_lumA.Rmd
- Edit line 15 "base.folder"
- Knit (takes about 10-15 minutes)
- If this is able to complete, then try new combinations of case_id and normal_id (lines 154-191)
- For Linux, you'll also need to install pandoc: https://github.com/rstudio/rmarkdown/blob/master/PANDOC.md


### Changelog:
<ul>
<li> 02/20/19: <ul>
<li>Consolidated many previous source() calls into core_functions.R. Also combined previous "DE_core_functions.R" & "drugs_core_functions.R" into core_functions.R</li>
<li>Consolidated many data files into new files: CCLE_OCTAD.RData, metadata.RData</li>
<li>Compound enrichment is no longer possible without using a workflow rmd (or taking the code from the EnrichmentAnalysis block).</li>
<li>All compound enrichment analysis is disabled by default due to taking a long time. Can be enabled by editing EnrichmentAnalysis chunk option "eval" to TRUE. Can further be modulated by changing "targets" variable in the same chunk. "target_type = 'ChemCluster'" should never be changed.</li>
<li>Updated core_functions.R to properly output drug_sensitivity_insilico_results.txt. It had been putting out a blank .txt.</li>
<li>Removed support for sRGES parallelization via doParallel as it was causing different issues for different hardware. Using compiler as stand-in until we can work this out. This has not significantly impacted overall duration.</li>
</ul>
<li> 02/25/19: <ul>
<li>Added random_gsea_score.RData to dataFolder. This is 17,000 random GSEAs for each enrichment type to increase padj resolution and speed. Overall run time now ~11min with all features enabled.</li>
<li>Changed breast_cancer_lumA.Rmd workflow example chunk 'EnrichmentAnalysis'. Now Enabled by default. Now runs loop of drug_enrichment() on all enrichment types, as above bullet enabled quick computation of all 5.</li>
<li>Changed core_functinos.R's drug_enrichment() to streamline drug enrichment analysis.</li>
<li>Made new auto-generating output folders for enrichment analysis subtypes.</li>
<li>Restricted output of enrichment analysis to top 50 per subtype. SEA still outputting ~1300 significant hits.</li>
<li> 02/26/19: <ul>
<li>Miscellaneous minor changes in workflow and core_functions.R</li>
</li>
</ul>

