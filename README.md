# octad_desktop
Creators: 
<ul>
<li>Bin Chen, PhD. MSU. contact: Bin.Chen@hc.msu.edu</li>
<li>William Zeng. MD Candidate. UCSF. contact: billy.zeng@ucsf.edu</li>
<li>Patrick Newbury, MD. contact: patrick.newbury@hc.msu.edu</li>
<li>Support: Anita Wen</li>
</ul>

Web version: octad.org

# Requisite data files can be downloaded here:
https://s3-us-west-2.amazonaws.com/chenlab-data-public/octad
In the parent folder for this pipeline, create a folder named "data". Place all files from this download into this data folder.
<ul>
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

# Citation
If you use our work, please cite us.

Pipeline for Open Cancer TherApeutic Discovery. Based on paper by Bin Chen, Phd using public data to repurpose drugs for Liver cancers.
https://www.nature.com/articles/ncomms16022
https://www.gastrojournal.org/article/S0016-5085(17)30264-0/abstract
http://www.dahshu.org/events/JournalClub/BinChen_Sep26_BigDataAnalytics.pdf


## Tasks:
- [x] Finish code updates as of 02/20/19
- [ ] Finalize manuscript
- [ ] Finalize octad.org
- [x] Precompute 10k random ssGSEA scores (17k!)
- [x] Deprecate Diff.Exp.R in favor of diffexp() in core_functions.R
- [x] Edit parameter chunk with Billy's optimization code
- [ ] Network analysis (Billy, Ben, Patrick)


