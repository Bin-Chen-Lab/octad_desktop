## Example Results

Download first for best viewing experience.

After completing the pipeline, these are the results files you will generate. Some results are optional and will not appear if not enabled in your workflow document.

1. all_lincs_score.csv
2. case_control_map.pdf
     - t-SNE visualization of cases & computed controls
3. case_normal_corMatrix.csv
4. case_normal_median_cor.csv
5. computedEmpGenes.csv
6. correlation graph.pdf
     - Visualization of computed controls' similarity to cases
7. res_geneinfo.csv
     - Pre-filtered version of dz_signature.csv
8. CellLineEval
     - Optional
     - Pipeline finds cell line most representative of your cases and compares sRGES vs IC50 and AUC50 from CTRPv2
     - drug_sensitivity_insilico_results.txt
     - See "...insilico_validation.html" file/s for visualization (first characters are selected cell line)
9. enrichment_analysis
     - Contains results of enrichment analysis
     - ChemCluster is enrichment based on compound physical structure clustering.
10. dz_sig_used.csv
     - Derived from dz_signature.csv
     - Genes in dz_signature.csv overlapping with LINCS L1000
     - These are the genes actually used to compute sRGES.
11. dz_signature.csv
     - All DE genes from DE_genes.csv filtered by user-determined hyperparameters.
12. highExpGenes.csv
     - Support file for gene enrichment analysis
13. lincs_reverse_expression.pdf
     - Based on dz_sig_used.csv
     - Case signature far left column, top few drugs are the other columns. The "perfect drug" would be the inverse of the case signature.
14. parameters.txt
     - log of user defined hyperparameters 
15. sRGES.csv
     - Results of pipeline. 
     - More negative sRGES is predicted to be more efficacious.
     - n is the number of tests conducted for this compound in L1000.
     - "mean" and "median" are not currently used.
16. session_info.txt
     - R config of the user
17. time_log.txt
     - Time at each step of the pipeline. Helpful for debugging.
