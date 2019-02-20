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
7. DE_genes.csv
     - pre-filtered version of dz_signature.csv
8. drug_sensitivity_insilico_results.txt
     - Optional
     - See "..._insilico_validation.html" file/s for visualization
9. dz_down_sig_genes_enriched.csv
10. dz_sig_used.csv
     - derived from dz_signature.csv
     - genes in dz_signature.csv overlapping with LINCS L1000
     - These are the genes actually used to compute sRGES.
11. dz_signature.csv
     - All DE genes from DE_genes.csv filtered by user-determined hyperparameters.
12. dz_up_sig_genes_enriched.csv
13. highExpGenes.csv
     - Support file for gene enrichment analysis
14. lincs_reverse_expression.pdf
15. parameters.txt
     - log of hyperparameters user defined
16. res_geneinfo.csv
17. sRGES.csv
     - results of pipeline. More negative sRGES is predicted to be more efficacious.
     - n is the number of tests conducted for this compound in L1000.
18. session_info.txt
     - R config of the user
19. time_log.txt
     - Time at each step of the pipeline. Helpful for debugging.
20. SKBR3auc_insilico_validation.html (& content folder)
     - Optional
     - Graph of (cell line most similar to samples in case_id vector) SKBR3 vs pharmacogx recomputed AUC 
     - See drug_sensitivity_insilico_results.txt for cor.tests.
21. SKBR3ic50_insilico_validation.html (& content folder)
     - See 20.
     - top line vs pharmacogx recomputed IC50.
22. enrichment_analysis
     - Optional
     - Folder containing results of user-determined enrichment analyses

