# Case 2: MYC amplified lung adenocarcinoma (goal: demonstrate its potential to screen compounds/targets for a cancer subgroup)
#target enrichment validation using dependency map
require(GSVA)
require(limma)

setwd("~/Documents/stanford/tumor_cell_line/pipeline/octad_case/")

############
#the enriched targets should be ensential for cell growth, reflected by gene effect score.
###########
load("dataset/cmpd_sets_chembl_targets.RData")
load(paste0("dataset/random_gsea_score_dist.RData"))
cmpdSets = cmpd_sets$cmpd.sets
names(cmpdSets) = cmpd_sets$cmpd.set.names
#############
ccle_sample_info = read.csv("dataset//sample_info.csv", stringsAsFactors = F)
#depedency scores from BROAD
gene_effect = read.csv("dataset//Achilles_gene_effect.csv", check.names = F, row.names = 1)
gene_names = sapply(colnames(gene_effect), function(x) unlist(strsplit(x, " "))[1])

ccle_cn = read.csv('dataset/CCLE_gene_cn.csv', check.names = F, row.names = 1)
ccle_cn_genes = sapply(colnames(ccle_cn), function(x) unlist(strsplit(x, " "))[1])

#sensitivity data
load("dataset/drug_cell_matrix.RData")
#########################

#######################
#identify MYC-amplified cell lines
table(ccle_sample_info$disease)
disease_name = "lung"
disease_celllines = as.character(ccle_sample_info$DepMap_ID[ccle_sample_info$disease %in% disease_name & ccle_sample_info$disease_sutype %in% "lung_NSC"])
gene = "MYC"

ccle_mutation_subset  = ccle_mutation[ccle_mutation$Hugo_Symbol %in% gene & 
                                        ccle_mutation$DepMap_ID %in% disease_celllines,]

#cutoff should be optmized
cutoff = 1
over_expressed_cells = disease_celllines[disease_celllines %in% rownames(ccle_cn[ccle_cn[, ccle_cn_genes == gene] > cutoff ,])] #
other_expressed_cells = disease_celllines[disease_celllines %in% rownames(ccle_cn[ccle_cn[, ccle_cn_genes == gene] < cutoff ,])] #

###################
#validate drug efficacy; looks like MYC specific
#MYC amplified cells
drug_pred = read.csv("dataset/octad_myc_lung_random_sRGES.csv") #prediction from octad_portal


drug_sensitivity_subset = drug_cell_matrix[,colnames(drug_cell_matrix) %in% as.character(ccle_sample_info[ccle_sample_info$DepMap_ID %in%  over_expressed_cells, "stripped_cell_line_name"])]
drug_sensitivity_subset_ave = data.frame(apply(drug_sensitivity_subset, 1, function(x) median(x, na.rm = T)))
colnames(drug_sensitivity_subset_ave) = "AUC"
drug_sensitivity_subset_ave$drug = toupper(rownames(drug_sensitivity_subset_ave))
drug_pred$pert_iname_upper = toupper(drug_pred$pert_iname)
drug_pred_subset = merge(drug_pred, drug_sensitivity_subset_ave, by.x="pert_iname_upper", by.y="drug")
plot(drug_pred_subset$sRGES, drug_pred_subset$AUC)
cor.test(drug_pred_subset$sRGES, drug_pred_subset$AUC)

#MYC non-amplified cells
drug_sensitivity_subset = drug_cell_matrix[,colnames(drug_cell_matrix) %in% as.character(ccle_sample_info[ccle_sample_info$DepMap_ID %in%  other_expressed_cells, "stripped_cell_line_name"])]
drug_sensitivity_subset_ave = data.frame(apply(drug_sensitivity_subset, 1, function(x) median(x, na.rm = T)))
colnames(drug_sensitivity_subset_ave) = "AUC"
drug_sensitivity_subset_ave$drug = toupper(rownames(drug_sensitivity_subset_ave))
drug_pred$pert_iname_upper = toupper(drug_pred$pert_iname)
drug_pred_subset = merge(drug_pred, drug_sensitivity_subset_ave, by.x="pert_iname_upper", by.y="drug")
plot(drug_pred_subset$sRGES, drug_pred_subset$AUC)
cor.test(drug_pred_subset$sRGES, drug_pred_subset$AUC)

################
#validate using drug efficacy data from MYC-amplied cell lines
drug_pred = drug_pred[drug_pred$n > 0, ]
rgess = matrix(-1*((drug_pred$sRGES)), ncol = 1)
rownames(rgess) = tolower(drug_pred$pert_iname)
gsea_results = gsva(rgess, cmpdSets, method = "ssgsea",  parallel.sz=8,  ssgsea.norm=F)
#convert to z scores such that scores from different conditions are comparable.
#caveats: random ssGSEA scores of a few targets are not normally distributed (likely these targets with a couple of ligands)
gsea_z_score = sapply(1:nrow(gsea_results), function(i){
  (gsea_results[i,1] - random_gsea_score_dist[["chembl_targets"]][["mean"]][i])/random_gsea_score_dist[["chembl_targets"]][["sd"]][i]
})

gsea_z_score = data.frame(target = names(gsea_z_score), z_score = gsea_z_score)


#############
##in target validation, seems like not MYC-specific

gene_effect_subset = gene_effect[over_expressed_cells, ] #other_expressed_cells
gene_effect_subset_ave = apply(gene_effect_subset, 2, function(x) median(x, na.rm = T))
#if (target_type == "sea_targets"){
  #gene_effect_score = data.frame(target = paste0(gene_names, "_HUMAN"), score = gene_effect_subset_ave)
#}else{
  gene_effect_score = data.frame(target = paste0(gene_names, ""), score = gene_effect_subset_ave)
#}

gsea_gene_effect = merge(gsea_z_score, gene_effect_score, by = "target")
gsea_gene_effect = gsea_gene_effect[gsea_gene_effect$z_score > -100, ] #may filter out based on z score


#apparently, essential targets have more positive scores (not neccessary to be with z score > 1.5)
plot(gsea_gene_effect$z_score, gsea_gene_effect$score)
cor.test(gsea_gene_effect$z_score, gsea_gene_effect$score, method = "spearman")

z_score_cutoff = 1.5
t.test(gsea_gene_effect$score[gsea_gene_effect$z_score < z_score_cutoff], gsea_gene_effect$score[gsea_gene_effect$z_score > z_score_cutoff])
t.test(gsea_gene_effect$z_score[gsea_gene_effect$score < -1], gsea_gene_effect$z_score[gsea_gene_effect$score > -1])


#pdf("~/Documents/stanford/tumor_cell_line/pipeline/doc/enrich_targets_LIHC.pdf")
#a = rbind(
#  data.frame(score = gsea_gene_effect$score[gsea_gene_effect$z_score > z_score_cutoff], group="Enriched Targets"),
#  data.frame(score = gsea_gene_effect$score[gsea_gene_effect$z_score < z_score_cutoff], group="Non-Enriched Targets"))
#createGGPlot(a, measureVar="score", groupVars = "group", main = "LIHC", ylab = "Gene Effect Score", xlab = "", "two.sided")
#dev.off()
