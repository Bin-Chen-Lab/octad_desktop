# HCC case (goal: reproduce the published work, and validate target enrichment method)
#target enrichment validation using dependency map
#run target enrichment of HCC predictions and correlate with crispr knock-out data
require(GSVA)
require(limma)

setwd("~/Documents/stanford/tumor_cell_line/pipeline/octad_case/")
################
#since the number of samples, the way to compute control samples, and the method to compute DE genes are different from the published work, it is critical to evaluate the new prediction.

sRGES_published = read.csv("dataset/LIHC_rges_ic50_normalized.csv")
sRGES_octad = read.csv("dataset/octad_sRGES.csv") #prediction from octad_portal

sRGES_compare = merge(sRGES_published, sRGES_octad, by = "pert_iname")

cor(sRGES_compare$standard_value, sRGES_compare$sRGES.x, method = "spearman") #0.6
cor(sRGES_compare$standard_value, sRGES_compare$sRGES.y, method = "spearman") #0.48
cor(sRGES_compare$sRGES.x, sRGES_compare$sRGES.y)
plot(sRGES_compare$sRGES.x, sRGES_compare$sRGES.y)

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
#############

drug_pred = sRGES_octad #read.csv("dataset/LIHC_sRGES_published.csv") #the published data seems the best so far
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
table(ccle_sample_info$disease)
disease_name = "liver"
disease_celllines = as.character(ccle_sample_info$DepMap_ID[ccle_sample_info$disease %in% disease_name])

gene_effect_subset = gene_effect[disease_celllines, ]
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

target_cmpd_count = sapply(1:length(cmpdSets), function(x) length(cmpdSets[[x]]))
gsea_gene_effect = merge(gsea_gene_effect, data.frame(target = names(cmpdSets), target_cmpd_count), by = "target")
gsea_gene_effect = gsea_gene_effect[gsea_gene_effect$target_cmpd_count > 5, ]
t.test(gsea_gene_effect$score[gsea_gene_effect$z_score < z_score_cutoff], gsea_gene_effect$score[gsea_gene_effect$z_score > z_score_cutoff])


#pdf("~/Documents/stanford/tumor_cell_line/pipeline/doc/enrich_targets_LIHC.pdf")
#a = rbind(
#  data.frame(score = gsea_gene_effect$score[gsea_gene_effect$z_score > z_score_cutoff], group="Enriched Targets"),
#  data.frame(score = gsea_gene_effect$score[gsea_gene_effect$z_score < z_score_cutoff], group="Non-Enriched Targets"))
#createGGPlot(a, measureVar="score", groupVars = "group", main = "LIHC", ylab = "Gene Effect Score", xlab = "", "two.sided")
#dev.off()
