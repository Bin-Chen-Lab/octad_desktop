#RGES
#calculates RGES by calculating reversal gene expression from LINCS1000

####DataFrames####
#lincs_signature : dataframe for lincs pertubation data 
  #row id : landmark genes 978 of them
    #need metadata to change them to Symbol
  #col id : experiment id
    #need metadata to parse them
#lincs_sig_info : experiment info for lincs pertubation data
  #row id corresponds to lincs_signature column id
#fda_drugs : fda drug data

####required input####
#outputFolder : folder to output the drug resutls
#pipelinedataFolder : root folder with lincs input data
#dz_signature : disease signature genes dataframe
  #required columns Symbol : this must be an UPPERCASE gene symbol needed to join with RGES landmark genes
  #must change this if they require another identifier

####Changes####
#switched plyr to dplyr


#an example of running RGES and summarizing RGES across multiple profiles.
#make sure change the workspace and code directory


#code_dir <- "../code/drug/"

library("dplyr")
library("ggplot2")


####parameters####
#store these parameters somewhere if you are doing multiple runs for comparison
#default parameters

#load dz_signature
#load output folder
#load LINCS drug gene expression profiles
if (choose_fda_drugs) {
  fda_drugs = read.csv("~/Documents/GitHub/OCTAD v180104/data/raw/repurposing_drugs_20170327.txt", comment.char = "!", header=T, sep="\t")
}

'~/FileZilla Downloads/Data for Web Portal Codes offline 062618/LINCS_RGES/lincs_probe_id_info.csv'

if (landmark == 1){
  #there are a few lincs dataframes in the datafolder but for some reason only this one has gene symbols...
  load(paste0(pipelineDataFolder,'/LINCS_RGES/lincs_signatures_cmpd_landmark_symbol.RData'))
}else{
  #don't bother
  load(paste0(dataFolder,"raw/lincs_signatures_cmpd_landmark_GSE92742.RData"))
}

#write paths
output_path <- paste0(outputFolder, "/all_lincs_score.csv")
sRGES_output_path <- paste0(outputFolder, "/sRGES.csv")
sRGES_output_path_drug <- paste0(outputFolder, "/sRGES_drug.csv")
dz_sig_output_path <- paste0(outputFolder, "/dz_sig_used.csv")


lincs_sig_info <- read.csv(paste0(pipelineDataFolder,"/LINCS_RGES/lincs_sig_info.csv"))


if (choose_fda_drugs){
  lincs_sig_info <- subset(lincs_sig_info, id %in% colnames(lincs_signatures) & tolower(pert_iname) %in% tolower(fda_drugs$pert_iname))
}else{
  lincs_sig_info <- lincs_sig_info %>% filter(id %in% colnames(lincs_signatures))
}

#remove duplicate instances
lincs_sig_info <- lincs_sig_info[!duplicated(lincs_sig_info$id),]
sig.ids <- lincs_sig_info$id


####compute RGES####
gene.list <- toupper(rownames(lincs_signatures))

dz_signature <- dz_signature %>% filter(Symbol %in% gene.list)
dz_genes_up <- dz_signature %>% filter(log2FoldChange>0) %>% arrange(desc(log2FoldChange))
dz_genes_down <- dz_signature %>% filter(log2FoldChange<0) %>% arrange(log2FoldChange)
###############8


#compute RGES
#caps gene selection to max gene size 
if (nrow(dz_genes_up) > max_gene_size){
  dz_genes_up <- dz_genes_up %>% head(max_gene_size)
}
if (nrow(dz_genes_down) > max_gene_size){
  #dz_genes_down <- data.frame(GeneID=dz_genes_down[1:max_gene_size,]) %>% left_join(dz_signature,by = 'GeneID')
  dz_genes_down <- dz_genes_down %>% head(max_gene_size)
}

write.csv(rbind(dz_genes_up, dz_genes_down),  dz_sig_output_path)

dz_cmap_scores <- NULL
#count <- 0

####slow loop####
for (exp_id in sig.ids) {
#count <- count + 1
 #print(count)
  
#cmap_exp_signature is a dataframe with columns ids: whatever id we're using e.g. Symbol
 #rank which is the reverse rank of the expression of the expression
cmap_exp_signature <- data.frame(gene.list,  
                                   rank(-1 * lincs_signatures[, as.character(exp_id)], 
                                        ties.method="random"))    
  colnames(cmap_exp_signature) <- c("ids","rank")
  #runs a function cmap_score_new from drugs_core_functions.R
  dz_cmap_scores <- c(dz_cmap_scores, 
                      cmap_score_new(dz_genes_up$Symbol,dz_genes_down$Symbol,
                                     cmap_exp_signature))
}



#random scores
N_PERMUTATIONS <- 10000 #default 100000
random_sig_ids <- sample(1:ncol(lincs_signatures),N_PERMUTATIONS,replace=T)
count <- 0
random_cmap_scores <- NULL
for (expr_id in random_sig_ids){
  count <- count + 1
  #print(count)
  cmap_exp_signature <- data.frame(gene.list,  rank(-1 * lincs_signatures[, as.character(exp_id)], ties.method="random"))    
  colnames(cmap_exp_signature) <- c("ids","rank")
  
  random_input_signature_genes <- sample(gene.list, (nrow(dz_genes_up)+nrow(dz_genes_down)))
  rand_dz_gene_up <- data.frame(GeneID=random_input_signature_genes[1:nrow(dz_genes_up)])
  rand_dz_gene_down <- data.frame(GeneID=random_input_signature_genes[(nrow(dz_genes_up)+1):length(random_input_signature_genes)])
  random_cmap_scores <- c(random_cmap_scores, cmap_score_new(rand_dz_gene_up,rand_dz_gene_down,cmap_exp_signature))
}

p <- sapply(dz_cmap_scores, function(score){
  sum(random_cmap_scores < score)/length(random_cmap_scores)
})

padj <- p.adjust(p, "fdr")
results <- data.frame(id = sig.ids, cmap_score = dz_cmap_scores, p, padj)

results <- merge(results, lincs_sig_info, by = "id")
results <- results[order(results$cmap_score),]
write.csv(results, output_path)

####################
#summarize RGES
lincs_drug_prediction <- read.csv(output_path)

#should use pert_dose > 0.01
lincs_drug_prediction_subset <- subset(lincs_drug_prediction,  pert_dose > 0 & pert_time %in% c(6, 24))
#pairs that share the same drug and cell id
lincs_drug_prediction_pairs <- merge(lincs_drug_prediction_subset, lincs_drug_prediction_subset, by=c("pert_iname", "cell_id")) 
#x is the reference
lincs_drug_prediction_pairs <- subset(lincs_drug_prediction_pairs, id.x != id.y & pert_time.x == 24 & pert_dose.x == 10) #, select <- c("cmap_score.x", "cmap_score.y", "pert_dose.y", "pert_time.y"))

#difference of RGES to the reference 
lincs_drug_prediction_pairs$cmap_diff <- lincs_drug_prediction_pairs$cmap_score.x - lincs_drug_prediction_pairs$cmap_score.y
lincs_drug_prediction_pairs$dose <- round(log(lincs_drug_prediction_pairs$pert_dose.y, 2), 1)

#estimate difference
lincs_drug_prediction_pairs$dose_bin <- ifelse(lincs_drug_prediction_pairs$pert_dose.y < 10, "low", "high")
diff <- tapply(lincs_drug_prediction_pairs$cmap_diff, paste(lincs_drug_prediction_pairs$dose_bin, lincs_drug_prediction_pairs$pert_time.y), mean)

#ignore weighting cell lines
if (weight_cell_line){
  lincs_cell_line_weight <- read.csv(paste0(outputFolder, "/lincs_cell_lines_cor.csv"))
  pred <- merge(lincs_drug_prediction, lincs_cell_line_weight, by ="cell_id")
}else{
  pred <- lincs_drug_prediction
  pred$cor <- 1
}
pred$RGES <- sapply(1:nrow(pred), function(id){getsRGES(pred$cmap_score[id], pred$cor[id], pred$pert_dose[id], pred$pert_time[id], diff, max(pred$cor))})

cmpd_freq <- table(pred$pert_iname)
pred <- subset(pred, pert_iname %in% names(cmpd_freq[cmpd_freq > 0]))

pred_merged <- pred %>% 
  group_by(pert_iname) %>% 
  dplyr::summarise(
    mean = mean(RGES),
    n = length(RGES),
    median = median(RGES),
    sd = sd(RGES))
pred_merged$sRGES <- pred_merged$mean
pred_merged <- pred_merged[order(pred_merged$sRGES), ]
write.csv(pred_merged, sRGES_output_path)

#limit to FDA approved drugs
if (choose_fda_drugs){
  pred_merged_drug <- merge(fda_drugs, pred_merged, by = "pert_iname")
  pred_merged_drug <- pred_merged_drug[order(pred_merged_drug$sRGES), ]
  write.csv(pred_merged_drug, sRGES_output_path_drug)
}
