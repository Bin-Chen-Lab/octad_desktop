cmap_score_new <- function(sig_up, sig_down, drug_signature) {
  #the old function does not support the input list with either all up genes or all down genes, this new function attempts to addess this.
  #we also modify the original CMap approach: whenever the sign of ks_up/ks_down, we substract the two scores such that the final scores would not enrich at 0.
  
  num_genes <- nrow(drug_signature)
  ks_up <- 0
  ks_down <- 0
  connectivity_score <- 0
  
  # I think we are re-ranking because the GeneID mapping changed the original rank range
  drug_signature[,"rank"] <- rank(drug_signature[,"rank"])
  
  # Merge the drug signature with the disease signature by GeneID. This becomes the V(j) from the algorithm description
  up_tags_rank <- merge(drug_signature, sig_up, by.x = "ids", by.y = 1)
  down_tags_rank <- merge(drug_signature, sig_down, by.x = "ids", by.y = 1)
  
  up_tags_position <- sort(up_tags_rank$rank)
  down_tags_position <- sort(down_tags_rank$rank)
  
  num_tags_up <- length(up_tags_position)
  num_tags_down <- length(down_tags_position)
  
  # 
  if(num_tags_up > 1) {
    a_up <- 0
    b_up <- 0
    
    a_up <- max(sapply(1:num_tags_up,function(j) {
      j/num_tags_up - up_tags_position[j]/num_genes
    }))
    b_up <- max(sapply(1:num_tags_up,function(j) {
      up_tags_position[j]/num_genes - (j-1)/num_tags_up
    }))
    
    if(a_up > b_up) {
      ks_up <- a_up
    } else {
      ks_up <- -b_up
    }
  }else{
    ks_up <- 0
  }
  
  if (num_tags_down > 1){
    
    a_down <- 0
    b_down <- 0
    
    a_down <- max(sapply(1:num_tags_down,function(j) {
      j/num_tags_down - down_tags_position[j]/num_genes
    }))
    b_down <- max(sapply(1:num_tags_down,function(j) {
      down_tags_position[j]/num_genes - (j-1)/num_tags_down
    }))
    
    if(a_down > b_down) {
      ks_down <- a_down
    } else {
      ks_down <- -b_down
    }
  }else{
    ks_down <- 0
  }
  
  if (ks_up == 0 & ks_down != 0){ #only down gene inputed
    connectivity_score <- -ks_down 
  }else if (ks_up !=0 & ks_down == 0){ #only up gene inputed
    connectivity_score <- ks_up
  }else if (sum(sign(c(ks_down,ks_up))) == 0) {
    connectivity_score <- ks_up - ks_down # different signs
  }else{
    connectivity_score <- ks_up - ks_down
  }
  
  return(connectivity_score)
}
#summary RGES function
getsRGES <- function(RGES, cor, pert_dose, pert_time, diff, max_cor){
  
  sRGES <- RGES
  pert_time <- ifelse(pert_time < 24, "short", "long")
  pert_dose <- ifelse(pert_dose < 10, "low", "high")
  if (pert_time == "short" & pert_dose == "low"){
    sRGES <- sRGES + diff[4]
  }
  if (pert_dose ==  "low" & pert_time == "long"){
    sRGES <- sRGES + diff[2]
  }
  if (pert_dose ==  "high" & pert_time == "short"){
    sRGES <- sRGES + diff[1]
  }
  return(sRGES * cor/max_cor) #
}

computeLINCSrges = function(dz_signature,choose_fda = F,parallel = F,maxGenes=900){
  require(dplyr)
  load(paste0(pipelineDataFolder,'/LINCS_RGES/lincs_signatures_cmpd_landmark_symbol.RData'))
  gene.list <- toupper(rownames(lincs_signatures))
  dz_sigUsed <- dz_signature %>% filter(toupper(Symbol) %in% gene.list)
  dz_genes_up <- dz_sigUsed %>% filter(log2FoldChange>0) %>% arrange(desc(log2FoldChange)) %>% head(maxGenes)
  dz_genes_down <- dz_sigUsed %>% filter(log2FoldChange<0) %>% arrange(log2FoldChange) %>% head(maxGenes)
  write.csv(rbind(dz_genes_up,dz_genes_down), paste0(outputFolder,'dz_sig_used.csv'))
  lincs_sig_info <- read.csv(paste0(pipelineDataFolder,"/LINCS_RGES/lincs_sig_info.csv"),
                             stringsAsFactors = F)
  lincs_sig_info <- lincs_sig_info %>% filter(id %in% colnames(lincs_signatures))
  if(choose_fda == T){
    fda_drugs = read.csv(paste0(pipelineDataFolder,"LINCS_RGES/repurposing_drugs_20170327.csv"),
                         comment.char = "!", header=T, sep="\t")
    lincs_sig_info <- lincs_sig_info %>% filter(id %in% colnames(lincs_signatures) & 
                                                  tolower(pert_iname) %in% tolower(fda_drugs$pert_iname))
  }
  lincs_sig_info <- lincs_sig_info[!duplicated(lincs_sig_info$id),]
  sig.ids <- lincs_sig_info$id
  if(parallel == T){
    require(doParallel)
    registerDoParallel(cores=4)
    dz_cmap_scores = foreach(exp_id = sig.ids,.combine = 'c')%dopar%{
      cmap_exp_signature <- data.frame(gene.list,  
                                       rank(-1 * lincs_signatures[, as.character(exp_id)], 
                                            ties.method="random"))    
      colnames(cmap_exp_signature) <- c("ids","rank")
      #runs a function cmap_score_new from drugs_core_functions.R
      cmap_score_new(dz_genes_up$Symbol,dz_genes_down$Symbol,
                     cmap_exp_signature)
      
    }    
  }else{
    dz_cmap_scores <- NULL
    for (exp_id in sig.ids){
      cmap_exp_signature <- data.frame(gene.list,  
                                       rank(-1 * lincs_signatures[, as.character(exp_id)], 
                                            ties.method="random"))    
      colnames(cmap_exp_signature) <- c("ids","rank")
      #runs a function cmap_score_new from drugs_core_functions.R
      dz_cmap_scores <- c(dz_cmap_scores, 
                          cmap_score_new(dz_genes_up$Symbol,dz_genes_down$Symbol,
                                         cmap_exp_signature))
    }
  }
  results <- data.frame(id = sig.ids, cmap_score = dz_cmap_scores)
  results <- merge(results, lincs_sig_info, by = "id")
  results <- results[order(results$cmap_score),] 
  return(results)
}



summarizeLincsRGES = function(lincs_rges,
                              weight_cell_line = F,
                              cell_lines = '',
                              choose_fda = F,
                              parallel = F){
  lincs_drug_prediction = lincs_rges
  lincs_drug_prediction_subset <- subset(lincs_drug_prediction,  pert_dose > 0 & pert_time %in% c(6, 24))
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
    lincs_cell_line_weight <- read.csv(paste0(pipelineDataFolder, "/lincs_cell_lines_cor.csv"))
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
      sRGES = mean(RGES),
      n = length(RGES),
      medRGES = median(RGES),
      sd = sd(RGES))
  pred_merged <- pred_merged[order(pred_merged$sRGES), ]
  return(pred_merged)
}

computeRandomLincsRGES = function(dz_signature,choose_fda = T,parallel = T,n_perm=10000){
  
  load(paste0(pipelineDataFolder,'/LINCS_RGES/lincs_signatures_cmpd_landmark_symbol.RData'))
  gene.list <- toupper(rownames(lincs_signatures))
  dz_sigUsed <- dz_signature %>% filter(toupper(Symbol) %in% gene.list)
  write.csv(dz_sigUsed, paste0(outputFolder,'dz_sig_used.csv'))
  dz_genes_up <- dz_sigUsed %>% filter(log2FoldChange>0) %>% arrange(desc(log2FoldChange))
  dz_genes_down <- dz_sigUsed %>% filter(log2FoldChange<0) %>% arrange(log2FoldChange)
  
  lincs_sig_info <- read.csv(paste0(pipelineDataFolder,"/LINCS_RGES/lincs_sig_info.csv"),
                             stringsAsFactors = F)
  if(choose_fda == T){
    fda_drugs = read.csv(paste0(pipelineDataFolder,"LINCS_RGES/repurposing_drugs_20170327.csv"),
                         comment.char = "!", header=T, sep="\t")
    lincs_sig_info <- lincs_sig_info %>% filter(id %in% colnames(lincs_signatures) & 
                                                  tolower(pert_iname) %in% tolower(fda_drugs$pert_iname))
  }
  lincs_sig_info <- lincs_sig_info[!duplicated(lincs_sig_info$id),]
  sig.ids <- lincs_sig_info$id
  
  
  if(parallel == T){
    require(doParallel)
    registerDoParallel(cores=4)
    N_PERMUTATIONS <- n_perm #default 100000
    random_sig_ids <- sample(colnames(lincs_signatures),N_PERMUTATIONS,replace=T)
    #random_cmap_scores <- NULL
    random_cmap_scores = foreach(exp_id = random_sig_ids,.combine = 'c')%dopar%{
      cmap_exp_signature <- data.frame(gene.list,  rank(-1 * lincs_signatures[, as.character(exp_id)], ties.method="random"))
      colnames(cmap_exp_signature) <- c("ids","rank")
      
      random_input_signature_genes <- sample(gene.list, (nrow(dz_genes_up)+nrow(dz_genes_down)))
      rand_dz_gene_up <- data.frame(GeneID=random_input_signature_genes[1:nrow(dz_genes_up)])
      rand_dz_gene_down <- data.frame(GeneID=random_input_signature_genes[(nrow(dz_genes_up)+1):length(random_input_signature_genes)])
      cmap_score_new(rand_dz_gene_up,rand_dz_gene_down,cmap_exp_signature)
    }
  }else{
    N_PERMUTATIONS <- n_perm #default 100000
    random_sig_ids <- sample(1:ncol(lincs_signatures),N_PERMUTATIONS,replace=T)
    random_cmap_scores <- NULL
    for (exp_id in random_sig_ids){
      #print(count)
      cmap_exp_signature <- data.frame(gene.list,  rank(-1 * lincs_signatures[, as.character(exp_id)], ties.method="random"))    
      colnames(cmap_exp_signature) <- c("ids","rank")
      
      random_input_signature_genes <- sample(gene.list, (nrow(dz_genes_up)+nrow(dz_genes_down)))
      rand_dz_gene_up <- data.frame(GeneID=random_input_signature_genes[1:nrow(dz_genes_up)])
      rand_dz_gene_down <- data.frame(GeneID=random_input_signature_genes[(nrow(dz_genes_up)+1):length(random_input_signature_genes)])
      random_cmap_scores <- c(random_cmap_scores, cmap_score_new(rand_dz_gene_up,rand_dz_gene_down,cmap_exp_signature))
    }
  }
  return(random_cmap_scores)
  }




drug_enrichment <- function(sRGES,target_type){
  require(GSVA)
  load(paste0(pipelineDataFolder,"cmpd_sets_", target_type, ".RData"))
  cmpdSets = cmpd_sets$cmpd.sets
  names(cmpdSets) = cmpd_sets$cmpd.set.names
  
  drug_pred = sRGES
  #create a random gene set
  random_times = 1000
  rgess = matrix(NA, nrow = nrow(drug_pred), ncol = random_times)
  for (i in 1:random_times){
    rgess[, i] = sample(drug_pred$sRGES, nrow(rgess))
  }
  rownames(rgess) = drug_pred$pert_iname
  rgess = cbind(rgess, drug_pred$sRGES)
  
  gsea_results = gsva(rgess, cmpdSets, method = "ssgsea",  parallel.sz=8)
  gsea_summary = data.frame(score = gsea_results[,101])
  
  #p_value: test bootstrapped random scores that are less than the actual score
  gsea_summary$p = apply(gsea_results, 1, function(x){
    sum(x[1:random_times] < x[random_times+1])/random_times
  })
  
  gsea_summary$target = row.names(gsea_summary)
  gsea_summary$padj = p.adjust(gsea_summary$p)
  gsea_summary = gsea_summary[,c('target','score','p','padj')]
  gsea_summary = gsea_summary[order(gsea_summary$p), ]
  return(gsea_summary)
}


visualizeLincsHits = function(rges,dz_sigUsed,drugs=''){
  require(pheatmap)
  require(dplyr)
  require("gplots")
  require("ggplot2")
  require("RColorBrewer")
  dz_sig = dz_sigUsed %>% select(gene,log2FoldChange)
  rges$RGES = rges$cmap_score
  rges_subset = rges[rges$pert_iname %in% drugs,]
  drug_cmap_score = aggregate(RGES ~ pert_iname, rges_subset, median)
  drug_instances_median  = merge(rges_subset, drug_cmap_score, by = c("pert_iname"))
  drug_instances_median$diff = abs(drug_instances_median$RGES.x - drug_instances_median$RGES.y)   
  #cmap_score.y is the median
  drug_instances_min_diff = aggregate(diff ~ pert_iname, drug_instances_median, min)
  drug_instances_select = merge(drug_instances_median, drug_instances_min_diff, by=c("pert_iname", "diff"))
  drug_instances_select = drug_instances_select[!duplicated(drug_instances_select$pert_iname), ]
  sig_id_selects = drug_instances_select$id
  load(paste0(pipelineDataFolder,"LINCS_RGES/lincs_signatures_cmpd_landmark_symbol.RData"))
  drug_dz_signature = merge(dz_sig, data.frame(gene = rownames(lincs_signatures), 
                                               lincs_signatures[, as.character(sig_id_selects)]),  by="gene", suffixes='')
  gene_ids = drug_dz_signature$gene
  drug_dz_signature_rank = drug_dz_signature[,-1]
  for (i in 1:ncol(drug_dz_signature_rank)){
    drug_dz_signature_rank[,i] = rank(-1 * drug_dz_signature_rank[,i] ) #highly expressed genes ranked on the top
  }
  gene_ids_rank <- gene_ids[order(drug_dz_signature_rank[,1])]
  drug_dz_signature_rank <- drug_dz_signature_rank[order(drug_dz_signature_rank[,1]),] #order by disease expression
  
  col_sorted = sort(cor(drug_dz_signature_rank, method="spearman")["log2FoldChange",-1])    
  drug_dz_signature_rank = drug_dz_signature_rank[,c("log2FoldChange", names(col_sorted))]
  
  drug_names = sapply(2:ncol(drug_dz_signature_rank), function(id){
    rges$pert_iname[paste("X",rges$id, sep="") == names(drug_dz_signature_rank)[id]][1]
  })
  dz = 'case'
  pdf(paste(outputFolder, "/lincs_reverse_expression.pdf", sep=""))
  #colPal <- bluered(100)
  colPal = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))
  par(mar=c(13, 6, 2, 0.5))
  axiscolor = sapply(c(dz, as.character(drug_names)), function(name){
    if (name == dz){
      "black"
    }else if (name %in% ""){
      "black"
    }else{
      "black"
    }
  })
  image(t(drug_dz_signature_rank), col=colPal,   axes=F, srt=45)
  axis(1,  at= seq(0,1,length.out=ncol( drug_dz_signature_rank ) ), labels= FALSE)
  text(x = seq(0,1,length.out=ncol( drug_dz_signature_rank ) ), c(-0.05),
       labels = c( dz,as.character(drug_names)),col=axiscolor, srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=0.6)
  dev.off()
}

