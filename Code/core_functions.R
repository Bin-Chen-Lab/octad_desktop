# Combined DE_core_functions.R & drugs_core_functions.R
# 02/19/19
# from drugs_core_functions.R

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
  load(paste0(dataFolder,'lincs_signatures_cmpd_landmark_symbol.RData'))
  gene.list <- toupper(rownames(lincs_signatures))
  dz_sigUsed <- dz_signature %>% filter(toupper(Symbol) %in% gene.list)
  dz_genes_up <- dz_sigUsed %>% filter(log2FoldChange>0) %>% arrange(desc(log2FoldChange)) %>% head(maxGenes)
  dz_genes_down <- dz_sigUsed %>% filter(log2FoldChange<0) %>% arrange(log2FoldChange) %>% head(maxGenes)
  write.csv(rbind(dz_genes_up,dz_genes_down), paste0(outputFolder,'dz_sig_used.csv'))
  lincs_sig_info <- read.csv(paste0(dataFolder,"lincs_sig_info.csv"),
                             stringsAsFactors = F)
  lincs_sig_info <- lincs_sig_info %>% filter(id %in% colnames(lincs_signatures))
  if(choose_fda == T){
    fda_drugs = read.csv(paste0(dataFolder,"repurposing_drugs_20170327.csv"),
                         comment.char = "!", header=T, sep="\t")
    lincs_sig_info <- lincs_sig_info %>% filter(id %in% colnames(lincs_signatures) & 
                                                  tolower(pert_iname) %in% tolower(fda_drugs$pert_iname))
  }else{
    lincs_sig_info <- lincs_sig_info[!duplicated(lincs_sig_info$id),]
  }
  sig.ids <- lincs_sig_info$id
  if(parallel == T){
    require(doParallel)
    registerDoParallel(cores=2)
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
  dz_cmap_scores <- data.frame(id = sig.ids, cmap_score = dz_cmap_scores)
  dz_cmap_scores$id = as.character(dz_cmap_scores$id)
  results <- left_join(dz_cmap_scores, lincs_sig_info, by = "id")
  
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
    lincs_cell_line_weight <- read.csv(paste0(dataFolder, "/lincs_cell_lines_cor.csv"))
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
  
  load(paste0(dataFolder,'lincs_signatures_cmpd_landmark_symbol.RData'))
  gene.list <- toupper(rownames(lincs_signatures))
  dz_sigUsed <- dz_signature %>% filter(toupper(Symbol) %in% gene.list)
  write.csv(dz_sigUsed, paste0(outputFolder,'dz_sig_used.csv'))
  dz_genes_up <- dz_sigUsed %>% filter(log2FoldChange>0) %>% arrange(desc(log2FoldChange))
  dz_genes_down <- dz_sigUsed %>% filter(log2FoldChange<0) %>% arrange(log2FoldChange)
  
  lincs_sig_info <- read.csv(paste0(dataFolder,"lincs_sig_info.csv"),
                             stringsAsFactors = F)
  if(choose_fda == T){
    fda_drugs = read.csv(paste0(dataFolder,"repurposing_drugs_20170327.csv"),
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
  
  enrichFolder.n <- paste0(enrichFolder,target_type,'/')
  if (!dir.exists(enrichFolder.n)) {
    dir.create(enrichFolder.n)
  }
  
  #load random scores
  load(paste0(dataFolder,"cmpd_sets_", target_type, ".RData"))
  cmpdSets = cmpd_sets$cmpd.sets
  names(cmpdSets) = cmpd_sets$cmpd.set.names
  
  drug_pred = sRGES
  
  rgess = matrix(-1*drug_pred$sRGES, ncol = 1)
  rownames(rgess) = drug_pred$pert_iname
  gsea_results = gsva(rgess, cmpdSets, method = "ssgsea",  parallel.sz=8)
  
  gsea_results = cbind(random_gsea_score[[target_type]], gsea_results)
  
  gsea_p = apply(gsea_results, 1, function(x){
    sum(x[1:ncol(random_gsea_score[[target_type]])] > x[ncol(random_gsea_score[[target_type]])+1])/ncol(random_gsea_score[[target_type]])
  })
  
  gsea_p = data.frame(target = names(gsea_p), p = gsea_p, padj = p.adjust(gsea_p))
  gsea_p = gsea_p[order(gsea_p$padj), ]
  # return(gsea_p)
  write.csv(gsea_p, paste0(enrichFolder.n, "/enriched_", target_type, ".csv"))
  top.out.num = nrow(gsea_p[which(gsea_p$padj<=0.05),])
  if (top.out.num == 0) {
    top.out.num = 1
  }
  if(top.out.num > 50){
    top.out.num <- 50
  }
  for (i in 1:top.out.num) {
    top_target = as.character(gsea_p$target[i])
    sRGES$rank = rank(sRGES$sRGES)
    target_drugs_score = sRGES$rank[sRGES$pert_iname %in% cmpdSets[[top_target]]]
    if (length(target_drugs_score) < 3) {
      next
    }
    pdf(paste0(enrichFolder.n, "/top_enriched_", top_target, "_", target_type, ".pdf"))
    barcodeplot(sRGES$sRGES, target_drugs_score, main = top_target, xlab = "sRGES")
    dev.off()
  }
  if (target_type == "ChemCluster"){
    clusternames <- as.character((gsea_p[which(gsea_p$padj<=0.05),])$target)
    topclusterlist <- cmpdSets[clusternames]
    cat(sapply(topclusterlist, toString), file = paste0(enrichFolder.n,"misc.csv"), sep="\n")
    clusterdf <- read.csv2(paste0(enrichFolder.n,"misc.csv"), header=FALSE)
    clusterdf$cluster <- clusternames
    clusterdf$pval <- (gsea_p[which(gsea_p$padj<=0.05),])$padj
    colnames(clusterdf)[1] <- "drugs.in.cluster"
    write.csv(clusterdf,file=paste0(enrichFolder.n,'drugstructureclusters.csv'),row.names = F)
  }
}



# drug_enrichment <- function(sRGES,target_type){
#   require(GSVA)
#   load(paste0(dataFolder,"cmpd_sets_", target_type, ".RData"))
#   cmpdSets = cmpd_sets$cmpd.sets
#   names(cmpdSets) = cmpd_sets$cmpd.set.names
#   
#   drug_pred = sRGES
#   
#   #FDA approved using repurposing_drugs_20170327.txt
#   
#   
#   
#   #create a random gene set
#   random_times = 1000
#   rgess = matrix(NA, nrow = nrow(drug_pred), ncol = random_times)
#   for (i in 1:random_times){
#     rgess[, i] = sample(drug_pred$sRGES, nrow(rgess))
#   }
#   rownames(rgess) = drug_pred$pert_iname
#   rgess = cbind(rgess, drug_pred$sRGES)
#   
#   gsea_results = gsva(rgess, cmpdSets, method = "ssgsea",  parallel.sz=8)
#   
#   gsea_p = apply(gsea_results, 1, function(x){
#     sum(x[1:random_times] > x[random_times+1])/random_times
#   })
#   
#   gsea_p = data.frame(target = names(gsea_p), p = gsea_p, padj = p.adjust(gsea_p))
#   gsea_p = gsea_p[order(gsea_p$p), ]
#   return(gsea_p)
# }

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
  load(paste0(dataFolder,"lincs_signatures_cmpd_landmark_symbol.RData"))
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

pick.out.cell.line <- function(expr.of.samples,expr.of.cell.lines,marker.gene){
  marker.gene           <- intersect(rownames(expr.of.samples),(marker.gene))  
  marker.gene           <- intersect(rownames(expr.of.cell.lines),(marker.gene)) 
  correlation.matrix    <- cor(expr.of.samples[marker.gene,],expr.of.cell.lines[marker.gene,],method='spearman')
  cell.line.median.cor  <- apply(correlation.matrix,2,median) %>% sort(decreasing = TRUE)
  best.cell.line        <- names(cell.line.median.cor)[1]
  p.value.vec           <- foreach(cell.line= setdiff(names(cell.line.median.cor),best.cell.line),.combine='c') %do% {
    v                     <- correlation.matrix[,cell.line]
    p.value               <- wilcox.test(correlation.matrix[,best.cell.line],v,alternative = 'greater',paired = TRUE)$p.value
  }
  names(p.value.vec) <- setdiff(names(cell.line.median.cor),best.cell.line)
  fdr.vec            <- p.adjust(p.value.vec,method='fdr')
  list(cell.line.median.cor=cell.line.median.cor,best.cell.line=best.cell.line,compare.fdr.vec=fdr.vec,correlation.matrix = correlation.matrix )
}


# from DE_core_functions.R

####Functions####
#computeRefTissue : returns the reference tissue given 'normal' and 'case' ids and 'expression' dataframe

computeRefTissue <- function(case_id = '', normal_id = '',
                             expSet=NULL,n_varGenes = 5000,
                             method='varGenes', #random 
                             #site_selection = 'top', 
                             #top site or any cor cutoff>quantile,
                             #all will select samples from 90th percentile
                             control_size = 50,
                             cor_cutoff='0%', #greater or equal than the cutoff 
                             output=T){
  #### required input ####
  # expSet: DF matrix with colnames consisting both the case and normal ids given
  # case_id: case to select correlated tissues
  # normal_id: all the normal tissues you want to check against the case
  # outputFolder: this needs to be set to output the intermediate files for some visualization
  #### options #### 
  #control_size : if you're selecting the max number of tissues it gives back
  #method: varGenes: will give the top varying ones, or random: random ones 
  #cor_cutoff : will give back tissue based on top_nth quantile of the correlation
  #must be a string in percents of 5s e.g. '5%' '10%' '25%' etc...
  #### Returns ####
  # GTEXid : don't worry about the name this returns the highest varying normals
  #### Outputs ####
  # /case_normal_corMatrix.csv : correlation matrix btw normal id and case id
  # /case_normal_median_cor.csv : median correlation btw normal id and case id
  
  
  require(dplyr)
  if(method == 'random'){
    GTEXid <- sample(normal_id,size = control_size)
    GTEXid
  }else if(method == 'varGenes'){
    expSet_normal <- expSet[,normal_id]
    expSet_case <- expSet[,case_id]
    #varGenes look at the top varying genes (IQR) within normal tissue expression and varies them to the case tissues
    iqr_gene <-apply(expSet_normal, 1, IQR) #get the IQR per gene
    varying_genes <-order(iqr_gene, decreasing=T)[1:min(n_varGenes,length(iqr_gene))]
    
    #get the correlation matrix for each normal id and each case id
    normal_dz_cor <-cor(expSet_normal[varying_genes, ], expSet_case[varying_genes, ], method = "spearman")
    normal_dz_cor_each <-apply(normal_dz_cor, 1, median) #getting the median correlation btw each normal tissue to the case overall
    normal_dz_cor_eachDF = data.frame(cor=sort(normal_dz_cor_each, decreasing=TRUE)) %>% 
      mutate(sample.id = row.names(.)) %>% select(sample.id,cor)
    cutoff = quantile(normal_dz_cor_eachDF$cor,probs=seq(0,1,0.05),na.rm=T)[cor_cutoff]
    GTEXid <- (normal_dz_cor_eachDF %>% 
                 arrange(desc(cor)) %>% 
                 filter(cor>=cutoff))$sample.id 
    GTEXid <- GTEXid[1:min(control_size,length(GTEXid))]
    
    if(output==T){
      tryCatch(write.csv(normal_dz_cor,file = paste0(outputFolder,'case_normal_corMatrix.csv')),
               error = function(c) "failed to write case normal cor matrix csv. Try checking if your outputFolder string is correct or exists")
      
      tryCatch(write.csv(normal_dz_cor_eachDF,row.names = F, paste0(outputFolder, "/case_normal_median_cor.csv")),
               error = function(c) "failed to write case normal median correlation csv. Try checking if your outputFolder string is correct or exists")
      
    }
    GTEXid
  }
}


computeCellLine <- function(case_id = '',
                            expSet = NULL,
                            LINCS_overlaps = T){
  require(plyr)
  require(dplyr)
  require(foreach)
  require(data.table)
  #helper function by Ke
  pick.out.cell.line <- function(expr.of.samples,expr.of.cell.lines,marker.gene){
    marker.gene           <- intersect(rownames(expr.of.samples),(marker.gene))  
    marker.gene           <- intersect(rownames(expr.of.cell.lines),(marker.gene)) 
    correlation.matrix    <- cor(expr.of.samples[marker.gene,],expr.of.cell.lines[marker.gene,],method='spearman')
    cell.line.median.cor  <- apply(correlation.matrix,2,median) %>% sort(decreasing = TRUE)
    best.cell.line        <- names(cell.line.median.cor)[1]
    p.value.vec           <- foreach(cell.line= setdiff(names(cell.line.median.cor),best.cell.line),.combine='c') %do% {
      v                     <- correlation.matrix[,cell.line]
      p.value               <- wilcox.test(correlation.matrix[,best.cell.line],v,alternative = 'greater',paired = TRUE)$p.value
    }
    names(p.value.vec) <- setdiff(names(cell.line.median.cor),best.cell.line)
    fdr.vec            <- p.adjust(p.value.vec,method='fdr')
    list(cell.line.median.cor=cell.line.median.cor,best.cell.line=best.cell.line,compare.fdr.vec=fdr.vec,correlation.matrix = correlation.matrix ) 
  }
  load(paste0(dataFolder,'CCLE_OCTAD.RData'))
  case_counts <- expSet[,case_id]
  if(LINCS_overlaps == T){
    
    CCLE.median                 <- apply(CCLE.overlaps,1,median)
  }else{
    CCLE.median = apply(CCLE.log2.read.count.matrix,1,median)
  }
  CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
  tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
  tmp.rank                    <- apply(tmp,2,rank)
  rank.mean                   <- apply(tmp.rank,1,mean)
  rank.sd                     <- apply(tmp.rank,1,sd)
  CCLE.rna.seq.marker.gene.1000                 <- names(sort(rank.sd,decreasing =TRUE))[1:1000]
  TCGA.vs.CCLE.polyA.expression.correlation.result  <- try(pick.out.cell.line(case_counts, CCLE.overlaps,CCLE.rna.seq.marker.gene.1000), 
                                                           outFile = paste0(outputFolder,"cellLineRankingError.txt"))
  topline <- data.frame(medcor = TCGA.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor) # could also do first 
  # 3 of TCGA.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor
  #topline <- rownames(topline)[1]
  # topline$cellLine = row.names(topline)
  # topline = topline %>% select(cellLine,medcor)
  # load(paste0(dataFolder,'CCLE_OCTAD.RData'))
  # topline = left_join(topline,CCLE_demoDat %>% select(CCLE.shortname,Age,Gender,Race,Site.Primary,Histology,Hist.Subtype1),
  #  by=c('cellLine'='CCLE.shortname'))
  return(topline)
}

# topLineEval <- function(topline = '',mysRGES){
#   require(plotly)
#   require(ggplot2)
#   require(dplyr)
#   load(paste0(dataFolder,'CCLE_OCTAD.RData'))
#   mysRGES$pert_iname <- tolower(mysRGES$pert_iname)
#   CTRPv2.sensprof$ic50_recomputed[!is.finite(CTRPv2.sensprof$ic50_recomputed)] <- 0
#   CTRPv2.ic50   <- dcast.data.table(CTRPv2.sensprof, drugid ~ cellid, value.var = "ic50_recomputed", fun.aggregate = mean)
#   CTRPv2.ic50 <- replace(CTRPv2.ic50, is.na(CTRPv2.ic50), 0)
#   colnames(CTRPv2.ic50) <- gsub("[^0-9A-Za-z///' ]","",colnames(CTRPv2.ic50))
#   colnames(CTRPv2.ic50) <- toupper(colnames(CTRPv2.ic50))
#   colnames(CTRPv2.ic50) <- gsub(" ","",colnames(CTRPv2.ic50))
#   CTRP.IC50   <- CTRPv2.ic50[,c('DRUGID',topline)]
#   
#   CTRPv2.sensprof$auc_recomputed[!is.finite(CTRPv2.sensprof$auc_recomputed)] <- 0
#   CTRPv2.auc   <- dcast.data.table(CTRPv2.sensprof, drugid ~ cellid, value.var = "auc_recomputed", fun.aggregate = mean)
#   CTRPv2.auc <- replace(CTRPv2.auc, is.na(CTRPv2.auc), 0)
#   colnames(CTRPv2.auc) <- gsub("[^0-9A-Za-z///' ]","",colnames(CTRPv2.auc))
#   colnames(CTRPv2.auc) <- toupper(colnames(CTRPv2.auc))
#   colnames(CTRPv2.auc) <- gsub(" ","",colnames(CTRPv2.auc))
#   CTRP.auc <- CTRPv2.auc[,c('DRUGID',topline)]
#   CTRP.auc <- as.data.frame(CTRP.auc)
#   CTRP.auc.n0 <- CTRP.auc[CTRP.auc[,2] != 0,]
#   
#   
#   deff <- CTRP.IC50
#   deffnd <- deff[!duplicated(deff$DRUGID),]
#   deffnd$DRUGID <- tolower(deffnd$DRUGID)
#   testdf <- merge(mysRGES, deffnd, by.x = "pert_iname", by.y = "DRUGID")
#   IC50.cortest <- cor.test(testdf$sRGES,log2(testdf[,8] +1))
#   ic50pval <- IC50.cortest$p.value
#   ic50rho  <- IC50.cortest$estimate
#   mylabel = c("p-value" = ic50pval,'Rho' = ic50rho)
#   # plot(testdf$sRGES, log2(testdf[,8] +1),xlab = "sRGES",ylab="log2(ic50+1)", main = "Top Line recomputed log2(ic50 + 1) vs sRGES")
#   
#   testdf$StronglyPredicted <- NA
#   testdf$StronglyPredicted <- ifelse(testdf$sRGES < -0.2,"Yes","No")
#   
#   StronglyPredicted <- testdf$StronglyPredicted
#   
#   Legend.title <- "Strongly <br>Predicted"
#   Legend.label1<- "No" 
#   Legend.label2<- "Yes"
#   Title <- "Top Line recomputed log2(ic50 + 1) vs sRGES"
#   xaxis <- "sRGES"
#   yaxis <- "log2(ic50+1)"
#   
#   p <- ggplot() +
#     geom_point(data = testdf, 
#                aes(x = testdf$sRGES, y = log2(testdf[,8] +1), 
#                    color = StronglyPredicted, 
#                    text = 
#                      paste('Drug: ', pert_iname, # These are the hover labels generated by p1
#                            '<br>sRGES: ', sRGES))) +
#     geom_smooth(data = testdf,          
#                 aes(x = testdf$sRGES, y = log2(testdf[,8] +1),
#                     label = ic50rho),
#                 method = 'lm',se = F, color = "black",size = 0.5) +
#     scale_color_discrete(
#       name = Legend.title) +
#     labs(x = xaxis, y = yaxis, title = Title) +
#     theme(legend.position = "right", legend.background = element_rect(fill="#F5F5F5"), legend.title = element_blank())
#   
#   p1 <- ggplotly(p, tooltip = c("text", "label")) %>%  layout(margin=list(l=15)) # change to white background
#   
#   ic50graph <- p1 %>% 
#     add_annotations( text=Legend.title,xref="paper",yref="paper",x=1.02,xanchor="left",y=0.8,yanchor="bottom",legendtitle=TRUE,showarrow=FALSE ) %>%
#     layout( legend=list(y=0.8, yanchor="top" ) )
#   
#   htmlwidgets::saveWidget(as_widget(ic50graph), paste0(outputFolder,topline,"ic50_insilico_validation.html"),selfcontained = F) 
#   # note: files are simply too large to set selfcontained = T. This just causes issues on linux machines.
#   
#   
#   # AUC
#   deff <- CTRP.auc.n0
#   deffnd <- deff[!duplicated(deff$DRUGID),]
#   deffnd$DRUGID <- tolower(deffnd$DRUGID)
#   testdf2 <- merge(mysRGES, deffnd, by.x = "pert_iname", by.y = "DRUGID")
#   AUC.cortest <- cor.test(testdf2$sRGES,testdf2[,8])
#   # plot(testdf2$sRGES, testdf2[,8],xlab = "sRGES",ylab=paste(topline,"AUC"), main = "Top Line recomputed AUC vs sRGES")
#   
#   testdf2$StronglyPredicted <- NA
#   testdf2$StronglyPredicted <- ifelse(testdf2$sRGES < -0.2,"Yes","No")
#   
#   StronglyPredicted <- testdf2$StronglyPredicted
#   
#   aucpval <- AUC.cortest$p.value
#   aucrho  <- AUC.cortest$estimate
#   
#   Legend.title <- "Strongly <br>Predicted"
#   Legend.label1<- "No"
#   Legend.label2<- "Yes"
#   Title <- "Top Line recomputed AUC vs sRGES"
#   xaxis <- "sRGES"
#   yaxis <- "AUC"
#   
#   
#   p <- ggplot() +
#     geom_point(data = testdf2, 
#                aes(x = testdf2$sRGES, y = testdf2[,8], 
#                    color = StronglyPredicted, 
#                    text = 
#                      paste('Drug: ', pert_iname, # These are the hover labels generated by p1
#                            '<br>sRGES: ', sRGES))) +
#     geom_smooth(data = testdf2,          
#                 aes(x = testdf2$sRGES, y = testdf2[,8],
#                     label = aucrho),
#                 method = 'lm',se = F, color = "black",size = 0.5) +
#     scale_color_discrete(
#       name = Legend.title) +
#     labs(x = xaxis, y = yaxis, title = Title) +
#     theme(legend.position = "right", legend.background = element_rect(fill="#F5F5F5"), legend.title = element_blank())
#   
#   p1 <- ggplotly(p, tooltip = c("text", "label")) %>%  layout(margin=list(l=15))
#   
#   aucgraph <- p1 %>% 
#     add_annotations( text=Legend.title,xref="paper",yref="paper",x=1.02,xanchor="left",y=0.8,yanchor="bottom",legendtitle=TRUE,showarrow=FALSE ) %>%
#     layout( legend=list(y=0.8, yanchor="top" ) )
#   
#   htmlwidgets::saveWidget(as_widget(aucgraph), paste0(outputFolder,topline,"auc_insilico_validation.html"),selfcontained = F) # test
#   
#   
#   # logging cortests
#   con <- file(paste0(outputFolder,"drug_sensitivity_insilico_results.txt"))
#   sink(con, append=TRUE)
#   sink(con, append=TRUE, type="message")
#   
#   print("AUC cortest")
#   print(AUC.cortest)
#   print("IC50 cortest")
#   print(IC50.cortest)
# 
#   sink() 
#   sink(type="message")
#   }


remLowExpr <- function(counts,counts_phenotype){
  x <-DGEList(counts = round(counts), group = counts_phenotype$sample_type )
  cpm_x <- cpm(x)
  #needs to be at least larger the than the size of the smallest set
  keep.exprs <- rowSums(cpm_x>1) >= min(table(counts_phenotype$sample_type)) 
  keep.exprs
}

#compute empirical control genes for RUVg
#consider deprecating since this is just first pass DE
compEmpContGenes <- function(counts, counts_phenotype, n_topGenes = 5000){
  set <- newSeqExpressionSet(round(counts),
                             phenoData = data.frame(counts_phenotype,row.names=counts_phenotype$sample))
  design <- model.matrix(~ sample_type, data = pData(set))
  y <- DGEList(counts=counts(set), group =  counts_phenotype$sample)
  y <- calcNormFactors(y, method="TMM") #upperquartile generate Inf in the LGG case
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit,2) #defaults to compare tumor to normal or tumor mutant to normal
  write.csv(lrt$table,paste0(outputFolder,'unnormalizedGenes.csv'))
  top <- topTags(lrt, n=nrow(set))$table
  #n_topGenes <- n_topGenes #5000: assume there are 5000 signficant genes
  
  #based on n_topGenes computing genes with low DE
  #the genes not computed significant is in the empirical set
  i = which(!(rownames(set) %in% rownames(top)[1:min(n_topGenes,dim(top)[1])]))
  empirical <- rownames(set)[i]
  empirical
}




diffExp <- function(case_id='',control_id='',expSet=dz_expr,
                    normalize_samples=T,
                    k=1,
                    n_topGenes=10000,
                    DE_method='edgeR'){
  require(dplyr)
  require(RUVSeq)
  require(edgeR)
  counts_phenotype <- rbind(data.frame(sample = case_id,sample_type = 'case'),
                            data.frame(sample = control_id, sample_type = 'control'))
  counts = expSet[,as.character(counts_phenotype$sample)]
  counts = 2^counts - 1 #unlog the counts it was log(2x + 1) in dz.expr.log2.readCounts
  counts_phenotype$sample = as.character(counts_phenotype$sample)
  counts_phenotype$sample_type = factor(counts_phenotype$sample_type, levels = c("control", "case"))
  highExpGenes <- remLowExpr(counts,counts_phenotype)
  counts = counts[highExpGenes,]
  set <- newSeqExpressionSet(round(counts),
                             phenoData = data.frame(counts_phenotype,row.names=counts_phenotype$sample))
  design <- model.matrix(~ sample_type, data = pData(set))
  y <- DGEList(counts=counts(set), group =  counts_phenotype$sample)
  y <- calcNormFactors(y, method="TMM") #upperquartile generate Inf in the LGG case
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit,2) #defaults to compare case control
  
  #counts dataframe === remove low counts ===> set === normalized ===> set1
  #if no empirical genes found it will just create a matrix without running RUVg
  if(normalize_samples == T){
    top <- topTags(lrt, n=nrow(set))$table
    i = which(!(rownames(set) %in% rownames(top)[1:min(n_topGenes,dim(top)[1])]))
    empirical <- rownames(set)[i]
    stopifnot(length(empirical)>0)
    write.csv(data.frame(empirical),file=paste0(outputFolder,"computedEmpGenes.csv"))
    set1 <- RUVg(set,empirical,k=k)
    
    print('computing DE via edgeR')
    
    #construct model matrix based on whether there was normalization ran
    if(normalize_samples == T){
      if(k==1){
        design <- model.matrix(~sample_type + W_1, data=pData(set1))
      }else if(k == 2){
        design <- model.matrix(~sample_type + W_1 + W_2, data = pData(set1))
      }else if (k == 3){
        design <- model.matrix(~sample_type + W_1 + W_2 + W_3, data = pData(set1))
      }
    }else{design <- model.matrix(~sample_type,data=pData(set1))}
    
    dgList <- DGEList(counts=counts(set1),group=set1$sample_type)
    dgList <- calcNormFactors(dgList, method="TMM") #using upperquartile seems to give issue for LGG
    dgList <- estimateGLMCommonDisp(dgList, design)
    dgList <- estimateGLMTagwiseDisp(dgList, design)
    fit <- glmFit(dgList, design)
    
    #see edgeRUsersGuide section on testing for DE genes for contrast
    lrt <- glmLRT(fit,2) 
    #second coefficient otherwise it'll default the W_1 term when normalize is on
  }
  res <- lrt$table
  colnames(res) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
  res$padj <- p.adjust(res$pvalue)
  res$identifier <- row.names(res)
  res = res %>% select(identifier,everything())
  return(res)
}

geneEnrich <- function(dz_signature, 
                       db.list=c( "KEGG_2016",  "GO_Biological_Process_2017", "GO_Cellular_Component_2017"),
                       suffix=''){
  #adapted from https://github.com/compbiomed/enrichR
  #edited for use by Billy Zeng, UCSF, Anita Wen UCD
  #' Perform functional enrichment on a set of genes.
  #' 
  #' This function interacts with Enrichr's REST API in order to perform functional enrichment of a single
  #' set of genes, for a set of specified databases which are already fronted by Enrichr.
  #' Databases are specified as seen in the web interface, with underscores for spaces
  #' (e.g. "WikiPathways_2016", "KEGG_2016", "GO_Biological_Process"). There's no way to query Enrichr
  #' to get these database names, so they can't be provided as options. You'll just have to guess. Sorry :/
  #' 
  #' @param gene.list a list of gene symbols
  #' @param databases a list of Enrichr-fronted databases, as mentioned above
  #' @param fdr.cutoff An FDR (adjusted p-value) threshold by which to limit the list of enriched pathways
  #' @keywords functional enrichment Enrichr
  #' @export
  
  enrichGeneList <- function (gene.list, databases=db.list, fdr.cutoff=NULL) {
    ######Step 1: Post gene list to EnrichR
    req.body <- list(list=paste(gene.list, collapse="\n"))
    post.req <- httr::POST("http://amp.pharm.mssm.edu/Enrichr/enrich", encode="multipart", body=I(req.body))
    
    #TODO: Real error handling..
    if (!grepl("success", httr::http_status(post.req)$category, ignore.case=T)) stop("Posting gene list to EnrichR failed")
    
    ######Step 2: Get results from posted gene list
    database.enrichments <- list()
    for (idx in 1:length(databases)) { 
      database <- databases[idx]
      get.req <- httr::GET(paste("http://amp.pharm.mssm.edu/Enrichr/enrich?backgroundType=", database, sep=""))
      if (!grepl("success", httr::http_status(get.req)$category, ignore.case=T)) stop("Retrieving results from EnrichR failed")
      
      response.content <- mungeResponseContent(httr::content(get.req)[[database]])
      
      if (length(response.content) > 1) {
        database.res <- data.table::rbindlist(response.content)
        database.res[, 1] <- rep(database, nrow(database.res))
        database.enrichments[[idx]] <- database.res[, paste("V", c(1, 2, 3, 7, 5,6), sep=""), with=F]
      }
    }
    
    query.results <- as.data.frame(data.table::rbindlist(database.enrichments))
    colnames(query.results) <- c("database", "category", "pval", "qval","Combined Score" ,"genes")
    
    if (!is.null(fdr.cutoff)) {
      query.results <- query.results[query.results$qval < fdr.cutoff, ]
    }
    
    query.results
  }
  
  
  
  #' Munge the Enrichr API response so it'll squeeze neatly (if untidily) into a dataframe.
  #' 
  #' The response from the Enrichr API is a list of lists, where each nested list item represents an enriched
  #' category. The 6th item of each category (i.e. response.content[[category.idx]][[6]]) corresponds to the 
  #' genes that overlapped with the gene set behind that category. This function bascically collapses that list of 
  #' genes into a single string. 
  #' 
  #' I'm sorry you ever had to look at this.
  #' 
  #' @param response.content result of calling httr::content on the GET request to the Enrichr API, after submitting a list for enrichment.
  #' 
  mungeResponseContent <- function (response.content) {
    munged.content <- response.content
    if (length(response.content) == 0) return(NA)
    
    for (idx in 1:length(response.content)) {
      munged.content[[idx]][[6]] <- paste(munged.content[[idx]][[6]], collapse=",")
    }
    
    munged.content
  }
  
  
  #munged.content not found
  
  
  databases = db.list
  fdr.cutoff = 0.25
  
  up.genes = dz_signature$Symbol[dz_signature$log2FoldChange > 0]
  if(length(up.genes)>0){
    (up.gene.res = enrichGeneList(up.genes[1:min(300, length(up.genes))], databases, fdr.cutoff))
    (up.gene.res = up.gene.res[order(up.gene.res$database, up.gene.res$p), ])
    write.csv(up.gene.res, paste0(outputFolder, "/dz_up_sig_genes_enriched", suffix,".csv"))  
  }
  dn.genes = dz_signature$Symbol[dz_signature$log2FoldChange < 0]
  if(length(dn.genes)>0){
    dn.gene.res = enrichGeneList(dn.genes[1:min(300, length(dn.genes))], databases, fdr.cutoff)
    dn.gene.res = dn.gene.res[order(dn.gene.res$database, dn.gene.res$p), ]
    write.csv(dn.gene.res, paste0(outputFolder, "/dz_down_sig_genes_enriched", suffix,".csv"))  
  }
  up.gene.res$dir = 'up'
  dn.gene.res$dir = 'dn'
  res = rbind(up.gene.res,dn.gene.res)
}

getptidDF = function(tcga_sample.id){
  require(tidyr)
  ptid = t(data.frame(tcga_sample.id %>% strsplit('-',useBytes = T)))
  ptDF = data.frame(cbind(tcga_sample.id,ptid))
  ptDF$ptid = paste(ptDF[,2],ptDF[,3],ptDF[,4],sep = "-")
  ptDF = ptDF %>% select(1,ptid,5)
  colnames(ptDF) = c('sample.id','ptid','tcga.sample.type')
  ptDF.case_control = ptDF %>% 
    group_by(ptid) %>% spread(tcga.sample.type,sample.id) %>% ungroup()
  return(ptDF.case_control)
}


# testcorTissue = function(case_id='',
#                          normal_id='',
#                          expSet=NULL,n_varGenes=10000,
#                          outRows=1,output=T){
#   require(dplyr)
#   expSet_normal <- expSet[,normal_id]
#   expSet_case <- as.matrix(expSet[,case_id])
#   iqr_gene <-apply(expSet_normal, 1, IQR) #get the IQR per gene
#   varying_genes <-order(iqr_gene, decreasing=T)[1:min(n_varGenes,length(iqr_gene))]
#   
#   #get the correlation matrix for each normal id and each case id
#   normal_dz_cor <-cor(expSet_normal[varying_genes, ], expSet_case[varying_genes, ], method = "spearman")
#   normal_dz_cor_each <-apply(normal_dz_cor, 1, median) #getting the median correlation btw each normal tissue to the case overall
#   normal_dz_cor_eachDF = data.frame(cor=sort(normal_dz_cor_each, decreasing=TRUE)) %>% 
#     mutate(sample.id = row.names(.)) %>% select(sample.id,cor)
#   phenoDF = read.csv('~/octad_desktop/Data/phenoDF_withage.csv',stringsAsFactors = F)    
#   normal_dz_cor_eachDF = normal_dz_cor_eachDF %>% left_join(phenoDF,by='sample.id')
#   normal_dz_cor_biopsy = normal_dz_cor_eachDF %>% 
#     group_by(sample.type,biopsy.site,cancer,data.source) %>% 
#     summarise(minCor = min(cor),medCor = median(cor),maxCor = max(cor)) %>% ungroup()
#   out = normal_dz_cor_biopsy %>% arrange(desc(medCor))
#   out$case_id = case_id
#   
#   if(output==T){
#     tryCatch(write.csv(normal_dz_cor,file = paste0(outputFolder,'case_normal_corMatrix.csv')),
#              error = function(c) "failed to write case normal cor matrix csv. Try checking if your outputFolder string is correct or exists")
#     
#     tryCatch(write.csv(normal_dz_cor_eachDF,row.names = F, paste0(outputFolder, "/case_normal_median_cor.csv")),
#              error = function(c) "failed to write case normal median correlation csv. Try checking if your outputFolder string is correct or exists")
#     tryCatch(write.csv(normal_dz_cor_biopsy,row.names = F, paste0(outputFolder, "/case_normal_biopsy_median_cor.csv")),
#              error = function(c) "failed to write case normal median correlation csv. Try checking if your outputFolder string is correct or exists")
#     
#   }
#   return(out[1:outRows,])
# }
# 
# cmap_score_new <- function(sig_up, sig_down, drug_signature) {
#   #the old function does not support the input list with either all up genes or all down genes, this new function attempts to addess this.
#   #we also modify the original CMap approach: whenever the sign of ks_up/ks_down, we substract the two scores such that the final scores would not enrich at 0.
#   
#   num_genes <- nrow(drug_signature)
#   ks_up <- 0
#   ks_down <- 0
#   connectivity_score <- 0
#   
#   # I think we are re-ranking because the GeneID mapping changed the original rank range
#   drug_signature[,"rank"] <- rank(drug_signature[,"rank"])
#   
#   # Merge the drug signature with the disease signature by GeneID. This becomes the V(j) from the algorithm description
#   up_tags_rank <- merge(drug_signature, sig_up, by.x = "ids", by.y = 1)
#   down_tags_rank <- merge(drug_signature, sig_down, by.x = "ids", by.y = 1)
#   
#   up_tags_position <- sort(up_tags_rank$rank)
#   down_tags_position <- sort(down_tags_rank$rank)
#   
#   num_tags_up <- length(up_tags_position)
#   num_tags_down <- length(down_tags_position)
#   
#   # 
#   if(num_tags_up > 1) {
#     a_up <- 0
#     b_up <- 0
#     
#     a_up <- max(sapply(1:num_tags_up,function(j) {
#       j/num_tags_up - up_tags_position[j]/num_genes
#     }))
#     b_up <- max(sapply(1:num_tags_up,function(j) {
#       up_tags_position[j]/num_genes - (j-1)/num_tags_up
#     }))
#     
#     if(a_up > b_up) {
#       ks_up <- a_up
#     } else {
#       ks_up <- -b_up
#     }
#   }else{
#     ks_up <- 0
#   }
#   
#   if (num_tags_down > 1){
#     
#     a_down <- 0
#     b_down <- 0
#     
#     a_down <- max(sapply(1:num_tags_down,function(j) {
#       j/num_tags_down - down_tags_position[j]/num_genes
#     }))
#     b_down <- max(sapply(1:num_tags_down,function(j) {
#       down_tags_position[j]/num_genes - (j-1)/num_tags_down
#     }))
#     
#     if(a_down > b_down) {
#       ks_down <- a_down
#     } else {
#       ks_down <- -b_down
#     }
#   }else{
#     ks_down <- 0
#   }
#   
#   if (ks_up == 0 & ks_down != 0){ #only down gene inputed
#     connectivity_score <- -ks_down 
#   }else if (ks_up !=0 & ks_down == 0){ #only up gene inputed
#     connectivity_score <- ks_up
#   }else if (sum(sign(c(ks_down,ks_up))) == 0) {
#     connectivity_score <- ks_up - ks_down # different signs
#   }else{
#     connectivity_score <- ks_up - ks_down
#   }
#   
#   return(connectivity_score)
# }
# #summary RGES function
# getsRGES <- function(RGES, cor, pert_dose, pert_time, diff, max_cor){
#   
#   sRGES <- RGES
#   pert_time <- ifelse(pert_time < 24, "short", "long")
#   pert_dose <- ifelse(pert_dose < 10, "low", "high")
#   if (pert_time == "short" & pert_dose == "low"){
#     sRGES <- sRGES + diff[4]
#   }
#   if (pert_dose ==  "low" & pert_time == "long"){
#     sRGES <- sRGES + diff[2]
#   }
#   if (pert_dose ==  "high" & pert_time == "short"){
#     sRGES <- sRGES + diff[1]
#   }
#   return(sRGES * cor/max_cor) #
# }
# 
# 
# drug_enrichment <- function(sRGES,target_type){
#   require(GSVA)
#   load(paste0(dataFolder,"cmpd_sets_", target_type, ".RData"))
#   cmpdSets = cmpd_sets$cmpd.sets
#   names(cmpdSets) = cmpd_sets$cmpd.set.names
#   
#   drug_pred = sRGES
#   #drug_pred = read.csv(paste0("~/Documents/stanford/tumor_cell_line/RGES_manuscript/release/data/LIHC/lincs_cancer_sRGES.csv"))
#   
#   #FDA approved using repurposing_drugs_20170327.txt
#   
#   
#   
#   #create a random gene set
#   random_times = 1000
#   rgess = matrix(NA, nrow = nrow(drug_pred), ncol = random_times)
#   for (i in 1:random_times){
#     rgess[, i] = sample(drug_pred$sRGES, nrow(rgess))
#   }
#   rownames(rgess) = drug_pred$pert_iname
#   rgess = cbind(rgess, drug_pred$sRGES)
#   
#   gsea_results = gsva(rgess, cmpdSets, method = "ssgsea",  parallel.sz=8)
#   
#   gsea_p = apply(gsea_results, 1, function(x){
#     sum(x[1:random_times] > x[random_times+1])/random_times
#   })
#   
#   gsea_p = data.frame(target = names(gsea_p), p = gsea_p, padj = p.adjust(gsea_p))
#   gsea_p = gsea_p[order(gsea_p$p), ]
#   return(gsea_p)
# }