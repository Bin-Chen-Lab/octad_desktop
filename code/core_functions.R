# Combined DE_core_functions.R & drugs_core_functions.R
# 02/19/19
# from drugs_core_functions.R
# last edit 06/29/19 by BZ

#### Functions ####
## Desktop Functions
# cmap_score_new
# getsRGES
# runsRGES
# drug_enrichment
# pick.out.cell.line
# computeRefTissue
# computeCellLine
# topLineEval
# remLowExpr
# compEmpContGenes
# diffExp
# geneEnrich
# visualize_drug_hits
## octad.org Functions
# id_mapping
# id_mapping_gene_ensembl
# queryGDC
# estimatePurity
# detectOutlier
# compute_tissue_cell_cor
# compute_tissue_lincs_cell_cor
# visualize_top_ref_tissue
# RDFHS prep
## deprecated Functions


####### cmap_score_new #######
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

####### getsRGES #######
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
  return(sRGES * cor/max_cor) 
}

####### runsRGES #######
runsRGES <- function(dz_signature,choose_fda_drugs = F,parallel = F,max_gene_size=900, landmark = 1,
                     cells=''){
  require("dplyr")
  require("ggplot2")
  require(data.table)
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
  #dataFolder : root folder with lincs input data
  #dz_signature : disease signature genes dataframe
  #required columns Symbol : this must be an UPPERCASE gene symbol needed to join with RGES landmark genes
  #must change this if they require another identifier
  
  ####parameters####
  #store these parameters somewhere if you are doing multiple runs for comparison
  
  # default parameters
  # landmark = 1
  # choose_fda_drugs = F
  # max_gene_size = 100
  # weight_cell_line = F
  #load dz_signature
  #load output folder
  #load LINCS drug gene expression profiles
  
  lincs_sig_info <- fread(paste0(dataFolder,"lincs_sig_info.csv"), stringsAsFactors = TRUE)
  if(cells!=''){
    lincs_sig_info$cell_id = toupper(lincs_sig_info$cell_id)
    lincs_sig_info = lincs_sig_info %>% filter(cell_id %in% cells)
  }
  
  if (landmark == 1){
    #there are a few lincs dataframes in the datafolder but for some reason only this one has gene symbols...
    load(paste0(dataFolder,'lincs_signatures_cmpd_landmark_symbol.RData'))
  }else{
    #don't bother
    load(paste0(dataFolder,"lincs_signatures_cmpd_landmark_GSE92742.RData"))
  }
  
  if (choose_fda_drugs) {
    fda_drugs = read.csv(paste0(dataFolder,"repurposing_drugs_20170327.csv"), stringsAsFactors = F)
    lincs_sig_info_FDA <- subset(lincs_sig_info, id %in% colnames(lincs_signatures) & tolower(pert_iname) %in% tolower(fda_drugs$pert_iname))
    FDAdf <- select(lincs_sig_info_FDA, pert_id, pert_iname)
    FDAdf <- unique(FDAdf[,1:2])
    write.csv(FDAdf,file = paste0(outputFolder,"FDA_approved_drugs.csv"),row.names = F)
    lincs_sig_info <- lincs_sig_info %>% filter(id %in% colnames(lincs_signatures))
  }else{
    lincs_sig_info <- lincs_sig_info %>% filter(id %in% colnames(lincs_signatures))
  }

  
  #write paths
  output_path <- paste0(outputFolder, "/all_",paste(cells,collapse='_'),"_lincs_score.csv")
  sRGES_output_path <- paste0(outputFolder, "/sRGES",paste(cells,collapse='_'),".csv")
  sRGES_output_path_drug <- paste0(outputFolder, "/sRGES_FDAapproveddrugs.csv")
  dz_sig_output_path <- paste0(outputFolder, "/dz_sig_used.csv")
  
  #remove duplicate instances
  lincs_sig_info <- lincs_sig_info[!duplicated(lincs_sig_info$id),]
  sig.ids <- lincs_sig_info$id
  
  ####compute RGES####
  gene.list <- toupper(rownames(lincs_signatures))
  
  dz_signature <- dz_signature %>% filter(Symbol %in% gene.list)
  dz_genes_up <- dz_signature %>% filter(log2FoldChange>0) %>% arrange(desc(log2FoldChange))
  dz_genes_down <- dz_signature %>% filter(log2FoldChange<0) %>% arrange(log2FoldChange)
  ###############
  
  
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
  
  
  
  
  if(parallel){
    require(doParallel)
    require(foreach)
    no_cores <- as.integer(detectCores() / 4)
    registerDoParallel(cores=no_cores)
    
    
    #Billy and Ke 
    #dz_cmap_scores <- NULL
    #count <- 0
    
    ####slow loop####
    dz_cmap_scores = foreach(exp_id = sig.ids,.combine = 'c') %dopar% {
      #count <- count + 1
      #print(count)
      # progress(match(exp_id, sig.ids))
      #cmap_exp_signature is a dataframe with columns ids: whatever id we're using e.g. Symbol
      #rank which is the reverse rank of the expression of the expression
      cmap_exp_signature <- data.frame(gene.list,  
                                       rank(-1 * lincs_signatures[, as.character(exp_id)], 
                                            ties.method="random"))    
      colnames(cmap_exp_signature) <- c("ids","rank")
      #runs a function cmap_score_new from drugs_core_functions.R
      cmap_score_new(dz_genes_up$Symbol,dz_genes_down$Symbol,
                     cmap_exp_signature)
    }
    
    #random scores
    N_PERMUTATIONS <- 10000 #default 100000
    random_sig_ids <- sample(colnames(lincs_signatures),N_PERMUTATIONS,replace=T)
    count <- 0
    #random_cmap_scores <- NULL
    
    random_cmap_scores = foreach(exp_id = random_sig_ids,.combine = 'c')%dopar%{
      count <- count + 1
      #print(count)
      cmap_exp_signature <- data.frame(gene.list,  rank(-1 * lincs_signatures[, as.character(exp_id)], ties.method="random"))    
      colnames(cmap_exp_signature) <- c("ids","rank")
      
      random_input_signature_genes <- sample(gene.list, (nrow(dz_genes_up)+nrow(dz_genes_down)))
      rand_dz_gene_up <- data.frame(GeneID=random_input_signature_genes[1:nrow(dz_genes_up)])
      rand_dz_gene_down <- data.frame(GeneID=random_input_signature_genes[(nrow(dz_genes_up)+1):length(random_input_signature_genes)])
      cmap_score_new(rand_dz_gene_up,rand_dz_gene_down,cmap_exp_signature)
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
  } # end if(parallel)
  
  
  
  
  if(!parallel){
    require(lme4)
    require(Rfast)
    dz_cmap_scores <- NULL
    #count <- 0
    
    cmap_exp_sig <- Rfast::colRanks(-1 * lincs_signatures, method = "max")  
    names.list <- list(rownames(lincs_signatures),colnames(lincs_signatures))
    dimnames(cmap_exp_sig) <- names.list
    
    ####slow loop####
    for (exp_id in sig.ids) {
      #count <- count + 1
      #print(count)
      
      #cmap_exp_signature is a dataframe with columns ids: whatever id we're using e.g. Symbol
      #rank which is the reverse rank of the expression of the expression
      # cmap_exp_signature <- data.frame(gene.list,
      #                                 Rfast::colRanks(-1 * lincs_signatures[, as.character(exp_id)],
      #                                          method = "min"))    
      #  colnames(cmap_exp_signature) <- c("ids","rank")
      #runs a function cmap_score_new from drugs_core_functions.R
      cmap_exp_signature <- data.frame(ids = gene.list, rank = cmap_exp_sig[,exp_id])
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
      cmap_exp_signature <- data.frame(gene.list,  
                                       rank(-1 * lincs_signatures[, as.character(exp_id)], ties.method="random"))    
      colnames(cmap_exp_signature) <- c("ids","rank")
      # cmap_exp_sig <- Rfast::colRanks(-1 * lincs_signatures, method = "min")    
      # names.list <- list(rownames(lincs_sig),colnames(lincs_sig))
      # dimnames(cmap_exp_signature) <- names.list
      
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
  }
}

####### drug_enrichment #######
drug_enrichment <- function(sRGES,target_type){
  require(GSVA)
  require(limma)
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
  gsea_results = gsva(rgess, cmpdSets, method = "ssgsea",  parallel.sz=8,ssgsea.norm = FALSE)
  
  gsea_results = cbind(random_gsea_score[[target_type]], gsea_results)
  gsea_summary = data.frame(score = gsea_results[,17000+1])
  
  gsea_p = apply(gsea_results, 1, function(x){
    sum(x[1:ncol(random_gsea_score[[target_type]])] > x[ncol(random_gsea_score[[target_type]])+1])/ncol(random_gsea_score[[target_type]])
  })
  
  gsea_p = data.frame(target = names(gsea_p),score = gsea_summary, p = gsea_p, padj = p.adjust(gsea_p, method = 'fdr'))
  gsea_p = gsea_p[order(gsea_p$padj), ]
  # return(gsea_p)
  write.csv(gsea_p, paste0(enrichFolder.n, "/enriched_", target_type, ".csv"),row.names = F)
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
    if(length(clusternames) !=0){
    topclusterlist <- cmpdSets[clusternames]
    cat(sapply(topclusterlist, toString), file = paste0(enrichFolder.n,"misc.csv"), sep="\n")
    clusterdf <- read.csv2(paste0(enrichFolder.n,"misc.csv"), header=FALSE)
    clusterdf$cluster <- clusternames
    clusterdf$pval <- (gsea_p[which(gsea_p$padj<=0.05),])$padj
    colnames(clusterdf)[1] <- "drugs.in.cluster"
    write.csv(clusterdf,file=paste0(enrichFolder.n,'drugstructureclusters.csv'),row.names = F)
    }
  }
}


####### pick.out.cell.line #######
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


####### computeRefTissue #######
#computeRefTissue : returns the reference tissue given 'normal' and 'case' ids and 'expression' dataframe
## required input ##
# expSet: DF matrix with colnames consisting both the case and normal ids given
# case_id: case to select correlated tissues
# normal_id: all the normal tissues you want to check against the case
# outputFolder: this needs to be set to output the intermediate files for some visualization
## options ##
#control_size : if you're selecting the max number of tissues it gives back
#method: varGenes: will give the top varying ones, or random: random ones 
#cor_cutoff : will give back tissue based on top_nth quantile of the correlation
#must be a string in percents of 5s e.g. '5%' '10%' '25%' etc...
## Returns ##
# GTEXid : don't worry about the name this returns the highest varying normals
## Outputs ##
# /case_normal_corMatrix.csv : correlation matrix btw normal id and case id
# /case_normal_median_cor.csv : median correlation btw normal id and case id

computeRefTissue <- function(case_id = '', normal_id = '',
                             expSet=NULL,n_varGenes = 5000,
                             method='varGenes', #random 
                             #site_selection = 'top', 
                             #top site or any cor cutoff>quantile,
                             #all will select samples from 90th percentile
                             control_size = 50,
                             cor_cutoff='0%', #greater or equal than the cutoff 
                             output=T){
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


####### computeCellLine #######
computeCellLine <- function(case_id = '',
                            expSet = NULL,
                            LINCS_overlaps = T,
                            returnDF = F){
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
  TCGA.vs.CCLE.polyA.expression.correlation.result  <- pick.out.cell.line(case_counts, CCLE.overlaps,CCLE.rna.seq.marker.gene.1000)
  correlation.dataframe <- TCGA.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor %>% as.data.frame()
  colnames(correlation.dataframe) <- "cor"
  cell.line.folder <- paste0(outputFolder,"CellLineEval/")
  if (!dir.exists(cell.line.folder)) {
    dir.create(cell.line.folder)
  }
  write.csv(correlation.dataframe, file = paste0(cell.line.folder,"CellLineCorrelations.csv"))
  
  topline <- data.frame(medcor = TCGA.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor) # could also do first 
  # 3 of TCGA.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor
  
  
  if(returnDF == T)
  {
    return(topline)
  }else{
    topline <- rownames(topline)[1]
    # topline$cellLine = row.names(topline)
    # topline = topline %>% select(cellLine,medcor)
    # load(paste0(dataFolder,'CCLE_OCTAD.RData'))
    # topline = left_join(topline,CCLE_demoDat %>% select(CCLE.shortname,Age,Gender,Race,Site.Primary,Histology,Hist.Subtype1),
    #  by=c('cellLine'='CCLE.shortname'))
    return(topline)
  }
  
}

####### topLineEval #######
topLineEval <- function(topline = c(''),
                        mysRGES=''){
  #topline can be single cell or multiple
  #topline = c('2004','22RV1')
  toplineName = paste(topline,collapse='_')
  require(plotly)
  require(ggplot2)
  require(dplyr)
  require(data.table)
  cell.line.folder <- paste0(outputFolder,"CellLineEval/")
  if (!dir.exists(cell.line.folder)) {
    dir.create(cell.line.folder)
  }
  load(paste0(dataFolder,'CCLE_OCTAD.RData'))
  if(mysRGES==''){
    mysRGES <- read.csv(paste0(outputFolder,'sRGES.csv'),stringsAsFactors = F)
  }else(print('using'))
  mysRGES$pert_iname <- toupper(mysRGES$pert_iname)
  CTRPv2.ic50 <- dcast.data.table(CTRPv2.sensprof, drugid ~ cellid, value.var = "ic50_recomputed", fun.aggregate = median)
  #CTRPv2.ic50 <- replace(CTRPv2.ic50, is.na(CTRPv2.ic50), 0)
  colnames(CTRPv2.ic50) <- gsub("[^0-9A-Za-z///' ]","",colnames(CTRPv2.ic50))
  colnames(CTRPv2.ic50) <- toupper(colnames(CTRPv2.ic50))
  colnames(CTRPv2.ic50) <- gsub(" ","",colnames(CTRPv2.ic50))
  CTRP.IC50   <- CTRPv2.ic50[,c('DRUGID',..topline)]
  CTRP.IC50.m = melt(CTRP.IC50,id.vars = 'DRUGID')
  CTRP.IC50.medianIC50 <- CTRP.IC50.m[,.(medIC50=median(value)), by=DRUGID]
  CTRP.IC50.medianIC50 <- CTRP.IC50.medianIC50[is.finite(CTRP.IC50.medianIC50$medIC50),]
  
  
  
  CTRPv2.auc   <- dcast.data.table(CTRPv2.sensprof, drugid ~ cellid, value.var = "auc_recomputed", fun.aggregate = median)
  #CTRPv2.auc <- replace(CTRPv2.auc, is.na(CTRPv2.auc), 0)
  colnames(CTRPv2.auc) <- gsub("[^0-9A-Za-z///' ]","",colnames(CTRPv2.auc))
  colnames(CTRPv2.auc) <- toupper(colnames(CTRPv2.auc))
  colnames(CTRPv2.auc) <- gsub(" ","",colnames(CTRPv2.auc))
  CTRP.auc   <- CTRPv2.auc[,c('DRUGID',..topline)]
  CTRP.auc.m = melt(CTRP.auc,id.vars = 'DRUGID')
  CTRP.auc.medianauc <- CTRP.auc.m[,.(medauc=median(value)), by=DRUGID]
  CTRP.auc.medianauc <- CTRP.auc.medianauc[is.finite(CTRP.auc.medianauc$medauc),]
  
  mysRGES$pert_iname = toupper(mysRGES$pert_iname)
  CTRP.IC50.medianIC50$DRUGID = toupper(CTRP.IC50.medianIC50$DRUGID)
  
  testdf <- merge(mysRGES, CTRP.IC50.medianIC50, by.x = "pert_iname", by.y = "DRUGID")
  IC50.cortest <- cor.test(testdf$sRGES,log10(testdf$medIC50))
  ic50pval <- IC50.cortest$p.value
  ic50rho  <- IC50.cortest$estimate
  mylabel = c("p-value" = ic50pval,'Rho' = ic50rho)

  testdf$StronglyPredicted <- NA
  testdf$StronglyPredicted <- ifelse(testdf$sRGES < -0.2,"Yes","No")
  
  StronglyPredicted <- testdf$StronglyPredicted
  
  Legend.title <- "Strongly <br>Predicted"
  Legend.label1<- "No" 
  Legend.label2<- "Yes"
  Title <- "Top Line recomputed log(ic50) vs sRGES"
  xaxis <- "sRGES"
  yaxis <- "log(ic50)"
  
  p <- ggplot() +
    geom_point(data = testdf, 
               aes(x = testdf$sRGES, y = log10(testdf$medIC50), 
                   color = StronglyPredicted, 
                   text = 
                     paste('Drug: ', pert_iname, # These are the hover labels generated by p1
                           '<br>sRGES: ', sRGES))) +
    geom_smooth(data = testdf,          
                aes(x = testdf$sRGES, y = log10(testdf$medIC50),
                    label = ic50rho),
                method = 'lm',se = F, color = "black",size = 0.5) +
    scale_color_discrete(
      name = Legend.title) +
    labs(x = xaxis, y = yaxis, title = Title) +
    theme(legend.position = "right", legend.background = element_rect(fill="#F5F5F5"), legend.title = element_blank())
  
  p1 <- ggplotly(p, tooltip = c("text", "label")) %>%  layout(margin=list(l=15)) # change to white background
  
  ic50graph <- p1 %>% 
    add_annotations( text=Legend.title,xref="paper",yref="paper",x=1.02,xanchor="left",y=0.8,yanchor="bottom",legendtitle=TRUE,showarrow=FALSE ) %>%
    layout( legend=list(y=0.8, yanchor="top" ) )
  
  htmlwidgets::saveWidget(as_widget(ic50graph), paste0(cell.line.folder,toplineName,"_ic50_insilico_validation.html"),selfcontained = F) 
  # note: files are simply too large to set selfcontained = T. This just causes issues on linux machines.
  
  
  # AUC
  testdf2 <- merge(mysRGES, CTRP.auc.medianauc, by.x = "pert_iname", by.y = "DRUGID")
  AUC.cortest <- cor.test(testdf2$sRGES,testdf2[,8])
  # plot(testdf2$sRGES, testdf2[,8],xlab = "sRGES",ylab=paste(topline,"AUC"), main = "Top Line recomputed AUC vs sRGES")
  
  testdf2$StronglyPredicted <- NA
  testdf2$StronglyPredicted <- ifelse(testdf2$sRGES < -0.2,"Yes","No")
  
  StronglyPredicted <- testdf2$StronglyPredicted
  
  aucpval <- AUC.cortest$p.value
  aucrho  <- AUC.cortest$estimate
  
  Legend.title <- "Strongly <br>Predicted"
  Legend.label1<- "No"
  Legend.label2<- "Yes"
  Title <- "Top Line recomputed AUC vs sRGES"
  xaxis <- "sRGES"
  yaxis <- "AUC"
  
  
  p <- ggplot() +
    geom_point(data = testdf2, 
               aes(x = testdf2$sRGES, y = testdf2[,8], 
                   color = StronglyPredicted, 
                   text = 
                     paste('Drug: ', pert_iname, # These are the hover labels generated by p1
                           '<br>sRGES: ', sRGES))) +
    geom_smooth(data = testdf2,          
                aes(x = testdf2$sRGES, y = testdf2[,8],
                    label = aucrho),
                method = 'lm',se = F, color = "black",size = 0.5) +
    scale_color_discrete(
      name = Legend.title) +
    labs(x = xaxis, y = yaxis, title = Title) +
    theme(legend.position = "right", legend.background = element_rect(fill="#F5F5F5"), legend.title = element_blank())
  
  p1 <- ggplotly(p, tooltip = c("text", "label")) %>%  layout(margin=list(l=15))
  
  aucgraph <- p1 %>% 
    add_annotations( text=Legend.title,xref="paper",yref="paper",x=1.02,xanchor="left",y=0.8,yanchor="bottom",legendtitle=TRUE,showarrow=FALSE ) %>%
    layout( legend=list(y=0.8, yanchor="top" ) )
  
  htmlwidgets::saveWidget(as_widget(aucgraph), paste0(cell.line.folder,toplineName,"_auc_insilico_validation.html"),selfcontained = F) # test
  
  
  # logging cortests
  con <- file(paste0(cell.line.folder,toplineName,"_drug_sensitivity_insilico_results.txt"))
  sink(con, append=TRUE)
  sink(con, append=TRUE, type="message")
  
  print("AUC cortest")
  print(AUC.cortest)
  print("IC50 cortest")
  print(IC50.cortest)
  
  sink() 
  sink(type="message")
  
}


####### remLowExpr #######
remLowExpr <- function(counts,counts_phenotype){
  x <-DGEList(counts = round(counts), group = counts_phenotype$sample_type )
  cpm_x <- cpm(x)
  #needs to be at least larger the than the size of the smallest set
  keep.exprs <- rowSums(cpm_x>1) >= min(table(counts_phenotype$sample_type)) 
  keep.exprs
}

####### compEmpContGenes #######
#compute empirical control genes for RUVg
#consider deprecating since this is just first pass DE
compEmpContGenes <- function(counts, coldata, n_topGenes = 5000){
  require(edgeR)
  require(RUVSeq)
  set <- newSeqExpressionSet(round(counts),
                             phenoData = data.frame(coldata,row.names= coldata$sample_id))
  
  design <- model.matrix(~ condition, data = pData(set))
  y <- DGEList(counts=counts(set), group =  coldata$condition)
  y <- calcNormFactors(y, method="TMM") #upperquartile generate Inf in the LGG case
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit) #defaults to compare tumor to normal or tumor mutant to normal
  
  top <- topTags(lrt, n=nrow(set))$table
  #n_topGenes <- n_topGenes #5000: assume there are 5000 signficant genes
  
  #based on n_topGenes computing genes with low DE
  #the genes not computed significant is in the empirical set
  empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:n_topGenes]))]
  empirical
}

####### diffExp #######
#modify v2 version, suport edgeR, Limma and deseq 2 now.
diffExp <- function(case_id='',control_id='',expSet=dz_expr,
                    normalize_samples=T,
                    k=1,
                    n_topGenes=10000,
                    DE_method='edgeR',
                    parallel_cores = 2){
  require(dplyr)
  require(RUVSeq)
  require(edgeR)
  counts_phenotype <- rbind(data.frame(sample = case_id,sample_type = 'case'),
                            data.frame(sample = control_id, sample_type = 'control'))
  counts_phenotype = counts_phenotype[counts_phenotype$sample %in% colnames(expSet),]
  
  counts = expSet[,as.character(counts_phenotype$sample)]
  counts = 2^counts - 1 #unlog the counts it was log(2x + 1) in dz.expr.log2.readCounts
  counts_phenotype$sample = as.character(counts_phenotype$sample)
  counts_phenotype$sample_type = factor(counts_phenotype$sample_type, levels = c("control", "case"))
  
  #remove lowly expressed transcripts
  highExpGenes <- remLowExpr(counts,counts_phenotype)
  counts = counts[highExpGenes,]
  
  set <- newSeqExpressionSet(round(counts),
                             phenoData = data.frame(counts_phenotype,row.names=counts_phenotype$sample))
  
  #normalize samples using RUVSeq
  if (normalize_samples == T){
    #compute empirical genes
    
    design <- model.matrix(~ sample_type, data = pData(set))
    y <- DGEList(counts=counts(set), group =  counts_phenotype$sample)
    y <- calcNormFactors(y, method="TMM") #upperquartile generate Inf in the LGG case
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit,2) #defaults to compare case control
    
    top <- topTags(lrt, n=nrow(set))$table
    i = which(!(rownames(set) %in% rownames(top)[1:min(n_topGenes,dim(top)[1])]))
    empirical <- rownames(set)[i]
    stopifnot(length(empirical)>0)
    write.csv(data.frame(empirical),file=paste0(outputFolder,"computedEmpGenes.csv"))
    set1 <- RUVg(set,empirical,k=k)
  }
  
  if(DE_method=='DESeq2'){
    library('DESeq2')
    print('computing DE via DESeq')
    row.names(counts_phenotype) = counts_phenotype$sample
    coldata = counts_phenotype
    
    if (normalize_samples == T){
      dds <- DESeqDataSetFromMatrix(countData = counts(set1),
                                    colData = pData(set1),
                                    design= ~ sample_type + W_1)
    }else{
      dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                    colData = coldata,
                                    design= ~ sample_type)
    }
    gc()
    
    if (parallel_cores > 1){
      dds <- DESeq(dds, parallel = T)
    }else{
      dds <- DESeq(dds)
    }
    
    #save(dds,file= paste0(outputFolder, "/dds", ".RData"))
    rnms <- resultsNames(dds)
    
    resRaw <- results(dds, contrast=c("sample_type","case","control"))
    res = data.frame(resRaw)
    res$identifier <- row.names(res)
    res = res %>% select(identifier,everything())
    return(res)
    
  }else if(DE_method=='edgeR'){
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
      dgList <- DGEList(counts=counts(set1),group=set1$sample_type)
      
    }else{
      design <- model.matrix(~ sample_type,data=pData(set))
      dgList <- DGEList(counts=counts(set),group=set$sample_type)
    }
    
    dgList <- calcNormFactors(dgList, method="TMM") #using upperquartile seems to give issue for LGG
    dgList <- estimateGLMCommonDisp(dgList, design)
    dgList <- estimateGLMTagwiseDisp(dgList, design)
    fit <- glmFit(dgList, design)
    
    #see edgeRUsersGuide section on testing for DE genes for contrast
    lrt <- glmLRT(fit,2) 
    #second coefficient otherwise it'll default the W_1 term when normalize is on
    
    res <- lrt$table
    colnames(res) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
    res$padj <- p.adjust(res$pvalue)
    res$identifier <- row.names(res)
    res = res %>% select(identifier,everything())
    return(res)
  }else if(DE_method=='limma'){
    #according to https://support.bioconductor.org/p/86461/, LIMMA + VOOM will not use normalized data
    print('computing DE via limma')
    require('Glimma')
    x <- counts
    nsamples <- ncol(x)
    lcpm <- cpm(x, log=TRUE)
    
    group <-counts_phenotype$sample_type
    
    ## ----design-----------------------------------------------------------------------------
    design <- model.matrix(~0 + group)
    colnames(design) <- gsub("group", "", colnames(design))
    
    contr.matrix <- makeContrasts(
      TumorvsNon = case - control, 
      levels = colnames(design))
    
    v <- voom(x, design, plot=F)
    vfit <- lmFit(v, design)
    vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
    efit <- eBayes(vfit)
    #plotSA(efit, main="Final model: Mean−variance trend")
    
    tfit <- treat(vfit, lfc=1)
    dt <- decideTests(tfit)
    summary(dt)
    
    tumorvsnormal <- topTreat(tfit, coef=1, n=Inf)
    tumorvsnormal <- tumorvsnormal[order(abs(tumorvsnormal$logFC), decreasing = T),]
    #tumorvsnormal.topgenes <- rownames(tumorvsnormal[1:50,])
    
    
    res <-tumorvsnormal
    colnames(res) <-c("log2FoldChange", "AveExpr", "t", "pvalue", "padj")
    res$identifier = row.names(res)
    return(res)
  }
  
}

# diffExp_old <- function(case_id='',control_id='',expSet=dz_expr,
#                     normalize_samples=T,
#                     k=1,
#                     n_topGenes=10000,
#                     DE_method='edgeR'){
#   require(dplyr)
#   require(RUVSeq)
#   require(edgeR)
#   counts_phenotype <- rbind(data.frame(sample = case_id,sample_type = 'case'),
#                             data.frame(sample = control_id, sample_type = 'control'))
#   counts = expSet[,as.character(counts_phenotype$sample)]
#   counts = 2^counts - 1 #unlog the counts it was log(2x + 1) in dz.expr.log2.readCounts
#   counts_phenotype$sample = as.character(counts_phenotype$sample)
#   counts_phenotype$sample_type = factor(counts_phenotype$sample_type, levels = c("control", "case"))
#   highExpGenes <- remLowExpr(counts,counts_phenotype)
#   counts = counts[highExpGenes,]
#   
#   
#   set <- newSeqExpressionSet(round(counts),
#                              phenoData = data.frame(counts_phenotype,row.names=counts_phenotype$sample))
#   design <- model.matrix(~ sample_type, data = pData(set))
#   y <- DGEList(counts=counts(set), group =  counts_phenotype$sample)
#   y <- calcNormFactors(y, method="TMM") #upperquartile generate Inf in the LGG case
#   y <- estimateGLMCommonDisp(y, design)
#   y <- estimateGLMTagwiseDisp(y, design)
#   fit <- glmFit(y, design)
#   lrt <- glmLRT(fit,2) #defaults to compare case control
#   
#   #counts dataframe === remove low counts ===> set === normalized ===> set1
#   #if no empirical genes found it will just create a matrix without running RUVg
#   if(normalize_samples == T){
#     top <- topTags(lrt, n=nrow(set))$table
#     i = which(!(rownames(set) %in% rownames(top)[1:min(n_topGenes,dim(top)[1])]))
#     empirical <- rownames(set)[i]
#     stopifnot(length(empirical)>0)
#     write.csv(data.frame(empirical),file=paste0(outputFolder,"computedEmpGenes.csv"))
#     set1 <- RUVg(set,empirical,k=k)
#     
#     print('computing DE via edgeR')
#     
#     #construct model matrix based on whether there was normalization ran
#     if(normalize_samples == T){
#       if(k==1){
#         design <- model.matrix(~sample_type + W_1, data=pData(set1))
#       }else if(k == 2){
#         design <- model.matrix(~sample_type + W_1 + W_2, data = pData(set1))
#       }else if (k == 3){
#         design <- model.matrix(~sample_type + W_1 + W_2 + W_3, data = pData(set1))
#       }
#     }else{design <- model.matrix(~sample_type,data=pData(set1))}
#     
#     dgList <- DGEList(counts=counts(set1),group=set1$sample_type)
#     dgList <- calcNormFactors(dgList, method="TMM") #using upperquartile seems to give issue for LGG
#     dgList <- estimateGLMCommonDisp(dgList, design)
#     dgList <- estimateGLMTagwiseDisp(dgList, design)
#     fit <- glmFit(dgList, design)
#     
#     #see edgeRUsersGuide section on testing for DE genes for contrast
#     lrt <- glmLRT(fit,2) 
#     #second coefficient otherwise it'll default the W_1 term when normalize is on
#   }
#   res <- lrt$table
#   colnames(res) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
#   res$padj <- p.adjust(res$pvalue)
#   res$identifier <- row.names(res)
#   res = res %>% select(identifier,everything())
#   return(res)
# }
####### diffExp #######
diffExp_v1 <- function(case_id='',control_id='',expSet=dz_expr,
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

# new one with limma & DESeq2
# note: normalize_samples only works with edgeR right now.
diffExp_v2 <- function(case_id='',control_id='',expSet=dz_expr,
                      normalize_samples=T,
                      k=1,
                      n_topGenes=10000,
                      DE_method='edgeR',
                      parallel_cores = 2){
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
  orderSeq = colnames(counts)
  orderSeq2 = counts_phenotype$sample
  counts_phenotype = counts_phenotype %>% filter(sample==orderSeq)
  
  if(DE_method=='DESeq2'){
    library('DESeq2')
    print('computing DE via DESeq')
    row.names(counts_phenotype) = counts_phenotype$sample
    coldata = counts_phenotype
    dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                  colData = coldata,
                                  design= ~ sample_type)
    gc()
    if (parallel_cores > 1){
      dds <- DESeq(dds, parallel = T)
    }else{
      dds <- DESeq(dds)
    }
    
    #save(dds,file= paste0(outputFolder, "/dds", ".RData"))
    rnms <- resultsNames(dds)
    
    resRaw <- results(dds, contrast=c("sample_type","case","control"))
    res = data.frame(resRaw)
    res$identifier <- row.names(res)
    res = res %>% select(identifier,everything())
    return(res)
    
  }else if(DE_method=='edgeR'){
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
  }else if(DE_method=='limma'){
    print('computing DE via limma')
    require('Glimma')
    x <- counts
    nsamples <- ncol(x)
    lcpm <- cpm(x, log=TRUE)
    
    group <-counts_phenotype$sample_type
    
    ## ----design-----------------------------------------------------------------------------
    design <- model.matrix(~0 + group)
    colnames(design) <- gsub("group", "", colnames(design))
    
    contr.matrix <- makeContrasts(
      TumorvsNon = case - control, 
      levels = colnames(design))
    
    v <- voom(x, design, plot=F)
    vfit <- lmFit(v, design)
    vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
    efit <- eBayes(vfit)
    #plotSA(efit, main="Final model: Mean−variance trend")
    
    tfit <- treat(vfit, lfc=1)
    dt <- decideTests(tfit)
    summary(dt)
    
    tumorvsnormal <- topTreat(tfit, coef=1, n=Inf)
    tumorvsnormal <- tumorvsnormal[order(abs(tumorvsnormal$logFC), decreasing = T),]
    #tumorvsnormal.topgenes <- rownames(tumorvsnormal[1:50,])
    
    
    res <-tumorvsnormal
    colnames(res) <-c("log2FoldChange", "AveExpr", "t", "pvalue", "padj")
    res$identifier = row.names(res)
    return(res)
  }
  
  
  
}


####### geneEnrich #######


#' geneEnrichOld <- function(dz_signature, 
#'                        db.list=c( "KEGG_2016",  "GO_Biological_Process_2017", "GO_Cellular_Component_2017"),
#'                        suffix=''){
#'   #adapted from https://github.com/compbiomed/enrichR
#'   #edited for use by Billy Zeng, UCSF, Anita Wen UCD
#'   #' Perform functional enrichment on a set of genes.
#'   #' 
#'   #' This function interacts with Enrichr's REST API in order to perform functional enrichment of a single
#'   #' set of genes, for a set of specified databases which are already fronted by Enrichr.
#'   #' Databases are specified as seen in the web interface, with underscores for spaces
#'   #' (e.g. "WikiPathways_2016", "KEGG_2016", "GO_Biological_Process"). There's no way to query Enrichr
#'   #' to get these database names, so they can't be provided as options. You'll just have to guess. Sorry :/
#'   #' 
#'   #' @param gene.list a list of gene symbols
#'   #' @param databases a list of Enrichr-fronted databases, as mentioned above
#'   #' @param fdr.cutoff An FDR (adjusted p-value) threshold by which to limit the list of enriched pathways
#'   #' @keywords functional enrichment Enrichr
#'   #' @export
#'   
#'   enrichGeneList <- function (gene.list, databases=db.list, fdr.cutoff=NULL) {
#'     ######Step 1: Post gene list to EnrichR
#'     req.body <- list(list=paste(gene.list, collapse="\n"))
#'     post.req <- httr::POST("http://amp.pharm.mssm.edu/Enrichr/enrich", encode="multipart", body=I(req.body))
#'     
#'     #TODO: Real error handling..
#'     if (!grepl("success", httr::http_status(post.req)$category, ignore.case=T)) stop("Posting gene list to EnrichR failed")
#'     
#'     ######Step 2: Get results from posted gene list
#'     database.enrichments <- list()
#'     for (idx in 1:length(databases)) { 
#'       database <- databases[idx]
#'       get.req <- httr::GET(paste("http://amp.pharm.mssm.edu/Enrichr/enrich?backgroundType=", database, sep=""))
#'       if (!grepl("success", httr::http_status(get.req)$category, ignore.case=T)) stop("Retrieving results from EnrichR failed")
#'       
#'       response.content <- mungeResponseContent(httr::content(get.req)[[database]])
#'       
#'       if (length(response.content) > 1) {
#'         database.res <- data.table::rbindlist(response.content)
#'         database.res[, 1] <- rep(database, nrow(database.res))
#'         database.enrichments[[idx]] <- database.res[, paste("V", c(1, 2, 3, 7, 5,6), sep=""), with=F]
#'       }
#'     }
#'     
#'     query.results <- as.data.frame(data.table::rbindlist(database.enrichments))
#'     colnames(query.results) <- c("database", "category", "pval", "qval","Combined Score" ,"genes")
#'     
#'     if (!is.null(fdr.cutoff)) {
#'       query.results <- query.results[query.results$qval < fdr.cutoff, ]
#'     }
#'     
#'     query.results
#'   }
#'   
#'   
#'   
#'   #' Munge the Enrichr API response so it'll squeeze neatly (if untidily) into a dataframe.
#'   #' 
#'   #' The response from the Enrichr API is a list of lists, where each nested list item represents an enriched
#'   #' category. The 6th item of each category (i.e. response.content[[category.idx]][[6]]) corresponds to the 
#'   #' genes that overlapped with the gene set behind that category. This function bascically collapses that list of 
#'   #' genes into a single string. 
#'   #' 
#'   #' I'm sorry you ever had to look at this.
#'   #' 
#'   #' @param response.content result of calling httr::content on the GET request to the Enrichr API, after submitting a list for enrichment.
#'   #' 
#'   mungeResponseContent <- function (response.content) {
#'     munged.content <- response.content
#'     if (length(response.content) == 0) return(NA)
#'     
#'     for (idx in 1:length(response.content)) {
#'       munged.content[[idx]][[6]] <- paste(munged.content[[idx]][[6]], collapse=",")
#'     }
#'     
#'     munged.content
#'   }
#'   
#'   
#'   #munged.content not found
#' 
#'   enrichFolder <- paste0(outputFolder,'enrichment_analysis/')
#'   if (!dir.exists(enrichFolder)) {
#'     dir.create(enrichFolder)
#'   }
#'   enrichFolder.n <- paste0(enrichFolder,'GeneEnrichment','/')
#'   if (!dir.exists(enrichFolder.n)) {
#'     dir.create(enrichFolder.n)
#'   }
#'   
#'   
#'   databases = db.list
#'   fdr.cutoff = 0.25
#'   
#'   up.genes = dz_signature$Symbol[dz_signature$log2FoldChange > 0]
#'   if(length(up.genes)>0){
#'     (up.gene.res = enrichGeneList(up.genes[1:min(300, length(up.genes))], databases, fdr.cutoff))
#'     (up.gene.res = up.gene.res[order(up.gene.res$database, up.gene.res$p), ])
#'     write.csv(up.gene.res, paste0(enrichFolder.n, "/dz_up_sig_genes_enriched", suffix,".csv"))  
#'   }
#'   dn.genes = dz_signature$Symbol[dz_signature$log2FoldChange < 0]
#'   if(length(dn.genes)>0){
#'     dn.gene.res = enrichGeneList(dn.genes[1:min(300, length(dn.genes))], databases, fdr.cutoff)
#'     dn.gene.res = dn.gene.res[order(dn.gene.res$database, dn.gene.res$p), ]
#'     write.csv(dn.gene.res, paste0(enrichFolder.n, "/dz_down_sig_genes_enriched", suffix,".csv"))  
#'   }
#'   up.gene.res$dir = 'up'
#'   dn.gene.res$dir = 'dn'
#'   res = rbind(up.gene.res,dn.gene.res)
#' }
#' 
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
  #' @keywords functional enrichment Enrichr
  #' @export
  
  enrichGeneList <- function (gene.list, databases=db.list, fdr.cutoff=NULL) {
    ######Step 1: Post gene list to EnrichR
    req.body <- list(list=paste(gene.list, collapse="\n"))
    post.req <- httr::POST("http://amp.pharm.mssm.edu/Enrichr/addList", 
                           encode="multipart", 
                           body=I(req.body))
    results = data.frame(matrix(ncol = 6,nrow=0))
    colnames(results) = c('database',
                          'category',
                          'pval_adj',
                          'Odds_Ratio',
                          'Combined Score',
                          'Genes')
    for(db in db.list){
      getEnrich = GET('http://amp.pharm.mssm.edu/Enrichr/export',
                      query=list(userListId = ids$userListId,
                                 filename='example',
                                 backgroundType='KEGG_2016'))
      query.results <- readr::read_tsv(content(getEnrich,type='text/tsv'))
      query.results$database = db
      query.results = query.results %>% select(database,
                                               category=Term,
                                               pval_adj=`Adjusted P-value`,
                                               Odds_Ratio = `Odds Ratio`,
                                               `Combined Score`,
                                              Genes)
      results = rbind(results,query.results)
    }
    results
  }
  
  enrichFolder <- paste0(outputFolder,'enrichment_analysis/')
  if (!dir.exists(enrichFolder)) {
    dir.create(enrichFolder)
  }
  enrichFolder.n <- paste0(enrichFolder,'GeneEnrichment','/')
  if (!dir.exists(enrichFolder.n)) {
    dir.create(enrichFolder.n)
  }
  
  
  databases = db.list

  up.genes = dz_signature$Symbol[dz_signature$log2FoldChange > 0]
  if(length(up.genes)>0){
    (up.gene.res = enrichGeneList(up.genes[1:min(300, length(up.genes))], databases))
    (up.gene.res = up.gene.res[order(up.gene.res$database, up.gene.res$pval_adj), ])
    write.csv(up.gene.res, paste0(enrichFolder.n, "/dz_up_sig_genes_enriched", suffix,".csv"))  
  }
  dn.genes = dz_signature$Symbol[dz_signature$log2FoldChange < 0]
  if(length(dn.genes)>0){
    dn.gene.res = enrichGeneList(dn.genes[1:min(300, length(dn.genes))], databases)
    dn.gene.res = dn.gene.res[order(dn.gene.res$database, dn.gene.res$pval_adj), ]
    write.csv(dn.gene.res, paste0(enrichFolder.n, "/dz_down_sig_genes_enriched", suffix,".csv"))  
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



####### visualize_drug_hits #######
visualize_drug_hits <- function(){
  require(pheatmap)
  require(dplyr)
  require("gplots")
  require("ggplot2")
  require("RColorBrewer")
  
  #CMAP score output
  
  lincs_drug_prediction = read.csv(paste(outputFolder, "/all_lincs_score.csv", sep=""),stringsAsFactors = F)
  dz_gene_ids = read.csv(paste(outputFolder, "/dz_sig_used.csv", sep=""),stringsAsFactors = F)
  
  dz_signature  = read.csv(paste0(outputFolder,'dz_signature.csv'),stringsAsFactors = F)
  dz_signature = dz_signature %>% filter(gene %in% dz_gene_ids$gene)
  
  dz_sig = dz_signature %>% select(gene,log2FoldChange)
  
  #drugs to visualize
  sRGES = read.csv(paste(outputFolder, "/sRGES.csv", sep=""))
  #choose top common drugs
  top_drugs = as.character(unique(sRGES$pert_iname[sRGES$n > 1 & !sRGES$pert_iname %in% sRGES$pert_iname[grep("BRD-|SA-", sRGES$pert_iname)]][1:20]))
  
  ##########
  #visualized gene reversed
  #only pick the signatures from close to the median
  #cell_lines = read.csv(paste("../table/", cancer, "_cell_lines.csv", sep=""))
  #lincs_drug_prediction_subset = subset(lincs_drug_prediction, cell_id %in% cell_lines$LINCS, select=c("id", "RGES", "pert_iname", "pert_dose", "pert_time"))
  lincs_drug_prediction$RGES = lincs_drug_prediction$cmap_score
  lincs_drug_prediction_subset = lincs_drug_prediction[lincs_drug_prediction$pert_iname %in% top_drugs,]
  ###selecting median still sounds weird... let's keep all signatures 
  drug_cmap_score = aggregate(RGES ~ pert_iname, lincs_drug_prediction_subset, median)
  drug_instances_median  = merge(lincs_drug_prediction_subset, drug_cmap_score, by = c("pert_iname"))
  drug_instances_median$diff = abs(drug_instances_median$RGES.x - drug_instances_median$RGES.y) #cmap_score.y is the median
  drug_instances_min_diff = aggregate(diff ~ pert_iname, drug_instances_median, min)
  drug_instances_select = merge(drug_instances_median, drug_instances_min_diff, by=c("pert_iname", "diff"))
  drug_instances_select = drug_instances_select[!duplicated(drug_instances_select$pert_iname), ]
  sig_id_selects = drug_instances_select$id
  
  #sig_id_selects = lincs_drug_prediction_subset$id
  landmark = 1
  if (landmark == 1){
    load(paste0(dataFolder,"lincs_signatures_cmpd_landmark_symbol.RData"))
  }else{
    load(paste0(dataFolder,"lincs_signatures_cmpd_landmark_GSE92742.RData"))
  }
  
  drug_dz_signature = merge(dz_sig, data.frame(gene = rownames(lincs_signatures), 
                                               lincs_signatures[, as.character(sig_id_selects)]),  by="gene", suffixes='')
  
  
  #########################
  ###
  #visualize the reversed gene expression
  #reorder drug dz signatures
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
    lincs_drug_prediction$pert_iname[paste("X",lincs_drug_prediction$id, sep="") == names(drug_dz_signature_rank)[id]][1]
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


####### id_mapping #######
id_mapping <-function(id, input_format = "hgnc_symbol", output_format = "ensembl_gene_id"){
  #ensembl_gene_id ensembl_transcript_id hgnc_symbol entrezgene
  library(biomaRt)
  
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  id_mapped <- getBM(attributes=c(input_format, output_format), filters =
                       input_format, values =id, mart = ensembl)
  return(id_mapped[1,2]) #only return the first hit
}

####### id_mapping_gene_ensembl #######
id_mapping_gene_ensembl <-function(id){
  #ensembl_gene_id ensembl_transcript_id hgnc_symbol entrezgene
  #mapping <- read.csv("~/Documents/GitHub/OCTAD/raw/gencode.v23.annotation.gene.probeMap.csv", stringsAsFactors = F)
  mapping <- read.csv(paste0(dataFolder,'raw/gencode.v23.annotation.gene.probeMap.csv'),stringsAsFactors = F)
  return(mapping$ensembl[mapping$gene == id][1])
}

####### queryGDC #######
queryGDC <- function(GENE="", PROJECT){
  library(curl)
  
  
  ## examples
  #GENE = id_mapping_gene_ensembl('IDH2')
  #PROJECT = 'TCGA-LGG'
  if (GENE!="") {
    all_gene<-curl_fetch_memory(paste0("https://api.gdc.cancer.gov/analysis/top_mutated_cases_by_gene?fields=submitter_id&pretty=true&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22genes.gene_id%22%2C%22value%22%3A%5B%22", GENE, "%22%5D%7D%7D%5D%7D&format=tsv&size=100000"))
    all_gene.df=read.table(text = rawToChar(all_gene$content), sep = '\t', header = TRUE)    
    all_gene.df$submitter_id=as.character(all_gene.df$submitter_id)    
  }
  
  ## retrieve all IDs with in PROJECT
  all_project<-curl_fetch_memory(paste0("https://api.gdc.cancer.gov/analysis/top_mutated_cases_by_gene?fields=submitter_id&pretty=true&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22",PROJECT,"%22%5D%7D%7D%5D%7D&format=tsv&size=100000"))
  all_project.df=read.table(text = rawToChar(all_project$content), sep = '\t', header = TRUE)    
  all_project.df$submitter_id=as.character(all_project.df$submitter_id)
  
  if(GENE!=""){
    ## combine into grouped matrix
    combined=as.data.frame(unique(append(unique(all_project.df$submitter_id),unique(all_gene.df$submitter_id))))
    colnames(combined)[1]="submitter_id"
  }else{combined=as.data.frame(unique(all_project.df$submitter_id))
  colnames(combined)[1]="submitter_id"}
  
  # add in gene info
  # is 1 
  if(GENE!=''){
    combined$gene=0
    combined[which(combined$submitter_id %in% all_gene.df$submitter_id),"gene"]=1
  }
  
  # add in project info
  combined$project=0
  combined[which(combined$submitter_id %in% all_project.df$submitter_id),"project"]=1  
  
  return(combined)
}

####### estimatePurity #######
estimatePurity  <- function(expr_matrix){
  library(estimate)
  
  #expr_matrix: with gene symbols as row names and samples as colnames
  #return purity score
  samples <- data.frame(NAME = rownames(expr_matrix), Description = NA,  expr_matrix)
  write("", file = "temp_samples.gct")
  write(paste(nrow(samples), "\t", ncol(samples) -2), file = "temp_samples.gct", append = T)
  write("", file = "temp_samples.gct", append = T)
  write.table(samples,  file = "temp_samples.gct", append = T, row.names=F, col.names=T, sep="\t", quote=F)
  
  in.file <- ("temp_samples.gct")
  out.file <- "temp_samples_output.gct"
  estimateScore(in.file, out.file)
  
  estimateScore <- read.delim(out.file, sep = "\t", skip = 3)
  
  return(as.numeric(estimateScore[3, -c(1,2)]))
}

####### detectOutlier #######
detectOutlier <- function(expSet,z_threshold=3,outlierPdf="/tissue_mds.pdf"){
  #refer to https://www.biostars.org/p/281767/
  x <-DGEList(counts = round(2^expSet - 1))
  x <- calcNormFactors(x, method = "TMM")
  lcpm <- cpm(x, log=TRUE)
  pca <- prcomp(t(lcpm))
  #plot(pca$x, pch = 20)
  pc1_z_score = as.numeric(scale(pca$x[,1]))
  pc2_z_score = as.numeric(scale(pca$x[,2]))
  outliers = rownames(pca$x)[abs(pc1_z_score) > z_threshold] #
  col.group = rep("black", length(pc1_z_score))
  col.group[rownames(pca$x) %in% outliers] = "red"
  pdf(paste0(outputFolder, outlierPdf))
  plot(pc1_z_score, pc2_z_score, xlab = "PC1", ylab = "PC2",col = col.group, pch = row.names(x$samples))
  dev.off()
  outliers
}

####### compute_tissue_cell_cor #######
compute_tissue_cell_cor <- function(dz_tissue_samples){
  load(paste0(dataFolder,"tissue_cell_line_cor.RData"))
  tumor_cell_cor <- tissue_cell_line_cor[, colnames(tissue_cell_line_cor) %in% dz_tissue_samples]
  #order based on median cor
  tumor_cell_cor_merged <- apply(tumor_cell_cor, 1, median)
  top_cell_lines <- names(head(sort(tumor_cell_cor_merged, decreasing = T), 15))
  tumor_cell_cor <- tumor_cell_cor[top_cell_lines, ]
  #tumor_cell_cor$tumor_type_name = factor(tumor_cell_cor$tumor_type_name, levels = tumor_cell_cor_merged$tumor_type_name)
  
  pdf(paste0(outputFolder, "/top_cell_lines.pdf"))
  par(mar=c(15,4.1,1.1,2.1))
  boxplot(t(tumor_cell_cor), las=2, cex.axis=1)
  dev.off()
}

####### compute_tissue_lincs_cell_cor #######
compute_tissue_lincs_cell_cor <- function(dz_tissue_samples){
  load(paste0(dataFolder,"tissue_cell_line_cor.RData"))
  ccle_mapping <- read.csv(paste0(dataFolder,"raw/ccle_lincs_mapping.csv"))
  
  tumor_cell_cor <- tissue_cell_line_cor[rownames(tissue_cell_line_cor) %in% ccle_mapping$CCLE.name, colnames(tissue_cell_line_cor) %in% dz_tissue_samples]
  
  #order based on median cor
  tumor_cell_cor_merged <- apply(tumor_cell_cor, 1, median)
  tumor_cell_cor_merged <- merge(data.frame(tumor_cell_cor_merged), ccle_mapping[, c("ccle_cell_line_name", "CCLE.name")], by.x = 0, by.y = "CCLE.name")
  names(tumor_cell_cor_merged) = c("CCLE_name", "cor", "cell_id")
  write.csv(tumor_cell_cor_merged, paste0(outputFolder, "/lincs_cell_lines_cor.csv"))
  
  tumor_cell_cor_merged <- tumor_cell_cor_merged[order(tumor_cell_cor_merged$cor, decreasing = T),]
  top_cell_lines <- tumor_cell_cor_merged$CCLE_name
  tumor_cell_cor <- tumor_cell_cor[top_cell_lines, ]
  #tumor_cell_cor$tumor_type_name = factor(tumor_cell_cor$tumor_type_name, levels = tumor_cell_cor_merged$tumor_type_name)
  
  
  pdf(paste0(outputFolder, "/lincs_cell_lines_cor.pdf"), width = 30)
  par(mar=c(12,4.1,4.1,2.1))
  boxplot(t(tumor_cell_cor), las=2, cex.axis=0.6)
  dev.off()
}

##### visualize_top_ref_tissue #####
visualize_top_ref_tissue <- function(){
  
  #tissue_ref_cor <- read.csv(paste0(outputFolder, "/GTEX_phenotype_cor.csv"))
  tissue_ref_cor <- GTEX_phenotype_cor
  reference_tissue_rank <- aggregate(cor ~ X_primary_site, GTEX_phenotype_cor, median)
  reference_tissue_rank <- reference_tissue_rank[order(reference_tissue_rank$cor, decreasing = T), ]
  
  top_refs <-  reference_tissue_rank[1:10, 1]
  
  tissue_ref_cor <- tissue_ref_cor[tissue_ref_cor$X_primary_site %in% top_refs, ]
  
  #order based on median cor
  tissue_ref_cor$ref <- factor(tissue_ref_cor$X_primary_site, levels = top_refs)
  
  #margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  pdf(paste0(outputFolder, "/top_reference_tissues.pdf"))
  par(mar=c(12,4.1,4.1,2.1))
  p <- ggplot(tissue_ref_cor, aes(ref, cor))
  print(p +   geom_boxplot(outlier.colour = "grey", notch=F, outlier.shape = NA) + geom_jitter() + theme_bw() +
          ylab("correlation") +
          xlab("") +  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
                            axis.text.y = element_text(size = 15), axis.title = element_text(size = 20), plot.margin = margin(l=35))
  )
  dev.off()
}

##### RDFHS prep #####
#prepare RDFHS matrix
prepare_rdfhs <- function(dataFolder){
  packages <- c("rhdf5")
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    print("Install required packages")
    source("https://bioconductor.org/biocLite.R")
    biocLite("rhdf5")
  }
  library("rhdf5")
  
  # load(paste0(dataFolder, "raw/treehouse/octad_readCounts.RData"))
  # load(paste0(dataFolder, "raw/treehouse/octad_tpm.RData"))
  load(paste0(dataFolder, "dz.expr.log2.readCounts.RData"))
  load(paste0(dataFolder, "dz.expr.log2.tpm.RData"))
  
  # rhdf5_file = paste0(dataFolder, "raw/treehouse/octad.h5")
  rhdf5_file = paste0(dataFolder, "octad.h5")
  
  if (file.exists(rhdf5_file)) {file.remove(rhdf5_file)}
  h5createFile(rhdf5_file)
  h5createGroup(rhdf5_file,"data/")
  h5createGroup(rhdf5_file,"meta/")
  h5createDataset(rhdf5_file, "data/count", c(nrow(dz.expr.log2.read.count), ncol(dz.expr.log2.read.count)), storage.mode = "double", chunk=c(200,200), level=7)
  h5createDataset(rhdf5_file, "data/tpm", c(nrow(dz.expr.log2.tpm), ncol(dz.expr.log2.tpm)), storage.mode = "double", chunk=c(200,200), level=7)
  
  h5write(dz.expr.log2.read.count, file=rhdf5_file, name="data/count")
  h5write(dz.expr.log2.tpm, file=rhdf5_file, name="data/tpm")
  
  h5write(rownames(dz.expr.log2.read.count), rhdf5_file ,"meta/transcripts")
  h5write(colnames(dz.expr.log2.read.count), rhdf5_file ,"meta/samples")
  
  H5close()
}


visualize_dz_sig_pathway <- function(res){
  geneEnrich(dz_signature = dz_signature)
  require(ggplot2)
  down.gene.res = read.csv(paste0(outputFolder,'enrichment_analysis/GeneEnrichment/dz_down_sig_genes_enriched.csv'),stringsAsFactors = F)
  up.gene.res = read.csv(paste0(outputFolder,'enrichment_analysis/GeneEnrichment/dz_up_sig_genes_enriched.csv'),stringsAsFactors = F)
  
  if (nrow(up.gene.res) > 0){
    up.gene.res$up_down = "up"
    up_down = up.gene.res
  }
  
  if (nrow(down.gene.res) > 0){
    down.gene.res$up_down = "down"
    up_down = down.gene.res
  }
  
  if (nrow(up.gene.res) > 0 & nrow(down.gene.res) > 0){
    up_down = rbind(up.gene.res, down.gene.res)
  }
  
  up_down$log_qval = -log(up_down$qval, 10)
  databases = unique(up_down$database)
  print (databases)
  for (database_val in databases){
    up_down_subset = subset(up_down, database == database_val)
    up_down_subset = up_down_subset[order(up_down_subset$pval), ]
    up_down_subset = up_down_subset[1:min(10, nrow(up_down_subset)), ]
    pdf(paste0(outputFolder,'enrichment_analysis/GeneEnrichment/', "dz_enriched_", database_val, ".pdf"))
    p <- ggplot(up_down_subset, aes(x=category, y = log_qval, fill= up_down)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      theme(panel.background = element_rect(fill = 'white', colour = 'white'), axis.text = element_text(colour = "black", size=10)) +
      xlab("") +
      ylab("-log(p value)") +
      scale_fill_manual(values=c("green", "red"))
    print(p)
    dev.off()
  }
}


visualize_dz_sig <- function(case_id, control_id){
  library(gplots)
  library(pheatmap)
  
  res <- read.csv(paste0(outputFolder,"dz_signature.csv"), stringsAsFactors = F)
  
  #only select top 100 genes for visualization
  res <- res[order(abs(res$log2FoldChange), decreasing = T), ]
  res <- res[1:min(200, nrow(res)), ]
  
  #print(head(res))
  #
  #load(paste0(dataFolder, 'raw/treehouse/octad_tpm.RData'))
  #expr <- dz.expr.log2.tpm[res$identifier, colnames(dz.expr.log2.tpm) %in% c(case_id, control_id)]
  
  library("rhdf5")
  rhdf5_file <- paste0(dataFolder, "octad.h5")
  transcripts = as.character(h5read(rhdf5_file, "meta/transcripts"))
  samples = as.character(h5read(rhdf5_file, "meta/samples"))
  expr = h5read(rhdf5_file, "data/tpm", index=list(which(transcripts %in% res$identifier), which(samples %in% c(case_id, control_id))))
  rownames(expr) = transcripts[which(transcripts %in% res$identifier)]
  colnames(expr) = samples[samples %in% c(case_id, control_id)]
  H5close()
  
  expr = 2^expr - 1
  #print(expr[1:3, 1:3])
  
  write.csv(expr, paste0(outputFolder, "enrichment_analysis/GeneEnrichment/dz_sig_expr_tpm.csv"))

  #for some reason, some genes' TPM is very low, they should be removed for visualization
  gene_mean = apply(expr, 1, mean)
  expr = expr[gene_mean > 1, ]

  df <- data.frame(row.names = c(case_id, control_id), type= c(rep("case", length(case_id)), rep("control", length(control_id))))
  colnames(df) <-c("type")
  mycol <- colorpanel(100,"blue","white","red")
  set.seed(0)
  # expr = t(scale(t(log(expr)))) # inconsistent but prettier. Sometimes produces NaNs, which kmeans (pheatmap) can't deal with.
  expr = t(scale(t(log2(expr +1)))) # uglier but reliable
  pheatmap(expr, cluster_rows=T, show_rownames=F,show_colnames=F,col=mycol,
           cluster_cols=F, annotation_col=df, file= paste0(outputFolder, "enrichment_analysis/GeneEnrichment/signature.pdf"))
  
}

####### deprecated functions #######
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
# visualizeLincsHits = function(rges,dz_sigUsed,drugs=''){
#   require(pheatmap)
#   require(dplyr)
#   require("gplots")
#   require("ggplot2")
#   require("RColorBrewer")
#   dz_sig = dz_sigUsed %>% select(gene,log2FoldChange)
#   rges$RGES = rges$cmap_score
#   rges_subset = rges[rges$pert_iname %in% drugs,]
#   drug_cmap_score = aggregate(RGES ~ pert_iname, rges_subset, median)
#   drug_instances_median  = merge(rges_subset, drug_cmap_score, by = c("pert_iname"))
#   drug_instances_median$diff = abs(drug_instances_median$RGES.x - drug_instances_median$RGES.y)   
#   #cmap_score.y is the median
#   drug_instances_min_diff = aggregate(diff ~ pert_iname, drug_instances_median, min)
#   drug_instances_select = merge(drug_instances_median, drug_instances_min_diff, by=c("pert_iname", "diff"))
#   drug_instances_select = drug_instances_select[!duplicated(drug_instances_select$pert_iname), ]
#   sig_id_selects = drug_instances_select$id
#   load(paste0(dataFolder,"lincs_signatures_cmpd_landmark_symbol.RData"))
#   drug_dz_signature = merge(dz_sig, data.frame(gene = rownames(lincs_signatures), 
#                                                lincs_signatures[, as.character(sig_id_selects)]),  by="gene", suffixes='')
#   gene_ids = drug_dz_signature$gene
#   drug_dz_signature_rank = drug_dz_signature[,-1]
#   for (i in 1:ncol(drug_dz_signature_rank)){
#     drug_dz_signature_rank[,i] = rank(-1 * drug_dz_signature_rank[,i] ) #highly expressed genes ranked on the top
#   }
#   gene_ids_rank <- gene_ids[order(drug_dz_signature_rank[,1])]
#   drug_dz_signature_rank <- drug_dz_signature_rank[order(drug_dz_signature_rank[,1]),] #order by disease expression
#   
#   col_sorted = sort(cor(drug_dz_signature_rank, method="spearman")["log2FoldChange",-1])    
#   drug_dz_signature_rank = drug_dz_signature_rank[,c("log2FoldChange", names(col_sorted))]
#   
#   drug_names = sapply(2:ncol(drug_dz_signature_rank), function(id){
#     rges$pert_iname[paste("X",rges$id, sep="") == names(drug_dz_signature_rank)[id]][1]
#   })
#   dz = 'case'
#   pdf(paste(outputFolder, "/lincs_reverse_expression.pdf", sep=""))
#   #colPal <- bluered(100)
#   colPal = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))
#   par(mar=c(13, 6, 2, 0.5))
#   axiscolor = sapply(c(dz, as.character(drug_names)), function(name){
#     if (name == dz){
#       "black"
#     }else if (name %in% ""){
#       "black"
#     }else{
#       "black"
#     }
#   })
#   image(t(drug_dz_signature_rank), col=colPal,   axes=F, srt=45)
#   axis(1,  at= seq(0,1,length.out=ncol( drug_dz_signature_rank ) ), labels= FALSE)
#   text(x = seq(0,1,length.out=ncol( drug_dz_signature_rank ) ), c(-0.05),
#        labels = c( dz,as.character(drug_names)),col=axiscolor, srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=0.6)
#   dev.off()
# }
# computeLINCSrges = function(dz_signature,choose_fda = F,parallel = F,maxGenes=900){
#   require(dplyr)
#   load(paste0(dataFolder,'lincs_signatures_cmpd_landmark_symbol.RData'))
#   gene.list <- toupper(rownames(lincs_signatures))
#   dz_sigUsed <- dz_signature %>% filter(toupper(Symbol) %in% gene.list)
#   dz_genes_up <- dz_sigUsed %>% filter(log2FoldChange>0) %>% arrange(desc(log2FoldChange)) %>% head(maxGenes)
#   dz_genes_down <- dz_sigUsed %>% filter(log2FoldChange<0) %>% arrange(log2FoldChange) %>% head(maxGenes)
#   write.csv(rbind(dz_genes_up,dz_genes_down), paste0(write.csv(rbind(dz_genes_up,dz_genes_down), paste0(outputFolder,'dz_sig_used.csv')),'dz_sig_used.csv'))
#   lincs_sig_info <- read.csv(paste0(dataFolder,"lincs_sig_info.csv"),
#                              stringsAsFactors = F)
#   lincs_sig_info <- lincs_sig_info %>% filter(id %in% colnames(lincs_signatures))
#   if(choose_fda == T){
#     fda_drugs = read.csv(paste0(dataFolder,"repurposing_drugs_20170327.csv"),
#                          comment.char = "!", header=T, sep="\t")
#     lincs_sig_info <- lincs_sig_info %>% filter(id %in% colnames(lincs_signatures) & 
#                                                   tolower(pert_iname) %in% tolower(fda_drugs$pert_iname))
#   }
#   lincs_sig_info <- lincs_sig_info[!duplicated(lincs_sig_info$id),]
#   sig.ids <- lincs_sig_info$id
#   if(parallel == T){
#     require(doParallel)
#     registerDoParallel(cores = (as.integer(detectCores()/4)))
#     dz_cmap_scores = foreach(exp_id = sig.ids,.combine = 'c')%dopar%{
#       cmap_exp_signature <- data.frame(gene.list,  
#                                        rank(-1 * lincs_signatures[, as.character(exp_id)], 
#                                             ties.method="random"))    
#       colnames(cmap_exp_signature) <- c("ids","rank")
#       #runs a function cmap_score_new from drugs_core_functions.R
#       cmap_score_new(dz_genes_up$Symbol,dz_genes_down$Symbol,
#                      cmap_exp_signature)
#       
#     }    
#   }else{
#     dz_cmap_scores <- NULL
#     for (exp_id in sig.ids){
#       cmap_exp_signature <- data.frame(gene.list,  
#                                        rank(-1 * lincs_signatures[, as.character(exp_id)], 
#                                             ties.method="random"))    
#       colnames(cmap_exp_signature) <- c("ids","rank")
#       #runs a function cmap_score_new from drugs_core_functions.R
#       dz_cmap_scores <- c(dz_cmap_scores, 
#                           cmap_score_new(dz_genes_up$Symbol,dz_genes_down$Symbol,
#                                          cmap_exp_signature))
#     }
#   }
#   #results <- data.frame(id = sig.ids, cmap_score = dz_cmap_scores)
#   results <- left_join(dz_cmap_scores, lincs_sig_info, by = "id")
#   results <- results[order(results$cmap_score),] 
#   return(results)
# }



# summarizeLincsRGES = function(lincs_rges,
#                               weight_cell_line = F,
#                               cell_lines = '',
#                               choose_fda = F,
#                               parallel = F){
#   lincs_drug_prediction = lincs_rges
#   lincs_drug_prediction_subset <- subset(lincs_drug_prediction,  pert_dose > 0 & pert_time %in% c(6, 24))
#   lincs_drug_prediction_pairs <- merge(lincs_drug_prediction_subset, lincs_drug_prediction_subset, by=c("pert_iname", "cell_id")) 
#   #x is the reference
#   lincs_drug_prediction_pairs <- subset(lincs_drug_prediction_pairs, id.x != id.y & pert_time.x == 24 & pert_dose.x == 10) #, select <- c("cmap_score.x", "cmap_score.y", "pert_dose.y", "pert_time.y"))
#   
#   #difference of RGES to the reference 
#   lincs_drug_prediction_pairs$cmap_diff <- lincs_drug_prediction_pairs$cmap_score.x - lincs_drug_prediction_pairs$cmap_score.y
#   lincs_drug_prediction_pairs$dose <- round(log(lincs_drug_prediction_pairs$pert_dose.y, 2), 1)
#   
#   #estimate difference
#   lincs_drug_prediction_pairs$dose_bin <- ifelse(lincs_drug_prediction_pairs$pert_dose.y < 10, "low", "high")
#   diff <- tapply(lincs_drug_prediction_pairs$cmap_diff, paste(lincs_drug_prediction_pairs$dose_bin, lincs_drug_prediction_pairs$pert_time.y), mean)
#   
#   #ignore weighting cell lines
#   if (weight_cell_line){
#     lincs_cell_line_weight <- read.csv(paste0(dataFolder, "/lincs_cell_lines_cor.csv"))
#     pred <- merge(lincs_drug_prediction, lincs_cell_line_weight, by ="cell_id")
#   }else{
#     pred <- lincs_drug_prediction
#     pred$cor <- 1
#   }
#   pred$RGES <- sapply(1:nrow(pred), function(id){getsRGES(pred$cmap_score[id], pred$cor[id], pred$pert_dose[id], pred$pert_time[id], diff, max(pred$cor))})
#   
#   cmpd_freq <- table(pred$pert_iname)
#   pred <- subset(pred, pert_iname %in% names(cmpd_freq[cmpd_freq > 0]))
#   
#   pred_merged <- pred %>% 
#     group_by(pert_iname) %>% 
#     dplyr::summarise(
#       sRGES = mean(RGES),
#       n = length(RGES),
#       medRGES = median(RGES),
#       sd = sd(RGES))
#   pred_merged <- pred_merged[order(pred_merged$sRGES), ]
#   return(pred_merged)
# }

# computeRandomLincsRGES = function(dz_signature,choose_fda = T,parallel = T,n_perm=10000){
#   
#   load(paste0(dataFolder,'lincs_signatures_cmpd_landmark_symbol.RData'))
#   gene.list <- toupper(rownames(lincs_signatures))
#   dz_sigUsed <- dz_signature %>% filter(toupper(Symbol) %in% gene.list)
#   write.csv(dz_sigUsed, paste0(outputFolder,'dz_sig_used.csv'))
#   dz_genes_up <- dz_sigUsed %>% filter(log2FoldChange>0) %>% arrange(desc(log2FoldChange))
#   dz_genes_down <- dz_sigUsed %>% filter(log2FoldChange<0) %>% arrange(log2FoldChange)
#   
#   lincs_sig_info <- read.csv(paste0(dataFolder,"lincs_sig_info.csv"),
#                              stringsAsFactors = F)
#   if(choose_fda == T){
#     fda_drugs = read.csv(paste0(dataFolder,"repurposing_drugs_20170327.csv"),
#                          comment.char = "!", header=T, sep="\t")
#     lincs_sig_info <- lincs_sig_info %>% filter(id %in% colnames(lincs_signatures) & 
#                                                   tolower(pert_iname) %in% tolower(fda_drugs$pert_iname))
#   }
#   lincs_sig_info <- lincs_sig_info[!duplicated(lincs_sig_info$id),]
#   sig.ids <- lincs_sig_info$id
#   
#   
#   if(parallel == T){
#     require(doParallel)
#     registerDoParallel(cores=4)
#     N_PERMUTATIONS <- n_perm #default 100000
#     random_sig_ids <- sample(colnames(lincs_signatures),N_PERMUTATIONS,replace=T)
#     #random_cmap_scores <- NULL
#     random_cmap_scores = foreach(exp_id = random_sig_ids,.combine = 'c')%dopar%{
#       cmap_exp_signature <- data.frame(gene.list,  rank(-1 * lincs_signatures[, as.character(exp_id)], ties.method="random"))
#       colnames(cmap_exp_signature) <- c("ids","rank")
#       
#       random_input_signature_genes <- sample(gene.list, (nrow(dz_genes_up)+nrow(dz_genes_down)))
#       rand_dz_gene_up <- data.frame(GeneID=random_input_signature_genes[1:nrow(dz_genes_up)])
#       rand_dz_gene_down <- data.frame(GeneID=random_input_signature_genes[(nrow(dz_genes_up)+1):length(random_input_signature_genes)])
#       cmap_score_new(rand_dz_gene_up,rand_dz_gene_down,cmap_exp_signature)
#     }
#   }else{
#     N_PERMUTATIONS <- n_perm #default 100000
#     random_sig_ids <- sample(1:ncol(lincs_signatures),N_PERMUTATIONS,replace=T)
#     random_cmap_scores <- NULL
#     for (exp_id in random_sig_ids){
#       #print(count)
#       cmap_exp_signature <- data.frame(gene.list,  rank(-1 * lincs_signatures[, as.character(exp_id)], ties.method="random"))    
#       colnames(cmap_exp_signature) <- c("ids","rank")
#       
#       random_input_signature_genes <- sample(gene.list, (nrow(dz_genes_up)+nrow(dz_genes_down)))
#       rand_dz_gene_up <- data.frame(GeneID=random_input_signature_genes[1:nrow(dz_genes_up)])
#       rand_dz_gene_down <- data.frame(GeneID=random_input_signature_genes[(nrow(dz_genes_up)+1):length(random_input_signature_genes)])
#       random_cmap_scores <- c(random_cmap_scores, cmap_score_new(rand_dz_gene_up,rand_dz_gene_down,cmap_exp_signature))
#     }
#   }
#   return(random_cmap_scores)
# }
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
