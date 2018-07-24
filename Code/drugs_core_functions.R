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


drug_enrichment <- function(sRGES,target_type){
  require(GSVA)
  load(paste0(pipelineDataFolder,"cmpd_sets_", target_type, ".RData"))
  cmpdSets = cmpd_sets$cmpd.sets
  names(cmpdSets) = cmpd_sets$cmpd.set.names
  
  drug_pred = sRGES
  #drug_pred = read.csv(paste0("~/Documents/stanford/tumor_cell_line/RGES_manuscript/release/data/LIHC/lincs_cancer_sRGES.csv"))
  
  #FDA approved using repurposing_drugs_20170327.txt
  
  
  
  #create a random gene set
  random_times = 1000
  rgess = matrix(NA, nrow = nrow(drug_pred), ncol = random_times)
  for (i in 1:random_times){
    rgess[, i] = sample(drug_pred$sRGES, nrow(rgess))
  }
  rownames(rgess) = drug_pred$pert_iname
  rgess = cbind(rgess, drug_pred$sRGES)
  
  gsea_results = gsva(rgess, cmpdSets, method = "ssgsea",  parallel.sz=8)
  
  gsea_p = apply(gsea_results, 1, function(x){
    sum(x[1:random_times] > x[random_times+1])/random_times
  })
  
  gsea_p = data.frame(target = names(gsea_p), p = gsea_p, padj = p.adjust(gsea_p))
  gsea_p = gsea_p[order(gsea_p$p), ]
  return(gsea_p)
}
