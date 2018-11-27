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

#dev on 9-23-18
#load('~/octad_desktop/Data/dz.expr.log2.readCounts.demo.RData')

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
  counts = 2^counts - 1 #unlog the counts it was log(2x + 1)
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
  res$id = row.names(res)
  res = res %>% select(id,everything())
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



testcorTissue = function(case_id='',
                         normal_id='',
                         expSet=NULL,n_varGenes=10000,
                         outRows=1,output=T){
  require(dplyr)
    expSet_normal <- expSet[,normal_id]
    expSet_case <- as.matrix(expSet[,case_id])
    iqr_gene <-apply(expSet_normal, 1, IQR) #get the IQR per gene
    varying_genes <-order(iqr_gene, decreasing=T)[1:min(n_varGenes,length(iqr_gene))]
    
    #get the correlation matrix for each normal id and each case id
    normal_dz_cor <-cor(expSet_normal[varying_genes, ], expSet_case[varying_genes, ], method = "spearman")
    normal_dz_cor_each <-apply(normal_dz_cor, 1, median) #getting the median correlation btw each normal tissue to the case overall
    normal_dz_cor_eachDF = data.frame(cor=sort(normal_dz_cor_each, decreasing=TRUE)) %>% 
      mutate(sample.id = row.names(.)) %>% select(sample.id,cor)
    phenoDF = read.csv('~/octad_desktop/Data/phenoDF_withage.csv',stringsAsFactors = F)    
    normal_dz_cor_eachDF = normal_dz_cor_eachDF %>% left_join(phenoDF,by='sample.id')
    normal_dz_cor_biopsy = normal_dz_cor_eachDF %>% 
      group_by(sample.type,biopsy.site,cancer,data.source) %>% 
      summarise(minCor = min(cor),medCor = median(cor),maxCor = max(cor)) %>% ungroup()
    out = normal_dz_cor_biopsy %>% arrange(desc(medCor))
    out$case_id = case_id
    
    if(output==T){
      tryCatch(write.csv(normal_dz_cor,file = paste0(outputFolder,'case_normal_corMatrix.csv')),
               error = function(c) "failed to write case normal cor matrix csv. Try checking if your outputFolder string is correct or exists")
      
      tryCatch(write.csv(normal_dz_cor_eachDF,row.names = F, paste0(outputFolder, "/case_normal_median_cor.csv")),
               error = function(c) "failed to write case normal median correlation csv. Try checking if your outputFolder string is correct or exists")
      tryCatch(write.csv(normal_dz_cor_biopsy,row.names = F, paste0(outputFolder, "/case_normal_biopsy_median_cor.csv")),
               error = function(c) "failed to write case normal median correlation csv. Try checking if your outputFolder string is correct or exists")
      
    }
    return(out[1:outRows,])
  }

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
#   load(paste0(pipelineDataFolder,"cmpd_sets_", target_type, ".RData"))
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

