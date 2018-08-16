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
compEmpContGenes <- function(counts, counts_phenotype, n_topGenes = 5000){
  set <- newSeqExpressionSet(round(counts),
                             phenoData = data.frame(counts_phenotype,row.names= counts_phenotype$sample))
  design <- model.matrix(~ sample_type, data = pData(set))
  y <- DGEList(counts=counts(set), group =  counts_phenotype$sample)
  y <- calcNormFactors(y, method="TMM") #upperquartile generate Inf in the LGG case
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit) #defaults to compare tumor to normal or tumor mutant to normal
  
  top <- topTags(lrt, n=nrow(set))$table
  #n_topGenes <- n_topGenes #5000: assume there are 5000 signficant genes
  
  #based on n_topGenes computing genes with low DE
  #the genes not computed significant is in the empirical set
  i = which(!(rownames(set) %in% rownames(top)[1:min(n_topGenes,dim(top)[1])]))
  empirical <- rownames(set)[i]
  empirical
}

geneEnrich <- function(dz_signature, 
                       db.list=c( "KEGG_2016",  "GO_Biological_Process_2017", "GO_Cellular_Component_2017")){
  enrichFullGeneList <- function(up.genes, dn.genes, databases=db.list, fdr.cutoff=NULL) {
    (up.gene.res <- enrichGeneList(up.genes, databases, fdr.cutoff))
    (up.gene.res$direction <- "UP")
    (dn.gene.res <- enrichGeneList(dn.genes, databases, fdr.cutoff))
    dn.gene.res$direction <- "DN"
    (rbind(up.gene.res, dn.gene.res))
  }
  #this enrichFullGeneList is not recalled back later.
  #is there a difference between this and enrichGeneList?
  #Also, changed all words with down to dn.
  #There is a difference in results betw dn.gene.res (8hits) and down.gene.res (7hits).
  #Not only in the num of hits but also the specific hits.
  
  
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
        database.enrichments[[idx]] <- database.res[, paste("V", c(1, 2, 3, 7, 6), sep=""), with=F]
      }
    }
    
    query.results <- as.data.frame(data.table::rbindlist(database.enrichments))
    colnames(query.results) <- c("database", "category", "pval", "qval", "genes")
    
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
  
  (up.genes = dz_signature$Symbol[dz_signature$log2FoldChange > 0])
  (up.gene.res = enrichGeneList(up.genes[1:min(300, length(up.genes))], databases, fdr.cutoff))
  (up.gene.res = up.gene.res[order(up.gene.res$database, up.gene.res$p), ])
  write.csv(up.gene.res, paste0(outputFolder, "/dz_up_sig_genes_enriched", ".csv"))
  
  
  dn.genes = dz_signature$Symbol[dz_signature$log2FoldChange < 0]
  dn.gene.res = enrichGeneList(dn.genes[1:min(300, length(dn.genes))], databases, fdr.cutoff)
  dn.gene.res = dn.gene.res[order(dn.gene.res$database, dn.gene.res$p), ]
  write.csv(dn.gene.res, paste0(outputFolder, "/dz_down_sig_genes_enriched", ".csv"))
}
  

