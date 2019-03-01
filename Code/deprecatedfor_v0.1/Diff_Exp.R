####ToDo####
#script to output before and after normalization
#ask Ke about the warnings()

####dataframes used####
#dz_expr : counts matrix for all the tissues

####Input Test####
# normalize_samples = F
# case_id <- colnames(dz_expr)[2:10]
# control_id <- colnames(dz_expr)[11:20]
# DE_method = 'edgeR'


####Input Variables Required####
#case_id : samples of case tissues in a character vector
#control_id : samples of control tissues in a character vector
#dz_expr : must have dataframe that contains both case_id and control_id in its
  #column names, #store it as "dz_expr" in order to run this script
#outputFolder : folder to place output files


####Input Parameters####
#must run the options in set_parameters.R

####Intermediate Variables####
  #dz_tissue : expr matrix of case id
  #ref_tissue : expr matrix of control id
#counts : expr matrix combined of dz_tissue and ref_tissue
  #rows genes
  #columns samples
#counts phenotype :
  # dataframe to annotate case vs. control 

####output####
#res either from limma or edgeR but either will contain:
  #1. rowname of the dz_expr will be the gene identifier this may differ from datasets
  #2. log2Foldchange
  #3. padj

####helper functions####

#computes the low expression and output expressions to keep
library(RUVSeq)
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

# dz_tissue <- dz_expr[,case_id]
# ref_tissue <- dz_expr[,control_id]

#unlog the counts it was log(2x + 1)
# dz_tissue <- 2^dz_tissue - 1
# ref_tissue <- 2^ref_tissue -1

counts <- cbind(dz_tissue, ref_tissue)

counts_phenotype <- rbind(data.frame(sample = case_id,sample_type = 'case'),
                          data.frame(sample = control_id, sample_type = 'control'))
#need to set the levels so that the design matrix won't be case as an intercept
counts_phenotype$sample_type <- factor(counts_phenotype$sample_type,levels=c('control','case'))

library(dplyr)


#counts dataframe === remove low counts ===> set === normalized ===> set1
#if no empirical genes found it will just create a matrix without running RUVg
if(normalize_samples == T){
  highExpGenes <- remLowExpr(counts,counts_phenotype)
  write.csv(data.frame(highExpGenes),file=paste0(outputFolder,"highExpGenes.csv"))
  counts <- counts[highExpGenes,]
  empiricalGenes <- compEmpContGenes(counts,counts_phenotype,n_topGenes = n_topGenes)
  stopifnot(length(empiricalGenes)>0)
  write.csv(data.frame(empiricalGenes),file=paste0(outputFolder,"computedEmpGenes.csv"))
  set <- newSeqExpressionSet(counts = round(counts),
                             phenoData = data.frame(counts_phenotype,row.names= counts_phenotype$sample))
  set1 <- RUVg(set,empiricalGenes,k=k)
  }else{set1 <- newSeqExpressionSet(counts = round(counts),
                                 phenoData = data.frame(counts_phenotype,row.names= counts_phenotype$sample))}

rm(counts,set) #no longer needs these dataframes

if (DE_method == "edgeR"){
  library(edgeR)
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
  res <- lrt$table
  colnames(res) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
  res$padj <- p.adjust(res$pvalue)
  write.csv(res,paste0(outputFolder,'DE_genes.csv')) 
} else if (DE_method == "limma"){
  print('computing DE via limma')
  library('Glimma')
  library(RColorBrewer)
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
  dgList <- calcNormFactors(dgList, method="TMM")
  x <- dgList
  nsamples <- ncol(x)
  col <- brewer.pal(nsamples, "Paired")
  lcpm <- cpm(x, log=TRUE)
  #plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
  #     main="", xlab="")
  #title(main="B. Filtered data", xlab="Log-cpm")
  #abline(v=0, lty=3)
  # for (i in 2:nsamples){
  #   den <- density(lcpm[,i])
  #   lines(den$x, den$y, col=col[i], lwd=2)
  # }
  
  group <-counts_phenotype$sample_type
  col.group <- group
  levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
  col.group <- as.character(col.group)
  #plotMDS(lcpm, labels=coldata$sample_type, col=col.group)
  
  # ----design
  #design <- model.matrix(~0 + group)
  colnames(design) <- gsub("group", "", colnames(design))
  
  contr.matrix <- makeContrasts(
    TumorvsNon = sample_typecase - (Intercept), 
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
  tumorvsnormal.topgenes <- rownames(tumorvsnormal[1:50,])
  'mycol <- colorpanel(1000,"blue","white","red")
  pdf( paste0(outputFolder, "/limma_sig.pdf"))
  heatmap.2(v$E[tumorvsnormal.topgenes,], scale="row",
  labRow=tumorvsnormal.topgenes, labCol=group, 
  col=mycol, trace="none", density.info="none", 
  margin=c(8,6), lhei=c(2,10), dendrogram="column")
  dev.off()'
  
  'df <- as.data.frame(coldata[,c("sample_type")])
  rownames(df) = (coldata$sample)
  colnames(df) = c("type")
  pheatmap(v$E[tumorvsnormal.topgenes,], cluster_rows=T, show_rownames=T,show_colnames=F, scale="row",col=mycol,
  cluster_cols=T, annotation_col=df, file= paste0(outputFolder, "/limma_sig.pdf"))
  '
  
  res <-tumorvsnormal
  colnames(res) <-c("log2FoldChange", "AveExpr", "t", "pvalue", "padj")
  write.csv(res,paste0(outputFolder,'DE_genes.csv'))
}

res$identifier <- row.names(res) #we don't know what row names and pheno they might have

