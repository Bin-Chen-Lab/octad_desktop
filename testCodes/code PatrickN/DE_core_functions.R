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
  

