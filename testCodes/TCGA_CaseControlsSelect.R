library(dplyr)
library(tidyr)

####filtering for case controls w/i liver####
PCA64 <- read.csv("FileZilla Downloads/Cancer Prediction/CancerPrediction_DAE_Data/Code - PCA/PCA Run 1 - Output/PCA64.csv",stringsAsFactors = F)
phenoDF = read.csv('FileZilla Downloads/Cancer Prediction/CancerPrediction_DAE_Data/integrated.OCTAD.basic.sample.meta.csv',stringsAsFactors = F)
phenoDF.liver = phenoDF %>% filter(cancer == 'Liver Hepatocellular Carcinoma',data.source=='TCGA')

ptid = t(data.frame(phenoDF.liver$sample.id %>% strsplit('-',useBytes = T)))
ptDF = data.frame(cbind(phenoDF.liver$sample.id,ptid))
ptDF$ptid = paste(ptDF$X2,ptDF$X3,ptDF$X4,sep = "-")
ptDF = ptDF %>% select(ptid,sample.id=X1,sample.type = X5)
ptDF.case_control = ptDF %>% 
  group_by(ptid) %>% spread(sample.type,sample.id) %>% ungroup() %>% 
  select(ptid,adjacent=`11`,primary=`01`) %>% filter(!is.na(adjacent))

