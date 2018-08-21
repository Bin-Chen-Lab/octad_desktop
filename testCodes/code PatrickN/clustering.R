#---------------
# package installation: do it just once

source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
install.packages("pheatmap")
install.packages("vegan")
install.packages("ape")
install.packages("rgl")

library(DESeq2)
library(pheatmap)
library(vegan)
library(ape)
library(rgl)


#----------------

# library(DESeq2)
# countData = read.table("myCounts.txt") # rows - genes, columns - samples
# conds=read.table("myConditions.txt") # rows - samples, columns - experimental conditions
# dds = DESeqDataSetFromMatrix(countData, colData=conds, design=~MyDesignFormula) # replace MyDesignFormula with something like ~factor1+factor2+factor1:factor2, using factors you have in the myConditions table  
# dds = DESeq(dds) # takes quite some time
# vst=varianceStabilizingTransformation(dds)  
# vsd=assay(vst)  # vsd is now the normalized log2-transformed data
# save(dds,vsd,conds,file="myData.RData") # save it for the future

#-----------------------
# next time, start from here

load("../melanoma/data/Skin Cutaneous Melanoma.RData")

# hierarchical clustering of samples and heatmap of sample similarities
library(pheatmap)
pheatmap(cor(log2.fpkm.matrix))  # yes it is that simple

# Principle coordiante analysis
library(vegan)
library(rgl)
library(ape)

# assembling table of conditions to lable PCoA plot:
# (in the chunk below, replace factor1 and factor2 with your actual factor names from myConditions table)
factor1=as.character(colData(dds)$factor1)
factor2=as.character(colData(dds)$factor2)
oneByTwo=paste(factor1,factor2,sep=".")
conditions=data.frame(cbind(factor1,factor2,oneByTwo))

# actual PCoA analysis
dds.pcoa=pcoa(vegdist(t(vsd),method="manhattan")/1000)
scores=dds.pcoa$vectors

# plotting
plot(scores[,1], scores[,2],col=as.numeric(as.factor(factor1)))
ordispider(scores,factor2,label=T)
ordiellipse(scores,factor2)

# interactive 3d plot - can rotate it by dragging mouse
radiusScale=2 # change it if the spheres come out too small or too big in the next one
plot3d(scores[,1], scores[,2],scores[,3],col=as.numeric(as.factor(factor1)),type="s",radius=radiusScale*as.numeric(as.factor(factor2)))

# formal permutation-based analysis of variance 
adonis(t(vsd)~factor1*factor2,data=conditions,method="manhattan")  




