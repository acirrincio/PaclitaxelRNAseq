library(openxlsx)
library(DESeq2)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
library(reshape2)
library(ggplot2)
library(viridis)
library(plyr)
library(dplyr)
library(DESeq2)
library(gplots)
library("pheatmap")
library("RColorBrewer")
library(AnnotationDbi)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(vidger)
library(tidyr)
library(apeglm)
library(devtools)
library(ComplexHeatmap)
library(dendextend)
library(openxlsx)
library(LSD)
library(biomaRt)
library(PoiClaClu)
library(ggpubr)
library(tidyverse)
library(rstatix)
library(ggpubr)



setwd("/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/DRG/Count_matrix_original_and_updated/")
#reading in count matrix
#countMatrixXXXXX<- readRDS("countMatrix.rds") ####this was for old count matrix 
#head(countMatrixXXXXX)
file = ("~/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/DRG/Count_matrix_original_and_updated/DRG_CountMatrix_lanesCollated_woPairedReads.txt")
countMatrix <- read.table(file)
View(countMatrix)
#read in sample info
setwd("/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/DRG/")


sampleinfo <- read.csv("sampleLookup.txt", sep = "\t", header = TRUE)
View(sampleinfo)
sampleinfo <- sampleinfo[65:128,] #getting rid of all the skin info
View(sampleinfo)
samplematch <- read.xlsx("sampleMatchup.xlsx",colNames = FALSE)
View(samplematch)



#trim off numbers to make columns match
colnames <- colnames(countMatrix)
View(colnames)
sampleinfo$col <- gsub('201926', 'X201926', sampleinfo$HIHG_ID)
sampleinfo$col2 <- gsub('-01', '', sampleinfo$col)
sampleinfo$id <- gsub('[0-9]+-','', sampleinfo$External_ID)
sampleinfo$id2 <- gsub('R-','', sampleinfo$id)
sampleinfo$id3 <- sampleinfo$id2
View(sampleinfo)

#fixing sample match

samplematch$id <- gsub('[0-9]+ - ','',samplematch$X1)
samplematch$id2 <- gsub('[(]','',samplematch$id)
samplematch$id3 <- gsub('[)]','',samplematch$id2)
samplematch$id3[34] <- c("DRG08")
#merging columns 
samplesMatched <- merge(samplematch,sampleinfo, by = 'id3' )
View(samplesMatched)

#subset because we only need two columns 
samplesMatched <- as.data.frame(samplesMatched)
samplesMatched <- samplesMatched %>% dplyr::select(col2,id3, X2,X3)
View(samplesMatched)
counts_origingal_columns <- countMatrix
View(counts_origingal_columns[1:100,])
View(countMatrix)

#going to change column names 
#samples matched is now the Col_data
coldata = samplesMatched
View(coldata)
colnames(coldata) = c("column","sample name","treatment","timepoint")
rownames(colData) = coldata$column


#Need to replace column '#' at top with DRG#
View(countMatrix)
countMatrix <- as.data.frame(countMatrix)
names(countMatrix)
names(countMatrix) <- coldata$`sample name`[match(names(countMatrix), coldata$column)]
names(countMatrix)
head(countMatrix)


#order of columns and order on the colData need to be the same 
#sort by  'column' column in coldata
colData <- arrange(colData,`sample name`)
View(colData)
coldata <- colData[,c("treatment","timepoint")]

#make new coldata without the old row#, keeping DRG#
colData$column = NULL
rownames(colData) = colData$`sample name`
View(colData)

#changing timepoint to number only....i.e. 'd.4' format to '4 

colData$timepoint <-gsub('d.','',colData$timepoint)

View(colData)

#Arrange columns in numerical order in countMatrix
colnames(countMatrix)

#countMatrix2 has the DRG names in order 

countMatrix2 <- countMatrix[,order(colnames(countMatrix))]
View(head(countMatrix2))
View(countMatrix2)

#we now have the count matrix with the samples ordered by DRG#
#now we can make sub-matriciesfor each day and colData information for our DEseq objects 
View(colData)
View(countMatrix2)


#need to remove the rows for # of mapped reads, multimapping, etc
row.names.remove <- c("N_multimapping", "N_noFeature","N_unmapped","N_ambiguous")
countMatrix2 <- countMatrix2[!(row.names(countMatrix2) %in% row.names.remove), ]
View(countMatrix2)


#converting the rownames to gene names 
###########################
###########################
###########################
###########################

#making a column name from the row names 

countMatrix2$ensembl_gene_id <- rownames(countMatrix2)


#now we have our huge dataFrame 
#first i will convert the ensembl names to gene symbols before proceeding 


##
##changing the gene names to gene symbols
##(here are the setps taken to get to the code below)
#ensembl <- useEnsembl(biomart = "genes")
#datasets <- listDatasets(ensembl)
#View(datasets)

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

#finding filters
#filters = listFilters(ensembl)
#View(filters)
#mgi_symbol

a <- getBM(attributes = c('mgi_symbol','ensembl_gene_id'),
           filters = 'ensembl_gene_id',
           values = countMatrix2$ensembl, 
           mart = ensembl)

zz <- join(a,countMatrix2,by ='ensembl_gene_id' )

View(zz)
#zz is the new RAW count matrix
#need to clearn up before saving 
#saving the count matrix to the general directory 

setwd('/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/DRG/')

#now doing cleanup for DEseq2
#getting rid of ensemblID column 
zz$ensembl_gene_id <- NULL
#making the rowname the gene symbol - dealing with unique names 
rownames(zz) <- make.names(zz$mgi_symbol,unique = TRUE)
#getting rid of mgi_symbol column
zz$mgi_symbol <- NULL



getwd()

#saving RAW matrix before proceeding with DEseq analysis and normalization

write.xlsx(zz,"RawCountMatrix_DRGs.xlsx", rowNames = TRUE)


#saving coldata sheet to make everything more straightforward in the future 

write.xlsx(colData, "drg_colData.xlsx")

#making colData for each day 
View(colData)
colData4 <- colData[1:16,]
colData7 <- colData[17:32,]
colData11 <- colData[33:48,]
colData23 <- colData[49:64,]

View(colData4)
View(colData7)
View(colData11)
View(colData23)

#making matrix for each day 
#all day
counts <- zz
counts4 <- counts[,1:16]
counts7 <- counts[,17:32]
counts11 <- counts[,33:48]
counts23 <- counts[,49:64]


#making DDS objects 
#going to filter out rows with <10 genes in the original dds only, then filter these out in the subsequent matrices

#all days
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData,design = ~ treatment)
dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#4 days
dds4 <- DESeqDataSetFromMatrix(countData = counts4, colData = colData4,design = ~ treatment)
dds4 <- DESeq(dds4)
dds4 <- dds4[keep,]
#7 days
dds7 <- DESeqDataSetFromMatrix(countData = counts7, colData = colData7,design = ~ treatment)
dds7 <- DESeq(dds7)
dds7 <- dds7[keep,]
#11 days
dds11 <- DESeqDataSetFromMatrix(countData = counts11, colData = colData11,design = ~ treatment)
dds11 <- DESeq(dds11)
dds11 <- dds11[keep,]
#23 days
dds23 <- DESeqDataSetFromMatrix(countData = counts23, colData = colData23,design = ~ treatment)
dds23 <- DESeq(dds23)
dds23 <- dds23[keep,]


#normalizing the datasets with rlog 

rld <- rlog(dds)
rld4 <- rlog(dds4)
rld7 <- rlog(dds7)
rld11 <- rlog(dds11)
rld23 <- rlog(dds23)

#normalizing the datasets with vst - this should be better with 
#sample sizes greater than 30 

vsd <- vst(dds, blind = TRUE)
vsd4 <- vst(dds4, blind = TRUE)
vsd7 <- vst(dds7, blind = TRUE)
vsd11 <- vst(dds11, blind = TRUE)
vsd23 <- vst(dds23, blind = TRUE)



##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
#######start here after loading environment RData file. 
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################


#checking number of rows in the transformed data just for kicks and giggles
nrow(rld)
nrow(rld4)
nrow(rld11)
nrow(rld23)

nrow(vsd)
nrow(vsd4)
nrow(vsd7)
nrow(vsd11)
nrow(vsd23)

##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
###### PCA [all treatments]
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################

#with rld transformation
plotPCA(rld,intgroup = c("timepoint","treatment"))
plotPCA(rld4,intgroup = c("treatment"),ntop = 42180)
plotPCA(rld7,intgroup = c("treatment"),ntop = 42180)
plotPCA(rld11,intgroup = c("treatment"),ntop = 42180)
plotPCA(rld23,intgroup = c("treatment"),ntop = 42180)

#with VSD transformation
plotPCA(vsd,intgroup = c("timepoint","treatment"))
plotPCA(vsd4,intgroup = c("treatment"),ntop = 25841) #+
        #geom_text(aes(label=name),vjust=2)
plotPCA(vsd7,intgroup = c("treatment"),ntop = 25841) #+
        #geom_text(aes(label=name),vjust=2)
plotPCA(vsd11,intgroup = c("treatment"),ntop = 25841) #+
        #geom_text(aes(label=name),vjust=2)
plotPCA(vsd23,intgroup = c("treatment"),ntop = 25841) #+
        #geom_text(aes(label=name),vjust=2)

dev.off()




##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
###### PCA (Paclitaxel vs Vehicle Only)
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################

#subsetting the vsds

vsd.sub <- vsd[,vsd$treatment %in% c("P&S","V&S")]
vsd4.sub <- vsd4[,vsd4$treatment %in% c("P&S","V&S")]
vsd7.sub <- vsd7[,vsd7$treatment %in% c("P&S","V&S")]
vsd11.sub <- vsd11[,vsd11$treatment %in% c("P&S","V&S")]
vsd23.sub <- vsd23[,vsd23$treatment %in% c("P&S","V&S")]


#with VSD transformation
plotPCA(vsd.sub,intgroup = c("timepoint","treatment"))
plotPCA(vsd4.sub,intgroup = c("treatment"),ntop = 25841) #+
        #geom_text(aes(label=name),vjust=2)
plotPCA(vsd7.sub,intgroup = c("treatment"),ntop = 25841) #+
        #geom_text(aes(label=name),vjust=5)
plotPCA(vsd11.sub,intgroup = c("treatment"),ntop = 25841) #+
        #geom_text(aes(label=name),vjust=2)
plotPCA(vsd23.sub,intgroup = c("treatment"),ntop = 25841) #+
        #geom_text(aes(label=name),vjust=2)

dev.off()




#######
#######
#######
#######
#######
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
#heatmap of sample-sample distances 
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
#all days 
#normalized
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(dds$treatment, dds$timepoint, sep="-")
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9,"Spectral")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         fontsize_row = 5,
         fontsize_col = 5,
         col=colors)
dev.off()


#also trying vsd transformation


sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(dds$treatment, dds$timepoint, sep="-")
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9,"Spectral")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         fontsize_row = 5,
         fontsize_col = 5,
         )
dev.off()


##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
############################################################################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
############################################################################################
##############################################
##############################################
##############################################
#sample distances (Paclitaxel vs. Vehicle ONLY)
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
View(colData)

#normal - vst transformed 

#subsetting the vsd object 

vsd.sub <- vsd[,vsd$treatment %in% c("P&S","V&S")]

sampleDists <- dist(t(assay(vsd.sub)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd.sub$treatment, vsd.sub$timepoint, sep="-")
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9,"Spectral")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         fontsize_row = 8,
         fontsize_col = 8)
dev.off()



##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################


#dispersion estimates
#all
plotDispEsts(dds)
dev.off()
#day4
plotDispEsts(dds4)
dev.off()
#day7
plotDispEsts(dds7)
dev.off()
#day11
plotDispEsts(dds11)
dev.off()
#day 23
plotDispEsts(dds23)
dev.off()


##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################


##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################


##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
#top 500 variably expressed genes
##############################################
##############################################
##############################################
##############################################




library("genefilter")
topVarGenes <- head(order(-rowVars(assay(vsd))),100)
mat <- assay(vsd)[ topVarGenes, ]
mat<- mat - rowMeans(mat)
colnames(mat) <- paste(dds$treatment, dds$timepoint, sep="-")
df <- as.data.frame(colData(dds)[,c("treatment","timepoint")])
pheatmap(mat, cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="Spectral")))(1000),
         cellwidth = 5,
         cellheight = 2.9,
         fontsize_row = 3,
         fontsize_col = 5, 
         
)
dev.off()

View(mat)

##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
#top 100 variably expressed genes (paclitaxel vs. vehicle)
##############################################
##############################################
##############################################
##############################################

dds.sub <- dds[ , dds$treatment %in% c("P&S", "V&S") ]




library("genefilter")
topVarGenes <- head(order(-rowVars(assay(vsd.sub))),100)
mat <- assay(vsd.sub)[ topVarGenes, ]
mat<- mat - rowMeans(mat)
colnames(mat) <- paste(dds.sub$treatment, dds.sub$timepoint, sep="-")
df <- as.data.frame(colData(dds.sub)[,c("treatment","timepoint")])
pheatmap(mat, cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="Spectral")))(1000),
         cellwidth = 5,
         cellheight = 2.9,
         fontsize_row = 3,
         fontsize_col = 5, 
         
)
dev.off()

View(mat)

##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################





##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
############  heatmap of the count matrix   ##################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("treatment","timepoint")])
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE)




##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################

#now volcano plot
#if significant (p-value<0.01), will have to have >1 fold change to be highlighted
#looking at log-fold change and p-values
#JUST NEED TO JUST CHANGE THE DDS 
#ALL COMBINATIONS OF TREATMENTS ARE ALREADY HERE
#
#PC vs VS
res <- results(dds4,contrast = c("treatment", "P&S", "V&S"))
res1 <- as.data.frame(res)
res1= mutate(res1, sig=ifelse(res1$pvalue<0.01, "FDR<0.01","Not Sig")) #so we can look at p-value
res1[which(abs(res1$log2FoldChange)<1),"sig"] = "Not Sig"
row.names(res1) = row.names((res))

#we will see on a better scale how sig. genes are 
#let's label genes on the volcano plot 
pdf("23days_volcano_PC_v_VS_pvalueThresh_01_lfcThresh_1.pdf")
EnhancedVolcano(res1,
                lab = rownames(res1),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.01,
                FCcutoff = 1.0,
                pointSize = 4.0,
                labSize = 4,
                legendLabSize = 12,
                legendIconSize = 4.0)
dev.off()
#write a CSV with the results 
write.csv(res1,"23days_volcano_PC_v_VS_pvalueThresh_01_lfcThresh_1_1.csv")
#



##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
#####need to make treatment contrast tables for STEM program that include the Log2FoldChanges 
#####looking at days 4,7,11,23
#####also making heatmaps 
#####these will be filtered by LFC and pValue cutoff 
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
#
#
#all possible pairwise combinations below 

#ps vs vs*
#ps vs vc
#pc vs ps*
#pc vs vc
#pc vs vs*
#vc vs vs

setwd("/Users/anthony/Desktop/RNAseq/BenSandra/Data&Figures_03242021/DRG/")

#we are only interested in these 3 for now:
#ps vs vs*
#pc vs ps*
#pc vs vs*

#first need to get the fold change lists from comparing treatments 
res4 <- results(dds4,contrast = c("treatment", "P&C", "V&S"))
res7 <- results(dds7,contrast = c("treatment", "P&C", "V&S"))
res11 <- results(dds11,contrast = c("treatment", "P&C", "V&S"))
res23 <- results(dds23,contrast = c("treatment", "P&C", "V&S"))
res4 <- as.data.frame(res4)
res7 <- as.data.frame(res7)
res11 <- as.data.frame(res11)
res23 <- as.data.frame(res23)

#remove all genes with p-values NA
res4 <- res4[!is.na(res4$pvalue), ] # Removes rows with pvalue = NA
res7 <- res7[!is.na(res7$pvalue), ] # Removes rows with pvalue = NA
res11 <- res11[!is.na(res11$pvalue), ] # Removes rows with pvalue = NA
res23 <- res23[!is.na(res23$pvalue), ] # Removes rows with pvalue = NA

# Removes rows with log2FoldChange = NA
res4 <- res4[!is.na(res4$log2FoldChange), ] 
res7 <- res7[!is.na(res7$log2FoldChange), ] 
res11 <- res11[!is.na(res11$log2FoldChange), ] 
res23 <- res23[!is.na(res23$log2FoldChange), ] 

#Keeps ONLY rows with p-values < 0.01
res4 <- res4[res4$pvalue < 0.01, ] 
res7 <- res7[res7$pvalue < 0.01, ] 
res11 <- res11[res11$pvalue < 0.01, ] 
res23 <- res23[res23$pvalue < 0.01, ] 


# Keeps rows with log2FoldChange  > 0.6 or < -0.6
res4.above <- res4[res4$log2FoldChange > 0.6, ] #day 4
res4.below <- res4[res4$log2FoldChange < -0.6, ]
res4.out <- rbind(res4.above, res4.below)

res7.above <- res7[res7$log2FoldChange > 0.6, ] #day 7
res7.below <- res7[res7$log2FoldChange < -0.6, ] 
res7.out <- rbind(res7.above, res7.below)

res11.above <- res11[res11$log2FoldChange > 0.6, ] #day 11
res11.below <- res11[res11$log2FoldChange < -0.6, ]
res11.out <- rbind(res11.above, res11.below)

res23.above <- res23[res23$log2FoldChange > 0.6, ] #day 23
res23.below <- res23[res23$log2FoldChange < -0.6, ]
res23.out <- rbind(res23.above, res23.below)

#make a cumulative list of genes that are significant from these lists 
list <- c(row.names(res4.out),row.names(res7.out),row.names(res11.out),row.names(res23.out))
length(list)
#remove duplicates from list
list <- unique(unlist(strsplit(list, " ")))
length(list)
#now we have a list of gene names that were significant from all of the time-points

#now will go back and get the log2fold changes for all of the genes in our list
res4 <- results(dds4,contrast = c("treatment", "P&C", "V&S"))
res7 <- results(dds7,contrast = c("treatment", "P&C", "V&S"))
res11 <- results(dds11,contrast = c("treatment", "P&C", "V&S"))
res23 <- results(dds23,contrast = c("treatment", "P&C", "V&S"))

res4 <- as.data.frame(res4)
res4 <- res4[list,]
res4 <- res4[,"log2FoldChange",drop = FALSE]
colnames(res4) = c("day4")

res7 <- as.data.frame(res7)
res7 <- res7[list,]
res7 <- res7[,"log2FoldChange",drop = FALSE]
colnames(res7) = c("day7")

res11 <- as.data.frame(res11)
res11 <- res11[list,]
res11 <- res11[,"log2FoldChange",drop = FALSE]
colnames(res11) = c("day11")

res23 <- as.data.frame(res23)
res23 <- res23[list,]
res23 <- res23[,"log2FoldChange",drop = FALSE]
colnames(res23) = c("day23")

#make final matrix
final <- cbind(res4,res7,res11,res23)
#remove rows with names containing NA
final <- final[grepl("^NA", rownames(final))==F,]
final <- as.data.frame(final)
nrow(final)
View(final)

write.xlsx(final,"Paclitaxel+CL vs Vehicle_Log2FC0-6_pvalue_01.xlsx", rowNames = TRUE)


#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#make heatmap of foldchanges (the filtered ones by LFC and p-value cutoff)
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################

pheatmap(final, cluster_rows= TRUE, show_rownames=TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="Spectral")))(1000),
         cluster_cols=FALSE,
         cellwidth = 100, 
         cellheight = 0.5,
         fontsize = 8,
         treeheight_row = 70,
         fontsize_row = 4.5,
         fontsize_col = 15,
         main = "Paclitaxel+CL vs Vehicle", 
         angle_col = "0",
         cutree_rows = 15,
         scale = 'row'
         
)

dev.off()

nrow(final)



##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
#going to make bar graphs of number of upregulated and downregulated genes 
#per day per contrast 
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################

#VC VS -
#PS VS -****
#PC VS -****
#PS PC -****
#PS VC -
#PC VC -

#   **** only interested in these contrasts for now 

#first need to get the fold change lists from comparing treatments 
res4 <- results(dds4,contrast = c("treatment", "P&C", "V&S"))
res7 <- results(dds7,contrast = c("treatment", "P&C", "V&S"))
res11 <- results(dds11,contrast = c("treatment", "P&C", "V&S"))
res23 <- results(dds23,contrast = c("treatment", "P&C", "V&S"))
res4 <- as.data.frame(res4)
res7 <- as.data.frame(res7)
res11 <- as.data.frame(res11)
res23 <- as.data.frame(res23)
#remove all genes with p-values NA
res4 <- res4[!is.na(res4$pvalue), ] # Removes rows with pvalue = NA
res7 <- res7[!is.na(res7$pvalue), ] # Removes rows with pvalue = NA
res11 <- res11[!is.na(res11$pvalue), ] # Removes rows with pvalue = NA
res23 <- res23[!is.na(res23$pvalue), ] # Removes rows with pvalue = NA
# Removes rows with log2FoldChange = NA
res4 <- res4[!is.na(res4$log2FoldChange), ] 
res7 <- res7[!is.na(res7$log2FoldChange), ] 
res11 <- res11[!is.na(res11$log2FoldChange), ] 
res23 <- res23[!is.na(res23$log2FoldChange), ] 
#Keeps ONLY rows with p-values < 0.01
res4 <- res4[res4$pvalue < 0.01, ] 
res7 <- res7[res7$pvalue < 0.01, ] 
res11 <- res11[res11$pvalue < 0.01, ] 
res23 <- res23[res23$pvalue < 0.01, ] 
# Keeps rows with log2FoldChange  > 0.6 or < -0.6
#separated into above and below 
res4.above <- res4[res4$log2FoldChange > 0.6, ] #day 4
res4.below <- res4[res4$log2FoldChange < -0.6, ]
res7.above <- res7[res7$log2FoldChange > 0.6, ] #day 7
res7.below <- res7[res7$log2FoldChange < -0.6, ] 
res11.above <- res11[res11$log2FoldChange > 0.6, ] #day 11
res11.below <- res11[res11$log2FoldChange < -0.6, ]
res23.above <- res23[res23$log2FoldChange > 0.6, ] #day 23
res23.below <- res23[res23$log2FoldChange < -0.6, ]
##going to make a count matrix of above and below genes 
Zz <- c(nrow(res4.above), #day 4
        nrow(res4.below),
        nrow(res7.above), #day 7
        nrow(res7.below),
        nrow(res11.above), #day 11
        nrow(res11.below),
        nrow(res23.above), #day 23
        nrow(res23.below))

Zz <- as.data.frame(Zz)

rownames(Zz) <- c("res4.above",
                  "res4.below",
                  "res7.above",
                  "res7.below",
                  "res11.above",
                  "res11.below",
                  "res23.above",
                  "res23.below")


View(Zz)


##############################################
##############################################
##############################################
##############################################
##############################################
# Plotting individual gene normalized counts
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################


setwd("/Users/anthony/Desktop/RNAseq/BenSandra/Data&Figures_03242021/DRG/")


AAA <- counts(dds, normalized = TRUE) #making normalized count matrix



View(AAA)


#P&S
#P&C
#V&S
#V&C


View(colData)

############################################matrix for V&C
VC4 <- AAA[,1:4] ################################4day column
#VC4 <- rowMeans(VC4)
VC4 <- as.data.frame(VC4)
colnames(VC4) <- c("4day","4day","4day","4day")
#view(VC4)
VC7 <- AAA[,21:24] ##############################7day column
#VC7 <- rowMeans(VC7)
VC7 <- as.data.frame(VC7)
colnames(VC7) <- c("7day","7day","7day","7day")
#view(VC7)
VC11 <- AAA[,37:40] ##############################11day column
#VC11 <- rowMeans(VC11)
VC11 <- as.data.frame(VC11)
colnames(VC11) <- c("11day","11day","11day","11day")
#view(VC11)
VC23 <- AAA[,53:56] ##############################23day column
#VC23 <- rowMeans(VC23)
VC23 <- as.data.frame(VC23)
colnames(VC23) <- c("23day","23day","23day","23day")
#view(VC23)
VC <- cbind(VC4,VC7,VC11,VC23) #here is our VC matrix
View(VC)

nrow(VC)
#write to .xlsx file 
write.xlsx(VC, 'Vehicle+CL.xlsx', rowNames= TRUE)

#keep this code for making tab delimited text files

############################################matrix for V&S
VS4 <- AAA[,5:8] ################################4day columns
#VS4 <- rowMeans(VS4)
VS4 <- as.data.frame(VS4)
colnames(VS4) <- c("4day","4day","4day","4day")
VS7 <- AAA[,17:20] ##############################7day columns
#VS7 <- rowMeans(VS7)
VS7 <- as.data.frame(VS7)
colnames(VS7) <- c("7day","7day","7day","7day")
VS11 <- AAA[,33:36] ##############################11day columns
#VS11 <- rowMeans(VS11)
VS11 <- as.data.frame(VS11)
colnames(VS11) <- c("11day","11day","11day","11day")
VS23 <- AAA[,49:52] ##############################23day columns
#VS23 <- rowMeans(VS23)
VS23 <- as.data.frame(VS23)
colnames(VS23) <- c("23day","23day","23day","23day")
VS <- cbind(VS4,VS7,VS11,VS23) #here is our VC matrix
View(VS)



#write to .xlsx file 
write.xlsx(VS, 'Vehicle+Saline.xlsx', rowNames= TRUE)


############################################matrix for P&C
PC4 <- AAA[,13:16] ################################4day column
#PC4 <- rowMeans(PC4)
PC4 <- as.data.frame(PC4)
colnames(PC4) <- c("4day","4day","4day","4day")
#view(PC4)
PC7 <- AAA[,29:32] ##############################7day column
#PC7 <- rowMeans(PC7)
PC7 <- as.data.frame(PC7)
colnames(PC7) <- c("7day","7day","7day","7day")
#View(VC7)
PC11 <- AAA[,45:48] ##############################11day column
#PC11 <- rowMeans(PC11)
PC11 <- as.data.frame(PC11)
colnames(PC11) <- c("11day","11day","11day","11day")
#view(VC11)
PC23 <- AAA[,61:64] ##############################23day column
#PC23 <- rowMeans(PC23)
PC23 <- as.data.frame(PC23)
colnames(PC23) <- c("23day","23day","23day","23day")
#View(VC23)
PC <- cbind(PC4,PC7,PC11,PC23) #here is our PC matrix

#check number of rows
nrow(PC)

#write to .xlsx file 
write.xlsx(PC,"Paclitaxel+CL.xlsx", rowNames = TRUE)

############################################matrix for P&S
PS4 <- AAA[,9:12] ################################4day column
#PS4 <- rowMeans(PS4)
PS4 <- as.data.frame(PS4)
colnames(PS4) <- c("4day","4day","4day","4day")
#view(PS4)
PS7 <- AAA[,25:28] ##############################7day column
#PS7 <- rowMeans(PS7)
PS7 <- as.data.frame(PS7)
colnames(PS7) <- c("7day","7day","7day","7day")
#view(PS7)
PS11 <- AAA[,41:44] ##############################11day column
#PS11 <- rowMeans(PS11)
PS11 <- as.data.frame(PS11)
colnames(PS11) <- c("11day","11day","11day","11day")
#view(PS11)
PS23 <- AAA[,57:60] ##############################23day column
#PS23 <- rowMeans(PS23)
PS23 <- as.data.frame(PS23)
colnames(PS23) <- c("23day","23day","23day","23day")
#view(PS23)
PS <- cbind(PS4,PS7,PS11,PS23) #here is our PC matrix

#check number of rows
nrow(PS)
#write to .xlsx file 
write.xlsx(PS,"Paclitaxel+Saline.xlsx", rowNames = TRUE)


View(PS)



##############################################
##############################################
##############################################
##############################################
##############################################
# Plotting individual gene normalized counts
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################

#
#
#
#
#
#
#

#
#plot an individual gene
#
#
#
#
#
#

goi <- c("Ankrd27")
tcounts <- t(counts(dds[goi, ], normalized=TRUE, replaced=FALSE)) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
#arranging variables 
tcounts$Timepoint <- factor(tcounts$timepoint, levels=c("4", "7", "11", "23"))
#all treatments
#tcounts$Treatment <- factor(tcounts$Treatment, levels=c("V&S", "P&S", "P&C"))
#PCTX and VEH only 
tcounts$treatment <- factor(tcounts$treatment, levels=c("V&S", "P&S"))

tcounts


#remove rows not included 
tcounts <- tcounts[!is.na(tcounts$treatment), ]
#convert timepoint into a factor
#tcounts$Timepoint <- factor(tcounts$Timepoint)
#statistical test
stat.test <- tcounts %>% group_by(timepoint) %>% t_test(expression~treatment, alternative = "less" ) 
stat.test <- stat.test %>% add_xy_position(x = "timepoint", dodge = 0.8)
stat.test
#plotting boxplot



g<- ggboxplot(tcounts, x= "timepoint", y= "expression", 
              add = "dotplot", title = goi, palette = 'npg',
              ylab = "Expression (Normalized gene count)", legend = "None",
              fill = "treatment",
              xlab = "Timepoint (days)") +
  theme(plot.title = element_text(hjust = 0.5))+
  grids(linetype = "solid")+
  stat_pvalue_manual(stat.test, label = "p",y.position = 3500, size = 8) +
  scale_y_continuous(limits=c(0,max(tcounts$expression)*1.3))

font <- 18
g <- ggpar(g, 
           font.main = font,
           font.x = font,
           font.y = font,
           font.xtickslab = font, 
           font.ytickslab = font)

g 
dev.off()


pdf(file = 
      "/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/pctx vs vehicle validation genes/drgs/Ankrd27.pdf", # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches
g
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()


####
####
####
####
####
####



#### making a catalog of ALL genes expression values 
####
####
####
####
####
####

goi <- res$row
stopifnot(all(goi %in% names(dds)))
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
        merge(colData(dds[goi,]), ., by="row.names") %>%
        gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

#arranging variables 
tcounts$timepoint <- factor(tcounts$timepoint, levels=c("d.4", "d.7", "d.11", "d.23"))
tcounts$treatment <- factor(tcounts$treatment, levels=c("V&S", "P&S", "P&C","V&C"))


pdf("multi-ggplot2-catalog-DRG.pdf")
for (i in goi) {
        p <- ggplot(filter(tcounts, gene==i), aes(timepoint, expression, fill=treatment)) + 
                geom_boxplot() + 
                facet_wrap(~gene, scales="free_y") + 
                labs(x="Timepoint", 
                     y="Expression (log normalized counts)", 
                     fill="Treatment" 
                )
                print(p)
}
dev.off()



#########






#write.table(example,".txt",sep="\t", row.names = TRUE)
#write.xlsx(example,".xlsx", row.names = TRUE)

# Save individual plotcounts per sample to variable
#plotCounts(x,....) x=dds, dds4, etc. 
d <- plotCounts(dds, gene = "Tubb3", intgroup = c("treatment", "timepoint"),returnData=TRUE, normalized = TRUE)
d <- t(d)
View(d)
barplot(d)
write.xlsx(d,"kif11.xlsx", row.names = TRUE)

##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################


##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
### making contrast tables with no filtering 
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################


#first need to get the fold change lists from comparing treatments 
res4 <- results(dds4,contrast = c("treatment", "P&C", "V&S"))
res7 <- results(dds7,contrast = c("treatment", "P&C", "V&S"))
res11 <- results(dds11,contrast = c("treatment", "P&C", "V&S"))
res23 <- results(dds23,contrast = c("treatment", "P&C", "V&S"))
res4 <- as.data.frame(res4)
res7 <- as.data.frame(res7)
res11 <- as.data.frame(res11)
res23 <- as.data.frame(res23)

res4 <- res4[,"log2FoldChange",drop = FALSE]
colnames(res4) = c("day4")
res7 <- res7[,"log2FoldChange",drop = FALSE]
colnames(res7) = c("day7")
res11 <- res11[,"log2FoldChange",drop = FALSE]
colnames(res11) = c("day11")
res23 <- res23[,"log2FoldChange",drop = FALSE]
colnames(res23) = c("day23")

View(res4)
View(res7)
View(res11)
View(res23)

#make final matrix
final <- cbind(res4,res7,res11,res23)
#remove rows with names containing NA
#final <- final[grepl("^NA", rownames(final))==F,]
final <- as.data.frame(final)
#final$gene <- rownames(final)
nrow(final)
View(final)

write.xlsx(final,"Paclitaxel+CL vs Vehicle_UNFILTERED.xlsx", rowNames = TRUE)




