#code for Skin 
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

#loading in the count matrix
getwd()
setwd("/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/SKIN/")
#now we have our huge dataFrame 
file = ("~/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/SKIN/RawMatrix/SKIN_CountMatrix_lanesCollated_woPairedReads.txt")

matrix <- read.table(file)
head(matrix)
View(matrix)

#read in sample info to match columns with sample details for proper 
#DESeq2 sorting/analysis 

setwd("/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/SKIN/")

sampleinfo <- read.csv("sampleLookup.txt", sep = "\t", header = TRUE)
View(sampleinfo)
sampleinfo <- sampleinfo[1:64,] #getting rid of all the DRG info
View(sampleinfo)
samplematch <- read.xlsx("sample_description.xlsx",colNames = TRUE)
View(samplematch)

#going to change days to single digit numbers....i.e. d.4 -> 4 

View(samplematch)
samplematch$Timepoint <-gsub('d.','',samplematch$Timepoint)


#trim off numbers to make columns match
colnames <- colnames(matrix)
View(colnames)
sampleinfo$col <- gsub('201926', 'X201926', sampleinfo$HIHG_ID)
sampleinfo$col2 <- gsub('-01', '', sampleinfo$col)
sampleinfo$id <- gsub('[0-9]+-','', sampleinfo$External_ID)
sampleinfo$id2 <- gsub('R-','', sampleinfo$id)
sampleinfo$id3 <- sampleinfo$id2
View(sampleinfo)
#fixing sample match

#samplematch$id <- gsub('[0-9]+ - ','',samplematch$X1)
#samplematch$id2 <- gsub('[(]','',samplematch$id)
#samplematch$id3 <- gsub('[)]','',samplematch$id2)
#samplematch$id3[34] <- c("DRG08")

samplematch$id3 <- samplematch$Sample
view(samplematch)
#merging columns by id3 column
samplesMatched <- merge(samplematch,sampleinfo, by = 'id3' )
View(samplesMatched)

#subset because we only need two columns 
samplesMatched <- as.data.frame(samplesMatched)
samplesMatched <- samplesMatched %>% dplyr::select(col2,id3,Treatment,Timepoint)
View(samplesMatched)
counts_origingal_columns <- matrix
View(counts_origingal_columns[1:100,])
View(matrix)

#going to change column names 
#samples matched is now the Col_data
coldata = samplesMatched
View(coldata)
colnames(coldata) = c("column","sample name","Treatment","Timepoint")
rownames(coldata) = coldata$column


#Need to replace column '#' at top with Skin#
View(matrix)
countMatrix <- as.data.frame(matrix)
names(countMatrix)
names(countMatrix) <- coldata$`sample name`[match(names(countMatrix), coldata$column)]
names(countMatrix)
head(countMatrix)


#order of columns and order on the colData need to be the same 
#sort by  'column' column in coldata
coldata <- arrange(coldata,`sample name`)
View(coldata)


#make new coldata without the old row#, keeping DRG#
#coldata$column = NULL
rownames(coldata) = coldata$`sample name`
View(coldata)

#Arrange columns in numerical order in countMatrix
#countMatrix2 has the SKIN names in order 

countMatrix2 <- countMatrix[,order(colnames(countMatrix))]
View(head(countMatrix2))
View(countMatrix2)

#make sure coldata now only has treatment and timepoint 
#making new column data with capital D
colData <- coldata[,c("Treatment","Timepoint")]
View(colData)
#we now have the count matrix with the samples ordered by DRG#
#now we can make sub-matriciesfor each day and colData information for our DEseq objects 
View(colData)
colData$sampleName <- rownames(colData)
View(countMatrix2)


#next i will convert the ensembl names to gene symbols before proceeding 
view(countMatrix2)
head(countMatrix2)
#need to remove the rows for # of mapped reads, multimapping, etc
row.names.remove <- c("N_multimapping", "N_noFeature","N_unmapped","N_ambiguous")
countMatrix2 <- countMatrix2[!(row.names(countMatrix2) %in% row.names.remove), ]
View(countMatrix2)
##
##changing the gene names to gene symbols
##(here are the setps taken to get to the code below)
#ensembl <- useEnsembl(biomart = "genes")
#datasets <- listDatasets(ensembl)
#View(datasets)
#making a column name from the row names 

matrix <- countMatrix2

matrix$ensembl_gene_id <- rownames(matrix)
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
#finding filters
#filters = listFilters(ensembl)
#View(filters)
#mgi_symbol
a <- getBM(attributes = c('mgi_symbol','ensembl_gene_id'),
         filters = 'ensembl_gene_id',
         values = matrix$ensembl, 
         mart = ensembl)
zz <- join(a,matrix,by ='ensembl_gene_id' )


#zz is the new RAW count matrix
#need to clean up  before saving RAW matrix


#need to set column name to sample name 

#now doing cleanup for DEseq2
#getting rid of ensemblID column 
zz$ensembl_gene_id <- NULL
#making the rowname the gene symbol - dealing with unique names 
rownames(zz) <- make.names(zz$mgi_symbol,unique = TRUE)
#getting rid of mgi_symbol column
zz$mgi_symbol <- NULL

head(zz)
#saving the RAW count matrix to the general directory 
setwd("/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/")
write.xlsx(zz,"RawCountMatrix_Skin.xlsx", rowNames = TRUE)

#saving the colData excel sheet for future analysis 
write.xlsx(coldata,"skin_colData.xlsx")

head(zz)
##############################################
##############################################
##############################################
#start of deseq2 analysis 
#the count matrix is called zz
##############################################
##############################################
##############################################

setwd('/Users/anthony/Desktop/RNAseq/BenSandra/')

#colData read in
experiment_description <- colData #loading in column/sample metadata
view(experiment_description)
#making metadata / colData for deseq objects with more appropriate variable names 


colData <- colData
day4colData <- colData[1:16,]
day7colData <- colData[17:32,]
day11colData <- colData[33:48,]
day23colData <- colData[49:64,]



#making separate tables to parse data 
allDays <- zz
day4 <- zz[,1:16]
day7 <- zz[,17:32]
day11 <- zz[,33:48]
day23 <- zz[,49:64]


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
#looking at RAW read counts 
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


#making a barchart for the counts per sample [all samples]
barplot(colSums(allDays), xlab = "sample", ylab = "number of reads",main = "All Samples")
dev.off() 
#making a barchart for the counts per sample [day 4 samples]
barplot(colSums(day4), xlab = "sample", ylab = "number of reads", main = "4 day Samples")
dev.off()
#making a barchart for the counts per sample [day 7 samples]
barplot(colSums(day7), xlab = "sample", ylab = "number of reads", main = "7 day Samples")
dev.off()
#making a barchart for the counts per sample [day 11 samples]
barplot(colSums(day11), xlab = "sample", ylab = "number of reads", main = "11 day Samples")
dev.off()
#making a barchart for the counts per sample [day 23 samples]
barplot(colSums(day23), xlab = "sample", ylab = "number of reads", main = "23 day Samples")
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
#making deseq objects for all days
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

#...All days
dds <- DESeqDataSetFromMatrix(countData = allDays, colData = colData, design = ~Treatment)
dds <- DESeq(dds)
#...4days
dds4 <- DESeqDataSetFromMatrix(countData = day4, colData = day4colData, design = ~Treatment)
dds4 <- DESeq(dds4)
#...7days
dds7 <- DESeqDataSetFromMatrix(countData = day7, colData = day7colData, design = ~Treatment)
dds7 <- DESeq(dds7)
#...11days
dds11 <- DESeqDataSetFromMatrix(countData = day11, colData = day11colData, design = ~Treatment)
dds11 <- DESeq(dds11)
#...23days
dds23 <- DESeqDataSetFromMatrix(countData = day23, colData = day23colData, design = ~Treatment)
dds23 <- DESeq(dds23)


#look for number of rows in deseq object
nrow(dds)
nrow(dds4)
nrow(dds7)
nrow(dds11)
nrow(dds23)

#52406 genes exactly in each 

#applying regularized log transformation of the count data ---- log2 transformation
rld = rlog(dds)
rld4 = rlog(dds4)
rld7 = rlog(dds7)
rld11 = rlog(dds11)
rld23 = rlog(dds23)


#applying vst transformation of the count data
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
#start here after LOADING ENVIRONMENT FILE
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
######heatmap of sample-sample distances 
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



#substituting names in vsd 

vsd$Treatment <- gsub('P&S', 'PCTX', vsd$Treatment)
vsd$Treatment <- gsub('P&C', 'PCTX+CL', vsd$Treatment)
vsd$Treatment <- gsub('V&S', 'VEH', vsd$Treatment)
vsd$Treatment <- gsub('V&C', 'CL', vsd$Treatment)

#now making sample distance matrix 

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Treatment, vsd$Timepoint, sep="-")
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9,"Spectral")) )(255)
p<-pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         fontsize_row = 8,
         fontsize_col = 8)


pdf(file = "/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/SKIN/sample distances/all_10_6_2021.pdf",
    height = 12,
    width = 12)
p
dev.off()


dev.off()


View(assay(vsd4))

sampleDists <- dist(t(assay(vsd4)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Treatment, vsd$Timepoint, sep="-")
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9,"Spectral")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         fontsize_row = 5,
         fontsize_col = 5)
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
##sample distances heatmap (only paclitaxel vs vehicle)
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


#normal - vst transformed 

#subsetting the vsd object 

vsd.sub <- vsd[,vsd$Treatment %in% c("P&S","V&S")]

#changing column names of P&S to PCTX and V&S to VEH???
#didn't do this but could 
array <- vsd.sub$Treatment
array[1] <- 'test'
arraynew <- vsd.sub$Treatment

sampleDists <- dist(t(assay(vsd.sub)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd.sub$Treatment, vsd.sub$Timepoint, sep="-")
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
#plotting PCA of interaction and individual groups
#printing & saving plots as PDF files 
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



#PCA with VSD transformation

pcaData <- plotPCA(vsd23,intgroup = c("Treatment"),ntop = nrow(vsd), returnData = TRUE) 
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca <- ggplot(pcaData, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  stat_ellipse(level=0.55)
  
pca

pdf(file = "/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/SKIN/PCA/pca23.pdf",
    height = 7,
    width = 7)
pca
dev.off()




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

vsd.sub <- vsd[,vsd$Treatment %in% c("P&S","V&S")]
vsd4.sub <- vsd4[,vsd4$Treatment %in% c("P&S","V&S")]
vsd7.sub <- vsd7[,vsd7$Treatment %in% c("P&S","V&S")]
vsd11.sub <- vsd11[,vsd11$Treatment %in% c("P&S","V&S")]
vsd23.sub <- vsd23[,vsd23$Treatment %in% c("P&S","V&S")]


#with VSD transformation
plotPCA(vsd.sub,intgroup = c("timepoint","treatment"))
plotPCA(vsd4.sub,intgroup = c("Treatment"),ntop = 25841)  #+
  #geom_text(aes(label=name),vjust=2) 
plotPCA(vsd7.sub,intgroup = c("Treatment"),ntop = 25841) #+
  #geom_text(aes(label=name),vjust=2)
plotPCA(vsd11.sub,intgroup = c("Treatment"),ntop = 25841) #+
  #geom_text(aes(label=name),vjust=2)
plotPCA(vsd23.sub,intgroup = c("Treatment"),ntop = 25841) #+
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


##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
#top 100 variably expressed genes
##############################################
##############################################
##############################################
##############################################

library("genefilter")

#substituting names in vsd 

vsd$Treatment <- gsub('P&S', 'PCTX', vsd$Treatment)
vsd$Treatment <- gsub('P&C', 'PCTX+CL', vsd$Treatment)
vsd$Treatment <- gsub('V&S', 'VEH', vsd$Treatment)
vsd$Treatment <- gsub('V&C', 'CL', vsd$Treatment)



topVarGenes <- head(order(-rowVars(assay(vsd))),500)
mat <- assay(vsd)[ topVarGenes, ]
mat<- mat - rowMeans(mat)
colnames(mat) <- paste(vsd$Treatment, vsd$Timepoint, sep="-")
df <- as.data.frame(colData(vsd)[,c("Treatment","Timepoint")])
p <- pheatmap(mat, cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "YlGnBu")))(1000),
         cellwidth = 5,
         cellheight = 1,
         fontsize_row = 3,
         fontsize_col = 5, 
         
)

pdf("/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/SKIN/skin_top250_var_genes.pdf",
    height = 10,
    width = 10)
p
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
#top 100 variably expressed genes (paclitaxel vs. vehicle)
##############################################
##############################################
##############################################
##############################################

#subsetting dds for just paclitaxel versus vehicle analysis 
dds.sub <- dds[ , dds$Treatment %in% c("P&S", "V&S") ]


library("genefilter")
topVarGenes <- head(order(-rowVars(assay(vsd.sub))),100)
mat <- assay(vsd.sub)[ topVarGenes, ]
mat<- mat - rowMeans(mat)
colnames(mat) <- paste(dds.sub$Treatment, dds.sub$Timepoint, sep="-")
df <- as.data.frame(colData(dds.sub)[,c("Treatment","Timepoint")])
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






#######################################
#######################################
#######################################
#######################################
#######################################
#dispersion estimates
#######################################
#######################################
#######################################


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




#####code for estimating size factors and showing what the different transformations
#####do to the data 


dds_new <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds_new, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation) 









##############################################
##############################################
##############################################
##############################################
##############################################
# Plotting individual gene normalized counts
# making excel matrices 
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################

getwd()

AAA <- counts(dds, normalized = TRUE) #making count matrix

View(AAA)
#P&S
#P&C
#V&S
#V&C
View(colData)

############################################matrix for V&C
VC4 <- AAA[,1:4] ################################4day column
VC4 <- as.data.frame(VC4)
colnames(VC4) <- c("4day","4day","4day","4day")

VC7 <- AAA[,21:24] ##############################7day column
VC7 <- as.data.frame(VC7)
colnames(VC7) <- c("7day","7day","7day","7day")

VC11 <- AAA[,37:40] ##############################11day column
VC11 <- as.data.frame(VC11)
colnames(VC11) <- c("11day","11day","11day","11day")
#view(VC11)
VC23 <- AAA[,53:56] ##############################23day column
VC23 <- as.data.frame(VC23)
colnames(VC23) <- c("23day","23day","23day","23day")

VC <- cbind(VC4,VC7,VC11,VC23) #here is our VC matrix
View(VC)
#check number of rows
nrow(VC)
#write to .xlsx file 
write.xlsx(VC,"Vehicle+CL.xlsx", rowNames = TRUE)


#keep this code for making tab delimited text files

############################################matrix for V&S
VS4 <- AAA[,5:8] ################################4day columns
VS4 <- as.data.frame(VS4)
colnames(VS4) <- c("4day","4day","4day","4day")

VS7 <- AAA[,17:20] ##############################7day columns
VS7 <- as.data.frame(VS7)
colnames(VS7) <- c("7day","7day","7day","7day")

VS11 <- AAA[,33:36] ##############################11day columns
VS11 <- as.data.frame(VS11)
colnames(VS11) <- c("11day","11day","11day","11day")

VS23 <- AAA[,49:52] ##############################23day columns
VS23 <- as.data.frame(VS23)
colnames(VS23) <- c("23day","23day","23day","23day")

VS <- cbind(VS4,VS7,VS11,VS23) #here is our Vehicle+Saline matrix

view(VS)
#write to .xlsx file 
write.xlsx(VS, 'Vehicle+Saline.xlsx', rowNames = TRUE)


############################################matrix for P&C
PC4 <- AAA[,13:16] ################################4day column
PC4 <- as.data.frame(PC4)
colnames(PC4) <- c("4day","4day","4day","4day")

PC7 <- AAA[,29:32] ##############################7day column
PC7 <- as.data.frame(PC7)
colnames(PC7) <- c("7day","7day","7day","7day")

PC11 <- AAA[,45:48] ##############################11day column
PC11 <- as.data.frame(PC11)
colnames(PC11) <- c("11day","11day","11day","11day")

PC23 <- AAA[,61:64] ##############################23day column
PC23 <- as.data.frame(PC23)
colnames(PC23) <- c("23day","23day","23day","23day")

PC <- cbind(PC4,PC7,PC11,PC23) #here is our PC matrix

#check number of rows
nrow(PC)
View(PC)
#write to .xlsx file 
write.xlsx(PC,"Paclitaxel+CL.xlsx", rowNames = TRUE)


############################################matrix for P&S
PS4 <- AAA[,9:12] ################################4day column
PS4 <- as.data.frame(PS4)
colnames(PS4) <- c("4day","4day","4day","4day")

PS7 <- AAA[,25:28] ##############################7day column
PS7 <- as.data.frame(PS7)
colnames(PS7) <- c("7day","7day","7day","7day")

PS11 <- AAA[,41:44] ##############################11day column
PS11 <- as.data.frame(PS11)
colnames(PS11) <- c("11day","11day","11day","11day")

PS23 <- AAA[,57:60] ##############################23day column
PS23 <- as.data.frame(PS23)
colnames(PS23) <- c("23day","23day","23day","23day")

#view(PS23)
PS <- cbind(PS4,PS7,PS11,PS23) #here is our PC matrix
 

#check number of rows
nrow(PS)


#pheatmap(VS[c("","Mmp9","Krt4","Krt8","Nox1","Nox4"),], cluster_cols = FALSE,scale ='row')

cview(PS)
#write to .xlsx file 
write.xlsx(PS,"Paclitaxel+Saline.xlsx", rowNames = TRUE)

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
#######
#making excel sheets that have the mean count per gene per day 
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

#vehicle and 

View(VC)
VC_4 <- rowMeans(VC[,0:4])
VC_7 <- rowMeans(VC[,5:8])
VC_11 <- rowMeans(VC[,9:12])
VC_23 <- rowMeans(VC[,13:16])
VC_means <- cbind(VC_4,VC_7,VC_11,VC_23)
View(VC_means)

#heatmap of means 

food <- VC_means[c("Mmp13","Mmp9","Mmp12","Mmp2","Timp4","Mfn2","Duox1"),]

pheatmap(food,scale = 'row',cluster_cols = FALSE)






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

goi <- c("Mmp7")
tcounts <- t((counts(dds[goi, ], normalized=TRUE, replaced=FALSE))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
#arranging variables 
tcounts$Timepoint <- factor(tcounts$Timepoint, levels=c("4", "7", "11", "23"))
#all treatments
tcounts$Treatment <- factor(tcounts$Treatment, levels=c("V&S","P&S","P&C","V&C"))
#PCTX and VEH only 
#tcounts$Treatment <- factor(tcounts$Treatment, levels=c("V&S","P&S"))

#remove rows not included 
tcounts <- tcounts[!is.na(tcounts$Treatment), ]
#convert timepoint into a factor
tcounts$Timepoint <- factor(tcounts$Timepoint)
#statistical test
stat.test <- tcounts %>% group_by(Timepoint) %>% t_test(expression~Treatment, alternative = "two.sided", ref.group = "V&S") 
stat.test <- stat.test %>% add_xy_position(x = "Timepoint", dodge = 0.8)
stat.test
#plotting boxplot

#a<-anova_test(tcounts, expression ~ Treatment + Timepoint * Treatment:Timepoint)
#View(a)


g<- ggboxplot(tcounts, x= "Timepoint", y= "expression", legend = "none",
               add = "dotplot", title = goi, palette = 'npg',
               ylab = "Expression (Normalized gene count)",
               fill = "Treatment",
               xlab = "Timepoint (days)") +
  theme(plot.title = element_text(hjust = 0.5))+
  grids(linetype = "solid")+
  stat_pvalue_manual(stat.test, label = "p",
                     y.position = 2500, size = 5, tip.length = 0.001, step.group.by = "Timepoint", step.increase = 0.2) +
  scale_y_continuous(limits=c(0,max(tcounts$expression)*1.5))

font <- 18
g <- ggpar(g, 
           font.main = font,
           font.x = font,
           font.y = font,
           font.xtickslab = font, 
           font.ytickslab = font,
           legend = "none")

g
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()


pdf(file = "/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/SKIN/gprofiler/gprofiler_gene plots/Ncf2.pdf", # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches

g
dev.off()
dev.off()
dev.off()
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
#####need to mak treatment contrast tables for STEM program that include the Log2FoldChanges 
#####looking at days 4,7,11,23
#####also making heatmaps 
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


#we are only interested in these 2 for now:
#ps vs vs*
#pc vs ps*

#first need to get the fold change lists from comparing treatments 
res4 <- results(dds4,contrast = c("Treatment", "P&C", "P&S"))
res7 <- results(dds7,contrast = c("Treatment", "P&C", "P&S"))
res11 <- results(dds11,contrast = c("Treatment", "P&C", "P&S"))
res23 <- results(dds23,contrast = c("Treatment", "P&C", "P&S"))
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


# Keeps rows with log2FoldChange  corresponding to > 1.5 or < -1.5 FoldChange 
res4.above <- res4[res4$log2FoldChange > log2(3/2), ] #day 4
res4.below <- res4[res4$log2FoldChange < log2(2/3), ]
res4.out <- rbind(res4.above, res4.below)

res7.above <- res7[res7$log2FoldChange > log2(3/2), ] #day 7
res7.below <- res7[res7$log2FoldChange < log2(2/3), ] 
res7.out <- rbind(res7.above, res7.below)

res11.above <- res11[res11$log2FoldChange > log2(3/2), ] #day 11
res11.below <- res11[res11$log2FoldChange < log2(2/3), ]
res11.out <- rbind(res11.above, res11.below)

res23.above <- res23[res23$log2FoldChange > log2(3/2), ] #day 23
res23.below <- res23[res23$log2FoldChange < log2(2/3), ]
res23.out <- rbind(res23.above, res23.below)

#make a cumulative list of genes that are significant from these lists 
list <- c(row.names(res4.out),row.names(res7.out),row.names(res11.out),row.names(res23.out))
length(list)
#remove duplicates from list
list <- unique(unlist(strsplit(list, " ")))
length(list)
#now we have a list of gene names that were significant from all of the time-points

#now will go back and get the log2fold changes for all of the genes in our list
res4 <- results(dds4,contrast = c("Treatment", "P&C", "P&S"))
res7 <- results(dds7,contrast = c("Treatment", "P&C", "P&S"))
res11 <- results(dds11,contrast = c("Treatment", "P&C", "P&S"))
res23 <- results(dds23,contrast = c("Treatment", "P&C", "P&S"))

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
setwd("/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/SKIN/Log2FoldChange Significant Genes/")
write.xlsx(final,"Paclitaxel+CL vs Paclitaxel_FC_1.5_pvalue_01.xlsx", rowNames = TRUE)


#look up specific gene to look at log2FC
x<-final[rownames(final)=='Mmp13',]

x$gene <- rownames(x)
x
rownames(x)
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
#make heatmap of foldchanges
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################

p<- pheatmap(final, cluster_rows= TRUE, show_rownames=TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="Spectral")))(10000),
         cluster_cols=FALSE,
         cellwidth = 42, 
         cellheight = 0.2,
         fontsize = 5,
         treeheight_row = 50,
         fontsize_row = 5,
         fontsize_col = 5,
         main = "Paclitaxel vs Vehicle", 
         angle_col = "0",
         cutree_rows = 15,
         scale = 'row'
         
)


pdf(file = "/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/SKIN/Log2FoldChange Significant Genes/paclitaxel vs vehicle.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches
p
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
#going to make bar graphs of number of upregulated and downregulated genes 
#per day per contrast 
#no filtering
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
#PC VS -
#PS PC -****
#PS VC -
#PC VC -

#   **** only interested in these contrasts for now 

#first need to get the fold change lists from comparing treatments 
res4 <- results(dds4,contrast = c("Treatment", "P&S", "V&S"))
res7 <- results(dds7,contrast = c("Treatment", "P&S", "V&S"))
res11 <- results(dds11,contrast = c("Treatment", "P&S", "V&S"))
res23 <- results(dds23,contrast = c("Treatment", "P&S", "V&S"))
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


setwd('/Users/anthony/Desktop/RNAseq/BenSandra/Data&Figures_03242021/')

counts(dds)

View(zz)

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
 
res4 <- results(dds4,contrast = c("Treatment", "P&S", "V&S"))
res7 <- results(dds7,contrast = c("Treatment", "P&S", "V&S"))
res11 <- results(dds11,contrast = c("Treatment", "P&S", "V&S"))
res23 <- results(dds23,contrast = c("Treatment", "P&S", "V&S"))
res4 <- as.data.frame(res4)
res7 <- as.data.frame(res7)
res11 <- as.data.frame(res11)
res23 <- as.data.frame(res23)

#only keeping LOG2FC column
#skip this block for unfiltered results 

View(res4)

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
#remember order is chronological from left to right 
final <- cbind(res4,res7,res11,res23)
#remove rows with names containing NA
#final <- final[grepl("^NA", rownames(final))==F,]
final <- as.data.frame(final)
nrow(final)
View(final)
write.xlsx(final,"Paclitaxel vs Vehicle_UNFILTERED.xlsx", rowNames = TRUE)

getwd()




###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
####going to make the gProfiler plots 
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
##first need to get significant gene list
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################


#first need to get the fold change lists from comparing treatments 
res4 <- results(dds4,contrast = c("Treatment", "P&S", "V&S"))
res7 <- results(dds7,contrast = c("Treatment", "P&S", "V&S"))
res11 <- results(dds11,contrast = c("Treatment", "P&S", "V&S"))
res23 <- results(dds23,contrast = c("Treatment", "P&S", "V&S"))
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
res4 <- results(dds4,contrast = c("Treatment", "P&S", "V&S"))
res7 <- results(dds7,contrast = c("Treatment", "P&S", "V&S"))
res11 <- results(dds11,contrast = c("Treatment", "P&S", "V&S"))
res23 <- results(dds23,contrast = c("Treatment", "P&S", "V&S"))

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






###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################

library(gprofiler2)
library(enrichplot)
#using rownames of the 'final' matrix to put into gprofiler 



gostres <- gost(rownames(final), organism = "mmusculus" , evcodes = TRUE)
goresults <- gostres$result
View(goresults)
goresults$adjustedP = -log10(goresults$p_value)


#going to select term names to keep from each category 

View(goresults)



#all biological process terms 
bp <- goresults %>% filter(grepl('GO:BP', source))
nrow(bp)
View(bp)

bp$term_name

list <- c(bp$term_name[8],bp$term_name[11],bp$term_name[15],bp$term_name[18],bp$term_name[20],bp$term_name[24],
          bp$term_name[25],bp$term_name[29],bp$term_name[32],bp$term_name[34],bp$term_name[33],
          bp$term_name[37],
          bp$term_name[38],
          bp$term_name[40],
          bp$term_name[43],
          bp$term_name[44],
          bp$term_name[45],
          bp$term_name[53],
          bp$term_name[68],
          bp$term_name[71],
          bp$term_name[74],
          bp$term_name[78],
          bp$term_name[80],
          bp$term_name[85],
          bp$term_name[91],
          bp$term_name[93],
          bp$term_name[101],
          bp$term_name[96],
          bp$term_name[106],
          bp$term_name[111],
          bp$term_name[119],
          bp$term_name[123],
          bp$term_name[128],
          bp$term_name[131],
          bp$term_name[130],
          bp$term_name[129],
          bp$term_name[132],
          bp$term_name[133],
          bp$term_name[136],
          bp$term_name[139],
          bp$term_name[142],
          bp$term_name[140],
          bp$term_name[179],
          bp$term_name[182],
          bp$term_name[225],
          bp$term_name[226],
          bp$term_name[236],
          bp$term_name[246],
          bp$term_name[241],
          bp$term_name[254],
          bp$term_name[257],
          bp$term_name[260],
          bp$term_name[262],
          bp$term_name[281],
          bp$term_name[300],
          bp$term_name[303],
          bp$term_name[315],
          bp$term_name[351],
          bp$term_name[352])

#filtering
          
bp <- bp[bp$term_name %in% list,] #only keeping the terms in the list 
bp <- bp[order(bp$adjustedP),] #reorder according to adjP
g <- ggbarplot(bp, "term_name", "adjustedP", fill = "term_name", title = "Biological process",
          legend = "none",lab.size = 0.1, xlab = '',ylab = "-log10(p-value)") +
  theme(axis.text.x = element_text(size = 10)) + scale_x_discrete(position = "top") + 
  coord_flip() + scale_x_discrete(position = "bottom")
pdf("/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/gprofiler/biologicalProcess.pdf",
    width = 11,
    height = 11)
g
dev.off()

write.xlsx(bp, file = "/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/gprofiler/bp.xlsx")





# Cellular component 
cc <- goresults %>% filter(grepl('GO:CC', source))
#filtering terms 
View(cc)

#sorting by adjusted p-value 
cc <- cc[order(cc$adjustedP),] #reorder according to adjP

g <- ggbarplot(cc, "term_name", "adjustedP", fill = "term_name", title = "Cellular component", legend = "none",
          ylab = "-log10(p-value)",
          xlab = '',
          lab.size = 0.1) +
  theme(axis.text.x = element_text(size = 10))+
  scale_x_discrete(position = "top") + 
  coord_flip() + scale_x_discrete(position = "bottom")+
  rotate_x_text(90)
g

pdf("/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/gprofiler/CellularComponent.pdf")
g
dev.off()

write.xlsx(cc, file = "/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/gprofiler/cc.xlsx")




#all KEGG
kegg <- goresults %>% filter(grepl('KEGG', source))
#reorder according to adjP
kegg<- kegg[order(kegg$adjustedP),]
g <- ggbarplot(kegg, "term_name", "adjustedP", fill = "term_name", title = "KEGG pathway",
          legend = "none",lab.size = 0.1, ylab = "-log10(p-value)", 
          palette = "ggplot2",
          xlab = '') + theme(axis.text.x = element_text(size = 10)) + scale_x_discrete(position = "top") + 
  coord_flip() + scale_x_discrete(position = "bottom")

pdf("/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/gprofiler/kegg.pdf")
g
dev.off()

write.xlsx(kegg, file = "/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/gprofiler/kegg.xlsx")




##REAC [reactome]
reac <- goresults %>% filter(grepl('REAC', source))
#reorder according to adjP
reac<- reac[order(reac$adjustedP),]
g <- ggbarplot(reac, "term_name", "adjustedP", fill = "term_name", title = "Reactome pathway",
               legend = "none",lab.size = 0.1, ylab = "-log10(p-value)", 
               palette = "ggplot2",
               xlab = '') + theme(axis.text.x = element_text(size = 10)) + scale_x_discrete(position = "top") + 
  coord_flip() + scale_x_discrete(position = "bottom")

reac

pdf("/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/gprofiler/reac_.pdf")
g
dev.off()

write.xlsx(reac, file = "/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/SKIN/gprofiler/REAC.xlsx")


##all HP
hp <- goresults %>% filter(grepl('HP', source))
#reorder according to adjP
hp<- hp[order(hp$adjustedP),]
View(hp)




#filtering list by terms
#first looking at terms 

hp$term_name
list <- c(hp$term_name[1],hp$term_name[17],hp$term_name[18],hp$term_name[22],hp$term_name[23])
#checking list 
list


#generating graph 
hp <- hp[hp$term_name %in% list,] #only keeping the terms in the list 
hp <- hp[order(hp_filtered$adjustedP),] #reorder according to adjP
g <- ggbarplot(hp, "term_name", "adjustedP", fill = "term_name", title = "Human process",
               legend = "none",lab.size = 0.1, ylab = "-log10(p-value)", 
               palette = "ggplot2",
               xlab = '') + theme(axis.text.x = element_text(size = 10)) + scale_x_discrete(position = "top") + 
  coord_flip() + scale_x_discrete(position = "bottom")


g #previewing graph

#saving graph 
pdf("/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/gprofiler/human process.pdf")
g
dev.off()

#saving xlsx sheet

write.xlsx(hp, file = "/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/SKIN/gprofiler/hp.xlsx")




#####all GO;MF (molecular function)
mf <- goresults %>% filter(grepl('GO:MF', source))

#reorder according to adjP
mf<- mf[order(mf$adjustedP),]


g <- ggbarplot(mf, "term_name", "adjustedP", fill = "term_name", title = "Molecular function",
               legend = "none",lab.size = 0.1, ylab = "-log10(p-value)", 
               palette = "ggplot2",
               xlab = '') + theme(axis.text.x = element_text(size = 10)) + scale_x_discrete(position = "top") + 
  coord_flip() + scale_x_discrete(position = "bottom")


g #previewing graph

#saving graph 
pdf("/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/gprofiler/molecular function.pdf")
g
dev.off()

#saving xlsx worksheet
write.xlsx(mf, file = "/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/SKIN/gprofiler/mf.xlsx")






####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####going to try timecourse analysis
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################

####################################
####################################
####################################
####################################
##this code below is already included in enviornment file 
####################################
####################################

full_model <- ~ Treatment + Timepoint + Treatment:Timepoint
reduced_model <- ~ Treatment + Timepoint

ddsTC <- DESeqDataSetFromMatrix(countData = allDays, colData = colData, design = ~ Treatment + Timepoint + Treatment:Timepoint)
dds_lrt_time <- DESeq(ddsTC, test="LRT", reduced = ~ Treatment + Timepoint)


#####
####
###

# Extract results for LRT
res_LRT <- results(dds_lrt_time)
# View results for LRT
res_LRT  


# Create a tibble for LRT results
res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

View(res_LRT_tb)



clusters <- degPatterns(cluster_rlog, metadata = meta, time="Timepint", col="Treatment")














