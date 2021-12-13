#PCA for all samples 


#loading in the count matrices (SKIN and DRG) and concatenating them 
getwd()
setwd("/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/PCA both tissues/")

#now we have our huge dataFrame 
skin = ("~/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/PCA both tissues/RawCountMatrix_Skin.xlsx")
drg = ("~/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/PCA both tissues/RawCountMatrix_DRGs.xlsx")

skin <- read.xlsx(skin)
drg <-read.xlsx(drg)

#look at matrices
head(skin)
head(drg)

#number of rows in each ...are they the same 
nrow(skin)
nrow(drg)

#need to concatenate the sample info columns to get sample info 
setwd("/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/PCA both tissues/")
skin_colData = ("~/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/PCA both tissues/skin_colData.xlsx")
drg_colData = ("~/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/PCA both tissues/drg_colData.xlsx")

skin_colData <- read.xlsx(skin_colData)
drg_colData <-read.xlsx(drg_colData)

#remove the 'column'column in the skin colData

skin_colData$column = NULL
View(skin_colData)
View(drg_colData)
#need to combine the colData matricies 

#first to combine the columns they need to have the same name 
#need to make sure the column names are the same
colnames(skin_colData)
colnames(drg_colData)

#will change the skin_colData rows to match drg_colData
colnames(skin_colData) <- c("sample.name","treatment","timepoint")
View(skin_colData)

colData <- rbind(skin_colData,drg_colData)
View(colData)

#saving combined colData for future reference
write.xlsx(colData, 'colData.xlsx')
View(colData)
##will concatenate the DRG and Skin matricies
#need to make sure they are in the same order as the colData
View(colData)
#going to add tissue column to colData
#first making a list of column names without the numbers
tissue <- sub("[[:digit:]][[:digit:]]", "", colData$sample.name)
#then going to make this a new column 
colData$tissue <- tissue
#skin goes first (left side), drgs(right side)
nrow(skin)
nrow(drg)
master <- cbind(skin,drg)
View(master)
View(skin)
#setting rownames in master matrix as gene names 
rownames(master) <- master$Var.1
View(master)
#get rid of var column 
master$Var.1<-NULL

#checking number of columns, should be 128
ncol(master)
colnames(master)

#need to get rid of extra column
master$Var.66 <-NULL

#now we will make our DESeq object 


##############################################
##############################################
##############################################
#start of deseq2 analysis 
#the count matrix is called zz
##############################################
##############################################
##############################################

View(colData)

#making column data for skin, drg, and all 

#all tissues 
all_colData <- colData
all_day4colData <- rbind(colData[1:16,],colData[65:80,])
all_day7colData <- rbind(colData[17:32,],colData[81:96,])
all_day11colData <- rbind(colData[33:48,],colData[97:112,])
all_day23colData <- rbind(colData[49:64,],colData[113:128,])

nrow(all_day4colData)
nrow(all_day7colData)
nrow(all_day11colData)
nrow(all_day23colData)

#making separate tables to parse data 
#all tissues 

all_allDays <- master
all_day4 <- cbind(master[,1:16],master[,65:80])
all_day7 <- cbind(master[,17:32],master[,81:96])
all_day11 <- cbind(master[,33:48],master[,97:112])
all_day23 <- cbind(master[,49:64],master[,113:128])


ncol(all_day4)
ncol(all_day7)
ncol(all_day11)
ncol(all_day23)

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
barplot(colSums(all_allDays), xlab = "sample", ylab = "number of reads",main = "All Samples")
dev.off() 
#making a barchart for the counts per sample [day 4 samples]
barplot(colSums(all_day4), xlab = "sample", ylab = "number of reads", main = "4 day Samples")
dev.off()
#making a barchart for the counts per sample [day 7 samples]
barplot(colSums(all_day7), xlab = "sample", ylab = "number of reads", main = "7 day Samples")
dev.off()
#making a barchart for the counts per sample [day 11 samples]
barplot(colSums(all_day11), xlab = "sample", ylab = "number of reads", main = "11 day Samples")
dev.off()
#making a barchart for the counts per sample [day 23 samples]
barplot(colSums(all_day23), xlab = "sample", ylab = "number of reads", main = "23 day Samples")
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
#making deseq objects for all days per tissue types
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
#all tissues 
##############################################
##############################################
##############################################

##############################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
##############################################
#...All days
all_dds <- DESeqDataSetFromMatrix(countData = all_allDays, colData = all_colData, design = ~treatment)
all_dds <- DESeq(all_dds)
#...4days
all_dds4 <- DESeqDataSetFromMatrix(countData = all_day4, colData = all_day4colData, design = ~treatment)
all_dds4 <- DESeq(all_dds4)
#...7days
all_dds7 <- DESeqDataSetFromMatrix(countData = all_day7, colData = all_day7colData, design = ~treatment)
all_dds7 <- DESeq(all_dds7)
#...11days
all_dds11 <- DESeqDataSetFromMatrix(countData = all_day11, colData = all_day11colData, design = ~treatment)
all_dds11 <- DESeq(all_dds11)
#...23days
all_dds23 <- DESeqDataSetFromMatrix(countData = all_day23, colData = all_day23colData, design = ~treatment)
all_dds23 <- DESeq(all_dds23)

#look for number of rows in deseq object
nrow(dds)
nrow(dds4)
nrow(dds7)
nrow(dds11)
nrow(dds23)

#applying regularized log transformation of the count data ---- log2 transformation
#rld = rlog(dds)
#rld4 = rlog(dds4)
#rld7 = rlog(dds7)
#rld11 = rlog(dds11)
#rld23 = rlog(dds23)

#applying vst transformation of the count data
#recc'd for n >30 , we have 128
#http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#the-variance-stabilizing-transformation-and-the-rlog
all_vsd <- vst(all_dds, blind = TRUE)
all_vsd4 <- vst(all_dds4, blind = TRUE)
all_vsd7 <- vst(all_dds7, blind = TRUE)
all_vsd11 <- vst(all_dds11, blind = TRUE)
all_vsd23 <- vst(all_dds23, blind = TRUE)



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
pcaData <- plotPCA(all_vsd11,intgroup = c("tissue"),ntop = nrow(vsd), returnData = TRUE) 

percentVar <- round(100 * attr(pcaData, "percentVar"))

pca <- ggplot(pcaData, aes(PC1, PC2, color=tissue)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  stat_ellipse(level=0.999)

pca

?pdf(file = "/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/PCA both tissues/pca.pdf",
    height = 7,
    width = 7)
pca
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
all_vsd$treatment <- gsub('P&S', 'PCTX', all_vsd$treatment)
all_vsd$treatment <- gsub('P&C', 'PCTX+CL', all_vsd$treatment)
all_vsd$treatment <- gsub('V&S', 'VEH', all_vsd$treatment)
all_vsd$treatment <- gsub('V&C', 'CL', all_vsd$treatment)

#now making sample distance matrix 
sampleDists <- dist(t(assay(all_vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(all_vsd$treatment, all_vsd$timepoint,all_vsd$tissue, sep="-")
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9,"Spectral")) )(255)
p<-pheatmap(sampleDistMatrix,
            clustering_distance_rows=sampleDists,
            clustering_distance_cols=sampleDists,
            col=colors,
            fontsize_row = 8,
            fontsize_col = 8)

pdf(file = "/Users/anthony/Desktop/RNAseq/BenSandra RNAseq/Code&Data&Figures/PCA both tissues/sd.pdf",
    height = 20,
    width = 20)
p
dev.off()

#DESeq2 sorting/analysis 