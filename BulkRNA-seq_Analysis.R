#Example of analysis of bulk RNA sequencing data
#Here, analysis for hESC-macrophages sorted from Biowires (human cardiac microtissues)

library(DESeq2)
library(ggplot2)
library(Matrix)
library(dplyr)
library(xlsx)
library(RColorBrewer)
library(pheatmap)
library(gplots)

#Set working directory
setwd("~/Desktop/Biowire_MFSort_BulkRNA-seq/Biowire analysis")

#CSV file with file names corresponding to sample names and conditions
sampleTable_HH <- read.csv("sampleTable_biowires.csv", sep= ",", header = T)

#read HT-seq counts
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_HH,
                                  directory = "~/Desktop/Human Cardiac Macrophages/Experiments - Biowire/Transcriptomics/2022 05 Biowire_MFSort_BulkRNA-seq/RawData_Biowire",
                                  design= ~condition)
dds
#dim: 26485 9 

dds <- dds[ rowSums(counts(dds)) > 1, ]
dds
#dim: 19206 9 

#set reference condition
dds$condition <- relevel(dds$condition, ref="MF")

#####################
#CMFBMF vs FBMF
#find DEGs
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "CM_FB_MF", "FB_MF"))

write.xlsx(res,"Biowires_CMFBMFvsFBMF_DEGs_NonFiltered.xlsx")

resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, pvalue < 0.05)
resSig

resOrdered <- res[order(res$log2FoldChange),]
up_and_down <- subset(resOrdered, padj<0.05 & abs(resOrdered$log2FoldChange)>=0.4)
down <- subset(res, padj<0.05 & log2FoldChange < -0.3)
up <- subset(res, padj<0.05 & log2FoldChange > 0.3)

write.xlsx(up,"Biowires_CMFBMFvsFBMF_DEGs_Filtered_Up_LCF0.4.xlsx")
write.xlsx(down,"Biowires_CMFBMFvsFBMF_DEGs_Filtered_Down_LFC0.4.xlsx")

#raw counts
rawCounts <- counts(dds)
#Normalization with DESeq2
dds<- estimateSizeFactors(dds)
sizeFactors(dds)
normalizedCounts <- counts(dds, normalized=TRUE)

write.xlsx(rawCounts,"Biowires_Raw_Counts.xlsx")
write.xlsx(normalizedCounts,"Biowires_Normalized_Counts.xlsx")

#PCA
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
assay(vsd)

plotPCA(vsd, intgroup="condition")
plotPCA(vsd, intgroup="id") 

pcaData <- plotPCA(vsdata, intgroup=c("condition", "id"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=id)) +
  geom_point(size=3, shape = 16) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

###VOLCANO PLOTS
library(EnhancedVolcano)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "CM_FB_MF", "MF"))
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',  pCutoff = 0.05,
                FCcutoff = 1)


#Scoring with in vivo signatures
library(singscore)
genes <- scan(text = ("COLEC12 SCN9A FXYD6 MAMDC2 PMP22 TTN HRH1 TGFBI EMB ITPR2 CCDC141 FOXO3 WLS MAN1A1 IGF1 CCL13 AHNAK ARHGAP5 KCTD12 GAS6 PPIC LRP1 SLC9A9 CD163L1 ZFHX3 GNG11 ZFP36L1 CYBB CFD FGFR1 ITSN1 EGFL7 KANK2 IGFBP4 TLR4 TMEM176B NINJ1 CD28 PER3 NBL1 DNM1 BCAT1 AKAP13 DPYSL3 SLC4A7 DUSP6 PDE4D NR1D2 C1QC LIMCH1 NORAD EPB41L1 EMP1 SYK ANTXR2 HOXB6 AP2A2 OSBPL1A IQGAP2 ARHGAP1 SPRED1 CCDC186 ANO6 CPM FTX PCDH12 TNS1 CALM2 TPM1 EPN2 GOLIM4 GNG12 TMEM176A POGZ PSAP RERE SGMS1 MYO5A MGAT1 GPR34 EPS15 HPGDS CD209 SNX2 GPX3 HOXB7 RB1CC1 CREBRF MYCBP2 ZCCHC24 TACC1 ARHGAP18 RNF213 CD14 CCND1 IFI16 RBMS1 SH3PXD2A"), what = "")

rankData <- rankGenes(dds)
scoredf <- simpleScore(rankData, upSet = genes)
scoredf

plotRankDensity(rankData[,3,drop = FALSE], upSet = genes, 
                isInteractive = FALSE)

scoredf$group <- c("CM+FB+MF", "CM+FB+MF", "CM+FB+MF", "FB+MF", "FB+MF", "FB+MF", "MF", "MF", "MF")
boxplot(scoredf$TotalScore~scoredf$group)




