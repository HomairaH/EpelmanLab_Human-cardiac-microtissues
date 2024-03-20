library(Seurat)
library(ggplot2)
library(dplyr)

#Example of analysis of human yolk sac data. Similar analysis was performed for human fetal macrophages and cardiomyocytes

#Load datasets
yolksac_4wk_1_rawdata <- Read10X("~/Documents/Fetal_RawData/YS/FCAImmP7504909/filtered")
yolksac_4wk_1 <- CreateSeuratObject(counts = yolksac_4wk_1_rawdata,  min.cells = 3, min.features = 200, project = "YolkSac_4wk_1")

yolksac_4wk_2_rawdata <- Read10X("~/Documents/Fetal_RawData/YS/FCAImmP7504910/filtered")
yolksac_4wk_2 <- CreateSeuratObject(counts = yolksac_4wk_2_rawdata,  min.cells = 3, min.features = 200, project = "YolkSac_4wk_2")

yolksac_4wk_3_rawdata <- Read10X("~/Documents/Fetal_RawData/YS/FCAImmP7504911/filtered")
yolksac_4wk_3 <- CreateSeuratObject(counts = yolksac_4wk_3_rawdata,  min.cells = 3, min.features = 200, project = "YolkSac_4wk_3")

yolksac_4wk_4_rawdata <- Read10X("~/Documents/Fetal_RawData/YS/FCAImmP7504912/filtered")
yolksac_4wk_4 <- CreateSeuratObject(counts = yolksac_4wk_4_rawdata,  min.cells = 3, min.features = 200, project = "YolkSac_4wk_4")

yolksac_4wk_5_rawdata <- Read10X("~/Documents/Fetal_RawData/YS/FCAImmP7504913/filtered")
yolksac_4wk_5 <- CreateSeuratObject(counts = yolksac_4wk_5_rawdata,  min.cells = 3, min.features = 200, project = "YolkSac_4wk_5")

#Merge all 5 datasets
YS <- merge(x = yolksac_4wk_1, y = c(yolksac_4wk_2, yolksac_4wk_3, yolksac_4wk_4, yolksac_4wk_5))
YS@assays

#Quality control and filtering
YS[["percent.mito"]] <- PercentageFeatureSet(YS, pattern = "^MT-")

VlnPlot(YS, features = c("nFeature_RNA", "nCount_RNA", "yolksac_4wk_1.percent.mito", "yolksac_4wk_1.percent.dag"), ncol = 4)

plot1 <- FeatureScatter(YS, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(YS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

YS <- subset(YS, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mito < 10)

#Normalization, dimensionality reduction, clustering and visualization
#Regressing out the percentage of transcripts mapping to mitochondrial genes and nCount_RNA
YS <- SCTransform(YS, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = TRUE, return.only.var.genes = F)
YS <- RunPCA(YS, verbose = FALSE)
ElbowPlot(YS)
YS <- FindNeighbors(YS, dims = 1:20, verbose = FALSE)
YS <- FindClusters(YS, verbose = FALSE, resolution = 0.4)
YS <- RunUMAP(YS, dims = 1:20, verbose = FALSE, n.neighbors = 30, min.dist = 0.3)
DimPlot(YS, label = TRUE) 
FeaturePlot(YS, features = c("C1QC", "TIMD4", "LYVE1", "FOLR2", "FCN1", "SPINK2"), order = T)

#Generate heatmap of the average expression of the top 30 DEGs in macrophages vs mono-MF vs progenitors
Idents(YS) <- "MF"
YS.markers <- FindAllMarkers(YS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
YS_sig <- YS.markers[which(YS.markers$p_val_adj < 0.01), ]
YS_sig <- YS_sig[!grepl("^MT-", rownames(YS_sig)), ]
YS_sig <- YS_sig[!grepl("^RP", rownames(YS_sig)), ]

top30 <- YS_sig %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
avg <- AverageExpression(YS, features = top30$gene, return.seurat = T)
DoHeatmap(avg, features=rownames(avg[["SCT"]]), draw.lines = FALSE) + theme(text = element_text(size = 4))

