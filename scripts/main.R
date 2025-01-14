library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
#Load the PMBC dataset
pbmc.data <- Read10X("data/filtered_gene_bc_matrices/hg19/") 

# Initialize the Seurat object with the raw (non-normalized data)
pbmc <- CreateSeuratObject(counts = pmbc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

#Examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"),1:30]

#The [[ operator can add columns to object metadata. This is a great place to stash QC stats. The percentage of MT genes in each cell
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Based on visualizing QC metrics, filter count data
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalizing the counts data. This line is the same as pbmc <- NormalizeData(pbmc)
pbmc <- NormalizeData(pbmc, normalize.method = "LogNormalize", scale.factor = 10000)

#Identifying highly variable features (genes)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc) ,10)

#plot variable features with and without lables
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling the data before PCA
all.genes <- rownames(pbmc) #getting a list of all the genes
pbmc <- ScaleData(pbmc, features = all.genes)

#Performing PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca") #Visualize top loadings
DimPlot(pbmc, reduction = "pca") + NoLegend() #Plot first 2 PC's
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) #Heatmap of first PC
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#Plotting an elbow plot. Appears that first 10 PC's capture most of the data
ElbowPlot(pbmc)

##Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10) #computes a similarity graph between cells based on their gene expression data, using a reduced number of dimensions (principal components)
pbmc <- FindClusters(pbmc, resolution = 0.5) #uses this graph to identify clusters of cells, which are groups of cells with similar expression profiles.

#look at cluster IDs of the first 5 cells
head(Idents(pbmc),5)

pmbc <- RunUMAP(pmbc, dims = 1:10)

#Running non-linear dimensional reduction (UMAP/tSNE)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

#Saving the object up to this point so that it can be easily loaded back in without reruning the above computational steps
saveRDS(pbmc, file = "output/pbmc_tutorial.rds")

##Finding differential expressed features (cluster biomarkers)

# find all markers of cluster 2 compared to all other clusters
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0,3))
head(cluster5.markers)

#Find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>% 
  group_by(cluster) %>% 
  filter(avg_log2FC > 1)

#Using ROC test for DE
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#Using a violin plot to plot the expression level of a gene in the different clusters
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

#plotting raw counts
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

#feature plot
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

#for each cluster, take the top 10 markers for each cluster (all markers if less than 10 in a cluster)
pbmc.markers %>% 
  group_by(cluster) %>% 
  filter(avg_log2FC > 1) %>% 
  slice_head(n = 10) %>% 
  ungroup() -> top10

#Heat map of the top 10 marker genes for each cluster
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

##Assigning cell type identity to clusters based on canonical markers
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size=10)))

ggsave(filename = "output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

saveRDS(pbmc, file = "output/pbmc3k_final.rds")

