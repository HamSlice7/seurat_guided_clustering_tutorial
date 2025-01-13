library(dplyr)
library(Seurat)
library(patchwork)

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




