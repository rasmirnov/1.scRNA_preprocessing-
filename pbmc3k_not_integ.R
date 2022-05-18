setwd("/mnt/tank/scratch/rasmirnov/code/Time_benchmarking/Datasets/pbmc3k")

#install.packages("BiocManager")
#BiocManager::install("MAST")

# STEP 1:Setup the Seurat Object
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

set.seed(1)

# STEP_1: Create a variable for seurat object
pbmc3k <- pbmc3k.SeuratData::pbmc3k
# ifnb <- ifnb.SeuratData::ifnb

# STEP 2: Standard Pre-processing workflow

# main steps: Quality Control (QC),filter cells (user-determined criteria), data normalization 

# 2.1_QC Metrics:

# The [[ operator add columns to our object --> store QC data
pbmc3k[["percent.mt"]] <- PercentageFeatureSet(pbmc3k, pattern = "^MT-")   # PercentageFeatureSet - mitochondrial QC metric (counts mit.genes&)
ifnb[["percent.mt"]] <- PercentageFeatureSet(ifnb, pattern = "^Mt-")

# Visualize QC metrics by a violin pplot
VlnPlot(pbmc3k, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
ggsave("VlnPlot_test.png", width = 5, height = 4)
# FeatureScatter needs to visualize feature-feature (but we can use it just for columns in object metadata)

plot1 <- FeatureScatter(pbmc3k, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc3k, feature1 = "nCount_RNA", feature2 ="nFeature_RNA")
plot1 + plot2
ggsave("Features_test.png", width = 7, height = 4)

#! We filtered 1)cells with RNA_features (unique) 200<x<2500 and x<5% mitochondirial reads
pbmc3k <- subset(pbmc3k, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


# 2.2_ Normalizing the data
#pbmc3k <- NormalizeData(pbmc3k, normalization.method = "LogNormalize", scale.factor = 10000)        # pbmc3k <- NormalizeData(pbmc3k) is the same (other parameters by default)
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("glmGamPoi")
pbmc3k <- SCTransform(pbmc3k, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE) 

# 2.3_Identification of highly variable Features
pbmc3k <- FindVariableFeatures(pbmc3k, 
                               selection.method = "vst", 
                               nfeatures = 2000)   #? vst read

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc3k), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc3k)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)          # repel google
plot1 + plot2
ggsave("VarFeatures_test.png", width = 10, height = 4)

# 2.4_Scaling the data (log scale transformation)
all.genes <- rownames(pbmc3k)
pbmc3k <- ScaleData(pbmc3k, features = all.genes)


# STEP 3: LINEAR DIMENSIONAL REDUCTION (PCA)
pbmc3k <- RunPCA(pbmc3k, features = VariableFeatures(object = pbmc3k))    # + check on the sceled data, so without features arg

# There are 3 key ways to visualize both cells and features that used in PCA:

# Firstly, let's check 5 different features of PCA
print(pbmc3k[["pca"]], dims = 1:5, nfeatures = 5)

# Way 1: VizDimLoadings
VizDimLoadings(pbmc3k, dim = 1:2, reduction = "pca")
ggsave("PCA_VizDim_test.png", width = 5, height = 6)

# STEP 4: DETERMINE DIMENSIONALITY OF THE DATASET (HOW MUCH PCs CHOOSE?)
ElbowPlot(pbmc3k)
ggsave("ElbowPlot_test.png", width = 4, height = 3)

# STEP 5: CLUSTER CELLS
pbmc3k <- FindNeighbors(pbmc3k, dims = 1:20)
pbmc3k <- FindClusters(pbmc3k, resolution = c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 5))

# Look at the cluster IDs of the first 5 cells
head(Idents(pbmc3k), 5)


# STEP 6: RUN NON-LINEAR DIMENSIONAL REDUCTION (UMAP/tSNE)
pbmc3k <- RunUMAP(pbmc3k, dims = 1:20)
DimPlot(pbmc3k, reduction = "umap")
ggsave("UMAP_test.png", width = 5, height = 4)

# STEP 7: Cell annotation function
analyze_object <- function(object, ident) {
  Idents(object) <- object[[ident]]
  if (length(levels(object)) == 1) {
    return(message(sprintf('%s: since only one cluster was identified, markers can not be found', ident)))
  }
  out_dir <- paste0('markers/', ident)
  dir.create(out_dir, recursive = T)
  whole.markers <- FindAllMarkers(object = object,
                                  assay='RNA',
                                  only.pos = TRUE,
                                  min.pct = 0.10,
                                  test.use = 'MAST')
  write.table(whole.markers, paste(out_dir, "markers.tsv", sep = '/'), sep="\t", quote=F, row.names=F)
}

idents <- grep('res', colnames(pbmc3k@meta.data), value = T)
sapply(c(idents, 'seurat_annotations'), function(ident) analyze_object(object = pbmc3k, ident = ident))



# SAVE
#save(pbmc3k, file = "/mnt/tank/scratch/rasmirnov/code/Time_benchmarking/Results/Standard_pipeline/pbmc3k_tutorial.RData")
save.image('pbmc3k_tutorial.RData')
