
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(data.table)
library(argparse)
library(tidyr)
library(glmGamPoi)                #?

set.seed(1)

parser <- ArgumentParser(description = 'Integration and processing of scrna-seq datasets')

parser$add_argument('-d',
                    '--data', 
                    type = 'character',
                    help = 'path to raw data')
parser$add_argument('-o',
                    '--obj',
                    type = 'character',
                    help = 'name of an object')

args <- parser$parse_args()

print(args)
####### в терминале!
# Rscript /mnt/tank/scratch/rasmirnov/code/Time_benchmarking/Scripts/Seurat_integration/unified_integ.R 
# -d /mnt/tank/scratch/rasmirnov/code/Time_benchmarking/Datasets/raw/retin -o retin 

# args$data <- '/mnt/tank/scratch/rasmirnov/code/Time_benchmarking/Datasets/raw/liver'
# args$obj <- 'liver'
# print(args)

setwd(args$data)

# setwd("/mnt/tank/scratch/rasmirnov/code/Time_benchmarking/Datasets/raw/oligo")
# mat <- fread("exprMatrix.tsv.gz")                             # через paste0 прилепить путь без setwd?
# meta <- data.frame(fread("meta_new.tsv"), row.names=1)

# STEP_1: Create a variable for seurat object: animal-mouse 
mat <- fread(paste0(args$data, '/', "exprMatrix.tsv.gz"))                           
meta <- data.frame(fread(paste0(args$data, '/', "meta_new.tsv")), row.names=1)
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)
### name of project
object <- CreateSeuratObject(counts = mat, project = args$obj, meta.data=meta)

# object <- CreateSeuratObject(counts = mat, project = 'oligo', meta.data=meta)
# object@meta.data %>% View()

#### 2.QC
# grep("^MT", rownames(object), value = T)
object@meta.data$Percent.Mitochond. <- NULL
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT")

###########! FOR EMBRYO DATASET ONLY
# object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Remove NA if there is someone
object <- subset(object, subset = sample != "NA")
object<- subset(object, subset = nCount_RNA > 0)
# object@meta.data %>% View()

#### 3.Split it based on sample:
object.list <- SplitObject(object, split.by = "sample")                 

### SCTransform
object.list <- lapply(X = object.list, 
                       FUN = SCTransform, 
                       method = "glmGamPoi",
                       vars.to.regress = "percent.mt",                
                       verbose = FALSE)

features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = 3000)
object.list <- PrepSCTIntegration(object.list = object.list, 
                                   anchor.features = features)

#### 4.Performing integration with SCTransform
immune.anchors <- FindIntegrationAnchors(object.list = object.list, 
                                         normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, 
                                     normalization.method = "SCT")


# Identification of highly variable Features
immune.combined.sct <- FindVariableFeatures(immune.combined.sct, 
                                            selection.method = "vst", 
                                            nfeatures = 2000) 

# Scaling the data (log scale transformation)
all.genes <- rownames(immune.combined.sct)
immune.combined.sct <- ScaleData(immune.combined.sct, 
                                 features = all.genes)

#### 5.PCA
immune.combined.sct <- RunPCA(immune.combined.sct, 
                              features = VariableFeatures(object = immune.combined.sct))
ElbowPlot(immune.combined.sct)

#### 6.CLUSTER CELLS
immune.combined.sct <- FindNeighbors(immune.combined.sct, dims = 1:20)
immune.combined.sct <- FindClusters(immune.combined.sct, 
                                    resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 5))           # 0.05 was removed

#### 7.UMAP
immune.combined.sct <- RunUMAP(immune.combined.sct, dims = 1:20)

#### 8. Cell annotation function 
analyze_object <- function(object, ident) {
  Idents(object) <- object[[ident]]
  if (length(levels(object)) == 1) {
    return(message(sprintf('%s: since only one cluster was identified, markers can not be found', ident)))
  }
  out_dir <- paste0('/mnt/tank/scratch/rasmirnov/code/Time_benchmarking/Datasets/processed/',
                  args$obj, '/', 'markers/', ident)
  dir.create(out_dir, recursive = T)
  whole.markers <- FindAllMarkers(object = object,
                                  assay='RNA',
                                  only.pos = TRUE,
                                  min.pct = 0.10,
                                  test.use = 'MAST')
  write.table(whole.markers, 
              paste(out_dir, "markers.tsv", sep = '/'), 
              sep="\t",
              quote=F, 
              row.names=F)
}

idents <- grep('_res', colnames(immune.combined.sct@meta.data), value = T)
sapply(c(idents, 'seurat_clusters'),       # or seurat_annotations
       function(ident) analyze_object(object = immune.combined.sct, ident = ident))


# name of an object 
# out_dir2 <- paste0('/mnt/tank/scratch/rasmirnov/code/Time_benchmarking/Datasets/processed/', args$obj)
# dir.create(out_dir2, recursive = T)

save.image(paste0('/mnt/tank/scratch/rasmirnov/code/Time_benchmarking/Datasets/processed/',
                  args$obj, '/', 
                  args$obj, '_tutorial.RData'))





