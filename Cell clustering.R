library(SingleCellExperiment)
library(SC3)
library(Seurat)
library(mclust)

################# SC3 ##########################
my.SC3 = function(count, k){
  sce = SingleCellExperiment(
    assays = list(
      counts = as.matrix(count),
      logcounts = log2(as.matrix(count) + 1)
    )
  )
  rowData(sce)$feature_symbol <- rownames(sce)
  sce = sc3(sce, ks = k, gene_filter=FALSE)
  sce
}

## apply SC3 on sc_celseq2_5cl_p1 dataset
count <- readRDS('sc_celseq2_5cl_p1.rds')

label <- sub('.*:','',colnames(count))

k <- length(unique(label))

result <- my.SC3(count=count, k=k)

ARI <- adjustedRandIndex(label, result@colData@listData$sc3_5_clusters)


############ Seurat ########################
## apply Seurat on sc_celseq2_5cl_p1 dataset
x.seurat <- CreateSeuratObject(count)
x.seurat <- NormalizeData(x.seurat)
x.seurat <- ScaleData(x.seurat)
x.seurat <- FindVariableFeatures(x.seurat, verbose = FALSE)

x.seurat <- RunPCA(x.seurat, features = VariableFeatures(object = x.seurat))  
x.seurat <- JackStraw(x.seurat, num.replicate = 100)
x.seurat <- ScoreJackStraw(x.seurat, dims = 1:20)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
x.seurat <- FindNeighbors(x.seurat, dims = 1:10)
x.seurat <- FindClusters(x.seurat, resolution = 0.5)

label <- sub('.*:','',colnames(count))

ARI <- adjustedRandIndex(as.factor(label), x.seurat$seurat_clusters)



