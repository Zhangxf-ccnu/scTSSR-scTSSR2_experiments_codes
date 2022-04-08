## DE ------------------------------------
# edgeR--
library(edgeR)

my.edgeR <- function(count, group){
  
  # round data if there are non-integer columns
  intcol<-vector("logical")
  for(i in 1:ncol(count)){
    intcol<-c(intcol,is.integer(count[,i]))
  }
  if (!all(TRUE==intcol)) {
    warning("WARNING! Non-integer expression levels. Data rounded")
    count<-round(count)
  }
  
  # remove out zero count genes in all samples
  #count.sel = count[rowSums(count)>0,]
  count.sel = count
  # creat main edgeR object
  d = DGEList(count.sel, group=group)
  
  # do default normalisation for edgeR
  d = calcNormFactors(d)
  # estimate the common and tagwise dispersion
  d = estimateCommonDisp(d)
  d = estimateTagwiseDisp(d)
  
  # determine differentially expressed genes (using exact test)
  dest = exactTest(d)
  
  dest.sort = topTags(dest, n = dim(dest)[1])
  
  # need to extract table data from de.com object,
  # then select only required columns, plus calculate adjust p-value
  res.edgeR <- data.frame(rownames(dest$table),round(dest$table$logFC,digits=2),signif(p.adjust(dest$table$PValue,method='BH'),digits=3),signif(dest$table$PValue,digits=3))
  
  res.edgeR.sort <- data.frame(rownames(dest.sort$table),round(dest.sort$table$logFC,digits=2),signif(p.adjust(dest.sort$table$PValue,method='BH'),digits=3),signif(dest.sort$table$PValue,digits=3))
  
  
  # name
  names(res.edgeR) <- c('GeneID','log2FC','p.adjust','pvalue')
  names(res.edgeR.sort) <- c('GeneID','log2FC','p.adjust','pvalue')
  
  out = list(res.edgeR=res.edgeR, res.edgeR.sort = res.edgeR.sort)
  out
}




log_normalization = function(x, percent=0.05, preprocess.only=FALSE){
  
  if (preprocess.only){
    n <- dim(x)[2]
    gene.exprs.count <- rowSums(x != 0)
    x = x[gene.exprs.count > n * percent, ]
    return(x)
  }
  else{
    n <- dim(x)[2]
    gene.exprs.count <- rowSums(x != 0)
    x = x[gene.exprs.count > n * percent, ]
    sf = colSums(x)/median(colSums(x))
    return(log(sweep(x, 2, sf, '/')+1))
    
  }
  
}



##### load datasets downloaded from GEO website
GSE75748_bulk_cell_type_ec <- read.csv("GSE75748_bulk_cell_type_ec.csv", row.names = 1)
GSE75748_sc_cell_type_ec <- read.csv("GSE75748_sc_cell_type_ec.csv", row.names = 1)


bulk_data <- as.matrix(GSE75748_bulk_cell_type_ec[,c(1:4, 8:9)])
bulk_group <- as.factor(c(rep("H1_ESC",4),rep("DEC",2)))
chu_sc_cell_type_se <- as.matrix(GSE75748_sc_cell_type_ec[,c(1:212,375:512)])
sc_group <- as.factor(c(rep("H1_ESC",212),rep("DEC",138)))

# preprocessing and normalization
sc_data <- log_normalization(chu_sc_cell_type_se, percent=0.05, preprocess.only=TRUE)

## select 2,000 highly variable genes
x.seurat <- CreateSeuratObject(sc_data)
x.seurat <- NormalizeData(x.seurat)
x.seurat <- ScaleData(x.seurat)
x.seurat <- FindVariableFeatures(x.seurat, verbose = FALSE)

filter_ID <- x.seurat@assays$RNA@var.features

bulk_data_filter <- bulk_data[filter_ID,]
sc_data_filter <- sc_data[filter_ID,]


## load preprocessed counts with selected 2,000 HVGs
bulk_data_filter <- readRDS('H1_DEC_bulk.rds')
sc_data_filter <- readRDS('H1_DEC_sc.rds')

### use edgeR to identify DE genes 
bulk_deg <- my.edgeR(bulk_data_filter, bulk_group)$res.edgeR.sort
raw_deg <- my.edgeR(sc_data_filter, sc_group)$res.edgeR.sort


### compute AUC 
bulk_pvalue <- bulk_deg$p.adjust
raw_pvalue <- raw_deg$p.adjust



n<-200  #### 400, 600, 800, 1000


bulk_tpr <- c()
raw_tpr <- c()


################

bulk_fpr <- c()
raw_fpr <- c()



degene <- bulk_deg$GeneID[1:n]
nodegene <- bulk_deg$GeneID[-c(1:n)]




for (i in 1:2000){
  
 
  raw_tpr[i]<- length(intersect(raw_deg[raw_deg$p.adjust<=raw_pvalue[i],]$GeneID,
                                degene))/length(degene)


  raw_fpr[i] <- length(intersect(raw_deg[raw_deg$p.adjust<=raw_pvalue[i],]$GeneID,
                                 nodegene))/length(nodegene)
 
}



raw_auc <- 0


for (i in 2:2000){
 
  raw_auc <- raw_auc + raw_tpr[i]*(raw_fpr[i]-raw_fpr[i-1])
  
}




#### compute Spearman correlation coefficient 
raw_padj <- raw_deg$p.adjust[match(bulk_deg$GeneID, raw_deg$GeneID)]


bulk_padj <- bulk_deg$p.adjust

cor.spearman <- cor(bulk_padj, raw_padj, method = c('spearman'))






