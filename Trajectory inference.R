library(TSCAN)
library(monocle)


my.TSCAN =  function(count, cellLabels){
  colnames(count) = c(1:ncol(count))
  procdata <- TSCAN::preprocess(count)
  lpsmclust <- TSCAN::exprmclust(procdata)
  
  lpsorder <- TSCAN::TSCANorder(lpsmclust, orderonly=F)
  
  Pseudotime = lpsorder$Pseudotime[match(colnames(count),lpsorder$sample_name)]
  cor.kendall = cor(Pseudotime, as.numeric(cellLabels), method = "kendall", use = "complete.obs")
  subpopulation <- data.frame(cell = colnames(count), sub = as.numeric(cellLabels)-1)
  POS <- orderscore(subpopulation, lpsorder)[1]
  out = list(cor.kendall=cor.kendall, POS=POS)
  out
}




my.monocle <- function(count, cellLabels){
  colnames(count) <- 1:ncol(count)
  geneNames <- rownames(count)
  rownames(count) <- 1:nrow(count)
  
  
  pd <- data.frame(timepoint = cellLabels)
  pd <- new("AnnotatedDataFrame", data=pd)
  fd <- data.frame(gene_short_name = geneNames)
  fd <- new("AnnotatedDataFrame", data=fd)
  
  dCellData <- newCellDataSet(count, phenoData = pd, featureData = fd, expressionFamily = uninormal())
  
  dCellData <- detectGenes(dCellData , min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(dCellData),
                                      num_cells_expressed >= 50))
  
  diff_test_res <- differentialGeneTest(dCellData[expressed_genes,],
                                        fullModelFormulaStr = "~timepoint",
                                        cores = 3)
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
  
  dCellData <- setOrderingFilter(dCellData, ordering_genes)
  
  dCellData <- reduceDimension(dCellData, max_components = 2,
                               method = 'DDRTree', norm_method='none')  
  dCellData <- orderCells(dCellData)
  
  cor.kendall = cor(dCellData@phenoData@data$Pseudotime, as.numeric(dCellData@phenoData@data$timepoint), 
                    method = "kendall", use = "complete.obs")
  
  lpsorder2 = data.frame(sample_name = colnames(count), State= dCellData@phenoData@data$State, 
                         Pseudotime = dCellData@phenoData@data$Pseudotime, rank = rank(dCellData@phenoData@data$Pseudotime))
  
  lpsorder_rank = dplyr::arrange(lpsorder2, rank)
  
  lpsorder_rank$Pseudotime = lpsorder_rank$rank
  
  lpsorder_rank = lpsorder_rank[-4]

  lpsorder_rank[1] <- lapply(lpsorder_rank[1], as.character)
  
  subpopulation <- data.frame(cell = colnames(count), sub = as.numeric(cellLabels)-1)
  
  POS <- TSCAN::orderscore(subpopulation, lpsorder_rank)[1]
 
  print(cor.kendall, POS)
  
}



### run TSCAN on Deng data 
deng <- readRDS('Deng.rds')
deng_cellLabels <- factor(colnames(deng),
                         levels=c('zygote', 'early 2-cell', 'mid 2-cell', 'late 2-cell',
                                  '4-cell', '8-cell', '16-cell', 'early blastocyst',
                                  'mid blastocyst', 'late blastocyst'))

deng.TSCAN <- my.TSCAN(deng, deng_cellLabels)

###  run Monocle2 on Deng data 
deng.monocle <- my.monocle(deng, deng_cellLabels)




