library(edgeR)

ReadCounts = read.csv('TCGA_BRCA_Gene_ReadCounts.txt', sep='\t', header=TRUE, row.names = 1)

dge <- DGEList(counts=data.matrix(ReadCounts))
dge <- calcNormFactors(dge)

hist(dge$samples$norm.factors)
plotMDS(dge, labels = NULL)
