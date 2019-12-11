library(preprocessCore)
library(randomcoloR)

data = read.csv('TCGA_BRCA_Gene_ReadCounts.txt', sep='\t', header=TRUE, row.names = 1)

normalized = normalize.quantiles(data.matrix(data))

write.csv(normalized, 'norm.csv')

sums = colSums(data)

data$TCGA.A1.A0SB.01
hist(colSums(data), main=NA, xlab='Number of Reads', ylab='Number of Samples')
hist(colSums(normalized))
