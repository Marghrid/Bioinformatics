library(ramify)
library(scales)

data = read.csv('TCGA_BRCA_Gene_ReadCounts.txt', sep='\t', header=TRUE, row.names = 1)
groups = as.numeric(substr(colnames(data), nchar(colnames(data))-1, nchar(colnames(data))-1))
colors <- c('darkcyan', 'red')

sample_n = 38
increment = 8000

list = c()
for (i in 1:sample_n) {
  list = append(list, c(length(unique(sample(rownames(data), size=i*increment, prob=data$TCGA.A1.A0SB.01, replace=TRUE))) / nrow(data)))
}

plot(1:length(list) * increment, unlist(list), type="l", col=alpha(colors[groups[which(colnames(data) == "TCGA.A1.A0SB.01")] + 1], 0.2), ylim=c(0,0.7), xlab="Number of Samples", ylab="Gene Coverage")

first = TRUE
for (person in colnames(data)) {
  if (first) {
    first = FALSE
    next
  }
  list = c()
  for (i in 1:sample_n) {
    list = append(list, c(length(unique(sample(rownames(data), size=i*increment, prob=data[[person]], replace=TRUE))) / nrow(data)))
  }
  lines(1:length(list) * increment, unlist(list), col=alpha(colors[groups[which(colnames(data) == person)] + 1], 0.2), xlab="Number of Samples", ylab="Gene Coverage")
}
