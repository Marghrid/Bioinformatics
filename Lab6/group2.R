library(randomcoloR)
library(ramify)

data = read.csv('TCGA_BRCA_Gene_ReadCounts.txt', sep='\t', header=TRUE)

sample_n = 38
increment = 8000

list = c()
for (i in 1:sample_n) {
  list = append(list, c(length(unique(sample(rownames(data), size=i*increment, prob=data$TCGA.A1.A0SB.01, replace=TRUE))) / nrow(data)))
}

plot(1:length(list), unlist(list), type="l", col='blue', ylim=c(0,0.7))

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
  lines(1:length(list), unlist(list), col=randomColor())
}
