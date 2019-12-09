library(randomcoloR)
library(ramify)

data = read.csv('TCGA_BRCA_Gene_ReadCounts.txt', sep='\t', header=TRUE)
stds = apply(data, 2, sd)
min = names(which.min(stds))
max = names(which.max(stds))

people = sample(colnames(data), size=5)
sample_n = 20000

sample_vector = rep(rownames(data), data[[min]])
total = length(unique(sample_vector))
list = c()

for (i in 1:sample_n) {
  list = append(list, c(length(unique(sample(sample_vector, size=i))) / total))
}

plot(1:length(list), unlist(list), type="l", col='blue')

sample_vector = rep(rownames(data), data[[max]])
total = length(unique(sample_vector))
list = c()

for (i in 1:sample_n) {
  list = append(list, c(length(unique(sample(sample_vector, size=i))) / total))
}

lines(1:length(list), unlist(list), col='red')
 
