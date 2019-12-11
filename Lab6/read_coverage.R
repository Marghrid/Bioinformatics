library(edgeR)


data = read.csv('outputs/counts_comparison.csv', sep=',', header=TRUE)
data
plot(data$est, data$target, xlab="Estimated Counts", ylab="Provided Counts")
