library(edgeR)
library(stringr)

ReadCounts = read.delim('TCGA_BRCA_Gene_ReadCounts.txt', header=TRUE, row.names = 1)
annotations = read.delim('TCGA_BRCA_ClinicalAnnotation.txt', header=TRUE, row.names=1)
groups = as.numeric(substr(colnames(ReadCounts), nchar(colnames(ReadCounts))-1, nchar(colnames(ReadCounts))-1))

rownames(annotations) = str_replace(str_replace(rownames(annotations), '-', '.'), '-', '.')

Tissue <- factor(groups)
Patient <- factor(substr(colnames(ReadCounts), 0, 12))
Gender <- annotations$Gender[Patient]
AgeDiag <- factor(annotations$Age.at.diagnosis..years.[Patient])
Menopause <- annotations$Menopausal.status[Patient]
PAM50 <- annotations$PAM50[Patient]
Vital <- annotations$Vital.status[Patient]
Procedure <- annotations$Surgical.procedure[Patient]
Histology <- annotations$Histology[Patient]
Anatomy <- annotations$Anatomic.subdividion[Patient]

design.matrix <- model.matrix(~Tissue)

dge <- DGEList(counts=data.matrix(ReadCounts), group = groups)

#par(mfrow=c(1,1))
# lcpm <- cpm(dge, log=TRUE)
# boxplot(lcpm, las=2, main="")
# title(main="A. Example: Unnormalised data",ylab="Log-cpm")
# 
dge <- calcNormFactors(dge)
# 
# lcpm <- cpm(dge, log=TRUE)
# boxplot(lcpm, las=2, main="")
# title(main="B. Example: Normalised data",ylab="Log-cpm")
# 
# cutoff <- 1
# drop <- which(apply(cpm(dge), 1, max) < cutoff)
# d <- dge[-drop,] 

# plotMD(cpm(dge, log=TRUE), column=1)
# abline(h=0, col="red", lty=2, lwd=2)

plotMDS(dge, pch=20, col = Anatomy, gene.selection = "common")
legend(-2.6,2.8,unique(Anatomy),col=1:length(Anatomy),pch=20)

v <- voom(dge,design.matrix,plot=TRUE)

linearfit = lmFit(v$E,design.matrix)

eBfit = eBayes(linearfit)

volcanoplot(eBfit,coef=2,style="B-statistic")

topTable(eBfit,coef=2)

PAM50
