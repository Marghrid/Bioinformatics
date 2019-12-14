library(edgeR)
library(stringr)
library(fitdistrplus)
library(ggplot2)

ReadCounts = read.delim('TCGA_BRCA_Gene_ReadCounts.txt', header=TRUE, row.names = 1)
annotations = read.delim('TCGA_BRCA_ClinicalAnnotation.txt', header=TRUE, row.names=1)
groups = as.numeric(substr(colnames(ReadCounts), nchar(colnames(ReadCounts))-1, nchar(colnames(ReadCounts))-1))

rownames(annotations) = str_replace(str_replace(rownames(annotations), '-', '.'), '-', '.')

# Tissue <- factor(groups)
Patient <- factor(substr(colnames(ReadCounts), 0, 12))
# Gender <- annotations$Gender[Patient]
# AgeDiag <- factor(annotations$Age.at.diagnosis..years.[Patient])
# Menopause <- annotations$Menopausal.status[Patient]
PAM50 <- addNA(annotations$PAM50[Patient])
# Vital <- annotations$Vital.status[Patient]
# Procedure <- annotations$Surgical.procedure[Patient]
# Histology <- annotations$Histology[Patient]
# Anatomy <- annotations$Anatomic.subdividion[Patient]
# Ethnicity <- annotations$Ethnicity[Patient]
# HER2 <- addNA(annotations$HER2[Patient])
# StageM <- annotations$Stage.M[Patient]
# StageT <- annotations$Stage.T[Patient]
# StageN <- annotations$Stage.N[Patient]
# Year <- annotations$Year.of.diagnosis[Patient]

data.frame(Sample=colnames(ReadCounts),Patient,Tissue,PAM50)

design <- model.matrix(~ PAM50)

dge <- DGEList(counts=data.matrix(ReadCounts), group = Tissue)

keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

v <- voom(dge,design,plot=TRUE)

linearfit = lmFit(v$E,design.matrix)

eBfit = eBayes(linearfit)

volcanoplot(eBfit,coef=2,style="B-statistic")

topTable(eBfit,coef=2)

