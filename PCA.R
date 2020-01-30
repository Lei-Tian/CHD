#! /usr/bin/env Rscript

library(DESeq2)
library(genefilter)
library(ggplot2)

tag <- "all"
res.path <- "./"
setwd(res.path)
reads.cnt.tbl <- read.table("A.out.txt",
                            stringsAsFactors=FALSE,
                            header=TRUE)

#families <- c("F1", "F1", "F1", "F1", "F1", "F1",
#  "F2", "F2", "F2", "F2", "F2", "F2")
#h339_WT	h341_TOF	h342_TOF	h343_TOF	h344_SRV	h345_SRV	h347_SLV	h348_SLV	h349_SLV	h384_WT	h385_WT	TK413_WT	TK416_WT	TK418_TOF	TK420_TOF
phenotypes <- c("CTR", "HD", "HD", "HD", "HD","HD","HD","HD","HD","CTR","CTR","CTR","CTR","HD","HD")
Group=phenotypes
#colData <- as.data.frame(cbind(families=families,
#            phenotypes=phenotypes))
colData <- as.data.frame(cbind(phenotypes=phenotypes))

rownames(reads.cnt.tbl) <- reads.cnt.tbl[ ,1]
reads.cnt.tbl <- reads.cnt.tbl[ , -1]

maxCounts=apply(reads.cnt.tbl,1,max)
reads.cnt.tbl=reads.cnt.tbl[which(maxCounts>=2),]

rownames(colData) <- names(reads.cnt.tbl)
print(colData)

cds <- DESeqDataSetFromMatrix(countData = reads.cnt.tbl,
                              colData = colData,
                              design = ~ phenotypes+1)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
vsd = varianceStabilizingTransformation(cds)
vsd.exp <- assay(vsd)
write.table(vsd.exp, file="vsd_exp1.txt", sep="\t", quote=FALSE)

cds <- DESeq(cds)
print(resultsNames(cds))

pca <- prcomp(t(vsd.exp))
write.table(cbind(pca$x,Group), file="pca1.txt", sep="\t", quote=FALSE)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
write.table(percentVar, file="percentVar1.txt", sep="\t", quote=FALSE)

