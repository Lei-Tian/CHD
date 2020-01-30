#! /usr/bin/env Rscript

library(DESeq2)
library(genefilter)
library(ggplot2)

source("/srv/gsfs0/projects/wu/tianlei/projects/CHD/scripts/SV/my_PCA.R")

exLim <- function(values, ratio) {
  lim.min <- min(values)
  lim.max <- max(values)
  lim.mid <- (lim.max-lim.min)/2 + lim.min
  lim.min.new <- lim.mid - (lim.max-lim.min)*(1+ratio)/2
  lim.max.new <- lim.mid + (lim.max-lim.min)*(1+ratio)/2
  lims <- list(min=lim.min.new, max=lim.max.new)
  lims
}

FC.cutoff <- 2.0

tag <- "SV"
res.path <- "./"
setwd(res.path)
reads.cnt.tbl <- read.table("A.out.txt",
                            stringsAsFactors=FALSE,
                            header=TRUE)

#families <- c("F1", "F1", "F1", "F1", "F1", "F1",
#  "F2", "F2", "F2", "F2", "F2", "F2")
phenotypes <- c("CTR", "SV","SV","SV","SV","SV","CTR","CTR","CTR","CTR")
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
write.table(vsd.exp, file="vsd_exp.txt", sep="\t", quote=FALSE)

cds <- DESeq(cds)
print(resultsNames(cds))

#ntop=500
#rv <- rowVars(vsd.exp)
#select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
#pca <- prcomp(t(vsd.exp[select,]))
pca <- prcomp(t(vsd.exp))
write.table(cbind(pca$x,Group), file="pca.txt", sep="\t", quote=FALSE)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
write.table(percentVar, file="percentVar.txt", sep="\t", quote=FALSE)

library(Rtsne) # Load package
set.seed(5)
tsne_out <- Rtsne(unique(t(vsd.exp)), perplexity=1, theta=0.0, dims=3) # Run TSNE
tsne.out.y <- as.data.frame(tsne_out$Y)
row.names(tsne.out.y)=rownames(colData)
names(tsne.out.y) <- c("tSNE1", "tSNE2", "tSNE3")
write.table(tsne.out.y, file="tSNE.txt", sep="\t", quote=FALSE)

print(colData(cds))
cds.1 <- DESeq(cds, test="LRT", reduced=~1)
print(resultsNames(cds.1))
res.LRT.1 <- results(cds.1)
write.table(res.LRT.1, file="res_LRT.txt", sep="\t", quote=FALSE)

#cds.2 <- DESeq(cds, test="Wald")#, reduced=~families)
#print(resultsNames(cds.2))
#res.LRT.2 <- results(cds.2)
#write.table(res.LRT.2, file="res_Wald_phenotypes.txt", sep="\t", quote=FALSE)

