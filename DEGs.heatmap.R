library(dendsort)
library("gplots")
library(RColorBrewer)
library(genefilter)

setwd("D:/SharedFolder/project/CHD")
RdBu = rev(brewer.pal(11, name="RdBu"))
a=read.table("vsd_exp.txt")
#ntop=100
#rv <- rowVars(a)
#select=order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
TOF=read.table("TOF/DEGs.txt")
SV=read.table("SV/DEGs.txt")
DEGs=unique(c(as.vector(TOF[,1]),as.vector(SV[,1])))

b=a[DEGs,]


pdf("DEGs.heatmap.pdf",width=6,height=9)
#for(dm in c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")){
#  for(hm in c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid")){
dm="euclidean"
hm="average"
    heatmap.2(as.matrix(b),dendrogram ="both",scale="row",distfun=function(x) dist(x,method=dm),
          hclustfun=function(x) hclust(x,method=hm), density.info="none",labRow=F,col=RdBu,
          trace="none",ColSideColors=c(rep("black",5),rep("darkcyan",5),rep("orange",5)),lhei=c(1,6))
#  }
#}
dev.off()
