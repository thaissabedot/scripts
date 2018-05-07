##RNA sequencing analysis - LGG primary x recurrent
setwd("/dados/camila/Recurrent.gliomas/")
##Create DGeList object
library(edgeR)
rnaseq.lgg.edger <- DGEList(counts = rnaseq.lgg[,-1], genes = rnaseq.lgg$gene_id)  ##20530 probes and 28 columns (samples)
setwd("/dados/camila/Recurrent.gliomas/RNASeqV2/Plots/LGG")
png(filename = "LGG_Librarysize(millions).png", res = 800, height = 4000, width = 4000)
barplot(rnaseq.lgg.edger$samples$lib.size*1e-6,names=1:28,ylab="Library size (millions)")
##Update gene symbols and remove outdated gene_ids
setwd("/dados/camila/Recurrent.gliomas/")
library(org.Hs.eg.db)
egSYMBOL <- toTable(org.Hs.egSYMBOL)
idfound.lgg <- rnaseq.lgg.edger$genes$genes %in% egSYMBOL$gene_id  ##length = 20530 probes
rnaseq.lgg.edger.clean <- rnaseq.lgg.edger[idfound.lgg,]  ##dim = 20440 probes
m.lgg <- match(rnaseq.lgg.edger.clean$genes$genes,egSYMBOL$gene_id)
rnaseq.lgg.edger.clean$genes$Symbol.Updated <- egSYMBOL$symbol[m.lgg]
##Filtering reads
lgg.o <- order(rowSums(rnaseq.lgg.edger.clean$counts),decreasing = TRUE)
rnaseq.lgg.edger.clean.o <- rnaseq.lgg.edger.clean[lgg.o,]
lgg.d <- duplicated(rnaseq.lgg.edger.clean.o$genes$Symbol.Updated)
rnaseq.lgg.edger.clean.o.d <- rnaseq.lgg.edger.clean.o[!lgg.d,]
nrow(rnaseq.lgg.edger.clean.o.d)  ##20440 probes
##Count lib size and getting data ready for normalization
rnaseq.lgg.edger.clean.o.d$samples$lib.size <- colSums(rnaseq.lgg.edger.clean.o.d$counts)
rownames(rnaseq.lgg.edger.clean.o.d$counts) <- rownames(rnaseq.lgg.edger.clean.o.d$genes) <- rnaseq.lgg.edger.clean.o.d$genes$Symbol.Updated
##Normalization (Stefano's quantile normalization)
install.packages("/dados/camila/Recurrent.gliomas/RNASeqV2/steFunctions_2014.10.12.tar.gz",repos = NULL,type = "source")
require(steFuctions)
require(EDASeq)
rawCounts.lgg <- rnaseq.lgg.edger.clean.o.d$counts
rawCounts.lgg <- floor(rawCounts.lgg)
tmp.lgg <- newSeqExpressionSet(rawCounts.lgg,featureData = rnaseq.lgg.edger.clean.o.d$genes)
tmp.lgg <- betweenLaneNormalization(tmp.lgg, which = "upper", offset = TRUE)
normCounts.lgg <- log(rawCounts.lgg + .1) + offst(tmp.lgg)
normCounts.lgg <- floor(exp(normCounts.lgg) - .1)
library(steFunctions)
tmp.lgg <- t(quantileNormalization(t(normCounts.lgg)))
normCounts.lgg <- floor(tmp.lgg)
rnaseq.lgg.edger.normCounts <- DGEList(counts = normCounts.lgg,genes = rnaseq.lgg.edger.clean.o.d$genes)
setwd("/dados/camila/Recurrent.gliomas/RNASeqV2/Plots/LGG")
png(filename = "LGG_normCounts_Librarysize(millions).png", res = 800, height = 4000, width = 4000)
barplot(rnaseq.lgg.edger.normCounts$samples$lib.size*1e-6,names=1:28,ylab="Library size (millions)")
png(filename = "MDSplot_LGG.png", res = 800, height = 4000, width = 4000)
plotMDS(rnaseq.lgg.edger.normCounts,pch = 20)
##Design matrix for paired samples
Patient.lgg <- factor(c("TCGA-DH-A669-01A-12R-A31N-07","TCGA-DH-A669-02A-11R-A31N-07","TCGA-TM-A7CF-01A-11R-A32Q-07","TCGA-TM-A7CF-02A-11R-A32Q-07","TCGA-TQ-A7RK-01A-11R-A33Z-07","TCGA-TQ-A7RK-02A-11R-A36H-07","TCGA-TQ-A7RV-01A-21R-A34F-07","TCGA-TQ-A7RV-02A-11R-A36H-07","TCGA-TQ-A8XE-01A-11R-A36H-07","TCGA-TQ-A8XE-02A-11R-A36H-07","TCGA-DU-5870-01A-11R-1708-07","TCGA-DU-5870-02A-12R-A36H-07","TCGA-DU-5872-01A-11R-1708-07","TCGA-DU-5872-02A-21R-A36H-07","TCGA-DU-6397-01A-11R-1708-07","TCGA-DU-6397-02A-12R-A36H-07","TCGA-DU-6404-01A-11R-1708-07","TCGA-DU-6404-02A-21R-A36H-07","TCGA-DU-6407-01A-13R-1708-07","TCGA-DU-6407-02A-12R-A36H-07","TCGA-FG-5963-01A-11R-1708-07","TCGA-FG-5963-02A-12R-A29R-07","TCGA-FG-5965-01B-11R-1896-07","TCGA-FG-5965-02A-11R-A29R-07","TCGA-DU-7304-01A-12R-2090-07","TCGA-DU-7304-02A-12R-A36H-07","TCGA-FG-A4MT-01A-11R-A26U-07","TCGA-FG-A4MT-02A-11R-A29R-07"))
Tumor.type.lgg <- factor(c("P","R","P","R","P","R","P","R","P","R","P","R","P","R","P","R","P","R","P","R","P","R","P","R","P","R","P","R"))
data.frame(Sample=colnames(rnaseq.lgg.edger.normCounts),Patient.lgg,Tumor.type.lgg)
design.lgg <- model.matrix(~ Patient.lgg + Tumor.type.lgg)
rownames(design.lgg) <- colnames(rnaseq.lgg.edger.normCounts)
rnaseq.lgg.edger.normCounts <- estimateGLMCommonDisp(rnaseq.lgg.edger.normCounts,design.lgg[,c(1,29)],verbose = TRUE)  ##Disp = 0.53972, BCV = 0.7347 
rnaseq.lgg.edger.normCounts <- estimateGLMTrendedDisp(rnaseq.lgg.edger.normCounts,design.lgg[,c(1,29)])
rnaseq.lgg.edger.normCounts <- estimateGLMTagwiseDisp(rnaseq.lgg.edger.normCounts,design.lgg[,c(1,29)])
setwd("/dados/camila/Recurrent.gliomas/RNASeqV2/Plots/LGG")
png(filename = "BCVplot_LGG.png", res = 800, height = 4000, width = 4000)
fit.lgg <- glmFit(rnaseq.lgg.edger.normCounts,design.lgg[,c(1,29)])
lrt.lgg <- glmLRT(fit.lgg)
topTags(lrt.lgg)
##Limma linear modeling
library(limma)
v.norm.lgg <- voom(rnaseq.lgg.edger.normCounts,design.lgg[,c(1,29)])
fit.limma.norm.lgg <- lmFit(v.norm.lgg,design.lgg[,c(1,29)])
fit.limma.norm.lgg <- eBayes(fit.limma.norm.lgg)
##Visualize results
y.lgg <- cpm(rnaseq.lgg.edger.normCounts,prior.count = 2,log = TRUE)
y.norm.lgg <- v.norm.lgg$E
##Plots
up.all.lgg <- y.lgg[lrt.lgg$table$logFC > 1.2 & lrt.lgg$table$PValue < 0.05,]  ##581 probes
up.all.lgg.order <- hclust(dist(up.all.lgg))  ##Cluster method: complete; distance: euclidean; number of objects: 581 
up.all.lgg <- up.all.lgg[rev(up.all.lgg.order$order),]
down.all.lgg <- y.lgg[lrt.lgg$table$logFC < -1.2 & lrt.lgg$table$PValue < 0.05,]  ##218 probes
down.all.lgg.order <- hclust(dist(down.all.lgg))  ##Cluster method: complete; distance: euclidean; number of objects: 218 
down.all.lgg <- down.all.lgg[down.all.lgg.order$order,]
notsig.all.lgg <- y.lgg[lrt.lgg$table$PValue > 0.05,]  ##18221 probes
rnaseqlog2.order.all.lgg <- rbind(down.all.lgg,up.all.lgg)  ##799 probes
primary.lgg <- rnaseqlog2.order.all.lgg[,c(1,3,5,7,9,11,13,15,17,19,21,23,25,27)]  ##799 probes and 14 columns (samples)
recurrent.lgg <- rnaseqlog2.order.all.lgg[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28)]  ##799 probes and 14 columns (samples)
primary.lgg.o <- hclust(dist(t(primary.lgg)))  ##Cluster method: complete; distance: euclidean; number of objects: 14
recurrent.lgg.o <- hclust(dist(t(recurrent.lgg)))  ##Cluster method: complete; distance: euclidean; number of objects: 14
topTable(fit.limma.norm.lgg,sort.by = "logFC",num = Inf,adjust.method = "BH")
lrt.all.g.lgg <- topTags(lrt.lgg,n = dim(lrt.lgg)[1])$table
range(lrt.all.g.lgg[lrt.all.g.lgg$FDR > 0,"FDR"])  # To determine the next lowest pvalue after 0. In this case this p-value is 0.0005464831
lrt.all.g.lgg[lrt.all.g.lgg$FDR == 0,"FDR"] <- 0.0005464831
lrt.all.g.lgg$type <- "All"
lrt.all.g.lgg$cutoffs <- "Notsig"
lrt.all.g.lgg[lrt.all.g.lgg$FDR < 0.05 & lrt.all.g.lgg$logFC > 1.2,"cutoffs"] <- "sigUP"
lrt.all.g.lgg[lrt.all.g.lgg$FDR < 0.05 & lrt.all.g.lgg$logFC < -1.2,"cutoffs"] <- "sigDOWN"
tmp2 <- lrt.all.g.lgg  ##For maintaining this object in ls()
library(ggplot2)
require(plyr) ##For subsetting in ggplot2
#Volcano (all primary x recurrent)
setwd("/dados/camila/Recurrent.gliomas/RNASeqV2/Plots/LGG/")
png(filename = "Volcanoplot_LGG.png", res = 800, height = 4000, width = 4000)
ggplot(lrt.all.g.lgg,aes(x = logFC,y = -log10(FDR),color = cutoffs)) +
  geom_point() +
  scale_color_manual(values = c("grey28", "royalblue4", "maroon4"),labels = c("Not Sig","Sig Down","Sig Up"))
rnaseqlog2.order.all.lgg <- cbind(primary.lgg[,primary.lgg.o$order],recurrent.lgg[,recurrent.lgg.o$order])  ##All primary x recurrent
rnaseqlog2.order.all.lgg <- rnaseqlog2.order.all.lgg[,c(13,6,7,10,11,1,3,2,12,8,14,4,5,9,24,19,20,21,27,28,15,25,16,26,17,23,18,22)]
#Heatmap (all primary x recurrent)
setwd("/dados/camila/Recurrent.gliomas/RNASeqV2/Plots/LGG/")
png(filename = "Heatmap_newcolor_LGG_primaryvsrecurrent.png", res = 800, height = 4000, width = 4000)
library(gplots)
library(RColorBrewer)
col <- colorRampPalette(c("mediumblue","lightgoldenrod1"))(5) 
heatmap.2(as.matrix(rnaseqlog2.order.all.lgg[rev(rownames(rnaseqlog2.order.all.lgg)),]),trace="none",col=col,labCol=substr(colnames(rnaseqlog2.order.all.lgg),9,16),cexCol=0.6,Colv=NA,Rowv=NA)