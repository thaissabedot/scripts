load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/GCIMPlow.90probes.Rda")
load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/low_high_450.rda")

a <- cbind(low[as.character(GCIMPlow.probes$probes),-c(1:4)],high[as.character(GCIMPlow.probes$probes),-c(1:4)])
aux <- hclust(dist(a))
p <- a[aux$order,]
heatmap.plus(as.matrix(p),
             Colv=NA,
             Rowv=NA,
             col=jet.colors(75),
             cexCol = 0.5,
             #labRow = NA,
             labCol = NA,
             scale = "none")

anno <- data.frame(probes=rownames(p),status=c(rep("hypo",11),rep("hyper",79)))
GCIMPlow.probes <- anno
save(GCIMPlow.probes,file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/GCIMPlow.90probes.Rda")

load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/GCIMPlow.90probes.Rda")

hg38_450k <- read.table("/media/data1/hg38_450k.bed")
colnames(hg38_450k) <- c("chrom","start","end","probeID")

hg38_450k <- merge(hg38_450k,GCIMPlow.probes,by.x="probeID",by.y="probes")
hg38_450k$chrom <- factor(hg38_450k$chrom)

h3k4me3.r <- dba.report(h3k4me3.a,th=1)
aux <- findOverlaps(makeGRangesFromDataFrame(hg38_450k),h3k4me3.r)

h3k27ac.r <- dba.report(h3k27ac.a,th=1)
aux <- findOverlaps(makeGRangesFromDataFrame(hg38_450k),h3k27ac.r) #36 e 40

ctcf.r <- dba.report(ctcf.a,th=1)
aux <- findOverlaps(makeGRangesFromDataFrame(hg38_450k),ctcf.r) #36 e 40

load("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/DiffBind/DiffBind_MUSIC_hmC.rda")
hmec.r <- dba.report(hmec.a,th=1)
newStyle <- mapSeqlevels(seqlevels(hmec.r),"UCSC")
hmec.r <- renameSeqlevels(hmec.r, newStyle)

aux <- findOverlaps(makeGRangesFromDataFrame(hg38_450k),hmec.r) #36 e 40
