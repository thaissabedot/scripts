hmec.r <- dba.report(hmec.a)
seqlevels(hmec.r) <- paste0("chr",seqlevels(hmec.r))
a <- hmec.r[mcols(hmec.r)[,4] < 0,]

aux <- findOverlaps(makeGRangesFromDataFrame(genes),a)

aux <- genes[unique(queryHits(aux)),]

ggplot(aux,aes(seqid,fill=RNAseq8.summary)) + geom_bar() + scale_fill_manual(values=c("gray48", "skyblue2", "orange", "black"),name="Legend") + scale_x_discrete(drop=FALSE)



a <- AnnotationTrack(start=31135000, end=31289000,chromosome="chrX", genome="hg19", name="AnnotationTrack")
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg19", chromosome = "chrX")

cpgIslands <- UcscTrack(genome = "hg19", chromosome = "chrX",track = "cpgIslandExt", from = 31135000, to = 31289000,trackType = "AnnotationTrack", start = "chromStart",end = "chromEnd", id = "name", shape = "box",fill = "#006400", name = "CpG Islands")


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
txTranscripts_v1 <- GeneRegionTrack(txdb_hg19, genome="hg19", chromosome="chrX", showId=TRUE, geneSymbol=TRUE, name="UCSC",from = 31135000, to = 31289000)

knownGenes <- UcscTrack(genome = "hg19", chromosome = "chrX",
track = "knownGene", from = 31135000, to = 31289000, trackType = "GeneRegionTrack",
rstarts = "exonStarts", rends = "exonEnds", gene = "name",
symbol = "name", transcript = "name", strand = "strand",
fill = "#8282d2", name = "UCSC Genes")

load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/all450.rda")
all450$Chromosome <- paste0("chr",all450$Chromosome)
colnames(all450) <- substr(colnames(all450),1,12)
load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/LGG-GBM.merged.data.FB20150930v3.Rdata")

aux.GR <- GRanges(seqnames = "chrX",ranges = IRanges(start = 31135000, end=31289000))
overlaps <- findOverlaps(makeGRangesFromDataFrame(all450,start.field="Genomic_Coor",end.field="Genomic_Coor"), aux.GR)
all450 <- all450[queryHits(overlaps),]


aux <- as.character(subset(pd, cartoon %in% c("GCIMP-low"))$case.id)
low <- all450[,colnames(all450) %in% aux]
#low <- apply(low,1,mean,na.rm=TRUE)
aux <- as.character(subset(pd, cartoon %in% c("GCIMP-high"))$case.id)
high <- all450[,colnames(all450) %in% aux]
#high <- apply(high,1,mean,na.rm=TRUE)
aux <- as.character(subset(pd, cartoon %in% c("IDHmut-codel"))$case.id)
codel <- all450[,colnames(all450) %in% aux]
#codel <- apply(codel,1,mean,na.rm=TRUE)
aux <- as.character(subset(pd, cartoon %in% c("Classic-like"))$case.id)
class <- all450[,colnames(all450) %in% aux]
#class <- apply(class,1,mean,na.rm=TRUE)
aux <- as.character(subset(pd, cartoon %in% c("Mesenchymal-like"))$case.id)
mes <- all450[,colnames(all450) %in% aux]
#mes <- apply(mes,1,mean,na.rm=TRUE)
aux <- as.character(subset(pd, cartoon %in% c("LGm6-GBM"))$case.id)
lgm6 <- all450[,colnames(all450) %in% aux]
#lgm6 <- apply(lgm6,1,mean,na.rm=TRUE)
aux <- as.character(subset(pd, cartoon %in% c("PA-like"))$case.id)
pa <- all450[,colnames(all450) %in% aux]
#pa <- apply(pa,1,mean,na.rm=TRUE)
all <- cbind(low,high,codel,class,mes,lgm6,pa)
all <- all[-7,]
all <- cbind(glial.cell[rownames(all),],all)
all <- cbind(all450[-7,c(3,4,4)],all)
all <- all[with(all,order(Genomic_Coor)),]

anno <- makeGRangesFromDataFrame(all450[,c(3,4)],start.field="Genomic_Coor",end.field="Genomic_Coor")
genome(anno) <- "hg19"
anno <- AnnotationTrack(anno,name="450k probes",col="black",fill="black")

meth <- makeGRangesFromDataFrame(all,TRUE,start.field="Genomic_Coor",end.field="Genomic_Coor.1")
genome(meth) <- "hg19"
#meth <- DataTrack(meth,name="DNA methylation",type="smooth",groups=c("GCIMP.low",rep("GCIMP.high",2),rep("Non-tumor",2)),col=c("firebrick","darkgreen","grey"))
meth <- DataTrack(meth,name="DNA methylation",type="smooth",groups=c(rep("Non-tumor brain",ncol(glial.cell)),rep("GCIMP-low",ncol(low)),rep("GCIMP-high",ncol(high)),rep("IDHmut-codel",ncol(codel)),rep("Classic-like",ncol(class)),rep("Mesenchymal-like",ncol(mes)),rep("LGm6-GBM",ncol(lgm6)),rep("PA-like",ncol(pa))),col=c("orange","firebrick","darkgreen","purple","darkblue","yellow","black","cyan"))
#meth <- DataTrack(meth,name="DNA methylation",type="smooth",groups=c("GCIMP.low",rep("GCIMP.high",2),rep("Non-tumor",2)))


plotTracks(list(itrack,gtrack,cpgIslands,txTranscripts_v1,anno,meth),legend=TRUE,from=31282400,to=31287500,span=0.15)

plotTracks(list(itrack,gtrack,cpgIslands,txTranscripts_v1,anno,meth),legend=TRUE,from=31284000,to=31289000)

all <- all[with(all,order(Genomic_Coor)),]
library(reshape)

all.m <- all[,-c(1,2,3)]
all.m$probeID <- rownames(all.m)
levels(all.m$probeID) <- rownames(all)
all.m <- melt(all.m,id="probeID")
all.m <- merge(all.m,pd[,c("case.id","cartoon")],by.x="variable",by.y="case.id",all.x=T)
all.m$cartoon <- as.character(all.m$cartoon)
all.m[is.na(all.m$cartoon),"cartoon"] <- "Non-tumor brain"
all.m$cartoon <- as.factor(all.m$cartoon)
all.m$cartoon <- factor(all.m$cartoon, levels = c("Non-tumor brain","GCIMP-low","GCIMP-high","IDHmut-codel","Classic-like","Mesenchymal-like","LGm6-GBM","PA-like"))

ggplot(all.m,aes(cartoon,value,fill=cartoon)) + 
  theme_classic() + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter() +
  facet_wrap( ~ probeID) +
  scale_fill_manual(values=c("darkgray","darkgreen","firebrick","purple","orange","yellow","darkblue","cyan")) 

