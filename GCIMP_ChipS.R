#"TCGA-DU-7010","TCGA-06-0129","TCGA-06-1805","TCGA-06-2570","TCGA-DU-6408","TCGA-DU-6542","TCGA-DU-7304","TCGA-DU-7007","TCGA-DU-7008","TCGA-06-0221"

load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/LGG-GBM.merged.data.FB20150930v3.Rdata") #Metadata

load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/LGG_GBM_nonTumorBrain_RNAseqNorm.rda") #RNAseq

CNV.focal <- read.table("/media/data1/thais.projects/GCIMP-low/CopyNumber_GISTIC/LGG-GBM.GISTIC2.focal.20150205.txt",header=T,sep="\t") #CNV - hg19

load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/IDHmut_min5_normalBrain_IDHwt.rda") #LGG.GBM.WGBS - hg19

h3k27ac.r <- dba.report(h3k27ac.a, th=1)
volcano <- as.data.frame(h3k27ac.r)
volcano$threshold <- "1"
volcano[volcano$FDR < 0.05 & volcano$Fold > 0, "threshold"] <- 2 #down no GCIMP-low
volcano[volcano$FDR < 0.05 & volcano$Fold < 0, "threshold"] <- 3 #up no GCIMP-low

ggplot(data=volcano,aes(x=Fold, y=-1*log10(FDR),color=threshold)) +
  geom_point() +
  xlab("Fold-change") + ylab("-1 * log10 of the Significance") +
  scale_color_manual(breaks=c("1","2","3"), # color scale (for points)
                     values=c("black", "green", "red"),
                     labels=c("Not Significant","Lost in GCIMP-low","Gained in GCIMP-low"),name="Legend")  +
  labs(title = "Volcano Plot - H3K27ac differential peaks") +
  theme(plot.title = element_text(hjust = 0.5))
table(volcano$threshold)
dim(volcano)

h3k4me3.r <- dba.report(h3k4me3.a, th=1)
volcano <- as.data.frame(h3k4me3.r)
volcano$threshold <- "1"
#volcano[volcano$FDR < 0.05 & volcano$Fold > 0, "threshold"] <- 2 #down no GCIMP-low
volcano[volcano$FDR < 0.003 & volcano$Fold < -2.5, "threshold"] <- 3 #up no GCIMP-low

ggplot(data=volcano,aes(x=Fold, y=-1*log10(FDR),color=threshold)) +
  geom_point() +
  xlab("Fold-change") + ylab("-1 * log10 of the Significance") +
  scale_color_manual(breaks=c("1","2","3"), # color scale (for points)
                     values=c("black", "green", "red"),
                     labels=c("Not Significant","Lost in GCIMP-low","Gained in GCIMP-low"),name="Legend")  +
  labs(title = "Volcano Plot - H3K4me3 differential peaks") +
  theme(plot.title = element_text(hjust = 0.5))
table(volcano$threshold)
dim(volcano)

ctcf.r <- dba.report(ctcf.a, th=1)
volcano <- as.data.frame(ctcf.r)
volcano$threshold <- "1"
volcano[volcano$FDR < 0.05 & volcano$Fold > 0, "threshold"] <- 2 #down no GCIMP-low
volcano[volcano$FDR < 0.05 & volcano$Fold < 0, "threshold"] <- 3 #up no GCIMP-low

ggplot(data=volcano,aes(x=Fold, y=-1*log10(FDR),color=threshold)) +
  geom_point() +
  xlab("Fold-change") + ylab("-1 * log10 of the Significance") +
  scale_color_manual(breaks=c("1","2","3"), # color scale (for points)
                     values=c("black", "green", "red"),
                     labels=c("Not Significant","Lost in GCIMP-low","Gained in GCIMP-low"),name="Legend")  +
  labs(title = "Volcano Plot - CTCF differential peaks") +
  theme(plot.title = element_text(hjust = 0.5))
table(volcano$threshold)
dim(volcano)

hmec.r <- dba.report(hmec.a, th=1)
volcano <- as.data.frame(hmec.r)
volcano$threshold <- "1"
volcano[volcano$FDR < 0.05 & volcano$Fold > 0, "threshold"] <- 2 #down no GCIMP-low
volcano[volcano$FDR < 0.05 & volcano$Fold < 0, "threshold"] <- 3 #up no GCIMP-low

ggplot(data=volcano,aes(x=Fold, y=-1*log10(FDR),color=threshold)) +
  geom_point() +
  xlab("Fold-change") + ylab("-1 * log10 of the Significance") +
  scale_color_manual(breaks=c("1","2","3"), # color scale (for points)
                     values=c("black", "green", "red"),
                     labels=c("Not Significant","Lost in GCIMP-low","Gained in GCIMP-low"),name="Legend")  +
  labs(title = "Volcano Plot - 5hmeC differential peaks") +
  theme(plot.title = element_text(hjust = 0.5))
table(volcano$threshold)
dim(volcano)


#chromHMM and homer
setwd("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/ChromHMM/GCIMP_PROJECT_HG38/Four_States_Results")
aux <- read.table("GCIMPhigh_s2_DU-6408_4_segments.bed")
aux <- aux[aux$V4 %in% "E2",]
write.table(aux,quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPhigh_DU-6408_CTCF.bed")
#findMotifsGenome.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPhigh_DU-6408_CTCF.bed hg38 /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/DU.6408 -size 100

aux <- read.table("GCIMPhigh_s3_DU-6542_4_segments.bed")
aux <- aux[aux$V4 %in% "E2",]
write.table(aux,quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPhigh_DU-6542_CTCF.bed")
#findMotifsGenome.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPhigh_DU-6542_CTCF.bed hg38 /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/DU.6542 -size 100

aux <- read.table("GCIMPhigh_s4_DU-7304_4_segments.bed")
aux <- aux[aux$V4 %in% "E2",]
write.table(aux,quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPhigh_DU-7304_CTCF.bed")
#findMotifsGenome.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPhigh_DU-7304_CTCF.bed hg38 /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/DU.7304 -size 100

aux <- read.table("GCIMPhigh_s9_DU-7007_4_segments.bed")
aux <- aux[aux$V4 %in% "E2",]
write.table(aux,quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPhigh_DU-7007_CTCF.bed")
#findMotifsGenome.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPhigh_DU-7007_CTCF.bed hg38 /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/DU.7007 -size 100

aux <- read.table("GCIMPhigh_s11_DU-7008_4_segments.bed")
aux <- aux[aux$V4 %in% "E2",]
write.table(aux,quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPhigh_DU-7008_CTCF.bed")
#findMotifsGenome.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPhigh_DU-7008_CTCF.bed hg38 /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/DU.7008 -size 100

aux <- read.table("GCIMPhigh_s12_06-0221_4_segments.bed")
aux <- aux[aux$V4 %in% "E2",]
write.table(aux,quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPhigh_06-0221_CTCF.bed")
#findMotifsGenome.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPhigh_06-0221_CTCF.bed hg38 /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/06.0221 -size 100

aux <- read.table("GCIMPlow_s5_DU-7010_4_segments.bed")
aux <- aux[aux$V4 %in% "E2",]
write.table(aux,quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPlow_DU-7010_CTCF.bed")
#findMotifsGenome.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPlow_DU-7010_CTCF.bed hg38 /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/DU.7010 -size 100

aux <- read.table("GCIMPlow_s6_06-0129_4_segments.bed")
aux <- aux[aux$V4 %in% "E2",]
write.table(aux,quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPlow_06-0129_CTCF.bed")
#findMotifsGenome.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPlow_06-0129_CTCF.bed hg38 /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/06-0129 -size 100

aux <- read.table("GCIMPlow_s7_06-1805_4_segments.bed")
aux <- aux[aux$V4 %in% "E2",]
write.table(aux,quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPlow_06-1805_CTCF.bed")
#findMotifsGenome.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPlow_06-1805_CTCF.bed hg38 /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/06.1805 -size 100

aux <- read.table("GCIMPlow_s8_06-2570_4_segments.bed")
aux <- aux[aux$V4 %in% "E2",]
write.table(aux,quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPlow_06-2570_CTCF.bed")
#findMotifsGenome.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/GCIMPlow_06-2570_CTCF.bed hg38 /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/chromHMM_4s/06.2570 -size 100


#Center CTCF
ctcf.r <- dba.report(ctcf.a)
ctcf.df <- as.data.frame(ctcf.r)
ctcf.g <- ctcf.df[ctcf.df$Fold < 0,]
write.table(ctcf.g[,c(1:3)],quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/DiffBind/CTCF_gained_hg38.bed")
ctcf.l <- ctcf.df[ctcf.df$Fold > 0,]
write.table(ctcf.l[,c(1:3)],quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/DiffBind/CTCF_lost_hg38.bed")

ctcf.g <- read.table("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/DiffBind/CTCF_gained_hg19.bed")
ctcf.l <- read.table("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/DiffBind/CTCF_lost_hg19.bed")
colnames(ctcf.l) <- c("chrom","start","end")
colnames(ctcf.g) <- c("chrom","start","end")
all.pr <- makeGRangesFromDataFrame(LGG.GBM.WGBS)

aux <- ctcf.g
aux$id <- seq(1,nrow(aux))
aux$centerCTCF <- floor((aux$end-aux$start)/2 + aux$start)
aux$start <- aux$centerCTCF - 500
aux$end <- aux$centerCTCF + 500
a <- findOverlaps(all.pr,makeGRangesFromDataFrame(aux))
new <- cbind(LGG.GBM.WGBS[queryHits(a),],aux[subjectHits(a),c("id","centerCTCF")])
#colnames(new)[12] <- "CTCF_id"
#obj <- data.frame(id=lost$id,position=floor((lost$end-lost$start)/2 + lost$start))
#new <- merge(new, aux[,c("id","centerCTCF")], by.x="CTCF_id",by.y="id",all.x=TRUE)
new$newStart <- new$start-new$centerCTCF


c <- melt(new[,c(4,6,8,10,11,13,15,17,19,21)],id=c("id","newStart"))

ggplot(c,aes(x=newStart,y=value,color=variable))  +
  geom_smooth(na.rm=TRUE,span=0.1,method="loess",se=F) + ylim(0,100)  +
  scale_x_continuous(breaks=arithSeq(-500,500))         +
  #xlim(-200,200) +
  scale_color_manual(values=c("darkgreen","firebrick","firebrick","gray","gray","#ffa500","#ffa500","cyan"),
                     labels=c("GCIMP-low","GCIMP-high1","GCIMP-high2","Non-tumor Brain1","Non-tumor Brain2","Mesenchymal-like1","Mesenchymal-like2","LGm6-GBM"),
                     name="Legend") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#+
# scale_x_continuous(breaks=arithSeq(-5000,5000))

arithSeq <- function(base, max){
  i <- base
  aux <- base
  while (i < (max)){
    i <-i + 25
    aux <- c(aux,i)
  }
  return(aux)
}

#Center enhancers
h3k27ac.r <- dba.report(h3k27ac.a)
h3k27ac.df <- as.data.frame(h3k27ac.r)
h3k27ac.g <- h3k27ac.df[h3k27ac.df$Fold < 0,]
write.table(h3k27ac.g[,c(1:3)],quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/DiffBind/H3K27ac_gained_hg38.bed")
h3k27ac.l <- h3k27ac.df[h3k27ac.df$Fold > 0,]
write.table(h3k27ac.l[,c(1:3)],quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/DiffBind/H3K27ac_lost_hg38.bed")

h3k27ac.g <- read.table("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/DiffBind/H3K27ac_gained_hg19.bed")
h3k27ac.l <- read.table("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/DiffBind/H3K27ac_lost_hg19.bed")
colnames(h3k27ac.l) <- c("chrom","start","end")
colnames(h3k27ac.g) <- c("chrom","start","end")

aux <- h3k27ac.l
aux$id <- seq(1,nrow(aux))
aux$centerCTCF <- floor((aux$end-aux$start)/2 + aux$start)
aux$start <- aux$centerCTCF - 1000
aux$end <- aux$centerCTCF + 1000
a <- findOverlaps(all.pr,makeGRangesFromDataFrame(aux))
new <- cbind(LGG.GBM.WGBS[queryHits(a),],aux[subjectHits(a),c("id","centerCTCF")])
#colnames(new)[12] <- "CTCF_id"
#obj <- data.frame(id=lost$id,position=floor((lost$end-lost$start)/2 + lost$start))
#new <- merge(new, aux[,c("id","centerCTCF")], by.x="CTCF_id",by.y="id",all.x=TRUE)
new$newStart <- new$start-new$centerCTCF


c <- melt(new[,c(4,6,8,10,11,13,15,17,19,21)],id=c("id","newStart"))

ggplot(c,aes(x=newStart,y=value,color=variable))  +
  geom_smooth(na.rm=TRUE,span=0.05,method="loess",se=F) + ylim(0,100)  +
  scale_x_continuous(breaks=arithSeq(-1000,1000))         +
  #xlim(-200,200) +
  scale_color_manual(values=c("darkgreen","firebrick","firebrick","gray","gray","#ffa500","#ffa500","cyan"),
                     labels=c("GCIMP-low","GCIMP-high1","GCIMP-high2","Non-tumor Brain1","Non-tumor Brain2","Mesenchymal-like1","Mesenchymal-like2","LGm6-GBM"),
                     name="Legend") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#+
# scale_x_continuous(breaks=arithSeq(-5000,5000))

arithSeq <- function(base, max){
  i <- base
  aux <- base
  while (i < (max)){
    i <-i + 50
    aux <- c(aux,i)
  }
  return(aux)
}

### Gene expression volcano
#symbol <- c("NKX6-1","TAL1","AR","FOXA2","LHX2","SOX9","SOX10","NFATC1","PDX1","NANOG") #H3K27ac gained know
#symbol <- c("HOXA10","CAMP","TBR1","ZNF189","TP53","ARNT","EBF1","IRF1") #H3K27ac gained de novo
#symbol <- c("FOSL1","FOSL2","ATF3","BATF","JUN","BACH2","NF1") #H3K27ac lost know
#symbol <- c("FOS","STAT5A","STAT5B","NF1","TBX20","SMAD2","SMAD3","ZBTB18","HMX1") #H3K27ac lost de novo
symbol <- c("FOS","JDP2","BATF","JUN","FOSL2","JUND") #H3K27ac lost de novo - TGAGTCAT motif
#symbol <- rownames(rna.seq.LGG.GBM.normalBrain[grep("^SOX+",rownames(rna.seq.LGG.GBM.normalBrain)),])
#symbol <- c("LIN9","RBBP4","MYBL1","MYBL2","MYB","FOXM1","CDK1","CCNB1")
symbol <- c("LOC284276")

normal.mrna <- rna.seq.LGG.GBM.normalBrain[rownames(rna.seq.LGG.GBM.normalBrain) %in% symbol,1:5]
normal.mrna <- normal.mrna+1
normal.mrna <- log2(normal.mrna)
aux <- as.character(subset(pd, cartoon %in% c("GCIMP-low"))$case.id)
low.mrna <- rna.seq.LGG.GBM.normalBrain[rownames(rna.seq.LGG.GBM.normalBrain) %in% symbol,colnames(rna.seq.LGG.GBM.normalBrain) %in% aux]
low.mrna <- low.mrna+1
low.mrna <- log2(low.mrna)
aux <- as.character(subset(pd, cartoon %in% c("GCIMP-high"))$case.id)
high.mrna <- rna.seq.LGG.GBM.normalBrain[rownames(rna.seq.LGG.GBM.normalBrain) %in% symbol,colnames(rna.seq.LGG.GBM.normalBrain) %in% aux]
high.mrna <- high.mrna+1
high.mrna <- log2(high.mrna)
aux <- as.character(subset(pd, cartoon %in% c("IDHmut-codel"))$case.id)
codel.mrna <- rna.seq.LGG.GBM.normalBrain[rownames(rna.seq.LGG.GBM.normalBrain) %in% symbol,colnames(rna.seq.LGG.GBM.normalBrain) %in% aux]
codel.mrna <- codel.mrna+1
codel.mrna <- log2(codel.mrna)
aux <- as.character(subset(pd, cartoon %in% c("Classic-like"))$case.id)
classic.mrna <- rna.seq.LGG.GBM.normalBrain[rownames(rna.seq.LGG.GBM.normalBrain) %in% symbol,colnames(rna.seq.LGG.GBM.normalBrain) %in% aux]
classic.mrna <- classic.mrna+1
classic.mrna <- log2(classic.mrna)
aux <- as.character(subset(pd, cartoon %in% c("Mesenchymal-like"))$case.id)
mes.mrna <- rna.seq.LGG.GBM.normalBrain[rownames(rna.seq.LGG.GBM.normalBrain) %in% symbol,colnames(rna.seq.LGG.GBM.normalBrain) %in% aux]
mes.mrna <- mes.mrna+1
mes.mrna <- log2(mes.mrna)
aux <- as.character(subset(pd, cartoon %in% c("LGm6-GBM"))$case.id)
lgm6.mrna <- rna.seq.LGG.GBM.normalBrain[rownames(rna.seq.LGG.GBM.normalBrain) %in% symbol,colnames(rna.seq.LGG.GBM.normalBrain) %in% aux]
lgm6.mrna <- lgm6.mrna+1
lgm6.mrna <- log2(lgm6.mrna)
aux <- as.character(subset(pd, cartoon %in% c("PA-like"))$case.id)
pa.mrna <- rna.seq.LGG.GBM.normalBrain[rownames(rna.seq.LGG.GBM.normalBrain) %in% symbol,colnames(rna.seq.LGG.GBM.normalBrain) %in% aux]
pa.mrna <- pa.mrna+1
pa.mrna <- log2(pa.mrna)

normal.mrna$symbol <- rownames(normal.mrna)
normal.m <- melt(normal.mrna,id="symbol")
low.mrna$symbol <- rownames(low.mrna)
low.m <- melt(low.mrna,id="symbol")
high.mrna$symbol <- rownames(high.mrna)
high.m <- melt(high.mrna,id="symbol")
codel.mrna$symbol <- rownames(codel.mrna)
codel.m <- melt(codel.mrna,id="symbol")
classic.mrna$symbol <- rownames(classic.mrna)
classic.m <- melt(classic.mrna,id="symbol")
mes.mrna$symbol <- rownames(mes.mrna)
mes.m <- melt(mes.mrna,id="symbol")
lgm6.mrna$symbol <- rownames(lgm6.mrna)
lgm6.m <- melt(lgm6.mrna,id="symbol")
pa.mrna$symbol <- rownames(pa.mrna)
pa.m <- melt(pa.mrna,id="symbol")

normal.m$Subtype <- "Non-tumor Brain"
low.m$Subtype <- "GCIMP-low"
high.m$Subtype <- "GCIMP-high"
codel.m$Subtype <- "IDHmut-codel"
classic.m$Subtype <- "Classic-like"
mes.m$Subtype <- "Mesenchymal-like"
lgm6.m$Subtype <- "LGm6-GBM"
pa.m$Subtype <- "PA-like"
all <- Reduce(rbind,list(normal.m,low.m,high.m,codel.m,classic.m, mes.m, lgm6.m,pa.m))
all$Subtype <- as.factor(all$Subtype)
all$Subtype <- factor(all$Subtype, levels = c("Non-tumor Brain","GCIMP-low","GCIMP-high","IDHmut-codel","Classic-like","Mesenchymal-like","LGm6-GBM","PA-like"))

all$symbol <- as.factor(all$symbol)
all$symbol <- factor(all$symbol, levels = symbol)

#plot gene expression
ggplot(all, aes(Subtype, value,fill = Subtype)) +
  geom_boxplot(outlier.shape=NA)  +
  geom_jitter() +
  facet_wrap( ~ symbol,scales="free") +
  theme(axis.text.x = element_blank(),strip.text = element_text(size=25),axis.text.y = element_text(size=20),legend.text=element_text(size=25),legend.title=element_text(size=25),legend.key.size=unit(2,"cm")) +
  scale_fill_manual(values = c("gray","darkgreen","firebrick","#9f20f0","#ffa500","#FFFF00","#006680","#00CCFF"))  #salvar 15 por 10 em pdf ou 2500 por 1500 em png


### H3k27ac gained - Peak annotation
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/DiffBind_H3K27ac_gained/TAAT_motif
#annotatePeaks.pl ../../../DiffBind/H3K27ac_gained_hg38.bed hg38 -m TAAT.motif > annotation.txt
#findMotifsGenome.pl ../../../DiffBind/H3K27ac_gained_hg38.bed hg38 . -find TAAT.motif > outputfile.txt
setwd("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/DiffBind_H3K27ac_gained/TAAT_motif")
a <- read.table("outputfile.txt",sep="\t",header=T,na.strings="")
peaks <- read.table("../../../DiffBind/H3K27ac_gained_hg38.bed") #402 peaks with TAAT signature
colnames(peaks) <- c("chrom","start","end")
rownames(peaks) <- NULL
peaks <- peaks[unique(a$PositionID),]
txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb_hg38)
symbols <- unlist(mapIds(org.Hs.eg.db, genes$gene_id, "SYMBOL", "ENTREZID", multiVals = "first"))
genes <- as.data.frame(genes)
genes$symbol <- as.matrix(symbols)
aux <- rna.seq.LGG.GBM.normalBrain[rna.seq.LGG.GBM.normalBrain$p.value.adj < 0.05 & rna.seq.LGG.GBM.normalBrain$DiffMean > 0,"entrezID"]
genes <- genes[genes$gene_id %in% aux,]
a <- distanceToNearest(makeGRangesFromDataFrame(genes),makeGRangesFromDataFrame(peaks))
rownames(peaks) <- NULL
a <- as.data.frame(a)
a <- a[a$distance < 20000,] 
aux <- genes[a$queryHits,]
write.table(aux$gene_id,quote=F,row.names=F,col.names = F,sep="\t",file="Nearest_Upgene.txt")

#Run find Motifs on Upregulated genes
symbol <- (rna.seq.LGG.GBM.normalBrain[rna.seq.LGG.GBM.normalBrain$p.value.adj < 0.05 & (rna.seq.LGG.GBM.normalBrain$DiffMean) > 3,"entrezID"])
write.table(symbol,quote=F,row.names=F,col.names = F,sep="\t",file="../../Upreg.genes/EntrezID.txt")

## DREAM-complex
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb_hg38)
symbols <- unlist(mapIds(org.Hs.eg.db, genes$gene_id, "SYMBOL", "ENTREZID", multiVals = "first"))
genes <- as.data.frame(genes)
genes$symbol <- as.matrix(symbols)

symbol <- c("LIN9","RBBP4","MYBL1","MYBL2","MYB","FOXM1","CDK1","CCNB1","EN1")
distanceToNearest(makeGRangesFromDataFrame(genes[genes$symbol %in% symbol,],TRUE),ctcf.r)


### Center H3K27ac peaks - Heatmap
setwd("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/DiffBind_H3K27ac_gained/TAAT_motif")
a <- read.table("outputfile.txt",sep="\t",header=T,na.strings="")
peaks <- read.table("../../../DiffBind/H3K27ac_gained_hg38.bed") #402 peaks with TAAT signature
colnames(peaks) <- c("chrom","start","end")
rownames(peaks) <- NULL
peaks <- peaks[unique(a$PositionID),]
#liftOver
write.table(peaks,file="peaks_TAAT_motif_hg38.bed",quote=F,row.names=F,sep="\t",col.names = F)
#liftOver peaks_TAAT_motif_hg38.bed /media/data1/biosoftwares/hg38ToHg19.over.chain peaks_TAAT_motif_hg19.bed unMapped.bed

peaks <- read.table("peaks_TAAT_motif_hg19.bed")
colnames(peaks) <- c("chrom","start","end")
aux <- peaks
aux$id <- seq(1,nrow(peaks))
aux$centerPeak <- floor((aux$end-aux$start)/2 + aux$start)
aux$start <- aux$centerPeak - 1000
aux$end <- aux$centerPeak + 1000

a <- findOverlaps(makeGRangesFromDataFrame(LGG.GBM.WGBS),makeGRangesFromDataFrame(aux)) #peak 255 doesnt have DNA methylation data on it
sub <- cbind(LGG.GBM.WGBS[queryHits(a),],aux[subjectHits(a),c("id","centerPeak")])
sub$newStart <- sub$start-sub$centerPeak
#order -> TCGA.06.0128, TCGA.16.1460, TCGA.19.1788, TCGA.14.1454, TCGA.14.1401, TCGA.14.3477, 59, 60

#aux2 <- subset(lost, p.valueADJ_lowXhigh_DNAmeth < 0.05 & DiffMean < 0)$id #83 CTCF binding sites e 1044 CpGs
#sub <- subset(new, id %in% aux2)

sub$id <- as.factor(sub$id)
teste1 <- NULL
#for(j in 1:nlevels(sub$id)){
for(j in 1:nlevels(sub$id)){
  #for(j in 1:2){
  print(j)
  aux <- sub[sub$id %in% levels(sub$id)[j],c("centerPeak","percentMeth.60","newStart","id")]
  aux <- na.omit(aux)
  i <- 1000
  while(i > -1001){
    a <- subset(aux, newStart > i - 50 & newStart < i+1)
    if(nrow(a) < 1)
      a <- cbind(i,nrow(a),NA,levels(sub$id)[j])
    else
      a <- cbind(i,nrow(a),mean(a[,"percentMeth.60"],na.rm=TRUE),levels(sub$id)[j])
    if(is.null(teste1)){
      teste1 <- a
    }
    else{
      teste1 <- rbind(teste1,a)
    }
    i<-i-50
  }
}

b <- teste1
teste1 <- as.data.frame(teste1)
colnames(teste1) <- c("seq","nCpGs","meth","Peak")
teste1$group <- NA
teste1$meth <- as.numeric(as.character(teste1$meth))
summary(teste1)
#teste1$seq <- as.numeric(as.character(teste1$seq))
#aux <- teste1[teste1$seq < 161 & teste1$seq > -159,]
#aux

summary(as.numeric(as.character(teste1$seq)))

arithSeq <- function(base, max){
  i <- base
  aux <- base
  while (i < (max)){
    i <-i + 50
    aux <- c(aux,i)
  }
  return(aux)
}

teste1$seq <- factor(teste1$seq, levels=arithSeq(-1000,1000))
teste1[teste1$meth <= 50 & !is.na(teste1$meth),"group"] <- "Unmethylated"
teste1[teste1$meth > 50 & !is.na(teste1$meth),"group"] <- "Methylated"
teste1[is.na(teste1$meth),"group"] <- "No.CpG"

library(reshape2)
aux <- dcast(teste1, Peak ~ seq, value.var = "meth")
aux <- aux[,18:26]
head(aux)
norm2 <- apply(aux,1,mean,na.rm=TRUE)
norm2 <- as.data.frame(norm2)
norm2$id <- seq(1,nrow(norm2))
norm2 <- norm2[with(norm2,order(norm2)),]
save(norm2,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/homer.motifs/DiffBind_H3K27ac_gained/TAAT_motif/ordem_low.rda")


library(ggplot2)
load("ordem_low.rda")
teste1$Peak <- factor(teste1$Peak, levels=norm2$id)
ggplot(teste1, aes(x=seq,y=factor(Peak))) +
  geom_tile(aes(fill=factor(group)))+
  scale_fill_manual(values=c("red","black","cyan")) +
  theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white"),axis.text.x = element_text(angle = 45, hjust = 1),axis.text.y = element_blank(),panel.spacing = unit(0, "lines"))


c <- melt(sub[,c(4,6,8,10,11,13,15,17,19,21)],id=c("id","newStart"))

ggplot(c,aes(x=newStart,y=value,color=variable))  +
  geom_smooth(na.rm=TRUE,span=0.2,method="loess") + ylim(0,100)  +
  scale_x_continuous(breaks=arithSeq(-1000,1000))         +
  #xlim(-200,200) +
  scale_color_manual(values=c("darkgreen","firebrick","firebrick","gray48","gray48","#ffa500","#ffa500","cyan"),
                     labels=c("GCIMP-low","GCIMP-high1","GCIMP-high2","Non-tumor Brain1","Non-tumor Brain2","Mesenchymal-like1","Mesenchymal-like2","LGm6-GBM"),
                     name="Legend") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#### Investigate where 5hmC losses are located
gencode <- readGFF("/home/thais/Downloads/gencode.v25.annotation.gtf")

ids <- read.table("/home/thais/Downloads/gencode.v25.metadata.EntrezGene")
colnames(ids) <- c("transcript_id","EntrezGene")
ids <- merge(ids, gencode[,c("transcript_id","gene_id")],by="transcript_id",all=TRUE)

ids <- ids[!duplicated(ids$EntrezGene),]

exon <- gencode[gencode$type %in% "transcript",]
hmec.r <- dba.report(hmec.a)
b <- hmec.r[mcols(hmec.r)[,4] > 0]
a <- distanceToNearest(b,makeGRangesFromDataFrame(exon))
a <- as.data.frame(a)
a <- a[a$distance == 0,]
table(exon[unique(a$subjectHits),"gene_type"])
exon <- exon[unique(a$subjectHits),]

b <- b[unique(a$queryHits),]
exon <- gencode[!(gencode$type %in% c("gene","transcript")),]
a <- distanceToNearest(b,makeGRangesFromDataFrame(exon))
a <- as.data.frame(a)
sum(a$distance == 0)
a <- a[a$distance != 0,]
b <- b[unique(a$queryHits),]


### Distribution of genes across chromosomes
hmec.r <- dba.report(hmec.a)
seqlevels(hmec.r) <- paste0("chr",seqlevels(hmec.r))

a <- hmec.r[mcols(hmec.r)[,4] > 0,]
aux <- findOverlaps(makeGRangesFromDataFrame(genes),a)
aux <- genes[unique(queryHits(aux)),]

ggplot(aux,aes(seqid,fill=RNAseqALL.summary)) + 
  geom_bar() + 
  scale_fill_manual(values=c("gray48", "skyblue2", "orange", "black"), name="Legend") + 
  scale_x_discrete(drop=FALSE)


## Copy number variation
library(TCGAbiolinks)
aux <- pd[pd$cartoon %in% c("GCIMP-low","GCIMP-high"),"case.id"]
setwd("/media/data2/CNV")
query.exp.hg38 <- GDCquery(project = c("TCGA-GBM","TCGA-LGG"), 
                           data.category = "Copy Number Variation", 
                           data.type = "Masked Copy Number Segment", 
                           barcode =  as.character(aux))
GDCdownload(query.exp.hg38)
cnv.data <- GDCprepare(query = query.exp.hg38,
                     save = TRUE, 
                     save.filename = "Masked_CNV.rda")
cnv.data <- cnv.data[substr(cnv.data$Sample,14,15) %in% c("01"),]
cnv.data <- as.data.frame(cnv.data)

save(cnv.data,file="/media/data2/CNV/Masked_CNV.rda")

#hg19 annotation
library(gaia)
aux <- read.table("/media/data2/CNV/genome.info.6.0_hg19.na31_minus_frequent_nan_probes_sorted_2.1.txt")
aux$V4 <- aux$V3 + 1
aux <- aux[,c(2,3,4,1)]
aux$V2 <- paste0("chr",aux$V2)
write.table(aux,row.names = F,col.names = F,quote=F,sep="\t",file="/media/data2/CNV/probes_hg19.txt")

probes <- read.table("/media/data2/CNV/snp6.na35.liftoverhg38.txt",header=T)
probes$chr <- as.character(probes$chr)
probes$chr[probes$chr %in% "X"] <- "23"
probes$chr[probes$chr %in% "Y"] <- "24"
probes$chr <- sapply(probes$chr,as.integer)
markerID <- apply(probes,1,function(x) paste0(x[2],":",x[3]))
probes <- probes[-which(duplicated(markerID)),]
probes <- probes[,c(1,2,3)]
colnames(probes) <- c("Probe.Name","Chromosome","Start")
markers_obj <- load_markers(as.data.frame(probes))

cnv.data <- as.data.frame(cnv.data)
cnv.data$Chromosome[cnv.data$Chromosome %in% "X"] <- "23"
cnv.data$Chromosome[cnv.data$Chromosome %in% "Y"] <- "24"
cnv.data$Chromosome <- sapply(cnv.data$Chromosome,as.integer)
cnv.data$Aberration <- NA
cnv.data[cnv.data$Segment_Mean < -0.3,"Aberration"] <- 0
cnv.data[cnv.data$Segment_Mean > 0.3,"Aberration"] <- 1
cnv.data <- cnv.data[!is.na(cnv.data$Aberration),]
cnv.data <- cnv.data[,-6]
colnames(cnv.data) <- c( "Sample.Name", "Chromosome", "Start", "End", "Num.of.Markers", "Aberration")
cnv_obj <- load_cnv(cnv.data, markers_obj, length(unique(cnv.data$Sample)))

setwd("/media/data2/aligned_hg38/hMEDIP/MACS_broad/")
lines <- data.frame(Sample=NULL,Width=NULL,N.of.peaks=NULL)
for(i in list.files(pattern="*woDup_peaks.broadPeak")){
  aux <- fread(i)
  aux <- data.frame(Sample="a",Width=mean(aux$V3-aux$V2),N.of.peaks=nrow(aux))
  lines <- rbind(lines,aux)
}
lines$Sample <- c("DU-6408","DU-6542","DU-7304","DU-7010","06-0129","06-1805","06-2570","DU-7007","DU-7008","06-0221")
lines$Subtype <- c(rep("GCIMP-high",3),rep("GCIMP-low",4),rep("GCIMP-high",3))
lines <- lines[c(4:7,1:3,8:10),]
lines$Sample <- factor(lines$Sample)
lines$Sample <- factor(lines$Sample, levels=levels(lines$Sample)[c(1,3,4,9,2,5,6,7,8,10)])

ggplot(lines,aes(Sample,Width,fill=Subtype)) + geom_col() + scale_fill_manual(values=c("firebrick","darkgreen")) + ggtitle("Width")


ggplot(genes,(aes(x=H3K4me3.Fold,y=-1*log10(H3K4me3.FDR),color=RNAseqALL.summary))) + geom_point()


ggplot(data=genes,aes(x=-1*H3K4me3.Fold,y=-1*log10(H3K4me3.FDR))) +
  geom_point(data=subset(genes,RNAseqALL.summary %in% c("No.info","No.difference")),aes(x=H3K4me3.Fold,y=-1*log10(H3K4me3.FDR), colour=RNAseqALL.summary),size=1) +
  geom_point(data=subset(genes,RNAseqALL.summary %in% c("Upregulated","Downregulated")),aes(x=H3K4me3.Fold,y=-1*log10(H3K4me3.FDR), colour=RNAseqALL.summary),size=1) +
  xlab("Gene expression (fold-change)") + ylab("-1 * log10 of the Significance") +
  labs(title = "Volcano Plot") +
  scale_color_manual(breaks=c("1","2","3","4","5","6"), # color scale (for points)
                     values=c("grey", "green", "red","darkgreen","firebrick","gray47"),
                     labels=c("Not Significant","Downregulated in GCIMP-low","Upregulated in GCIMP-low","Downregulated in GCIMP-low (near CTCF)","Upregulated in GCIMP-low (near CTCF)","Not Significant (near CTCF)"),name="Legend")



ggplot(data=genes,aes(x=RNAseq.hg38.FC,y=-1*log10(RNAseq.hg38.pv.adj),color=RNAseq.hg38.summary)) + geom_point()  + theme_bw() +
  xlab("Gene expression fold-change") + ylab("-1 * log10 of the Significance (FDR)") + ggtitle("Differential express genes",subtitle="G-CIMP-low vs G-CIMP-high") +
  scale_color_manual(values=c("skyblue2", "darkgrey","darkgrey", "orange"),name="Gene Exp Summary")

##### LSD1
a <- read.csv("/media/data1/thais.projects/GCIMP-low/gene_orthologs_mm10_human_May032018.txt")
s <- read.table("~/Downloads/GSM1941292.txt",skip=6)
colnames(s) <- c("chr","start","end","refseq","score","strand","symbol")
sum(unique(as.character(s$symbol)) %in% as.character(a$Gene.name))
s <- merge(s,a,by.x="symbol",by.y="Gene.name",all.x=T)
s <- merge(s,genes[,c("gene_name","RNAseq.hg38.pv.adj","RNAseq.hg38.FC","RNAseqALL.summary","RNAseq.hg38.summary")],by.x="Human.gene.name",by.y="gene_name",all.x=T)
