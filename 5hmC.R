require(rtracklayer)
gencode <- readGFF("/media/data2/TCGA_RNAseq/GDC_gencodeV22/gencode.v22.annotation.gtf")
gencode <- gencode[grep("+_PAR_+",gencode$gene_id,invert=T),]
aux <- strsplit(gencode$gene_id, "\\." )
aux <- unlist(lapply(aux, "[[", 1 ))
gencode$gene_id <- aux #Ensembl id list
genes <- gencode[gencode$type %in% "gene",]

tss <- promoters(makeGRangesFromDataFrame(genes), upstream = 0, downstream = 1)

hmec.r <- dba.report(hmec.a,th=1)
a <- suppressWarnings(distanceToNearest(hmec.r,tss,ignore.strand=T))
cat.5hmC <- as.data.frame(hmec.r)
rownames(cat.5hmC) <- NULL
cat.5hmC$geneID <- NA
cat.5hmC$distance.to.TSS <- NA
cat.5hmC[unique(queryHits(a)),"distance.to.TSS"] <- mcols(a)[,1]
a <- a[mcols(a)[,1] < 2000,]
cat.5hmC$type <- NA
cat.5hmC[unique(queryHits(a)),"type"] <- "TSS" #no duplicates
aux <- as.character(genes[subjectHits(a),"gene_id"])
cat.5hmC[queryHits(a),"geneID"] <- aux
cat.5hmC$id <- seq(1,nrow(cat.5hmC))

aux <- cat.5hmC[is.na(cat.5hmC$type),]
a <- suppressWarnings(findOverlaps(makeGRangesFromDataFrame(aux),makeGRangesFromDataFrame(genes),ignore.strand=T))
b <- data.frame(id = as.character(aux[queryHits(a),"id"]), geneID = as.character(genes[subjectHits(a),"gene_id"]))
#b <- aggregate(geneID ~ id, data = b, c)
aux <- aux[unique(queryHits(a)),]
genes <- gencode[gencode$type %in% "exon",]
a <- suppressWarnings(findOverlaps(makeGRangesFromDataFrame(aux),makeGRangesFromDataFrame(genes),ignore.strand=T))
aux$type <- "intron"
aux[unique(queryHits(a)),"type"] <- "exon"

id <- aux[aux$type %in% "exon","id"]
cat.5hmC[cat.5hmC$id %in% id,"type"] <- "exon"

id <- aux[aux$type %in% "intron","id"]
cat.5hmC[cat.5hmC$id %in% id,"type"] <- "intron"
cat.5hmC[is.na(cat.5hmC$type),"type"] <- "intergenic"
cat.5hmC$status <- NA
cat.5hmC[cat.5hmC$Fold < 0 & cat.5hmC$FDR < 0.01, "status"] <- "gain"
#cat.5hmC[cat.5hmC$Fold > 0 & cat.5hmC$FDR < 0.01, "status"] <- "loss"
cat.5hmC[cat.5hmC$Fold > 2.5 & cat.5hmC$FDR < 0.00001, "status"] <- "loss"
cat.5hmC[is.na(cat.5hmC$status),"status"] <- "No difference"

ggplot(cat.5hmC,aes(type,fill=type)) + 
  geom_bar()+
  facet_wrap(~ status,scales="free")

load("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/DiffBind/enhancer_annot_v22.rda")
cat.5hmC$overlap.enhancer <- "No"
cat.5hmC$enhancer.status <- NA
aux <- enhancer[enhancer$type %in% "intergenic",]
a <- suppressWarnings(findOverlaps(makeGRangesFromDataFrame(cat.5hmC),makeGRangesFromDataFrame(aux),ignore.strand=T))
cat.5hmC[(queryHits(a)),"overlap.enhancer"] <- "Yes"
cat.5hmC[(queryHits(a)),"enhancer.status"] <- as.character(aux[(subjectHits(a)),"status"])

aux <- cat.5hmC[cat.5hmC$overlap.enhancer]
ggplot(cat.5hmC,aes(overlap.enhancer,fill=overlap.enhancer)) + 
  geom_bar()+
  facet_wrap(~ status,scales="free")

#save(cat.5hmC,file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/cat_5hmc_v2.rda")

library(GEOquery)
a <- getGEO('GSE73895',GSEMatrix=TRUE)
meta <- (pData(phenoData(a[[1]])))
meta$meth <- rep(c("5mC","5hmC"),30)
meta$sample <- paste0("tissue_",rep(1:30,rep(2,30)))
array.5hmC <- exprs(a$GSE73895_series_matrix.txt.gz)
identical(colnames(array.5hmC),as.character(meta$geo_accession)) #TRUE
colnames(array.5hmC) <- as.character(meta$title)

probe.high5hmC <- read.table("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/probes_high_5hmC.txt",header=T)
aux <- meta[meta$meth %in% "5hmC","title"]
aux <- array.5hmC[,as.character(aux)]
aux <- as.data.frame(aux)

aux$mean <- apply(aux[,1:30],1,mean,na.rm=T)
aux$sd <- apply(aux[,1:30],1,sd,na.rm=T)
aux <- aux[aux$sd > 0.15,]

aux <- aux[as.character(probe.high5hmC$probeID),]


heatmap.plus(as.matrix(aux),
             col=col,
             #col=c("blue","yellow","orange","black","green","red"),
             cexCol = 0.5,
             labRow = NA,
             #labCol = NA,
             scale = "none")


load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/low_high_450.rda")

low.p <- low[as.character(probe.high5hmC$probeID),-c(1:4)]
aux <- hclust(dist(t(low.p)))
low.p <- low.p[,aux$order]
high.p <- high[as.character(probe.high5hmC$probeID),-c(1:4)]
aux <- hclust(dist(t(high.p)))
high.p <- high.p[,aux$order]
aux <- cbind(low.p,high.p)

a <- na.omit(aux)
p <- setdiff(rownames(aux),rownames(a))
low.p <-  hclust(dist(a))
aux <- rbind(aux[p,],a[low.p$order,])

low.p <- low[as.character(probe.high5hmC$probeID),-c(1:4)]

heatmap.plus(as.matrix(aux),
             col=col,
             #col=c("blue","yellow","orange","black","green","red"),
             cexCol = 0.5,
             Rowv=NA,
             Colv=NA,
             labRow = NA,
             ColSideColors = cbind(c(rep("darkgreen",ncol(low.p)),(rep("firebrick",ncol(high.p)))),c(rep("darkgreen",ncol(low.p)),(rep("firebrick",ncol(high.p))))),
             labCol = NA,
             scale = "none")

load("/media/data1/450k_hg38_GDC.rda")

probe.high5hmC <- hg38.450k[as.character(hg38.450k$Composite.Element.REF) %in% as.character(probe.high5hmC$probeID),]
probe.high5hmC <- probe.high5hmC[!(probe.high5hmC$Chromosome %in% "*"),]
hmec.r <- dba.report(hmec.a,th=1)
aux <- findOverlaps(makeGRangesFromDataFrame(probe.high5hmC),(hmec.r))
aux <- cat.5hmC[subjectHits(aux),]

# Heatmap de metilação e 5hmC - loss
load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/hg38_IDHmut_min5_normalBrain_IDHwt_CpGtype_mean.rda")
hmec.r <- dba.report(hmec.a, bCounts = T, th=0.01)
#hmec.r <- hmec.r[mcols(hmec.r)[,4] > 0 & mcols(hmec.r)[,6] < 0.01,]
hmec.r <- hmec.r[mcols(hmec.r)[,4] > 2.5 & mcols(hmec.r)[,6] < 0.00001,]
#hmec.r <- hmec.r[1:1000,]
#all.GR <- makeGRangesFromDataFrame(LGG.GBM.WGBS)
aux <- findOverlaps(hmec.r,all.ps)
aux <- as.data.frame(aux)
hmC <- hmec.r
hmC <- as.data.frame(hmC)
hmC$meanLow <- NA
hmC$meanHigh1 <- NA
hmC$meanHigh2 <- NA
hmC$meanNormal <- NA
for(i in 1:nrow(hmC)){
  a <- aux[aux$queryHits == i,]
  a <- LGG.GBM.WGBS[a$subjectHits,]
  a <- apply(a[,c("percentMeth.0128","percentMeth.1460","percentMeth.1788","meanNormal")],2,mean,na.rm=T)
  hmC[i,c("meanLow","meanHigh1","meanHigh2","meanNormal")] <- a
}

#a <- (as.matrix(hmC[,c("Conc_G.CIMP.low","meanLow")]))
#heatpairs(as.matrix(a))

ha2 = HeatmapAnnotation(df = data.frame(Sample = c("GCIMP-low","GCIMP-high","GCIMP-high","Non-tumor Brain","Non-tumor Brain")),
                        col = list(Sample = c("Non-tumor Brain" = "gray","GCIMP-low" = "darkgreen","GCIMP-high"="firebrick")))

ha1 = HeatmapAnnotation(df = data.frame(Sample = c(rep("GCIMP-low",4),rep("GCIMP-high",6))),
                        col = list(Sample = c("GCIMP-low" = "darkgreen","GCIMP-high"="firebrick")))


colors <- c("lightpink1","lightpink2","lightpink3","mediumpurple2","mediumpurple3","mediumpurple4")
#brewer.pal(6, "PuOr")

#order by hc
Heatmap(log2(hmC[,c(21:12)]), name = "ChIP-seq", col=colors, cluster_columns = FALSE,cluster_rows = TRUE,show_row_names=FALSE,show_column_names = TRUE,top_annotation = ha1) +
  Heatmap(hmC[,c(22:25)], name = "WGBS", col=jet.colors(75), cluster_columns = FALSE,show_row_names=FALSE,show_column_names = F,top_annotation = ha2)

#order by p-value
aux <- hmC[with(hmC, order(-Fold, -FDR)),]
Heatmap(log2(aux[,c(21:12)]), name = "ChIP-seq", col=colors, cluster_columns = FALSE,cluster_rows = FALSE,show_row_names=FALSE,show_column_names = TRUE,top_annotation = ha1) +
  Heatmap(aux[,c(22:25)], name = "WGBS", col=jet.colors(75), cluster_columns = FALSE,show_row_names=FALSE,show_column_names = F,top_annotation = ha2)


#Overlap of intronic loss of 5hmC with gene expression
a <- cat.5hmC[cat.5hmC$type %in% "intron" & cat.5hmC$status %in% "loss",]
aux <- findOverlaps(makeGRangesFromDataFrame(a),makeGRangesFromDataFrame(genes))
a <- genes[unique(subjectHits(aux)),]
table(a$RNAseqALL.summary)


#################### PAPER
## Figura 1
#Panel A - Density plot (Distribution of peaks across CpG types)
hmec.r <- dba.report(hmec.a,th=1)
island <- read.table("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/CpGislands_hg38_april17.txt")
island <- island[,2:4]
colnames(island) <- c("chrom","start","end")
b <- cat.5hmC #objeto que você tem as probes a serem classificadas
#b <- as.data.frame(b)
b$CpG_type <- "Open-Sea"
aux <- findOverlaps(makeGRangesFromDataFrame(cat.5hmC),makeGRangesFromDataFrame(island),type="within")
b[queryHits(aux),"CpG_type"] <- "CpG-island"
shore <- island[,c(1,2)]
shore$startS <- shore$start - 2000
shore <- shore[,c(1,3,2)]
aux <- island[,c(1,3)]
aux$end2 <- aux$end + 2000
colnames(shore) <- c("chrom","start","end")
colnames(aux) <- c("chrom","start","end")
shore <- rbind(shore,aux)
aux <- findOverlaps(makeGRangesFromDataFrame(cat.5hmC),makeGRangesFromDataFrame(shore),type="within")
b[queryHits(aux),"CpG_type"] <- ifelse(b[queryHits(aux),"CpG_type"] %in% "CpG-island","CpG-island","CpG-shore")

library(reshape)

aux <- b[,c(7,16)]
colnames(aux) <- c("Peaks","Location")
aux$Subtype <- "G-CIMP-high"

b.m <- b[,c(8,16)]
colnames(b.m) <- c("Peaks","Location")
b.m$Subtype <- "G-CIMP-low"
b.m <- rbind(b.m,aux)

ggplot(b.m, aes(Peaks, colour=Location, fill=Location, color=Location)) + 
  geom_density(alpha=0.6) +
  facet_wrap(~ Subtype) +
  scale_x_continuous(breaks = c(0:15)) +
  scale_fill_manual(values = c("orange","darkslateblue","turquoise")) +
  scale_color_manual(values = c("orange","darkslateblue","turquoise")) +
  theme_bw()

ggplot(b.m, aes(Location, Peaks,fill = Subtype,color=Location)) + 
  #geom_jitter() +
  #geom_boxplot(aes(fill = Location)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  #facet_wrap(~ Subtype) +
  scale_fill_manual(values = c("firebrick","darkgreen")) +
  scale_color_manual(values = c("orange","darkslateblue","turquoise")) +
  theme_bw()


#Panel B -l Scatter plot (5hmC and 450k)
load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/tudao.27.450.rda")
colnames(LGG.GBM) <- substr(colnames(LGG.GBM),1,12)
LGG.GBM <- LGG.GBM[,colnames(LGG.GBM) %in% c("TCGA-DU-7010","TCGA-06-0129","TCGA-06-1805","TCGA-06-2570","TCGA-DU-6408","TCGA-DU-6542","TCGA-DU-7304","TCGA-DU-7007","TCGA-DU-7008","TCGA-06-0221")]
a <- read.table("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/jhu-usc.edu_GBM.HumanMethylation450.2.lvl-3.TCGA-06-0221-01A-01D-A45W-05.gdc_hg38.txt",header=T,sep="\t")
LGG.GBM <- merge(LGG.GBM,a[,c(1,2)],by.x=0,by.y="Composite.Element.REF")
rownames(LGG.GBM) <- as.character(LGG.GBM$Row.names)
colnames(LGG.GBM)[ncol(LGG.GBM)] <- "TCGA-06-0221"
LGG.GBM <- LGG.GBM[,-1]
LGG.GBM <- LGG.GBM[,c(1,2,3,4,6,10,5,7,8,9)]
#save(LGG.GBM, file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/gcimp.10samples.merged.450k27.rda")
a <- read.table("/media/data1/hg38_450k.bed")
load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/GCIMPlow.90probes.Rda")
a <- a[as.character(a$V4) %in% as.character(GCIMPlow.probes$probes),]
colnames(a) <- c("chrom","start","end","probeID")
hmec.r <- dba.report(hmec.a, bCounts=T, th=1)
aux <- findOverlaps(hmec.r,makeGRangesFromDataFrame(a))
b <- hmec.r[queryHits(aux),]
b <- as.data.frame(b)
colnames(b)[12:21] <- gsub("^s[0-9]{1,}_+","TCGA.",colnames(b)[12:21])
colnames(b)[12:21] <- gsub("[.]","-",colnames(b)[12:21])
LGG.GBM$meanLow <- apply(LGG.GBM[,c(7:10)],1,mean,na.rm=T)
LGG.GBM$meanHigh <- apply(LGG.GBM[,c(1:6)],1,mean,na.rm=T)
LGG.GBM <- LGG.GBM[as.character(a$probeID),]
LGG.GBM <- LGG.GBM[subjectHits(aux),]
aux <- cbind(LGG.GBM[,c("meanLow","meanHigh")],b[,c("Conc_G.CIMP.high","Conc_G.CIMP.low")])
LGG.GBM <- LGG.GBM[,c(10:1)]
b <- b[,colnames(LGG.GBM)]

b <- log2(b+1)
rm <- rowMeans(b)
sx <- apply(b, 1, sd)
zz <- sweep(b,1,rm)
zz <- sweep(zz, 1, sx, "/")
zzz <- t(scale(t(b)))

rm <- rowMeans(LGG.GBM)
sx <- apply(LGG.GBM, 1, sd)
zz <- sweep(LGG.GBM,1,rm)
zz <- sweep(zz, 1, sx, "/")
xxx <- t(scale(t(LGG.GBM)))

#library(RColorBrewer)
colors <- c("lightpink1","lightpink2","lightpink3","mediumpurple2","mediumpurple3","mediumpurple4")

Heatmap(LGG.GBM, name = "DNA methylation", col=jet.colors(75), cluster_columns = FALSE,show_row_names=FALSE,show_column_names = T,top_annotation = ha1,cluster_rows=T)+
#  Heatmap(zzz, name = "ChIP-seq", col=brewer.pal(11, "PuOr"), cluster_columns = FALSE,cluster_rows = F,show_row_names=FALSE,show_column_names = TRUE,top_annotation = ha1) 
Heatmap(zzz, name = "ChIP-seq", col=colors, cluster_columns = FALSE,show_row_names=FALSE,show_column_names = TRUE,top_annotation = ha1) 


aux <- cbind(LGG.GBM,b)
a <- cbind(melt(aux[,1:10]),melt(aux[,11:20]))
colnames(a) <- c("sample","meth","a","hmec")
a$subtype <- c(rep("G-CIMP-low",36),rep("G-CIMP-high",54))

ggplot(a,aes(meth,hmec,color=subtype)) + 
  geom_point(size=2) +
  scale_color_manual(values=c("firebrick","darkgreen")) +
  theme_bw()

## Volcano
load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/cat_5hmc_v3.rda")
aux <- findOverlaps(makeGRangesFromDataFrame(cat.5hmC),makeGRangesFromDataFrame(hg38.450k))
hmec.r$probes.450k <- "No"
hmec.r[queryHits(aux),"probes.450k"] <- "Yes"
hg38.450k <- hg38.450k[as.character(hg38.450k$Composite.Element.REF) %in% as.character(GCIMPlow.probes$probes),]
aux <- findOverlaps(makeGRangesFromDataFrame(hmec.r),makeGRangesFromDataFrame(hg38.450k))
hmec.r[queryHits(aux),"probes.450k"] <- "Define.GCIMP"
hmec.r$thr.cmb <- "No difference"
hmec.r[hmec.r$thr %in% "Gains","thr.cmb"] <- "gain"
hmec.r[hmec.r$thr %in% "Losses","thr.cmb"] <- "loss"
hmec.r[hmec.r$probes.450k %in% "Define.GCIMP","thr.cmb"] <- "Define.GCIMP"
#hmec.r[hmec.r$probes.450k %in% "Yes","thr.cmb"] <- "450k"
#hmec.r[hmec.r$thr %in% "Gains" & hmec.r$probes.450k %in% "Yes","thr.cmb"] <- "Gains (450k)"
#hmec.r[hmec.r$thr %in% "Gains" & hmec.r$probes.450k %in% "Define.GCIMP","thr.cmb"] <- "Gains (Low)"
#hmec.r[hmec.r$thr %in% "Losses" & hmec.r$probes.450k %in% "Yes","thr.cmb"] <- "Losses (450k)"
#hmec.r[hmec.r$thr %in% "Losses" & hmec.r$probes.450k %in% "Define.GCIMP","thr.cmb"] <- "Losses (Low)"

ggplot(data=cat.5hmC,aes(x=Fold, y=-1*log10(FDR))) +
  #geom_point(size=2)+
  geom_point(data=subset(hmec.r,status %in% c("gain","loss","No difference")),aes(x=Fold, y=-1*log10(FDR),color=thr.cmb),size=2) +
  geom_point(data=subset(hmec.r,status %in% c("Define.GCIMP")),aes(x=Fold, y=-1*log10(FDR),color=thr.cmb),size=2) +
  #geom_point(data=subset(hmec.r,thr.cmb %in% c("Gains","Losses","No.diff","450k")),aes(x=Fold, y=-1*log10(FDR),color=thr.cmb),size=2) +
  #geom_point(data=subset(hmec.r,thr.cmb %in% c("Gains (450k)","Losses (450k)","Define.GCIMP")),aes(x=Fold, y=-1*log10(FDR),color=thr.cmb),size=2) +
  xlab("5hmC (fold-change)") + ylab("-1 * log10 of the Significance") +
  labs(title = "5hmC (hMeDIP-Seq)") +
  geom_hline(yintercept = -1*log10(0.05),linetype="dashed") + 
  scale_color_manual(values=c("gray99","yellow","dodgerblue2","gray48"),
                     labels=c("450k CpG probes that define GCIMP-low","Decreased affinity in GCIMP-low","Increased affinity in GCIMP-low","Not Significant"),name="Legend") +
 # scale_color_manual(values=c("gray99","black","dodgerblue2","lightblue","yellow","wheat","gray48"), labels=c("Not Significant (Overlaps 450k array probes)","450k array probes that define G-CIMP-low","Increased affinity in GCIMP-low","Increased affinity in GCIMP-low (Overlaps 450k array probes)","Decreased affinity in GCIMP-low","Decreased affinity in GCIMP-low (Overlaps 450k array probes)","Not Significant"),name="Legend") +
  theme_bw()


ggplot(data=subset(hmec.r,thr %in% "No.diff"),aes(thr,fill=probes.450k)) + geom_bar(aes(y = (..count..)/sum(..count..)))


ggplot(data=volcano,aes(x=DiffMean, y=-1*log10(p.value.adj))) +
  geom_point(data=subset(volcano,threshold %in% c("1","2","3")),aes(x=DiffMean, y=-1*log10(p.value.adj), colour=threshold),size=1) +
  geom_point(data=subset(volcano,threshold %in% c("4","5","6")),aes(x=DiffMean, y=-1*log10(p.value.adj), colour=threshold),size=2) 

####### Figure 2
### Panel A - Enriched Heatmap
#### 5hmC
hmec.r <- dba.report(hmec.a)
a <- hmec.r[elementMetadata(hmec.r)[,4] > 0 & elementMetadata(hmec.r)[,6] < 0.01,] # <0 gain .. >0 loss
a <- as.data.frame(a)
peaks.GR <- makeGRangesFromDataFrame(a[1:2000,])

a <- fread("/media/data2/aligned_hg38/hMEDIP/bedGraph_unzip/01_3185HF_s2_DU-6408_hMeDIP_hs_i86.noDup.sorted.q30.bedGraph",skip=1)
colnames(a) <- c("chrom","start","end","coverage")
a$coverage <- log2(a$coverage)
a <- makeGRangesFromDataFrame(a,TRUE)

DU.6408 = normalizeToMatrix(a, peaks.GR, value_column = "coverage", extend = 200, mean_mode = "w0", w = 20)
rm(a); gc();

a <- fread("/media/data2/aligned_hg38/hMEDIP/bedGraph_unzip/02_3185HF_s3_DU-6542_hMeDIP_hs_i87.noDup.sorted.q30.bedGraph",skip=1)
colnames(a) <- c("chrom","start","end","coverage")
a$coverage <- log2(a$coverage)
a <- makeGRangesFromDataFrame(a,TRUE)

DU.6542 = normalizeToMatrix(a, peaks.GR, value_column = "coverage", extend = 200, mean_mode = "w0", w = 20)
rm(a); gc();

a <- fread("/media/data2/aligned_hg38/hMEDIP/bedGraph_unzip/03_3185HF_s4_DU-7304_hMeDIP_hs_i88.noDup.sorted.q30.bedGraph",skip=1)
colnames(a) <- c("chrom","start","end","coverage")
a$coverage <- log2(a$coverage)
a <- makeGRangesFromDataFrame(a,TRUE)

DU.7304 = normalizeToMatrix(a, peaks.GR, value_column = "coverage", extend = 200, mean_mode = "w0", w = 20)
rm(a); gc();

a <- fread("/media/data2/aligned_hg38/hMEDIP/bedGraph_unzip/04_3185HF_s5_DU-7010_hMeDIP_hs_i89.noDup.sorted.q30.bedGraph",skip=1)
colnames(a) <- c("chrom","start","end","coverage")
a$coverage <- log2(a$coverage)
a <- makeGRangesFromDataFrame(a,TRUE)

DU.7010 = normalizeToMatrix(a, peaks.GR, value_column = "coverage", extend = 200, mean_mode = "w0", w = 20)
rm(a); gc();

a <- fread("/media/data2/aligned_hg38/hMEDIP/bedGraph_unzip/05_3185HF_s6_06-0129_hMeDIP_hs_i90.noDup.sorted.q30.bedGraph",skip=1)
colnames(a) <- c("chrom","start","end","coverage")
a$coverage <- log2(a$coverage)
a <- makeGRangesFromDataFrame(a,TRUE)

X6.0129 = normalizeToMatrix(a, peaks.GR, value_column = "coverage", extend = 200, mean_mode = "w0", w = 20)
rm(a); gc();

a <- fread("/media/data2/aligned_hg38/hMEDIP/bedGraph_unzip/06_3185HF_s7_06-1805_hMeDIP_hs_i91.noDup.sorted.q30.bedGraph",skip=1)
colnames(a) <- c("chrom","start","end","coverage")
a$coverage <- log2(a$coverage)
a <- makeGRangesFromDataFrame(a,TRUE)

X6.1805 = normalizeToMatrix(a, peaks.GR, value_column = "coverage", extend = 200, mean_mode = "w0", w = 20)
rm(a); gc();

a <- fread("/media/data2/aligned_hg38/hMEDIP/bedGraph_unzip/07_3185HF_s8_06-2570_hMeDIP_hs_i92.noDup.sorted.q30.bedGraph",skip=1)
colnames(a) <- c("chrom","start","end","coverage")
a$coverage <- log2(a$coverage)
a <- makeGRangesFromDataFrame(a,TRUE)

X6.2570 = normalizeToMatrix(a, peaks.GR, value_column = "coverage", extend = 200, mean_mode = "w0", w = 20)
rm(a); gc();

a <- fread("/media/data2/aligned_hg38/hMEDIP/bedGraph_unzip/08_3185HF_s9_DU-7007_hMeDIP_hs_i93.noDup.sorted.q30.bedGraph",skip=1)
colnames(a) <- c("chrom","start","end","coverage")
a$coverage <- log2(a$coverage)
a <- makeGRangesFromDataFrame(a,TRUE)

DU.7007 = normalizeToMatrix(a, peaks.GR, value_column = "coverage", extend = 200, mean_mode = "w0", w = 20)
rm(a); gc();

a <- fread("/media/data2/aligned_hg38/hMEDIP/bedGraph_unzip/09_merged.noDup.sorted.q30.bedGraph",skip=1)
colnames(a) <- c("chrom","start","end","coverage")
a$coverage <- log2(a$coverage)
a <- makeGRangesFromDataFrame(a,TRUE)

DU.7008 = normalizeToMatrix(a, peaks.GR, value_column = "coverage", extend = 200, mean_mode = "w0", w = 20)
rm(a); gc();

a <- fread("/media/data2/aligned_hg38/hMEDIP/bedGraph_unzip/10_merged.noDup.sorted.q30.bedGraph",skip=1)
colnames(a) <- c("chrom","start","end","coverage")
a$coverage <- log2(a$coverage)
a <- makeGRangesFromDataFrame(a,TRUE)

X6.0221 = normalizeToMatrix(a, peaks.GR, value_column = "coverage", extend = 200, mean_mode = "w0", w = 20)
rm(a); gc()

low = getSignalsFromList(list(X6.0129,X6.1805,X6.2570,DU.7010))
high = getSignalsFromList(list(X6.0221,DU.7008,DU.7007,DU.7304,DU.6542,DU.6408))

all.ps <- makeGRangesFromDataFrame(LGG.GBM.WGBS,TRUE)
meth.low = normalizeToMatrix(all.ps, peaks.GR, value_column = "percentMeth.0128", mean_mode = "absolute",
                             extend = 200, w = 20, empty_value = NA, smooth = TRUE)
meth.1460 = normalizeToMatrix(all.ps, peaks.GR, value_column = "percentMeth.1460", mean_mode = "absolute",
                              extend = 200, w = 20, empty_value = NA, smooth = TRUE)
meth.1788 = normalizeToMatrix(all.ps, peaks.GR, value_column = "percentMeth.1788", mean_mode = "absolute",
                              extend = 200, w = 20, empty_value = NA, smooth = TRUE)

meth.high = getSignalsFromList(list(meth.1460,meth.1788))

meth_col_fun = colorRamp2(seq(0,100,10), jet.colors(11))
col_fun = circlize::colorRamp2(seq(0,5,by=0.25), jet.colors(21))

ht_list <- 
  EnrichedHeatmap(low,column_title = "GCIMP-low - 5hmC", name = "GCIMP-low - 5hmC", #split = partition,
                  #top_annotation = HeatmapAnnotation(lines = anno_enriched(ylim=c(0,3))), 
                  top_annotation_height = unit(2, "cm") ,col = col_fun)  +
  EnrichedHeatmap(high,column_title = "GCIMP-high - 5hmC", name = "GCIMP-high - 5hmC",
                  # top_annotation = HeatmapAnnotation(lines = anno_enriched(ylim=c(0,3))), 
                  top_annotation_height = unit(2, "cm") ,col = col_fun) + 
  EnrichedHeatmap(meth.low, col = meth_col_fun, column_title = "GCIMP-low - DNA meth", name = "GCIMP-low - DNA meth", 
                  #top_annotation = HeatmapAnnotation(lines = anno_enriched(ylim=c(0,80))), 
                  top_annotation_height = unit(2, "cm")) +
  EnrichedHeatmap(meth.high, col = meth_col_fun, column_title = "GCIMP-high - DNA meth", name = "GCIMP-high - DNA meth",
                  #top_annotation = HeatmapAnnotation(lines = anno_enriched(ylim=c(0,80))), 
                  top_annotation_height = unit(2, "cm"))

draw(ht_list, main_heatmap = "GCIMP-low - 5hmC")

### ChromHMM - losses
hmec.r <- dba.report(hmec.a)
aux <- hmec.r[elementMetadata(hmec.r)[,4] > 0 & elementMetadata(hmec.r)[,6] < 0.01,]
hmec.f <- as.data.frame(aux)

hmec.f <- cat.5hmC[cat.5hmC$status %in% "loss",]
aux <- makeGRangesFromDataFrame(hmec.f)

setwd("/media/data1/thais.projects/GCIMP-low/chromHMM_Costello/samples/hg38")
for(i in list.files()){
  a <- import.bed(i)
  A <- findOverlaps(aux, a, select="first")
  tmp <- as.data.frame(a)
  tmp <- tmp[A,]
  hmec.f <- cbind(hmec.f,  tmp[,6])
  colnames(hmec.f)[ncol(hmec.f)] <- gsub("_5_segments_hg38.bed","",i)
}

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  E1 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "palevioletred1", col = NA))
  },
  E2 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "lightyellow1", col = NA))
  },
  E3 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "lightskyblue", col = NA))
  },
  E4 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "palegreen", col = NA))
  },
  E5 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "mediumpurple1", col = NA))
  },
  No = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "white", col = NA))
  }
)
col = c("E1" = "palevioletred1", "E2" = "lightyellow1", "E3" = "lightskyblue", "E4" = "palegreen", "E5" = "mediumpurple1", "No" = "white")

aux <- cat.5hmC[cat.5hmC$status %in% "loss",17:24]

#aux[] <- lapply(aux, as.numeric)
#aux <- na.omit(aux)
#a <- hclust(dist(aux[,5:8]))


aux[] <- lapply(aux, as.character)
aux[is.na(aux)] <- "No"
aux[] <- lapply(aux, as.factor)
hc.m.gain <- as.hclust(agnes(daisy(data.frame(aux), metric="gower"), method="ward"))
aux <- aux[hc.m.gain$order,]

#save(hmec.f,file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/hmec.f.order_loss_0.01.rda")
hmec.f <- aux

oncoPrint(as.matrix(hmec.f[,1:4]), get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Non-tumor brain", 
          row_title = "1,765 differentially bound peaks with decreased affinity in G-CIMP-low",column_order = NULL, show_pct = F,
          row_order = NULL, top_annotation = NULL,
          show_row_names = FALSE, show_column_names = TRUE,
          heatmap_legend_param = list(title = "ChromHMM", at = c("E1", "E2", "E3", "E4", "E5", "No"), 
                                      labels = c("Active Promoter", "Heterochromatin", "Poised Enhancer", "Active Enhancer", "Low H3K27ac","Not Available"))) +
  oncoPrint(as.matrix(hmec.f[,5:8]), get_type = function(x) strsplit(x, ";")[[1]],
            alter_fun = alter_fun, col = col, 
            column_title = "Glioblastoma IDH wildtype", column_order = NULL, show_pct = F,
            row_order = NULL, top_annotation = NULL,
            show_row_names = FALSE, show_column_names = TRUE,
            heatmap_legend_param = list(title = "ChromHMM", at = c("E1", "E2", "E3", "E4", "E5", "No"), 
                                        labels = c("Active Promoter", "Heterochromatin", "Poised Enhancer", "Active Enhancer", "Low H3K27ac","Not Available")))

Heatmap(log2(hmC[rownames(aux),c(21:12)]), name = "ChIP-seq", col=jet.colors(75), cluster_columns = FALSE,cluster_rows = FALSE,show_row_names=FALSE,show_column_names = TRUE,top_annotation = ha1) +
  Heatmap(hmC[,c(22:25)], name = "WGBS", col=jet.colors(75), cluster_columns = FALSE,show_row_names=FALSE,show_column_names = F,top_annotation = ha2) 

#### Forest plot
load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/cat_5hmc.rda")
aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference"),]
aux[!(aux$type %in% "TSS"),"type"] <- "No.TSS"
#aux <- aux[68200:68700,]
table(aux$status,aux$type,aux$CpG_type)
aux <- table(aux$status,aux$type,aux$CpG_type)
#dimnames(aux) <- list(Status = c("loss", "no.difference"),
#                      location = c("No.TSS", "TSS"),
#                      CpG.type = c("island","shore","op"))
mantelhaen.test(aux)
#max 2147483647

#odds ratio
aux <- makeGRangesFromDataFrame(cat.5hmC)
library(rtracklayer)
setwd("/media/data1/thais.projects/GCIMP-low/chromHMM_Costello/samples/hg38")
for(i in list.files()){
  a <- import.bed(i)
  A <- findOverlaps(aux, a, select="first")
  tmp <- as.data.frame(a)
  tmp <- tmp[A,]
  cat.5hmC <- cbind(cat.5hmC,  tmp[,6])
  colnames(cat.5hmC)[ncol(cat.5hmC)] <- gsub("_5_segments_hg38.bed","",i)
}

cat.5hmC$chromHMM.gbm <- "Unclassified"
cat.5hmC[rowSums(cat.5hmC[,17:24] == "E1",na.rm=T) > 2,"chromHMM.gbm"] <- "Active_promoter"
cat.5hmC[rowSums(cat.5hmC[,21:24] == "E2",na.rm=T) > 2,"chromHMM.gbm"] <- "Heterochromatin"
cat.5hmC[rowSums(cat.5hmC[,21:24] == "E3",na.rm=T) > 2,"chromHMM.gbm"] <- "Poised_enhancer"
cat.5hmC[rowSums(cat.5hmC[,21:24] == "E4",na.rm=T) > 2,"chromHMM.gbm"] <- "Active_enhancer"
cat.5hmC[rowSums(cat.5hmC[,21:24] == "E5",na.rm=T) > 2,"chromHMM.gbm"] <- "Low_H3K27ac"

cat.5hmC$chromHMM.normal <- "Unclassified"
cat.5hmC[rowSums(cat.5hmC[,17:20] == "E1",na.rm=T) > 2,"chromHMM.normal"] <- "Active_promoter"
cat.5hmC[rowSums(cat.5hmC[,17:20] == "E2",na.rm=T) > 2,"chromHMM.normal"] <- "Heterochromatin"
cat.5hmC[rowSums(cat.5hmC[,17:20] == "E3",na.rm=T) > 2,"chromHMM.normal"] <- "Poised_enhancer"
cat.5hmC[rowSums(cat.5hmC[,17:20] == "E4",na.rm=T) > 2,"chromHMM.normal"] <- "Active_enhancer"
cat.5hmC[rowSums(cat.5hmC[,17:20] == "E5",na.rm=T) > 2,"chromHMM.normal"] <- "Low_H3K27ac"

### Odds-ratio - loss
library(fmsb)
aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-island",]
aux[!(aux$type %in% "TSS"),"type"] <- "ZNo.TSS"
a <- table(aux$type,aux$status)
a <- fisher.test(a)
#a <- oddsratio(as.numeric(a[1,2]), as.numeric(a[2,2]), as.numeric(a[1,1]), as.numeric(a[2,1]), conf.level=0.95) 

df <- data.frame(Location="Promoter",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-island")

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-shore",]
aux[!(aux$type %in% "TSS"),"type"] <- "ZNo.TSS"
a <- table(aux$type,aux$status)
a <- fisher.test(a)
#a <- oddsratio(as.numeric(a[1,2]), as.numeric(a[2,2]), as.numeric(a[1,1]), as.numeric(a[2,1]), conf.level=0.95) 

aux <- data.frame(Location="Promoter",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-shore")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "Open-Sea",]
aux[!(aux$type %in% "TSS"),"type"] <- "ZNo.TSS"
a <- table(aux$type,aux$status)
a <- fisher.test(a)
#a <- oddsratio(as.numeric(a[1,2]), as.numeric(a[2,2]), as.numeric(a[1,1]), as.numeric(a[2,1]), conf.level=0.95) 

aux <- data.frame(Location="Promoter",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="Open-sea")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-island",]
aux[!(aux$type %in% "exon"),"type"] <- "No.exon"
a <- table(aux$type,aux$status)
a <- fisher.test(a)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Exon",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-island")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-shore",]
aux[!(aux$type %in% "exon"),"type"] <- "No.exon"
a <- table(aux$type,aux$status)
a <- fisher.test(a)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Exon",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-shore")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "Open-Sea",]
aux[!(aux$type %in% "exon"),"type"] <- "No.exon"
a <- table(aux$type,aux$status)
a <- fisher.test(a)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Exon",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="Open-sea")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-island",]
aux[!(aux$type %in% "intron"),"type"] <- "No.intron"
a <- table(aux$type,aux$status)
a <- fisher.test(a)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Intron",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-island")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-shore",]
aux[!(aux$type %in% "intron"),"type"] <- "No.intron"
a <- table(aux$type,aux$status)
a <- fisher.test(a)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Intron",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-shore")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "Open-Sea",]
aux[!(aux$type %in% "intron"),"type"] <- "No.intron"
a <- table(aux$type,aux$status)
a <- fisher.test(a)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Intron",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="Open-sea")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-island",]
aux[!(aux$type %in% "intergenic"),"type"] <- "No.intergenic"
a <- table(aux$type,aux$status)
a <- fisher.test(a)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Intergenic",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-island")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-shore",]
aux[!(aux$type %in% "intergenic"),"type"] <- "No.intergenic"
a <- table(aux$type,aux$status)
a <- fisher.test(a)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Intergenic",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-shore")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "Open-Sea",]
aux[!(aux$type %in% "intergenic"),"type"] <- "No.intergenic"
a <- table(aux$type,aux$status)
a <- fisher.test(a)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Intergenic",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="Open-sea")
df <- rbind(df,aux)

#levels(df$Location) <- as.character(df$Location)

#ggplot(data=df, aes(x=Location, y=Mean, ymin=Lower, ymax=Upper, color = CpG.type, label = P.value)) +
ggplot(data=df, aes(x=Location, y=Mean, ymin=Lower, ymax=Upper, color = CpG.type)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Location") + ylab("Mean (95% CI)") +
  theme_bw()  + # use a white background
  #scale_y_continuous(breaks = seq(from=ifelse(min(df$Mean)>1,1,floor(min(df$Mean))),to=ceiling(max(df$Mean)),by=1)) +
  #geom_text(vjust=0) +
  scale_color_manual(values = c("orange","darkslateblue","turquoise")) 


aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-island",]
aux[!(aux$chromHMM.gbm %in% "Active_enhancer"),"chromHMM.gbm"] <- "Not.enhancer"
#a <- table(aux$status,aux$chromHMM.gbm)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

a <- table(aux$chromHMM.gbm,aux$status)
a <- fisher.test(a)

df <- data.frame(Location="Active_enhancer",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-island")

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-shore",]
aux[!(aux$chromHMM.gbm %in% "Active_enhancer"),"chromHMM.gbm"] <- "Not.enhancer"
#a <- table(aux$status,aux$chromHMM.gbm)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

a <- table(aux$chromHMM.gbm,aux$status)
a <- fisher.test(a)

aux <- data.frame(Location="Active_enhancer",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-shore")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "Open-Sea",]
aux[!(aux$chromHMM.gbm %in% "Active_enhancer"),"chromHMM.gbm"] <- "Not.enhancer"
#a <- table(aux$status,aux$chromHMM.gbm)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

a <- table(aux$chromHMM.gbm,aux$status)
a <- fisher.test(a)

aux <- data.frame(Location="Active_enhancer",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="Open-sea")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-island",]
aux[!(aux$chromHMM.gbm %in% "Active_promoter"),"chromHMM.gbm"] <- "Not.promoter"
#a <- table(aux$status,aux$chromHMM.gbm)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

a <- table(aux$chromHMM.gbm,aux$status)
a <- fisher.test(a)

aux <- data.frame(Location="Active_promoter",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-island")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-shore",]
aux[!(aux$chromHMM.gbm %in% "Active_promoter"),"chromHMM.gbm"] <- "Not.promoter"
#a <- table(aux$status,aux$chromHMM.gbm)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

a <- table(aux$chromHMM.gbm,aux$status)
a <- fisher.test(a)

aux <- data.frame(Location="Active_promoter",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-shore")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "Open-Sea",]
aux[!(aux$chromHMM.gbm %in% "Active_promoter"),"chromHMM.gbm"] <- "Not.promoter"
#a <- table(aux$status,aux$chromHMM.gbm)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

a <- table(aux$chromHMM.gbm,aux$status)
a <- fisher.test(a)

aux <- data.frame(Location="Active_promoter",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="Open-sea")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-island",]
aux[!(aux$chromHMM.gbm %in% "Heterochromatin"),"chromHMM.gbm"] <- "Not.Heterochromatin"
#a <- table(aux$status,aux$chromHMM.gbm)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

a <- table(aux$chromHMM.gbm,aux$status)
a <- fisher.test(a)

aux <- data.frame(Location="Heterochromatin",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-island")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-shore",]
aux[!(aux$chromHMM.gbm %in% "Heterochromatin"),"chromHMM.gbm"] <- "Not.Heterochromatin"
#a <- table(aux$status,aux$chromHMM.gbm)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

a <- table(aux$chromHMM.gbm,aux$status)
a <- fisher.test(a)

aux <- data.frame(Location="Heterochromatin",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-shore")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "Open-Sea",]
aux[!(aux$chromHMM.gbm %in% "Heterochromatin"),"chromHMM.gbm"] <- "Not.Heterochromatin"
#a <- table(aux$status,aux$chromHMM.gbm)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95)

a <- table(aux$chromHMM.gbm,aux$status)
a <- fisher.test(a)

aux <- data.frame(Location="Heterochromatin",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="Open-sea")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-island",]
aux[!(aux$chromHMM.gbm %in% "Low_H3K27ac"),"chromHMM.gbm"] <- "Not.Low_H3K27ac"
#a <- table(aux$status,aux$chromHMM.gbm)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

a <- table(aux$chromHMM.gbm,aux$status)
a <- fisher.test(a)

aux <- data.frame(Location="Low_H3K27ac",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-island")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-shore",]
aux[!(aux$chromHMM.gbm %in% "Low_H3K27ac"),"chromHMM.gbm"] <- "Not.Low_H3K27ac"
#a <- table(aux$status,aux$chromHMM.gbm)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

a <- table(aux$chromHMM.gbm,aux$status)
a <- fisher.test(a)

aux <- data.frame(Location="Low_H3K27ac",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-shore")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "Open-Sea",]
aux[!(aux$chromHMM.gbm %in% "Low_H3K27ac"),"chromHMM.gbm"] <- "Not.Low_H3K27ac"
#a <- table(aux$status,aux$chromHMM.gbm)
#a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

a <- table(aux$chromHMM.gbm,aux$status)
a <- fisher.test(a)

aux <- data.frame(Location="Low_H3K27ac",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="Open-sea")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-island",]
aux[!(aux$chromHMM.gbm %in% "Poised_enhancer"),"chromHMM.gbm"] <- "ZNot.Poised_enhancer"
#a <- table(aux$status,aux$chromHMM.gbm)
#a <- oddsratio(as.numeric(a[1,2]), as.numeric(a[2,2]), as.numeric(a[1,1]), as.numeric(a[2,1]), conf.level=0.95) 

a <- table(aux$chromHMM.gbm,aux$status)
a <- fisher.test(a)

aux <- data.frame(Location="Poised_enhancer",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-island")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "CpG-shore",]
aux[!(aux$chromHMM.gbm %in% "Poised_enhancer"),"chromHMM.gbm"] <- "ZNot.Poised_enhancer"
#a <- table(aux$status,aux$chromHMM.gbm)
#a <- oddsratio(as.numeric(a[1,2]), as.numeric(a[2,2]), as.numeric(a[1,1]), as.numeric(a[2,1]), conf.level=0.95) 

a <- table(aux$chromHMM.gbm,aux$status)
a <- fisher.test(a)

aux <- data.frame(Location="Poised_enhancer",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-shore")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("loss","No difference") & cat.5hmC$CpG_type %in% "Open-Sea",]
aux[!(aux$chromHMM.gbm %in% "Poised_enhancer"),"chromHMM.gbm"] <- "ZNot.Poised_enhancer"
#a <- table(aux$status,aux$chromHMM.gbm)
#a <- oddsratio(as.numeric(a[1,2]), as.numeric(a[2,2]), as.numeric(a[1,1]), as.numeric(a[2,1]), conf.level=0.95) 

a <- table(aux$chromHMM.gbm,aux$status)
a <- fisher.test(a)

aux <- data.frame(Location="Poised_enhancer",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="Open-sea")
df <- rbind(df,aux)

df$pvalue.0.05 <- df$P.value < 0.05
df$pvalue.0.005 <- df$P.value < 0.005
df$pvalue.0.0005 <- df$P.value < 0.0005
df[df$Upper > 12,"Upper"] <- 12

ggplot(data=df, aes(x=Location, y=Mean, ymin=Lower, ymax=Upper, color = CpG.type)) +#, label = P.value)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Location") + ylab("Mean (95% CI)") +
  theme_bw()  + # use a white background
  #scale_y_continuous(breaks = seq(from=ifelse(min(df$Mean)>1,1,floor(min(df$Mean))),to=ceiling(max(df$Mean)),by=1)) +
  #geom_text(vjust=0) +
  scale_color_manual(values = c("orange","darkslateblue","turquoise")) 

# pie chart
pie(as.numeric(table(cat.5hmC$type,cat.5hmC$status)[,2]), labels = names(table(cat.5hmC$type,cat.5hmC$status)[,2]), main="Genomic Location")

pie(as.numeric(table(cat.5hmC$chromHMM.gbm,cat.5hmC$status)[,2]), labels = names(table(cat.5hmC$chromHMM.gbm,cat.5hmC$status)[,2]), main="Genomic Location",col=c("palegreen","palevioletred1","lightyellow1","mediumpurple1","lightskyblue","white"))

col = c("E1" = "palevioletred1", "E2" = "lightyellow1", "E3" = "lightskyblue", "E4" = "palegreen", "E5" = "mediumpurple1", "No" = "white")


df <- data.frame(value=as.numeric(table(cat.5hmC$chromHMM.gbm,cat.5hmC$status)[,2]), labels = names(table(cat.5hmC$chromHMM.gbm,cat.5hmC$status)[,2]))
df$value <- (df$value/sum(df$value))*100
library(scales)

# GREAT analysis
write.table(cat.5hmC[cat.5hmC$status %in% "loss" & cat.5hmC$seqnames %in% paste0("chr",c(1:22,"X","Y")),1:3],quote=F,row.names = F,col.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/loss_peaks_hg38.bed")

write.table(cat.5hmC[cat.5hmC$status %in% "No difference" & cat.5hmC$seqnames %in% paste0("chr",c(1:22,"X","Y")),1:3],quote=F,row.names = F,col.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/nodiff_peaks_hg38.bed")

# Homer analysis
write.table(cat.5hmC[cat.5hmC$status %in% "loss" & cat.5hmC$seqnames %in% paste0("chr",c(1:22,"X","Y")) & cat.5hmC$chromHMM.gbm %in% "Active_enhancer",1:3],quote=F,row.names = F,col.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/loss_peaks_act_enhancer_hg382.bed")

write.table(cat.5hmC[cat.5hmC$status %in% "No difference" & cat.5hmC$seqnames %in% paste0("chr",c(1:22,"X","Y")) & cat.5hmC$chromHMM.gbm %in% "Active_enhancer",1:3],quote=F,row.names = F,col.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/nodiff_peaks_act_enhancer_hg382.bed")

write.table(cat.5hmC[cat.5hmC$status %in% "gain" & cat.5hmC$seqnames %in% paste0("chr",c(1:22,"X","Y")) & cat.5hmC$chromHMM.gbm %in% "Active_enhancer",1:3],quote=F,row.names = F,col.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/gain_peaks_act_enhancer_hg38.bed")

write.table(cat.5hmC[cat.5hmC$status %in% "gain" & cat.5hmC$seqnames %in% paste0("chr",c(1:22,"X","Y")) & cat.5hmC$chromHMM.gbm %in% "Active_promoter",1:3],quote=F,row.names = F,col.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/gain_peaks_act_promoter_hg38.bed")



### Gains
# Fig. 4
# gene NR5A2
load("/media/data2/TCGA_RNAseq/gene_annotation_v22.rda")
genes[genes$gene_name %in% "NR5A2",]
aux <- promoters(makeGRangesFromDataFrame(genes[genes$gene_name %in% "NR5A2",]), upstream = 2000, downstream = 1000)
a <- findOverlaps(aux,all.ps)
aux <- all.ps[subjectHits(a),]
aux <- as.data.frame(aux)

ha2 = HeatmapAnnotation(df = data.frame(Sample = c("GCIMP-low","GCIMP-high","GCIMP-high","Non-tumor Brain","Non-tumor Brain")),
                        col = list(Sample = c("Non-tumor Brain" = "gray","GCIMP-low" = "darkgreen","GCIMP-high"="firebrick")))

Heatmap(aux[,c(6,8,10,12,13)], name = "WGBS", col=jet.colors(75), cluster_columns = FALSE,show_row_names=FALSE,show_column_names = F,top_annotation = ha2) 


#pie
pie(as.numeric(table(cat.5hmC$chromHMM.gbm,cat.5hmC$status)[,1]), labels = names(table(cat.5hmC$chromHMM.gbm,cat.5hmC$status)[,1]), main="Genomic Location",col=c("palegreen","palevioletred1","lightyellow1","mediumpurple1","lightskyblue","white"))

pie(as.numeric(table(cat.5hmC$type,cat.5hmC$status)[,1]), labels = names(table(cat.5hmC$type,cat.5hmC$status)[,1]), main="Genomic Location")

#NR5A2 expression
load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/LGG_GBM_nonTumorBrain_RNAseqNorm.rda")
symbol <- "NR5A2"

normal.mrna <- rna.seq.LGG.GBM.normalBrain[rownames(rna.seq.LGG.GBM.normalBrain) %in% symbol,1:5]
normal.mrna <- apply(normal.mrna,1,mean,na.rm=TRUE)
normal.mrna <- as.data.frame(normal.mrna)
colnames(normal.mrna) <- "value"
normal.mrna$type <- "Non-tumor Brain"
aux <- as.character(subset(pd, clustM.supervised2 %in% c("GCIMP-low"))$case.id)
low.mrna <- rna.seq.LGG.GBM.normalBrain[rownames(rna.seq.LGG.GBM.normalBrain) %in% symbol,colnames(rna.seq.LGG.GBM.normalBrain) %in% aux]
low.mrna <- apply(low.mrna,1,mean,na.rm=TRUE)
low.mrna <- as.data.frame(low.mrna)
colnames(low.mrna) <- "value"
low.mrna$type <- "GCIMP-low"
aux <- as.character(subset(pd, clustM.supervised2 %in% c("GCIMP-high"))$case.id)
high.mrna <- rna.seq.LGG.GBM.normalBrain[rownames(rna.seq.LGG.GBM.normalBrain) %in% symbol,colnames(rna.seq.LGG.GBM.normalBrain) %in% aux]
high.mrna <- apply(high.mrna,1,mean,na.rm=TRUE)
high.mrna <- as.data.frame(high.mrna)
colnames(high.mrna) <- "value"
high.mrna$type <- "GCIMP-high"
aux <- as.character(subset(pd, clustM.supervised2 %in% c("IDHmut-codel"))$case.id)
codel.mrna <- rna.seq.LGG.GBM.normalBrain[rownames(rna.seq.LGG.GBM.normalBrain) %in% symbol,colnames(rna.seq.LGG.GBM.normalBrain) %in% aux]
codel.mrna <- apply(codel.mrna,1,mean,na.rm=TRUE)
codel.mrna <- as.data.frame(codel.mrna)
colnames(codel.mrna) <- "value"
codel.mrna$type <- "IDHmut-codel"
aux <- as.character(subset(pd, clustM.supervised2 %in% c("Classic-like"))$case.id)
classic.mrna <- rna.seq.LGG.GBM.normalBrain[rownames(rna.seq.LGG.GBM.normalBrain) %in% symbol,colnames(rna.seq.LGG.GBM.normalBrain) %in% aux]
classic.mrna <- apply(classic.mrna,1,mean,na.rm=TRUE)
classic.mrna <- as.data.frame(classic.mrna)
colnames(classic.mrna) <- "value"
classic.mrna$type <- "Classic-like"
aux <- as.character(subset(pd, clustM.supervised2 %in% c("Mesenchymal-like"))$case.id)
mes.mrna <- rna.seq.LGG.GBM.normalBrain[rownames(rna.seq.LGG.GBM.normalBrain) %in% symbol,colnames(rna.seq.LGG.GBM.normalBrain) %in% aux]
mes.mrna <- apply(mes.mrna,1,mean,na.rm=TRUE)
mes.mrna <- as.data.frame(mes.mrna)
colnames(mes.mrna) <- "value"
mes.mrna$type <- "Mesenchymal-like"
aux <- as.character(subset(pd, clustM.supervised2 %in% c("LGm6-GBM"))$case.id)
lgm6.mrna <- rna.seq.LGG.GBM.normalBrain[rownames(rna.seq.LGG.GBM.normalBrain) %in% symbol,colnames(rna.seq.LGG.GBM.normalBrain) %in% aux]
lgm6.mrna <- apply(lgm6.mrna,1,mean,na.rm=TRUE)
lgm6.mrna <- as.data.frame(lgm6.mrna)
colnames(lgm6.mrna) <- "value"
lgm6.mrna$type <- "LGm6-GBM"
aux <- as.character(subset(pd, clustM.supervised2 %in% c("PA-like"))$case.id)
pa.mrna <- rna.seq.LGG.GBM.normalBrain[rownames(rna.seq.LGG.GBM.normalBrain) %in% symbol,colnames(rna.seq.LGG.GBM.normalBrain) %in% aux]
pa.mrna <- apply(pa.mrna,1,mean,na.rm=TRUE)
pa.mrna <- as.data.frame(pa.mrna)
colnames(pa.mrna) <- "value"
pa.mrna$type <- "PA-like"

aux <- rbind(normal.mrna,low.mrna,high.mrna,codel.mrna,classic.mrna,mes.mrna,lgm6.mrna,pa.mrna)
aux$type <- as.factor(aux$type)
aux$type <- factor(aux$type, levels = c("Non-tumor Brain","GCIMP-low","GCIMP-high","IDHmut-codel","Classic-like","Mesenchymal-like","LGm6-GBM","PA-like"))

dodge <- position_dodge(width=0.9)

ggplot(aux, aes(type, value,fill=type)) +
  geom_col(position = dodge) +
  scale_fill_manual(values = c("gray","darkgreen","firebrick","#9f20f0","#ffa500","#FFFF00","#006680","#00CCFF"),name="Glioma subtypes") +
  ylab("RNAseq values (log2 transformed") +
  xlab("Glioma subtypes") +
  theme(axis.text.x = element_blank(),strip.text = element_text(size=15),axis.text.y = element_text(size=10),legend.text=element_text(size=10),legend.title=element_text(size=10)) +
  labs(title = "PRRX1") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = dodge, width = 0.25)


### Statistical test: RNAseq data normalized by Michele Ceccarelli
load("/media/data2/TCGA_RNAseq/Michele_norm/GcimpHighAndLow_20170712_full.RData")
load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/LGG-GBM.merged.data.FB20150930v3.Rdata") #pdata from Sep 2015 (LGG-GBM project). Created by Floris Barthel

RNASeq.NORM <- RNASeq.NORM[,substr(colnames(RNASeq.NORM),14,15) %in% "01"]

aux <- as.character(subset(pd, cartoon %in% c("GCIMP-low"))$case.id)
low.mrna <- RNASeq.NORM[,substr(colnames(RNASeq.NORM),1,12) %in% aux]
aux <- as.character(subset(pd, cartoon %in% c("GCIMP-high"))$case.id)
high.mrna <- RNASeq.NORM[,substr(colnames(RNASeq.NORM),1,12) %in% aux]

RNASeq.NORM.pv <- cbind(low.mrna,high.mrna) #16 and 233
RNASeq.NORM.pv <- RNASeq.NORM.pv+1
RNASeq.NORM.pv <- log2(RNASeq.NORM.pv)

a <- apply(RNASeq.NORM.pv,1,function(probe) {
  zz <- t.test(as.matrix(probe[1:16]),as.matrix(probe[17:249]),na.action=na.omit)
  z <- zz$p.value
  return(z)
})

w.p.values.adj <- p.adjust(a, method = "fdr")
RNASeq.NORM.pv <- as.data.frame(RNASeq.NORM.pv)
RNASeq.NORM.pv$p.v <- a
RNASeq.NORM.pv$p.v.adj <- w.p.values.adj
RNASeq.NORM.pv$meanLow <- apply(RNASeq.NORM.pv[,1:16],1,mean,na.rm=T)
RNASeq.NORM.pv$meanHigh <- apply(RNASeq.NORM.pv[,17:249],1,mean,na.rm=T)
RNASeq.NORM.pv$FC <- RNASeq.NORM.pv$meanLow - RNASeq.NORM.pv$meanHigh
RNASeq.NORM.pv$status <- "No.diff"
RNASeq.NORM.pv[RNASeq.NORM.pv$p.v.adj < 0.05 & RNASeq.NORM.pv$FC < 0,"status"] <- "Downregulated"
RNASeq.NORM.pv[RNASeq.NORM.pv$p.v.adj < 0.05 & RNASeq.NORM.pv$FC > 0,"status"] <- "Upregulated"

save(RNASeq.NORM.pv,file="/media/data2/TCGA_RNAseq/Michele_norm/RNASeq_norm_test.rda")

# Heatmap de metilação e 5hmC - gains
load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/hg38_IDHmut_min5_normalBrain_IDHwt_CpGtype_mean.rda")
hmec.r <- dba.report(hmec.a, bCounts = T, th=0.01)
hmec.r <- hmec.r[mcols(hmec.r)[,4] < 0 & mcols(hmec.r)[,6] < 0.01,]#hmec.r <- hmec.r[1:1000,]
#all.ps <- makeGRangesFromDataFrame(LGG.GBM.WGBS)
aux <- findOverlaps(hmec.r,all.ps)
aux <- as.data.frame(aux)
hmC <- hmec.r
hmC <- as.data.frame(hmC)
hmC$meanLow <- NA
hmC$meanHigh1 <- NA
hmC$meanHigh2 <- NA
hmC$meanNormal <- NA
for(i in 1:nrow(hmC)){
  a <- aux[aux$queryHits == i,]
  a <- LGG.GBM.WGBS[a$subjectHits,]
  a <- apply(a[,c("percentMeth.0128","percentMeth.1460","percentMeth.1788","meanNormal")],2,mean,na.rm=T)
  hmC[i,c("meanLow","meanHigh1","meanHigh2","meanNormal")] <- a
}

a <- (as.matrix(hmC[,c("Conc_G.CIMP.low","meanLow")]))
heatpairs(as.matrix(a))

ha2 = HeatmapAnnotation(df = data.frame(Sample = c("GCIMP-low","GCIMP-high","GCIMP-high","Non-tumor Brain","Non-tumor Brain")),
                        col = list(Sample = c("Non-tumor Brain" = "gray","GCIMP-low" = "darkgreen","GCIMP-high"="firebrick")))

ha1 = HeatmapAnnotation(df = data.frame(Sample = c(rep("GCIMP-low",4),rep("GCIMP-high",6))),
                        col = list(Sample = c("GCIMP-low" = "darkgreen","GCIMP-high"="firebrick")))


colors <- c("lightpink1","lightpink2","lightpink3","mediumpurple2","mediumpurple3","mediumpurple4")
#brewer.pal(6, "PuOr")

aux <- hmC[with(hmC, order(-Fold, -FDR)),]
Heatmap(log2(hmC[,c(21:12)]), name = "ChIP-seq", col=colors, cluster_columns = FALSE,cluster_rows = FALSE,show_row_names=FALSE,show_column_names = TRUE,top_annotation = ha1) +
  Heatmap(hmC[,c(22:25)], name = "WGBS", col=jet.colors(75), cluster_columns = FALSE,show_row_names=FALSE,show_column_names = F,top_annotation = ha2) 


### Odds-ratio - gain
library(fmsb)
aux <- cat.5hmC[cat.5hmC$status %in% c("gain","No difference") & cat.5hmC$CpG_type %in% "CpG-island",]
aux[!(aux$type %in% "TSS"),"type"] <- "No.TSS"
a <- table(aux$status,aux$type)
a <- oddsratio(as.numeric(a[1,2]), as.numeric(a[2,2]), as.numeric(a[1,1]), as.numeric(a[2,1]), conf.level=0.95) 

df <- data.frame(Location="Promoter",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-island")

aux <- cat.5hmC[cat.5hmC$status %in% c("gain","No difference") & cat.5hmC$CpG_type %in% "CpG-shore",]
aux[!(aux$type %in% "TSS"),"type"] <- "No.TSS"
a <- table(aux$status,aux$type)
a <- oddsratio(as.numeric(a[1,2]), as.numeric(a[2,2]), as.numeric(a[1,1]), as.numeric(a[2,1]), conf.level=0.95) 

aux <- data.frame(Location="Promoter",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-shore")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("gain","No difference") & cat.5hmC$CpG_type %in% "Open-Sea",]
aux[!(aux$type %in% "TSS"),"type"] <- "No.TSS"
a <- table(aux$status,aux$type)
a <- oddsratio(as.numeric(a[1,2]), as.numeric(a[2,2]), as.numeric(a[1,1]), as.numeric(a[2,1]), conf.level=0.95) 

aux <- data.frame(Location="Promoter",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="Open-sea")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("gain","No difference") & cat.5hmC$CpG_type %in% "CpG-island",]
aux[!(aux$type %in% "exon"),"type"] <- "No.exon"
a <- table(aux$status,aux$type)
a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Exon",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-island")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("gain","No difference") & cat.5hmC$CpG_type %in% "CpG-shore",]
aux[!(aux$type %in% "exon"),"type"] <- "No.exon"
a <- table(aux$status,aux$type)
a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Exon",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-shore")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("gain","No difference") & cat.5hmC$CpG_type %in% "Open-Sea",]
aux[!(aux$type %in% "exon"),"type"] <- "No.exon"
a <- table(aux$status,aux$type)
a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Exon",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="Open-sea")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("gain","No difference") & cat.5hmC$CpG_type %in% "CpG-island",]
aux[!(aux$type %in% "intron"),"type"] <- "No.intron"
a <- table(aux$status,aux$type)
a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Intron",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-island")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("gain","No difference") & cat.5hmC$CpG_type %in% "CpG-shore",]
aux[!(aux$type %in% "intron"),"type"] <- "No.intron"
a <- table(aux$status,aux$type)
a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Intron",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-shore")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("gain","No difference") & cat.5hmC$CpG_type %in% "Open-Sea",]
aux[!(aux$type %in% "intron"),"type"] <- "No.intron"
a <- table(aux$status,aux$type)
a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Intron",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="Open-sea")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("gain","No difference") & cat.5hmC$CpG_type %in% "CpG-island",]
aux[!(aux$type %in% "intergenic"),"type"] <- "No.intergenic"
a <- table(aux$status,aux$type)
a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Intergenic",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-island")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("gain","No difference") & cat.5hmC$CpG_type %in% "CpG-shore",]
aux[!(aux$type %in% "intergenic"),"type"] <- "No.intergenic"
a <- table(aux$status,aux$type)
a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Intergenic",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="CpG-shore")
df <- rbind(df,aux)

aux <- cat.5hmC[cat.5hmC$status %in% c("gain","No difference") & cat.5hmC$CpG_type %in% "Open-Sea",]
aux[!(aux$type %in% "intergenic"),"type"] <- "No.intergenic"
a <- table(aux$status,aux$type)
a <- oddsratio(as.numeric(a[1,1]), as.numeric(a[1,2]), as.numeric(a[2,1]), as.numeric(a[2,2]), conf.level=0.95) 

aux <- data.frame(Location="Intergenic",Mean=a$estimate, Lower=a$conf.int[1],Upper=a$conf.int[2],P.value=a$p.value, CpG.type="Open-sea")
df <- rbind(df,aux)

#levels(df$Location) <- as.character(df$Location)

ggplot(data=df, aes(x=Location, y=Mean, ymin=Lower, ymax=Upper, color = CpG.type, label = P.value)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Location") + ylab("Mean (95% CI)") +
  theme_bw()  + # use a white background
  #scale_y_continuous(breaks = seq(from=ifelse(min(df$Mean)>1,1,floor(min(df$Mean))),to=ceiling(max(df$Mean)),by=1)) +
  geom_text(vjust=0) +
  scale_color_manual(values = c("orange","darkslateblue","turquoise")) +
  scale_y_continuous(breaks=c(0:35))


a <- read.delim("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/homer_gain_act_promoter/NR5A2_peaks.txt")
#findMotifsGenome.pl gain_peaks_act_promoter_hg38.bed hg38 homer_gain_act_promoter/ -find homer_gain_act_promoter/knownResults/known1.motif > NR5A2_peaks.txt
aux <- read.table("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/gain_peaks_act_promoter_hg38.bed")
aux <- aux[a$PositionID,]
colnames(aux) <- c("chrom","start","end")
a <- promoters(makeGRangesFromDataFrame(genes))
a <- distanceToNearest(a,makeGRangesFromDataFrame(aux))
a <- a[mcols(a)[,1] < 2000,]
a <- genes[queryHits(a),]

# GREAT analysis
write.table(cat.5hmC[cat.5hmC$status %in% "gain" & cat.5hmC$seqnames %in% paste0("chr",c(1:22,"X","Y")),1:3],quote=F,row.names = F,col.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/gain_peaks_hg38.bed")

write.table(cat.5hmC[cat.5hmC$status %in% "No difference" & cat.5hmC$seqnames %in% paste0("chr",c(1:22,"X","Y")),1:3],quote=F,row.names = F,col.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/nodiff_peaks_hg38.bed")

#### Figura 2B - Distance to TSS
aux <- genes[genes$gene_type %in% "protein_coding",]
a <- cat.5hmC[cat.5hmC$status %in% "loss" & cat.5hmC$seqnames %in% paste0("chr",c(1:22,"X","Y","M")),]
a$seqnames <- factor(a$seqnames)
tss <- promoters(makeGRangesFromDataFrame(aux), upstream = 0, downstream = 1)
aux <- precede(makeGRangesFromDataFrame(a),tss)
tss <- tss[aux,]
aux <- distance(makeGRangesFromDataFrame(a),tss)
a$precede <- aux
aux <- genes[genes$gene_type %in% "protein_coding",]
tss <- promoters(makeGRangesFromDataFrame(aux), upstream = 0, downstream = 1)
aux <- follow(makeGRangesFromDataFrame(a),tss)
tss <- tss[aux,]
aux <- distance(makeGRangesFromDataFrame(a),tss)
a$follow <- aux

aux <- genes[genes$gene_type %in% "protein_coding",]
tss <- promoters(makeGRangesFromDataFrame(aux), upstream = 0, downstream = 1)
aux <- findOverlaps(makeGRangesFromDataFrame(a),tss)
a$overlap <- FALSE
a[queryHits(aux),"overlap"] <- TRUE
a$distance <- NA
a$distance <- ifelse(a$precede < a$follow,-1*a$precede,a$follow)
a[a$overlap,"distance"] <- 0

ggplot(a,aes(distance)) + geom_histogram(breaks = c(-550000,-500000,-50000,-5000,0,5000,50000,500000)) 
table(a$distance < -500000)/nrow(a) #?
ggplot(a,aes(distance)) + geom_bar(aes(y = (..count..)/sum(..count..))) #+ scale_y_continuous(labels = scales::percent) 



###### Hi-C data - /media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/HiC
library(data.table)
a <- fread("GSM2176968_GB182_A005_r1_3.hicup.HiC_Interactions.txt")
a <- a[,-1]
colnames(a) <- c("chrom1","start1","strand1","chrom2","start2","strand2")
a <- a[a$chrom1 %in% paste0("chr",c(1:22,"X","Y")),]
a <- a[a$chrom2 %in% paste0("chr",c(1:22,"X","Y")),]
a$chrom1 <- factor(a$chrom1, levels = paste0("chr",c(1:22,"X","Y")))
a$chrom2 <- factor(a$chrom2, levels = paste0("chr",c(1:22,"X","Y")))

aux <- a[a$chrom1 %in% c("chr21","chrY"),]
ggplot(aux, aes(start1,start2)) + geom_point() + facet_grid(chrom1 ~ chrom2)


## Gene
aux <- cat.5hmC[cat.5hmC$status %in% "loss" & cat.5hmC$type %in% "intron" & cat.5hmC$CpG_type %in% "CpG-island" ,]
aux <- cat.5hmC[cat.5hmC$status %in% "loss" & cat.5hmC$type %in% "intron" & cat.5hmC$CpG_type %in% "CpG-shore" ,]
aux <- cat.5hmC[cat.5hmC$status %in% "loss" & cat.5hmC$type %in% "intron" & cat.5hmC$CpG_type %in% "Open-Sea" ,]
#aux <- cat.5hmC[cat.5hmC$status %in% "loss" & cat.5hmC$type %in% "TSS" ,]
aux <- as.character(aux$geneID)
aux <- genes[as.character(genes$gene_id) %in% unique(aux),]
table(aux$RNAseq.hg38.summary)/nrow(aux)
table(genes$RNAseq.hg38.summary)/nrow(genes)

## Cistrome - ATF3
a <- read.table("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/Gargiulo_2013_gbm_ATF3_peaks.bed")
colnames(a) <- c("chrom","start","end","peakID","score")

hmec.r <- cat.5hmC[cat.5hmC$status %in% "loss",]
aux <- findOverlaps(makeGRangesFromDataFrame(hmec.r),makeGRangesFromDataFrame(a))
hmec.r <- hmec.r[queryHits(aux),] #16%

ids <- list(hmec.r[hmec.r$type %in% "TSS","geneID"])

## Cistrome - JUN
a <- read.table("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/Yan_2013_colon_JUN.bed")
colnames(a) <- c("chrom","start","end","peakID","score")

hmec.r <- cat.5hmC[cat.5hmC$status %in% "loss",]
aux <- findOverlaps(makeGRangesFromDataFrame(hmec.r),makeGRangesFromDataFrame(a))
hmec.r <- hmec.r[queryHits(aux),] #11%

ids <- append(ids,list(hmec.r[hmec.r$type %in% "TSS","geneID"]))

## Cistrome - FOSL1
a <- read.table("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/Gertz_2013_colon_FOSL1.bed")
colnames(a) <- c("chrom","start","end","peakID","score")

hmec.r <- cat.5hmC[cat.5hmC$status %in% "loss",]
aux <- findOverlaps(makeGRangesFromDataFrame(hmec.r),makeGRangesFromDataFrame(a))
hmec.r <- hmec.r[queryHits(aux),] #22%
ids <- append(ids,list(hmec.r[hmec.r$type %in% "TSS","geneID"]))

prom.enhan <- hmec.r[hmec.r$chromHMM.gbm %in% "Active_promoter" & hmec.r$type %in% "TSS",]

tss <- promoters(makeGRangesFromDataFrame(genes), upstream = 0, downstream = 1)

a <- suppressWarnings(distanceToNearest(makeGRangesFromDataFrame(prom.enhan),tss,ignore.strand=T))
#range(mcols(a)[,1]) - #less than 2000
aux <- genes[subjectHits(a),]
write.table(as.character(aux[aux$RNAseq.hg38.summary %in% "Downregulated","gene_id"]),quote=F,row.names = F,col.names = F,sep=" ",file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/FOSL1_TSS_downregulatedGenes2.txt")
write.table(as.character(aux[aux$RNAseq.hg38.summary %in% "Upregulated","gene_id"]),quote=F,row.names = F,col.names = F,sep=" ",file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/FOSL1_TSS_upregulatedGenes2.txt")


aux <- cat.5hmC[cat.5hmC$status %in% "loss" & cat.5hmC$type %in% "TSS","geneID"]
aux <- genes[genes$gene_id %in% aux,]
write.table(as.character(aux[aux$RNAseq.hg38.summary %in% "Upregulated","gene_id"]),quote=F,row.names = F,col.names = F,sep=" ",file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/Loss_5hmC_upregulatedGenes_TSS.txt")


a <-
  matrix(c(16, 1, 2058, 988),
         nrow = 2,
         dimnames =
           list(c("TSS", "Non.TSS"),
                c("Lost", "No.diff")))

a <-
  matrix(c(16, 2058, 1, 988),
         nrow = 2,
         dimnames =
           list(c("Lost", "No.diff"),
                c("TSS", "No.TSS")))
### HiC
setwd("/media/data2/HiC/")
a <- fread("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/HiC/GSM2176962_AA86_A002_r1_3.hicup.HiC_Interactions.txt")
a$end <- a$V3+1
a$score <- 0
a$end2 <- a$V6+1
write.table(a[,c(2,3,8,1,9,4)],quote=F,row.names=F,col.names=F,sep="\t",file="GSM2176962_AA86_first_hg19.bed")
write.table(a[,c(5,6,10,1,9,7)],quote=F,row.names=F,col.names=F,sep="\t",file="GSM2176962_AA86_sec_hg19.bed")

b <- fread("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/HiC/GSM2176966_GB176_A006_r1_3.hicup.HiC_Interactions.txt")
b$end <- b$V3+1
b$score <- 0
b$end2 <- b$V6+1
write.table(b[,c(2,3,8,1,9,4)],quote=F,row.names=F,col.names=F,sep="\t",file="GSM2176966_GB176_first_hg19.bed")
write.table(b[,c(5,6,10,1,9,7)],quote=F,row.names=F,col.names=F,sep="\t",file="GSM2176966_GB176_sec_hg19.bed")

c <- fread("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/HiC/GSM2176967_GB180_A012_r1_3.hicup.HiC_Interactions.txt")
c$end <- c$V3+1
c$score <- 0
c$end2 <- c$V6+1
write.table(c[,c(2,3,8,1,9,4)],quote=F,row.names=F,col.names=F,sep="\t",file="GSM2176967_GB180_first_hg19.bed")
write.table(c[,c(5,6,10,1,9,7)],quote=F,row.names=F,col.names=F,sep="\t",file="GSM2176967_GB180_sec_hg19.bed")

d <- fread("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/HiC/GSM2176968_GB182_A005_r1_3.hicup.HiC_Interactions.txt")
d$end <- d$V3+1
d$score <- 0
d$end2 <- d$V6+1
write.table(d[,c(2,3,8,1,9,4)],quote=F,row.names=F,col.names=F,sep="\t",file="GSM2176968_GB182_first_hg19.bed")
write.table(d[,c(5,6,10,1,9,7)],quote=F,row.names=F,col.names=F,sep="\t",file="GSM2176968_GB182_sec_hg19.bed")

e <- fread("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/HiC/GSM2176969_GB183_A019_r1_3.hicup.HiC_Interactions.txt")
e$end <- e$V3+1
e$score <- 0
e$end2 <- e$V6+1
write.table(e[,c(2,3,8,1,9,4)],quote=F,row.names=F,col.names=F,sep="\t",file="GSM2176969_GB183_first_hg19.bed")
write.table(e[,c(5,6,10,1,9,7)],quote=F,row.names=F,col.names=F,sep="\t",file="GSM2176969_GB183_sec_hg19.bed")

f <- fread("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/HiC/GSM2176970_GB238_A007_r1_3.hicup.HiC_Interactions.txt")
f$end <- f$V3+1
f$score <- 0
f$end2 <- f$V6+1
write.table(f[,c(2,3,8,1,9,4)],quote=F,row.names=F,col.names=F,sep="\t",file="GSM2176970_GB238_first_hg19.bed")
write.table(f[,c(5,6,10,1,9,7)],quote=F,row.names=F,col.names=F,sep="\t",file="GSM2176970_GB238_sec_hg19.bed")

#liftOver

a <- fread("GSM2176962_AA86_first_hg38.bed")
aux <- fread("GSM2176962_AA86_sec_hg38.bed")
a <- merge(a,aux,by="V4")
a <- a[,c(1,2,3,6,7,8,11)]
colnames(a) <- c("id","chrom","start","strand","int.chrom","int.start","int.strand")

b <- fread("GSM2176966_GB176_first_hg38.bed")
aux <- fread("GSM2176966_GB176_sec_hg38.bed")
b <- merge(b,aux,by="V4")
b <- b[,c(1,2,3,6,7,8,11)]
colnames(b) <- c("id","chrom","start","strand","int.chrom","int.start","int.strand")

c <- fread("GSM2176967_GB180_first_hg38.bed")
aux <- fread("GSM2176967_GB180_sec_hg38.bed")
c <- merge(c,aux,by="V4")
c <- c[,c(1,2,3,6,7,8,11)]
colnames(c) <- c("id","chrom","start","strand","int.chrom","int.start","int.strand")

d <- fread("GSM2176968_GB182_first_hg38.bed")
aux <- fread("GSM2176968_GB182_sec_hg38.bed")
d <- merge(d,aux,by="V4")
d <- d[,c(1,2,3,6,7,8,11)]
colnames(d) <- c("id","chrom","start","strand","int.chrom","int.start","int.strand")

e <- fread("GSM2176969_GB183_first_hg38.bed")
aux <- fread("GSM2176969_GB183_sec_hg38.bed")
e <- merge(e,aux,by="V4")
e <- e[,c(1,2,3,6,7,8,11)]
colnames(e) <- c("id","chrom","start","strand","int.chrom","int.start","int.strand")

f <- fread("GSM2176970_GB238_first_hg38.bed")
aux <- fread("GSM2176970_GB238_first_hg38.bed")
f <- merge(f,aux,by="V4")
f <- f[,c(1,2,3,6,7,8,11)]
colnames(f) <- c("id","chrom","start","strand","int.chrom","int.start","int.strand")

HiC.common <- Reduce(subsetByOverlaps, list(makeGRangesFromDataFrame(a,end="start",TRUE), makeGRangesFromDataFrame(b,end="start",TRUE), makeGRangesFromDataFrame(c,end="start",TRUE),makeGRangesFromDataFrame(d,end="start",TRUE),makeGRangesFromDataFrame(e,end="start",TRUE),makeGRangesFromDataFrame(f,end="start",TRUE)))

HiC.common <- HiC.common[seqnames(HiC.common) %in% paste0("chr",c(1:22,"X","Y"))]
HiC.common <- HiC.common[mcols(HiC.common)[,2] %in% paste0("chr",c(1:22,"X","Y"))]

save(HiC.common,file="/media/data2/HiC/HiC_commomRegions.rda")


### Enhancer and HiC
loss <- cat.5hmC[cat.5hmC$status %in% "loss",]
rownames(loss) <- NULL
aux <- findOverlaps(makeGRangesFromDataFrame(loss),HiC.common)
loss <- loss[unique(queryHits(aux)),]
enh.act <- enhancer[enhancer$status %in% "gain",]
#aux <- findOverlaps(makeGRangesFromDataFrame(loss),makeGRangesFromDataFrame(enh.act))
loss <- loss[loss$chromHMM.gbm %in% "Active_enhancer" & loss$chromHMM.normal != "Active_enhancer",]
aux <- as.data.frame(aux)
aux <- aux[aux$queryHits %in% rownames(loss),]
HiC.sub <- HiC.common[aux$subjectHits,]
HiC.sub <- as.data.frame(HiC.sub)
HiC.int <- makeGRangesFromDataFrame(HiC.sub,seqnames.field = "int.chrom",start.field = "int.start",end.field = "int.start",strand.field = "int.strand")
aux <- findOverlaps(makeGRangesFromDataFrame(genes),HiC.int)
g <- genes[unique(queryHits(aux)),]

write.table(cat.5hmC[cat.5hmC$chromHMM.gbm %in% "Active_enhancer" & cat.5hmC$status %in% "loss",1:3],row.names = F,col.names = F,quote=F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/5hmC_enhancer_loss.bed")

### Average profile plot
txdb <- makeTxDbFromGFF("/media/data2/TCGA_RNAseq/GDC_gencodeV22/gencode.v22.annotation.gtf")
hmec.r <- dba.report(hmec.a, th=1)
mcols(hmec.r) <- mcols(hmec.r[,2])

plotAvgProf2(a, TxDb=txdb, upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
files <- list.files("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/macs2_calls_2.0",pattern = "*woDup_peaks.narrowPeak")
tagMatrixList <- list()
for(i in files){
  a <- read.table(paste("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/macs2_calls_2.0",i,sep="/"))
  a <- a[,1:5]
  colnames(a) <- c("chrom","start","end","name","score")
  a <- makeGRangesFromDataFrame(a,TRUE)
  tagMatrix <- getTagMatrix(a, windows=promoter)
  tagMatrixList <- append(tagMatrixList,list(tagMatrix))
}
names(tagMatrixList) <- c("TCGA-DU-6408","TCGA-DU-6542","TCGA-DU-7304","TCGA-DU-7010","TCGA-06-0129","TCGA-06-1805","TCGA-06-2570","TCGA-DU-7007","TCGA-DU-7008","TCGA-06-0221")
plotAvgProf(tagMatrixList, xlim=c(-2000, 2000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)

##DiffBind
hmec.r <- dba.report(hmec.a, th=1)
hmec.r <- as.data.frame(hmec.r)
hmec.r$high <- 2^(hmec.r$Conc_G.CIMP.high)
hmec.r$low <- 2^(hmec.r$Conc_G.CIMP.low)
a <- makeGRangesFromDataFrame(hmec.r[,c(1:3,12)],TRUE)
tagMatrix <- getTagMatrix(a, windows=promoter)
tagMatrixList <- append(tagMatrixList,list(tagMatrix))
a <- makeGRangesFromDataFrame(hmec.r[,c(1:3)],TRUE)
tagMatrix <- getTagMatrix(a, windows=promoter)
tagMatrixList <- append(tagMatrixList,list(tagMatrix))
names(tagMatrixList) <- c("GCIMP-high","GCIMP-low")
plotAvgProf(tagMatrix, xlim=c(-2000, 2000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")


tss <- promoters(makeGRangesFromDataFrame(genes), upstream = 2000, downstream = 2000)
hmec.r <- dba.report(hmec.a)
aux <- findOverlaps(hmec.r,tss)
hmec.tss <- hmec.r[unique(queryHits(aux)),]
hmec.tss <- as.data.frame(hmec.tss)
ggplot(hmec.tss,aes())

########## Average plot HOMER
### All TSS
tss <- promoters(makeGRangesFromDataFrame(genes), upstream = 0, downstream = 1)
tss <- as.data.frame(tss)
tss$id <- rownames(tss)
write.table(tss[,c(1,2,3,6,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss.bed",sep="\t")

#bed2pos.pl tss.bed > tss_pos.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss_pos.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/output_hMEDIP.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss_pos.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K27ac_06-0221_woDup/ GCIMPhigh_H3K27ac_DU-6408_woDup/ GCIMPhigh_H3K27ac_DU-6542_woDup/ GCIMPhigh_H3K27ac_DU-7007_woDup/ GCIMPhigh_H3K27ac_DU-7008_woDup/ GCIMPhigh_H3K27ac_DU-7304_woDup/ GCIMPlow_H3K27ac_06-0129_woDup/ GCIMPlow_H3K27ac_06-1805_woDup/ GCIMPlow_H3K27ac_06-2570_woDup/ GCIMPlow_H3K27ac_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/output_H3K27ac.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K4me3/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss_pos.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K4me3_06-0221_woDup/ GCIMPhigh_H3K4me3_DU-6408_woDup/ GCIMPhigh_H3K4me3_DU-6542_woDup/ GCIMPhigh_H3K4me3_DU-7007_woDup/ GCIMPhigh_H3K4me3_DU-7008_woDup/ GCIMPhigh_H3K4me3_DU-7304_woDup/ GCIMPlow_H3K4me3_06-0129_woDup/ GCIMPlow_H3K4me3_06-1805_woDup/ GCIMPlow_H3K4me3_06-2570_woDup/ GCIMPlow_H3K4me3_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/output_H3K4me3.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/CTCF/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss_pos.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_CTCF_06-0221_woDup/ GCIMPhigh_CTCF_DU-6408_woDup/ GCIMPhigh_CTCF_DU-6542_woDup/ GCIMPhigh_CTCF_DU-7007_woDup/ GCIMPhigh_CTCF_DU-7008_woDup/ GCIMPhigh_CTCF_DU-7304_woDup/ GCIMPlow_CTCF_06-0129_woDup/ GCIMPlow_CTCF_06-1805_woDup/ GCIMPlow_CTCF_06-2570_woDup/ GCIMPlow_CTCF_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/output_CTCF.txt

library(reshape2)
library(ggplot2)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/output_hMEDIP.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
aa$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
aa <- aa[aa$var %in% "cov",]
aa$data <- "5hmC"

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/output_H3K27ac.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "H3K27ac"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/output_CTCF.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
aa <- aa[aa$var %in% "cov",]
bb$data <- "CTCF"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/output_H3K4me3.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
aa <- aa[aa$var %in% "cov",]
bb$data <- "H3K4me3"

aa <- rbind(aa,bb)

ggplot(bb, aes(x = Position, y = value, color = type)) + geom_smooth(span=0.1,method="loess",aes(fill=type),alpha=0.1) + facet_wrap(~data,scales="free") + scale_color_manual(values=c("firebrick","darkgreen"),                                                                                        labels=c("GCIMP-high","GCIMP-low"), name="Legend") +
  scale_fill_manual(values=c("firebrick","darkgreen"))



### WGBS - Real CTCF
#cd /media/data1/thais.projects/GCIMP-low/CTCF_Flavahan
#bed2pos.pl CTCF_predicted_hg38.bed > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/WGBS/Real.CTCF/CTCF_predicted_hg38.txt
#cd /media/data1/thais.projects/GCIMP-low/GBM.WGBS/Methylation_-_Bisulfite_Sequencing/JHU_USC__IlluminaHiSeq_WGBS/Level_3/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/WGBS/Real.CTCF/CTCF_predicted_hg38.txt hg38 -size 4000 -hist 10 -d TCGA.06.0128.01A/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/WGBS/Real.CTCF/output_WGBS_All_CTCF.txt


### 5hmC - Gains
write.table(cat.5hmC[cat.5hmC$status %in% "gain",c(1,2,3,15,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/all_gains.bed",sep="\t")

#bed2pos.pl all_gains.bed > all_gains.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/all_gains.bed hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_hMEDIP.txt

### 5hmC - Losses
#All
write.table(cat.5hmC[cat.5hmC$status %in% "loss",c(1,2,3,15,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/all_loss.bed",sep="\t")

#bed2pos.pl all_loss.bed > all_loss.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/all_loss.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/output_hMEDIP.txt

#TSS
write.table(cat.5hmC[cat.5hmC$status %in% "loss" & cat.5hmC$type %in% "TSS",c(1,2,3,15,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/tss_loss.bed",sep="\t")

#bed2pos.pl tss_loss.bed > tss_loss.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/tss_loss.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/output_hMEDIP_tss.txt

#intron
write.table(cat.5hmC[cat.5hmC$status %in% "loss" & cat.5hmC$type %in% "intron",c(1,2,3,15,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/intron_loss.bed",sep="\t")

#bed2pos.pl intron_loss.bed > intron_loss.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/intron_loss.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/output_hMEDIP_intron.txt

#exon
write.table(cat.5hmC[cat.5hmC$status %in% "loss" & cat.5hmC$type %in% "exon",c(1,2,3,15,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/exon_loss.bed",sep="\t")

#bed2pos.pl exon_loss.bed > exon_loss.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/exon_loss.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/output_hMEDIP_exon.txt

#intergenic
write.table(cat.5hmC[cat.5hmC$status %in% "loss" & cat.5hmC$type %in% "intergenic",c(1,2,3,15,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/intergenic_loss.bed",sep="\t")

#bed2pos.pl intergenic_loss.bed > intergenic_loss.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/intergenic_loss.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/output_hMEDIP_intergenic.txt

#Loss - Real TSS
aux <- cat.5hmC[cat.5hmC$status %in% "loss" & cat.5hmC$type %in% "TSS","geneID"]
g <- genes[genes$gene_id %in% unique(as.character(aux)),]
tss <- promoters(makeGRangesFromDataFrame(g), upstream = 0, downstream = 1)
tss <- as.data.frame(tss)
tss$id <- rownames(tss)
write.table(tss[,c(1,2,3,6,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/tss_loss_bygene.bed",sep="\t")

#No difference - All
write.table(cat.5hmC[cat.5hmC$status %in% "No difference" & cat.5hmC$FDR > 0.1,c(1,2,3,15,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/no.diff/all_nodiff.bed",sep="\t")

#bed2pos.pl all_nodiff.bed > all_nodiff.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/no.diff/all_nodiff.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/no.diff/output_hMEDIP_nodiff.txt


#bed2pos.pl tss_loss_bygene.bed > tss_loss_bygene.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/tss_loss_bygene.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/output_hMEDIP_TSS_gene.txt


### 5hmC - Gains
#All
write.table(cat.5hmC[cat.5hmC$status %in% "gain",c(1,2,3,15,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/all_gains.bed",sep="\t")

#bed2pos.pl all_gains.bed > all_gains.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/all_gains.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_hMEDIP.txt

#TSS
write.table(cat.5hmC[cat.5hmC$status %in% "gain" & cat.5hmC$type %in% "TSS",c(1,2,3,15,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/tss_gain.bed",sep="\t")

#bed2pos.pl tss_gain.bed > tss_gain.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/tss_gain.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_hMEDIP_tss.txt

#intron
write.table(cat.5hmC[cat.5hmC$status %in% "gain" & cat.5hmC$type %in% "intron",c(1,2,3,15,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/intron_gain.bed",sep="\t")

#bed2pos.pl intron_gain.bed > intron_gain.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/intron_gain.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_hMEDIP_intron.txt

#exon
write.table(cat.5hmC[cat.5hmC$status %in% "gain" & cat.5hmC$type %in% "exon",c(1,2,3,15,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/exon_gain.bed",sep="\t")

#bed2pos.pl exon_gain.bed > exon_gain.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/exon_gain.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_hMEDIP_exon.txt

#intergenic
write.table(cat.5hmC[cat.5hmC$status %in% "gain" & cat.5hmC$type %in% "intergenic",c(1,2,3,15,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/intergenic_gain.bed",sep="\t")

#bed2pos.pl intergenic_gain.bed > intergenic_gain.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/intergenic_gain.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_hMEDIP_intergenic.txt

#Gain - Real TSS
aux <- cat.5hmC[cat.5hmC$status %in% "gain" & cat.5hmC$type %in% "TSS","geneID"]
g <- genes[genes$gene_id %in% unique(as.character(aux)),]
tss <- promoters(makeGRangesFromDataFrame(g), upstream = 0, downstream = 1)
tss <- as.data.frame(tss)
tss$id <- rownames(tss)
write.table(tss[,c(1,2,3,6,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/tss_gain_bygene.bed",sep="\t")

#bed2pos.pl tss_gain_bygene.bed > tss_gain_bygene.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/tss_gain_bygene.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_hMEDIP_TSS_gene.txt


#WGBS
#cd /media/data1/thais.projects/GCIMP-low/GBM.WGBS/Methylation_-_Bisulfite_Sequencing/JHU_USC__IlluminaHiSeq_WGBS/Level_3/homer/min5reads
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/WGBS/Real.CTCF/CTCF_predicted_hg38.txt hg38 -size 2000 -hist 1 -d TCGA.06.0128.01A/ TCGA.14.1454.01A/ TCGA.14.3477.01A/ TCGA.14.1401.01A/ TCGA.16.1460.01A/ TCGA.19.1788.01A/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/WGBS/Real.CTCF/output_WGBS_RealCTCF_min5reads.txt

#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/all_gains.txt hg38 -size 2000 -hist 1 -d TCGA.06.0128.01A/ TCGA.14.1454.01A/ TCGA.14.3477.01A/ TCGA.14.1401.01A/ TCGA.16.1460.01A/ TCGA.19.1788.01A/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_WGBS_5hmC_gains_min5reads.txt

#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/all_loss.txt hg38 -size 2000 -hist 1 -d TCGA.06.0128.01A/ TCGA.14.1454.01A/ TCGA.14.3477.01A/ TCGA.14.1401.01A/ TCGA.16.1460.01A/ TCGA.19.1788.01A/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/output_WGBS_5hmC_loss_min5reads.txt

library(reshape2)
#a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/WGBS/Real.CTCF/output_WGBS_RealCTCF_min5reads.txt")
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_WGBS_5hmC_gains_min5reads.txt")
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/output_WGBS_5hmC_loss_min5reads.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
#aa$type <- c(rep("low", 3*201), rep("high", 2*3*201))
aa$type <- c(rep("low", 3*2001),rep("mesenchymal", 3*2001), rep("LGm6-GBM", 3*2001),rep("mesenchymal", 3*2001),rep("high", 2*3*2001))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 2001), rep("tag5", 2001), rep("tag3", 2001)), 6))
aa <- aa[aa$var %in% "cov",]
aa$data <- "WGBS.5hmC.gains"

ggplot(aa, aes(x = Position, y = value, color = variable))  + 
  geom_smooth(aes(fill=type),alpha=0.1,span=0.1,method="loess") + scale_color_manual(values=c("darkgreen","yellow","blue","yellow","firebrick","firebrick")) + scale_x_continuous(breaks=arithSeq(-1000,1000))

arithSeq <- function(base, max){
  i <- base
  aux <- base
  while (i < (max)){
    i <-i + 40
    aux <- c(aux,i)
  }
  return(aux)
}
                                                                                      
#labels=c("GCIMP-high","GCIMP-low"), name="Legend") +
#scale_fill_manual(values=c("firebrick","darkgreen"))

#5hmC - Loss
library(reshape2)
library(ggplot2)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/output_hMEDIP.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
aa$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
aa <- aa[aa$var %in% "cov",]
aa$data <- "5hmC.all"

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/output_hMEDIP_exon.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "5hmC.Exon"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/output_hMEDIP_intron.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "5hmC.Intron"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/output_hMEDIP_tss.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "5hmC.TSS"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/output_hMEDIP_TSS_gene.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "5hmC.TSS.byGene"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/output_hMEDIP_intergenic.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "5hmC.Intergenic"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/no.diff/output_hMEDIP_nodiff.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "5hmC.NoDiff"

aa <- rbind(aa,bb)

ggplot(aa, aes(x = Position, y = value, color = type)) + geom_smooth(span=0.1,method="loess",aes(fill=type),alpha=0.1) + facet_wrap(~data,scales="free") + scale_color_manual(values=c("firebrick","darkgreen"),                                                                                        labels=c("GCIMP-high","GCIMP-low"), name="Legend") +
  scale_fill_manual(values=c("firebrick","darkgreen"))

#5hmC - Gain
library(reshape2)
library(ggplot2)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_hMEDIP.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
aa$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
aa <- aa[aa$var %in% "cov",]
aa$data <- "5hmC.all"

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_hMEDIP_exon.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "5hmC.Exon"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_hMEDIP_intron.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "5hmC.Intron"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_hMEDIP_tss.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "5hmC.TSS"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_hMEDIP_TSS_gene.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "5hmC.TSS.byGene"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_hMEDIP_intergenic.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "5hmC.Intergenic"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/no.diff/output_hMEDIP_nodiff.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "5hmC.NoDiff"

aa <- rbind(aa,bb)

ggplot(aa, aes(x = Position, y = value, color = type)) + geom_smooth(span=0.1,method="loess",aes(fill=type),alpha=0.1) + facet_wrap(~data,scales="free") + scale_color_manual(values=c("firebrick","darkgreen"),                                                                                        labels=c("GCIMP-high","GCIMP-low"), name="Legend") +
  scale_fill_manual(values=c("firebrick","darkgreen"))

### H3K27ac
h3k27ac.r <- dba.report(h3k27ac.a)
h3k27ac.r <- h3k27ac.r[mcols(h3k27ac.r)[,4] > 0,]
h3k27ac.r <- as.data.frame(h3k27ac.r)
h3k27ac.r$id <- 1:nrow(h3k27ac.r)

write.table(h3k27ac.r[,c(1,2,3,12,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/loss/all_loss.bed",sep="\t")

#bed2pos.pl all_loss.bed > all_loss.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/loss/all_loss.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K27ac_06-0221_woDup/ GCIMPhigh_H3K27ac_DU-6408_woDup/ GCIMPhigh_H3K27ac_DU-6542_woDup/ GCIMPhigh_H3K27ac_DU-7007_woDup/ GCIMPhigh_H3K27ac_DU-7008_woDup/ GCIMPhigh_H3K27ac_DU-7304_woDup/ GCIMPlow_H3K27ac_06-0129_woDup/ GCIMPlow_H3K27ac_06-1805_woDup/ GCIMPlow_H3K27ac_06-2570_woDup/ GCIMPlow_H3K27ac_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/loss/output_H3k27ac_all.txt

h3k27ac.r <- dba.report(h3k27ac.a)
h3k27ac.r <- h3k27ac.r[mcols(h3k27ac.r)[,4] < 0,]
h3k27ac.r <- as.data.frame(h3k27ac.r)
h3k27ac.r$id <- 1:nrow(h3k27ac.r)

write.table(h3k27ac.r[,c(1,2,3,12,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/all_gain.bed",sep="\t")

#bed2pos.pl all_gain.bed > all_gain.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/all_gain.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K27ac_06-0221_woDup/ GCIMPhigh_H3K27ac_DU-6408_woDup/ GCIMPhigh_H3K27ac_DU-6542_woDup/ GCIMPhigh_H3K27ac_DU-7007_woDup/ GCIMPhigh_H3K27ac_DU-7008_woDup/ GCIMPhigh_H3K27ac_DU-7304_woDup/ GCIMPlow_H3K27ac_06-0129_woDup/ GCIMPlow_H3K27ac_06-1805_woDup/ GCIMPlow_H3K27ac_06-2570_woDup/ GCIMPlow_H3K27ac_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/output_H3k27ac_all.txt

###H3K4me3
h3k4me3.r <- dba.report(h3k4me3.a)
h3k4me3.r <- h3k4me3.r[mcols(h3k4me3.r)[,4] > 0,]
h3k4me3.r <- as.data.frame(h3k4me3.r)
h3k4me3.r$id <- 1:nrow(h3k4me3.r)

write.table(h3k4me3.r[,c(1,2,3,12,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K4me3/loss/all_loss.bed",sep="\t")

#bed2pos.pl all_loss.bed > all_loss.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K4me3/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K4me3/loss/all_loss.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K4me3_06-0221_woDup/ GCIMPhigh_H3K4me3_DU-6408_woDup/ GCIMPhigh_H3K4me3_DU-6542_woDup/ GCIMPhigh_H3K4me3_DU-7007_woDup/ GCIMPhigh_H3K4me3_DU-7008_woDup/ GCIMPhigh_H3K4me3_DU-7304_woDup/ GCIMPlow_H3K4me3_06-0129_woDup/ GCIMPlow_H3K4me3_06-1805_woDup/ GCIMPlow_H3K4me3_06-2570_woDup/ GCIMPlow_H3K4me3_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K4me3/loss/output_H3K4me3_all.txt

h3k4me3.r <- dba.report(h3k4me3.a)
h3k4me3.r <- h3k4me3.r[mcols(h3k4me3.r)[,4] < 0,]
h3k4me3.r <- as.data.frame(h3k4me3.r)
h3k4me3.r$id <- 1:nrow(h3k4me3.r)

write.table(h3k4me3.r[,c(1,2,3,12,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K4me3/gain/all_gain.bed",sep="\t")

#bed2pos.pl all_gain.bed > all_gain.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K4me3/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K4me3/gain/all_gain.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K4me3_06-0221_woDup/ GCIMPhigh_H3K4me3_DU-6408_woDup/ GCIMPhigh_H3K4me3_DU-6542_woDup/ GCIMPhigh_H3K4me3_DU-7007_woDup/ GCIMPhigh_H3K4me3_DU-7008_woDup/ GCIMPhigh_H3K4me3_DU-7304_woDup/ GCIMPlow_H3K4me3_06-0129_woDup/ GCIMPlow_H3K4me3_06-1805_woDup/ GCIMPlow_H3K4me3_06-2570_woDup/ GCIMPlow_H3K4me3_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K4me3/gain/output_H3K4me3_all.txt

###CTCF
ctcf.r <- dba.report(ctcf.a)
ctcf.r <- ctcf.r[mcols(ctcf.r)[,4] > 0,]
ctcf.r <- as.data.frame(ctcf.r)
ctcf.r$id <- 1:nrow(ctcf.r)

write.table(ctcf.r[,c(1,2,3,12,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/CTCF/loss/all_loss.bed",sep="\t")

#bed2pos.pl all_loss.bed > all_loss.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/CTCF/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/CTCF/loss/all_loss.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_CTCF_06-0221_woDup/ GCIMPhigh_CTCF_DU-6408_woDup/ GCIMPhigh_CTCF_DU-6542_woDup/ GCIMPhigh_CTCF_DU-7007_woDup/ GCIMPhigh_CTCF_DU-7008_woDup/ GCIMPhigh_CTCF_DU-7304_woDup/ GCIMPlow_CTCF_06-0129_woDup/ GCIMPlow_CTCF_06-1805_woDup/ GCIMPlow_CTCF_06-2570_woDup/ GCIMPlow_CTCF_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/CTCF/loss/output_CTCF_all.txt

ctcf.r <- dba.report(ctcf.a)
ctcf.r <- ctcf.r[mcols(ctcf.r)[,4] < 0,]
ctcf.r <- as.data.frame(ctcf.r)
ctcf.r$id <- 1:nrow(ctcf.r)

write.table(ctcf.r[,c(1,2,3,12,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/CTCF/gain/all_gain.bed",sep="\t")

#bed2pos.pl all_gain.bed > all_gain.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/CTCF/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/CTCF/gain/all_gain.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_CTCF_06-0221_woDup/ GCIMPhigh_CTCF_DU-6408_woDup/ GCIMPhigh_CTCF_DU-6542_woDup/ GCIMPhigh_CTCF_DU-7007_woDup/ GCIMPhigh_CTCF_DU-7008_woDup/ GCIMPhigh_CTCF_DU-7304_woDup/ GCIMPlow_CTCF_06-0129_woDup/ GCIMPlow_CTCF_06-1805_woDup/ GCIMPlow_CTCF_06-2570_woDup/ GCIMPlow_CTCF_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/CTCF/gain/output_CTCF_all.txt

## Figura com os gains e losses de cada experimento
library(reshape2)
library(ggplot2)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_hMEDIP.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
aa$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
aa <- aa[aa$var %in% "cov",]
aa$data <- "5hmC.gain"

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/loss/output_hMEDIP.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "5hmC.loss"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/output_H3k27ac_all.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "H3K27ac.gain"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/loss/output_H3k27ac_all.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "H3K27ac.loss"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K4me3/gain/output_H3K4me3_all.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "H3K4me3.gain"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K4me3/loss/output_H3K4me3_all.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "H3K4me3.loss"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/CTCF/gain/output_CTCF_all.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "CTCF.gain"

aa <- rbind(aa,bb)

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/CTCF/loss/output_CTCF_all.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "CTCF.loss"

aa <- rbind(aa,bb)

ggplot(aa, aes(x = Position, y = value, color = type)) + geom_smooth(span=0.1,method="loess",aes(fill=type),alpha=0.1) + facet_wrap(~data,scales="free") + scale_color_manual(values=c("firebrick","darkgreen"),                                                                                        labels=c("GCIMP-high","GCIMP-low"), name="Legend") +
  scale_fill_manual(values=c("firebrick","darkgreen"))

### ELMER CpG (n=886)
a <- read.csv("/media/data1/Random/ELMER/getPair.hyper.pairs.significant.csv") #886 rows - 669 unique CpGs
load("/media/data1/450k_hg38_GDC.rda")
coor <- hg38.450k[hg38.450k$Composite.Element.REF %in% a$Probe,]
coor$random <- 1
coor$strand <- "*"

write.table(coor[,c(2,3,4,1,5,6)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/ELMER/CpGs/all_CpGs_669.bed",sep="\t")

#bed2pos.pl all_CpGs_669.bed > all_CpGs_669.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/ELMER/CpGs/all_CpGs_669.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/ELMER/CpGs/output_CpGs_all.txt

### ELMER CpG (n=886)
a <- read.csv("/media/data1/Random/ELMER/getPair.hyper.pairs.significant.csv") #886 rows - 669 unique CpGs and 465 unique genes
g <- genes[genes$gene_id %in% a$GeneID,]
tss <- promoters(makeGRangesFromDataFrame(g), upstream = 0, downstream = 1)
tss <- as.data.frame(tss)
tss$id <- rownames(tss)
write.table(tss[,c(1,2,3,6,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/ELMER/TSS/tss_ELMER_465.bed",sep="\t")

#bed2pos.pl tss_ELMER_465.bed > tss_ELMER_465.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/hMEDIP/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/ELMER/TSS/tss_ELMER_465.bed hg38 -size 4000 -hist 10 -d GCIMPhigh_hMEDIP_06-0221_woDup/ GCIMPhigh_hMEDIP_DU-6408_woDup/ GCIMPhigh_hMEDIP_DU-6542_woDup/ GCIMPhigh_hMEDIP_DU-7007_woDup/ GCIMPhigh_hMEDIP_DU-7008_woDup/ GCIMPhigh_hMEDIP_DU-7304_woDup/ GCIMPlow_hMEDIP_06-0129_woDup/ GCIMPlow_hMEDIP_06-1805_woDup/ GCIMPlow_hMEDIP_06-2570_woDup/ GCIMPlow_hMEDIP_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/ELMER/TSS/output_TSS_all.txt

library(reshape2)
library(ggplot2)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/ELMER/CpGs/output_CpGs_all.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
aa$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
aa <- aa[aa$var %in% "cov",]
aa$data <- "5hmC.CpGs"

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/ELMER/TSS/output_TSS_all.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "5hmC.TSS"

aa <- rbind(aa,bb)

ggplot(aa, aes(x = Position, y = value, color = type)) + geom_smooth(span=0.1,method="loess",aes(fill=type),alpha=0.1) + facet_wrap(~data,scales="free") + scale_color_manual(values=c("firebrick","darkgreen"),                                                                                        labels=c("GCIMP-high","GCIMP-low"), name="Legend") +
  scale_fill_manual(values=c("firebrick","darkgreen"))


#####Levels of H3K27ac and H3K4me3 at regions with gain of 5hmC n=808

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/all_gains.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K27ac_06-0221_woDup/ GCIMPhigh_H3K27ac_DU-6408_woDup/ GCIMPhigh_H3K27ac_DU-6542_woDup/ GCIMPhigh_H3K27ac_DU-7007_woDup/ GCIMPhigh_H3K27ac_DU-7008_woDup/ GCIMPhigh_H3K27ac_DU-7304_woDup/ GCIMPlow_H3K27ac_06-0129_woDup/ GCIMPlow_H3K27ac_06-1805_woDup/ GCIMPlow_H3K27ac_06-2570_woDup/ GCIMPlow_H3K27ac_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_H3k27ac_on_gains5hmC.txt


#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K4me3/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/all_gains.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K4me3_06-0221_woDup/ GCIMPhigh_H3K4me3_DU-6408_woDup/ GCIMPhigh_H3K4me3_DU-6542_woDup/ GCIMPhigh_H3K4me3_DU-7007_woDup/ GCIMPhigh_H3K4me3_DU-7008_woDup/ GCIMPhigh_H3K4me3_DU-7304_woDup/ GCIMPlow_H3K4me3_06-0129_woDup/ GCIMPlow_H3K4me3_06-1805_woDup/ GCIMPlow_H3K4me3_06-2570_woDup/ GCIMPlow_H3K4me3_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_H3K4me3_on_gains5hmC.txt

library(reshape2)
library(ggplot2)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_H3k27ac_on_gains5hmC.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
aa$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
aa <- aa[aa$var %in% "cov",]
aa$data <- "H3K27ac"

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_H3K4me3_on_gains5hmC.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "H3K4me3"

aa <- rbind(aa,bb)

ggplot(aa, aes(x = Position, y = value, color = type)) + geom_smooth(span=0.1,method="loess",aes(fill=type),alpha=0.1) + facet_wrap(~data,scales="free") + scale_color_manual(values=c("firebrick","darkgreen"),                                                                                        labels=c("GCIMP-high","GCIMP-low"), name="Legend") +
  scale_fill_manual(values=c("firebrick","darkgreen"))


#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/tss_gain_bygene.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K27ac_06-0221_woDup/ GCIMPhigh_H3K27ac_DU-6408_woDup/ GCIMPhigh_H3K27ac_DU-6542_woDup/ GCIMPhigh_H3K27ac_DU-7007_woDup/ GCIMPhigh_H3K27ac_DU-7008_woDup/ GCIMPhigh_H3K27ac_DU-7304_woDup/ GCIMPlow_H3K27ac_06-0129_woDup/ GCIMPlow_H3K27ac_06-1805_woDup/ GCIMPlow_H3K27ac_06-2570_woDup/ GCIMPlow_H3K27ac_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_H3k27ac_on_TSS_byGene.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K4me3/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/tss_gain_bygene.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K4me3_06-0221_woDup/ GCIMPhigh_H3K4me3_DU-6408_woDup/ GCIMPhigh_H3K4me3_DU-6542_woDup/ GCIMPhigh_H3K4me3_DU-7007_woDup/ GCIMPhigh_H3K4me3_DU-7008_woDup/ GCIMPhigh_H3K4me3_DU-7304_woDup/ GCIMPlow_H3K4me3_06-0129_woDup/ GCIMPlow_H3K4me3_06-1805_woDup/ GCIMPlow_H3K4me3_06-2570_woDup/ GCIMPlow_H3K4me3_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_H3k4me3_on_TSS_byGene.txt


library(reshape2)
library(ggplot2)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_H3k27ac_on_TSS_byGene.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
aa$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
aa <- aa[aa$var %in% "cov",]
aa$data <- "H3K27ac"

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/5hmC/gain/output_H3k4me3_on_TSS_byGene.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "H3K4me3"

aa <- rbind(aa,bb)

ggplot(aa, aes(x = Position, y = value, color = type)) + geom_smooth(span=0.1,method="loess",aes(fill=type),alpha=0.1) + facet_wrap(~data,scales="free") + scale_color_manual(values=c("firebrick","darkgreen"),                                                                                        labels=c("GCIMP-high","GCIMP-low"), name="Legend") +
  scale_fill_manual(values=c("firebrick","darkgreen"))




#cd /media/data1/thais.projects/GCIMP-low/GBM.WGBS/Methylation_-_Bisulfite_Sequencing/JHU_USC__IlluminaHiSeq_WGBS/Level_3/homer/min5reads/
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/CTCF/gain/all_gain.txt hg38 -size 4000 -hist 1 -d TCGA.06.0128.01A/ TCGA.14.1454.01A/ TCGA.14.3477.01A/ TCGA.14.1401.01A/ TCGA.16.1460.01A/ TCGA.19.1788.01A/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/WGBS/CTCF/gain/output_WGBS_CTCF.txt


library(reshape2)
library(ggplot2)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/WGBS/CTCF/gain/output_WGBS_CTCF.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
#aa$type <- c(rep("low", 3*201), rep("high", 2*3*201))
aa$type <- c(rep("low", 3*4001),rep("mesenchymal", 3*4001), rep("LGm6-GBM", 3*4001),rep("mesenchymal", 3*4001),rep("high", 2*3*4001))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 4001), rep("tag5", 4001), rep("tag3", 4001)), 6))
aa <- aa[aa$var %in% "cov",]
aa$data <- "WGBS.CTCF.gains"

ggplot(aa, aes(x = Position, y = value, color = type))  +  
  geom_smooth(aes(fill=type),alpha=0.1,span=0.05,method="loess") +
 scale_color_manual(values=c("firebrick","blue","darkgreen","yellow")) + 
  scale_x_continuous(breaks=arithSeq(-2000,2000)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
facet_wrap(~ variable)

ctcf.r <- dba.report(ctcf.a,th=0.05)
real.ctcf <- read.table("/media/data1/thais.projects/GCIMP-low/CTCF_Flavahan/CTCF_predicted_hg38.bed")
colnames(real.ctcf) <- c("chrom","start","end","id","score","strand")
aux <- findOverlaps(makeGRangesFromDataFrame(real.ctcf),ctcf.r)
real.ctcf <- real.ctcf[unique(queryHits(aux)),]
ctcf.r <- as.data.frame(ctcf.r)
ctcf.r$id <- rownames(ctcf.r)
write.table(ctcf.r[,c(1,2,3,12,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/CTCF/all/all_62k_ctcf.bed",sep="\t")


a <- read.table("/media/data1/thais.projects/GCIMP-low/CTCF_Flavahan/FromHuy/lost_in_idh1mut_hg38.bed")
colnames(a) <- c("chrom","start","end","id","score","strand")
real.ctcf <- read.table("/media/data1/thais.projects/GCIMP-low/CTCF_Flavahan/CTCF_predicted_hg38.bed")
colnames(real.ctcf) <- c("chrom","start","end","id","score","strand")
aux <- findOverlaps(makeGRangesFromDataFrame(real.ctcf),makeGRangesFromDataFrame(a))

#bed2pos.pl all_62k_ctcf.bed > all_62k_ctcf.txt
#cd /media/data1/thais.projects/GCIMP-low/GBM.WGBS/Methylation_-_Bisulfite_Sequencing/JHU_USC__IlluminaHiSeq_WGBS/Level_3/homer/min5reads/
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/CTCF/all/all_62k_ctcf.txt hg38 -size 4000 -hist 1 -d TCGA.06.0128.01A/ TCGA.14.1454.01A/ TCGA.14.3477.01A/ TCGA.14.1401.01A/ TCGA.16.1460.01A/ TCGA.19.1788.01A/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/WGBS/CTCF/all/output_WGBS_CTCF.txt

library(reshape2)
library(ggplot2)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/WGBS/CTCF/all/output_WGBS_CTCF.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
#aa$type <- c(rep("low", 3*201), rep("high", 2*3*201))
aa$type <- c(rep("low", 3*4001),rep("mesenchymal", 3*4001), rep("LGm6-GBM", 3*4001),rep("mesenchymal", 3*4001),rep("high", 2*3*4001))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 4001), rep("tag5", 4001), rep("tag3", 4001)), 6))
aa <- aa[aa$var %in% "cov",]
aa$data <- "WGBS.CTCF.all"

ggplot(aa, aes(x = Position, y = value, color = type))  +  
  geom_smooth(aes(fill=type),alpha=0.1,span=0.01,method="loess") +
  scale_color_manual(values=c("firebrick","blue","darkgreen","yellow")) + 
  scale_x_continuous(breaks=arithSeq(-2000,2000)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


g <- genes[genes$H3K27ac.summary %in% "Gained" & genes$H3K4me3.summary %in% "Gained" & genes$RNAseq.hg38.summary %in% "Upregulated",]
tss <- promoters(makeGRangesFromDataFrame(g), upstream = 0, downstream = 1)
tss <- as.data.frame(tss)
tss$id <- rownames(tss)
write.table(tss[,c(1,2,3,6,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss_gained_of_Ac_me3_upr.bed",sep="\t")

#bed2pos.pl tss_gained_of_Ac_me3_upr.bed > tss_gained_of_Ac_me3_upr.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss_gained_of_Ac_me3_upr.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K27ac_06-0221_woDup/ GCIMPhigh_H3K27ac_DU-6408_woDup/ GCIMPhigh_H3K27ac_DU-6542_woDup/ GCIMPhigh_H3K27ac_DU-7007_woDup/ GCIMPhigh_H3K27ac_DU-7008_woDup/ GCIMPhigh_H3K27ac_DU-7304_woDup/ GCIMPlow_H3K27ac_06-0129_woDup/ GCIMPlow_H3K27ac_06-1805_woDup/ GCIMPlow_H3K27ac_06-2570_woDup/ GCIMPlow_H3K27ac_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/output_H3k27ac_on_TSS_with_gain_Ac_me3_upr.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K4me3/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss_gained_of_Ac_me3_upr.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K4me3_06-0221_woDup/ GCIMPhigh_H3K4me3_DU-6408_woDup/ GCIMPhigh_H3K4me3_DU-6542_woDup/ GCIMPhigh_H3K4me3_DU-7007_woDup/ GCIMPhigh_H3K4me3_DU-7008_woDup/ GCIMPhigh_H3K4me3_DU-7304_woDup/ GCIMPlow_H3K4me3_06-0129_woDup/ GCIMPlow_H3K4me3_06-1805_woDup/ GCIMPlow_H3K4me3_06-2570_woDup/ GCIMPlow_H3K4me3_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/output_H3k4me3_on_TSS_with_gain_Ac_me3_upr.txt


library(reshape2)
library(ggplot2)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/output_H3k4me3_on_TSS_with_gain_Ac_me3_upr.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
aa$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
aa <- aa[aa$var %in% "cov",]
aa$data <- "H3K4me3"

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/output_H3k27ac_on_TSS_with_gain_Ac_me3_upr.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "H3K27ac"

aa <- rbind(aa,bb)

ggplot(aa, aes(x = Position, y = value, color = variable)) + 
  geom_density() +
  #geom_smooth(span=0.5,method="loess",aes(fill=type),alpha=0.1) + 
  facet_wrap(~data,scales="free") + scale_color_manual(values=c("firebrick","darkgreen"),                                                                                        labels=c("GCIMP-high","GCIMP-low"), name="Legend") +
  scale_fill_manual(values=c("firebrick","darkgreen"))


g <- genes[genes$H3K27ac.summary %in% "Gained" & genes$H3K4me3.summary %in% "Gained" & genes$RNAseq.hg38.summary %in% "Upregulated",]
tss <- promoters(makeGRangesFromDataFrame(g), upstream = 0, downstream = 1)
tss <- as.data.frame(tss)
tss$id <- rownames(tss)
write.table(tss[,c(1,2,3,6,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss_gained_of_Ac_me3_upr.bed",sep="\t")

#bed2pos.pl tss_gained_of_Ac_me3_upr.bed > tss_gained_of_Ac_me3_upr.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss_gained_of_Ac_me3_upr.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K27ac_06-0221_woDup/ GCIMPhigh_H3K27ac_DU-6408_woDup/ GCIMPhigh_H3K27ac_DU-6542_woDup/ GCIMPhigh_H3K27ac_DU-7007_woDup/ GCIMPhigh_H3K27ac_DU-7008_woDup/ GCIMPhigh_H3K27ac_DU-7304_woDup/ GCIMPlow_H3K27ac_06-0129_woDup/ GCIMPlow_H3K27ac_06-1805_woDup/ GCIMPlow_H3K27ac_06-2570_woDup/ GCIMPlow_H3K27ac_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/output_H3k27ac_on_TSS_with_gain_Ac_me3_upr.txt


### H3K27ac Gain and Loss
g <- enhancer[enhancer$status %in% "gain" & enhancer$type %in% "TSS",]
write.table(enhancer[,c(1,2,3,14,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/TSS_gain.bed",sep="\t")

#bed2pos.pl TSS_gain.bed > TSS_gain.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/TSS_gain.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K27ac_06-0221_woDup/ GCIMPhigh_H3K27ac_DU-6408_woDup/ GCIMPhigh_H3K27ac_DU-6542_woDup/ GCIMPhigh_H3K27ac_DU-7007_woDup/ GCIMPhigh_H3K27ac_DU-7008_woDup/ GCIMPhigh_H3K27ac_DU-7304_woDup/ GCIMPlow_H3K27ac_06-0129_woDup/ GCIMPlow_H3K27ac_06-1805_woDup/ GCIMPlow_H3K27ac_06-2570_woDup/ GCIMPlow_H3K27ac_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/output_gain_of_H3k27ac_on_TSS.txt

### H3K27ac Gain and Loss
g <- enhancer[enhancer$status %in% "gain" & enhancer$type %in% "TSS",]
g <- genes[genes$gene_id %in% g$geneID,]
write.table(enhancer[,c(1,2,3,14,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/TSS_gain.bed",sep="\t")

#bed2pos.pl TSS_gain.bed > TSS_gain.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/TSS_gain.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K27ac_06-0221_woDup/ GCIMPhigh_H3K27ac_DU-6408_woDup/ GCIMPhigh_H3K27ac_DU-6542_woDup/ GCIMPhigh_H3K27ac_DU-7007_woDup/ GCIMPhigh_H3K27ac_DU-7008_woDup/ GCIMPhigh_H3K27ac_DU-7304_woDup/ GCIMPlow_H3K27ac_06-0129_woDup/ GCIMPlow_H3K27ac_06-1805_woDup/ GCIMPlow_H3K27ac_06-2570_woDup/ GCIMPlow_H3K27ac_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/output_gain_of_H3k27ac_on_TSS.txt

#Butterfly plot - H3K27ac and H3K4me3
# plot p-value
bt2 <- genes[,c(8,15,16,18,19,28,30,31,32)]
bt2$meth <- 0
bt2$meth <- bt2$WGBS.low - bt2$WGBS.high
#bt2[bt2$WGBS.low < 0.25 & !is.na(bt2$WGBS.low),"meth"] <- 0.1
#bt2[bt2$WGBS.low > 0.25 & bt2$WGBS.low < 0.75 &!is.na(bt2$WGBS.low),"meth"] <- 1
#bt2[bt2$WGBS.low > 0.75 & !is.na(bt2$WGBS.low),"meth"] <- 3
#bt2$meth <- bt2$meth + 1
#bt2[is.na(bt2$meth),"meth"] <- 0.1
bt2 <- na.omit(bt2)
bt2$p.ac <- log10(bt2$H3K27ac.FDR)
bt2$p.ac2 <- bt2$p.ac
bt2[bt2$H3K27ac.Fold < 0, "p.ac2"] <- -1 * bt2[bt2$H3K27ac.Fold < 0, "p.ac"]
bt2$p.me3 <- log10(bt2$H3K4me3.FDR)
bt2$p.me <- bt2$p.me3
bt2[bt2$H3K4me3.Fold < 0, "p.me"] <- -1 * bt2[bt2$H3K4me3.Fold < 0, "p.me3"]

bt2$threshold <- "1"

f <- subset(bt2, p.ac2 > -1*log10(0.05)   )
bt2[rownames(f),"threshold"] <- "6"
g <- subset(bt2, p.ac2 < log10(0.05)  )
bt2[rownames(g),"threshold"] <- "6"
h <- subset(bt2, p.me > -1*log10(0.05)  )
bt2[rownames(h),"threshold"] <- "7"
i <- subset(bt2, p.me < log10(0.05)  )
bt2[rownames(i),"threshold"] <- "7"
b <- subset(bt2,p.ac2 < log10(0.05) &  p.me < log10(0.05)  ) 
bt2[rownames(b),"threshold"] <- "2"
c <- subset(bt2,  p.ac2 < log10(0.05) &  p.me > -1*log10(0.05)  )
bt2[rownames(c),"threshold"] <- "3"
d <- subset(bt2, p.ac2 > -1*log10(0.05) &  p.me > -1*log10(0.05)  )
bt2[rownames(d),"threshold"] <- "4"
e <- subset(bt2, p.ac2 > -1*log10(0.05) &  p.me < log10(0.05)  )
bt2[rownames(e),"threshold"] <- "5"
table(bt2$threshold)
bt2$geneID <- rownames(bt2)


ggplot(bt2, aes(x = p.me, y =p.ac2, color = threshold, shape=RNAseq.hg38.summary, size=meth)) + 
  geom_point() + 
  scale_color_manual(values = c("black","darkturquoise", "purple", "darkolivegreen2", "deeppink"),
                     labels=c("Not significant","Depleted in both H3K27ac and H3K4me3","Enriched in both H3K27ac and H3K4me3","Enriched/Depleted only in H3K27ac","Enriched/Depleted only in H3K4me3")) +
  labs(colour = "Significant (padj<0.05)", title = "H3K4me3 and H3K27ac - TSS", x = "log10 (FDR) (H3K4me3)", y = "log10 (FDR) (H3K27ac)") +
 # facet_wrap(~ CNV.GCIMPlow) + 
  theme(axis.text.x = element_text(face="bold", size = 15), 
        axis.title.x = element_text(face="bold", size = 15), 
        axis.title.y = element_text(face="bold", size = 15),
        axis.text.y = element_text( face="bold", size = 15),
        plot.title = element_text(face = "bold",size =20))    +
  theme_bw() + scale_x_continuous(breaks=floor(min(bt2$p.me)):ceiling(max(bt2$p.me))) +
  scale_y_continuous(breaks=floor(min(bt2$p.ac2)):ceiling(max(bt2$p.ac2))) +
  geom_label_repel(data=bt2[(bt2$p.ac2 < -2 & bt2$p.me < -1.4) | (bt2$p.ac2 > 4 & bt2$p.me > 3),],aes(label=gene_name),box.padding = unit(1, "lines"),point.padding = unit(1, "lines")) 

####### Average profile plot
h3k27ac.r <- dba.report(h3k27ac.a)
h3k27ac.r <- h3k27ac.r[mcols(h3k27ac.r)[,4] < 0,]
h3k27ac.r <- as.data.frame(h3k27ac.r)
h3k27ac.r$id <- 1:nrow(h3k27ac.r)

write.table(h3k27ac.r[,c(1,2,3,12,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/all_gain.bed",sep="\t")


#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/all_gain.txt hg38 -size 250000 -hist 10 -d GCIMPhigh_H3K27ac_06-0221_woDup/ GCIMPhigh_H3K27ac_DU-6408_woDup/ GCIMPhigh_H3K27ac_DU-6542_woDup/ GCIMPhigh_H3K27ac_DU-7007_woDup/ GCIMPhigh_H3K27ac_DU-7008_woDup/ GCIMPhigh_H3K27ac_DU-7304_woDup/ GCIMPlow_H3K27ac_06-0129_woDup/ GCIMPlow_H3K27ac_06-1805_woDup/ GCIMPlow_H3K27ac_06-2570_woDup/ GCIMPlow_H3K27ac_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/output_H3k27ac_gain.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GBM.IDHwt/H3K27ac
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/all_gain.txt hg38 -size 250000 -hist 10 -d GSM1306369/ GSM1306371/ GSM1306373/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GBM.IDHwt/H3K27ac/output_H3k27ac_gain.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/all_gain.txt hg38 -size 250000 -hist 10 -d GSM1112810/ GSM916035/ GSM1112812/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files/output_H3k27ac_gain.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/control
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/all_gain.txt hg38 -size 250000 -hist 10 -d blood/files/GSM1841087/ pancreas/files/GSM1013129/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/control/output_H3k27ac_gain.txt

library(reshape2)
library(ggplot2)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/output_H3k27ac_gain.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
aa$type <- c(rep("GCIMP-high", 6*3*5001), rep("GCIMP-low", 4*3*5001))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 5001), rep("tag5", 5001), rep("tag3", 5001)), 10))
aa <- aa[aa$var %in% "cov",]
aa$data <- "H3K27ac"

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GBM.IDHwt/H3K27ac/output_H3k27ac_gain.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("GBM IDHwt", 3*3*5001))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 5001), rep("tag5", 5001), rep("tag3", 5001)), 3))
bb <- bb[bb$var %in% "cov",]
bb$data <- "H3K27ac"

aa <- rbind(aa,bb)

c <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files/output_H3k27ac_gain.txt")
colnames(c)[1] <- "Position"
cc <-(melt(c,id="Position"))
cc$type <- c(rep("Non-tumor brain", 3*3*5001))
cc$type <- as.factor(cc$type)
cc$var <- as.factor(rep(c(rep("cov", 5001), rep("tag5", 5001), rep("tag3", 5001)), 3))
cc <- cc[cc$var %in% "cov",]
cc$data <- "H3K27ac"

aa <- rbind(aa,cc)

d <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/control/output_H3k27ac_gain.txt")
colnames(d)[1] <- "Position"
dd <-(melt(d,id="Position"))
dd$type <- c(rep("Control", 2*3*5001))
dd$type <- as.factor(dd$type)
dd$var <- as.factor(rep(c(rep("cov", 5001), rep("tag5", 5001), rep("tag3", 5001)), 2))
dd <- dd[dd$var %in% "cov",]
dd$data <- "H3K27ac"

aa <- rbind(aa,dd)

ggplot(aa, aes(x = Position, y = value, color = type)) + 
  #geom_density() +
  geom_smooth(span=0.5,method="loess",size=2) + 
 # facet_wrap(~type) + 
  scale_color_manual(values=c("firebrick","grey10","darkgrey","lightslateblue","khaki3"),                                                                                        labels=c("GCIMP-high","GCIMP-low","GBM IDHwt","Non-tumor brain","Control (blood and pancreas)"), name="Legend") #+
  #scale_fill_manual(values=c("firebrick","darkgreen","darkgrey","purple","black"))


h3k27ac.r <- dba.report(h3k27ac.a)
h3k27ac.r <- h3k27ac.r[mcols(h3k27ac.r)[,4] < 0,]
h3k27ac.r <- as.data.frame(h3k27ac.r)
h3k27ac.r$id <- 1:nrow(h3k27ac.r)

write.table(h3k27ac.r[,c(1,2,3,12,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/all_gain.bed",sep="\t")


#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/loss/all_loss.txt hg38 -size 50000 -hist 10 -d GCIMPhigh_H3K27ac_06-0221_woDup/ GCIMPhigh_H3K27ac_DU-6408_woDup/ GCIMPhigh_H3K27ac_DU-6542_woDup/ GCIMPhigh_H3K27ac_DU-7007_woDup/ GCIMPhigh_H3K27ac_DU-7008_woDup/ GCIMPhigh_H3K27ac_DU-7304_woDup/ GCIMPlow_H3K27ac_06-0129_woDup/ GCIMPlow_H3K27ac_06-1805_woDup/ GCIMPlow_H3K27ac_06-2570_woDup/ GCIMPlow_H3K27ac_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/loss/output_H3k27ac_loss.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GBM.IDHwt/H3K27ac
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/loss/all_loss.txt hg38 -size 50000 -hist 10 -d GSM1306369/ GSM1306371/ GSM1306373/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GBM.IDHwt/H3K27ac/output_H3k27ac_loss.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/loss/all_loss.txt hg38 -size 50000 -hist 10 -d GSM1112810/ GSM916035/ GSM1112812/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files/output_H3k27ac_loss.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/control
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/loss/all_loss.txt hg38 -size 50000 -hist 10 -d blood/files/GSM1841087/ pancreas/files/GSM1013129/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/control/output_H3k27ac_loss.txt

library(reshape2)
library(ggplot2)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/loss/output_H3k27ac_loss.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
aa$type <- c(rep("GCIMP-high", 6*3*5001), rep("GCIMP-low", 4*3*5001))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 5001), rep("tag5", 5001), rep("tag3", 5001)), 10))
aa <- aa[aa$var %in% "cov",]
aa$data <- "H3K27ac"

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GBM.IDHwt/H3K27ac/output_H3k27ac_loss.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("GBM IDHwt", 3*3*5001))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 5001), rep("tag5", 5001), rep("tag3", 5001)), 3))
bb <- bb[bb$var %in% "cov",]
bb$data <- "H3K27ac"

aa <- rbind(aa,bb)

c <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files/output_H3k27ac_loss.txt")
colnames(c)[1] <- "Position"
cc <-(melt(c,id="Position"))
cc$type <- c(rep("Non-tumor brain", 3*3*5001))
cc$type <- as.factor(cc$type)
cc$var <- as.factor(rep(c(rep("cov", 5001), rep("tag5", 5001), rep("tag3", 5001)), 3))
cc <- cc[cc$var %in% "cov",]
cc$data <- "H3K27ac"

aa <- rbind(aa,cc)

d <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/control/output_H3k27ac_loss.txt")
colnames(d)[1] <- "Position"
dd <-(melt(d,id="Position"))
dd$type <- c(rep("Control (blood and pancreas)", 2*3*5001))
dd$type <- as.factor(dd$type)
dd$var <- as.factor(rep(c(rep("cov", 5001), rep("tag5", 5001), rep("tag3", 5001)), 2))
dd <- dd[dd$var %in% "cov",]
dd$data <- "H3K27ac"

aa <- rbind(aa,dd)

ggplot(aa, aes(x = Position, y = value, color = type)) + 
  #geom_density() +
  geom_smooth(span=0.5,method="loess",size=2) + 
  #facet_wrap(~type) + 
  scale_color_manual(values=c("firebrick","darkgreen","grey10","lightslateblue","khaki3"),                                                                                        labels=c("GCIMP-high","GCIMP-low","GBM IDHwt","Non-tumor brain","Control (blood and pancreas)"), name="Legend") #+
#scale_fill_manual(values=c("firebrick","darkgreen","darkgrey","purple","black"))

####TSS

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss_pos.txt hg38 -size 4000 -hist 10 -ghist -d GCIMPhigh_H3K27ac_06-0221_woDup/ GCIMPhigh_H3K27ac_DU-6408_woDup/ GCIMPhigh_H3K27ac_DU-6542_woDup/ GCIMPhigh_H3K27ac_DU-7007_woDup/ GCIMPhigh_H3K27ac_DU-7008_woDup/ GCIMPhigh_H3K27ac_DU-7304_woDup/ GCIMPlow_H3K27ac_06-0129_woDup/ GCIMPlow_H3K27ac_06-1805_woDup/ GCIMPlow_H3K27ac_06-2570_woDup/ GCIMPlow_H3K27ac_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/tss/output_H3k27ac_tss.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GBM.IDHwt/H3K27ac
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss_pos.txt hg38 -size 4000 -hist 10 -ghist -d GSM1306369/ GSM1306371/ GSM1306373/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GBM.IDHwt/H3K27ac/output_H3k27ac_tss.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss_pos.txt hg38 -size 4000 -hist 10 -ghist -d GSM1112810/ GSM916035/ GSM1112812/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files/output_H3k27ac_tss.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/control
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss_pos.txt hg38 -size 4000 -hist 10 -ghist -d blood/files/GSM1841087/ pancreas/files/GSM1013129/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/control/output_H3k27ac_tss.txt

library(reshape2)
library(ggplot2)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/tss/output_H3k27ac_tss.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
aa$type <- c(rep("GCIMP-high", 6*3*401), rep("GCIMP-low", 4*3*401))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
aa <- aa[aa$var %in% "cov",]
aa$data <- "H3K27ac"

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GBM.IDHwt/H3K27ac/output_H3k27ac_tss.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("GBM IDHwt", 3*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 3))
bb <- bb[bb$var %in% "cov",]
bb$data <- "H3K27ac"

aa <- rbind(aa,bb)

c <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files/output_H3k27ac_tss.txt")
colnames(c)[1] <- "Position"
cc <-(melt(c,id="Position"))
cc$type <- c(rep("Non-tumor brain", 3*3*401))
cc$type <- as.factor(cc$type)
cc$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 3))
cc <- cc[cc$var %in% "cov",]
cc$data <- "H3K27ac"

aa <- rbind(aa,cc)

d <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/control/output_H3k27ac_tss.txt")
colnames(d)[1] <- "Position"
dd <-(melt(d,id="Position"))
dd$type <- c(rep("Control (blood and pancreas)", 2*3*401))
dd$type <- as.factor(dd$type)
dd$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 2))
dd <- dd[dd$var %in% "cov",]
dd$data <- "H3K27ac"

aa <- rbind(aa,dd)

ggplot(aa, aes(x = Position, y = value, color = type)) + 
  #geom_density() +
  geom_smooth(span=0.1,method="loess",size=2) + 
  #facet_wrap(~type) + 
  scale_color_manual(values=c("firebrick","darkgreen","grey10","lightslateblue","khaki3"),                                                                                        labels=c("GCIMP-high","GCIMP-low","GBM IDHwt","Non-tumor brain","Control (blood and pancreas)"), name="Legend") + theme_gdocs(base_size=16) + ggtitle(label = "Genome-wide H3K27ac across tissues", subtitle = "Centered on TSS") + labs(x = "Position relative to TSS (bp)", y = "ChIP-Fragment Depth (per bp per peak)", colour = "Tissue Type")
#scale_fill_manual(values=c("firebrick","darkgreen","darkgrey","purple","black"))

aux <- genes[genes$H3K27ac.summary %in% "Gained" & genes$H3K4me3.summary %in% "Gained",]
tss <- promoters(makeGRangesFromDataFrame(aux), upstream = 0, downstream = 1)
tss <- as.data.frame(tss)
tss$id <- paste0("A",1:nrow(tss))
write.table(tss[,c(1,2,3,6,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/gain_TSS.bed",sep="\t")

aux <- enhancer[enhancer$type %in% "TSS" & enhancer$status %in% "loss",]
aux <- genes[genes$gene_id %in% aux$geneID,]
tss <- promoters(makeGRangesFromDataFrame(aux), upstream = 0, downstream = 1)
tss <- as.data.frame(tss)
tss$id <- paste0("A",1:nrow(tss))
write.table(tss[,c(1,2,3,6,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/gain_ac_TSS.bed",sep="\t")

#bed2pos.pl loss_ac_TSS.bed > loss_ac_TSS.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/loss/loss_ac_TSS.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K27ac_06-0221_woDup/ GCIMPhigh_H3K27ac_DU-6408_woDup/ GCIMPhigh_H3K27ac_DU-6542_woDup/ GCIMPhigh_H3K27ac_DU-7007_woDup/ GCIMPhigh_H3K27ac_DU-7008_woDup/ GCIMPhigh_H3K27ac_DU-7304_woDup/ GCIMPlow_H3K27ac_06-0129_woDup/ GCIMPlow_H3K27ac_06-1805_woDup/ GCIMPlow_H3K27ac_06-2570_woDup/ GCIMPlow_H3K27ac_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/tss/output_H3k27ac_tss_loss_ac.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GBM.IDHwt/H3K27ac
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/loss/loss_ac_TSS.txt hg38 -size 4000 -hist 10 -d GSM1306369/ GSM1306371/ GSM1306373/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GBM.IDHwt/H3K27ac/output_H3k27ac_tss_loss_ac.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/loss/loss_ac_TSS.txt hg38 -size 4000 -hist 10 -d GSM1112810/ GSM916035/ GSM1112812/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files/output_H3k27ac_tss_loss_ac.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/control
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/loss/loss_ac_TSS.txt hg38 -size 4000 -hist 10 -d blood/files/GSM1841087/ pancreas/files/GSM1013129/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/control/output_H3k27ac_tss_loss_ac.txt

library(reshape2)
library(ggplot2)
library(ggthemes)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/tss/output_H3k27ac_tss_loss_ac.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
aa$type <- c(rep("GCIMP-high", 6*3*401), rep("GCIMP-low", 4*3*401))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
aa <- aa[aa$var %in% "cov",]
aa$data <- "H3K27ac"

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GBM.IDHwt/H3K27ac/output_H3k27ac_tss_loss_ac.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("GBM IDHwt", 3*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 3))
bb <- bb[bb$var %in% "cov",]
bb$data <- "H3K27ac"

aa <- rbind(aa,bb)

c <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files/output_H3k27ac_tss_loss_ac.txt")
colnames(c)[1] <- "Position"
cc <-(melt(c,id="Position"))
cc$type <- c(rep("Non-tumor brain", 3*3*401))
cc$type <- as.factor(cc$type)
cc$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 3))
cc <- cc[cc$var %in% "cov",]
cc$data <- "H3K27ac"

aa <- rbind(aa,cc)

d <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/control/output_H3k27ac_tss_loss_ac.txt")
colnames(d)[1] <- "Position"
dd <-(melt(d,id="Position"))
dd$type <- c(rep("Control (blood and pancreas)", 2*3*401))
dd$type <- as.factor(dd$type)
dd$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 2))
dd <- dd[dd$var %in% "cov",]
dd$data <- "H3K27ac"

aa <- rbind(aa,dd)

ggplot(aa, aes(x = Position, y = value, color = type)) + 
  #geom_density() +
  geom_smooth(span=0.1,method="loess",size=2) + 
  #facet_wrap(~type) + 
  scale_color_manual(values=c("firebrick","darkgreen","grey10","lightslateblue","khaki3"),                                                                                        labels=c("GCIMP-high","GCIMP-low","GBM IDHwt","Non-tumor brain","Control (blood and pancreas)"), name="Legend") +
   theme_gdocs(base_size=16) + ggtitle(label = "Loss of H3K27ac across tissues", subtitle = "Centered on TSS") + labs(x = "Position relative to TSS (bp)", y = "ChIP-Fragment Depth (per bp per peak)", colour = "Tissue Type")
#scale_fill_manual(values=c("firebrick","darkgreen","darkgrey","purple","black"))

write.table(as.character(aux$geneID),quote=F,sep="\t",row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/geneID_TSS.txt")


write.table(as.character(aux$gene_entrezID),quote=F,sep="\t",row.names = F,col.names = T,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K27ac/gain/geneID_up_gained_ac_me3.txt")

#### Profile plot H3K27ac normal and IDHwt
library(reshape2)
library(ggplot2)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files/output_H3k27ac_tss.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
aa$type <- c(rep("Non-tumor brain", 3*3*401))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 3))
aa <- aa[aa$var %in% "cov",]
aa$data <- "H3K27ac"

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GBM.IDHwt/H3K27ac/output_H3k27ac_tss.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("GBM IDHwt", 3*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 3))
bb <- bb[bb$var %in% "cov",]
bb$data <- "H3K27ac"

aa <- rbind(aa,bb)
ggplot(aa, aes(x = Position, y = value, color = type)) + geom_smooth(span=0.1,method="loess",aes(fill=type),alpha=0.1,size=2) + scale_color_manual(values=c("lightslateblue","grey10"),                                                                                        labels=c("Non-tumor brain","GBM IDHwt"), name="Legend") +
  ggtitle(label = "Genome-wide H3K27ac across tissues", subtitle = "Centered on TSS") + labs(x = "Position relative to TSS (bp)", y = "ChIP-Fragment Depth (per bp per peak)", colour = "Tissue Type") + theme_gdocs(base_size=16)


####TSS

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K4me3/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss_pos.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K4me3_06-0221_woDup/ GCIMPhigh_H3K4me3_DU-6408_woDup/ GCIMPhigh_H3K4me3_DU-6542_woDup/ GCIMPhigh_H3K4me3_DU-7007_woDup/ GCIMPhigh_H3K4me3_DU-7008_woDup/ GCIMPhigh_H3K4me3_DU-7304_woDup/ GCIMPlow_H3K4me3_06-0129_woDup/ GCIMPlow_H3K4me3_06-1805_woDup/ GCIMPlow_H3K4me3_06-2570_woDup/ GCIMPlow_H3K4me3_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K4me3/tss/output_H3K4me3_tss.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GBM.IDHwt/H3K4me3
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss_pos.txt hg38 -size 4000 -hist 10 -d GSM1866117/ GSM1866119/ GSM1866127/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GBM.IDHwt/H3K4me3/output_H3K4me3_tss.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/H3K4me3
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/TSS/tss_pos.txt hg38 -size 4000 -hist 10 -d GSM669992/ GSM773012/ GSM772996/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/H3K4me3/output_H3K4me3_tss.txt

library(reshape2)
library(ggplot2)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/H3K4me3/tss/output_H3K4me3_tss.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
aa$type <- c(rep("GCIMP-high", 6*3*401), rep("GCIMP-low", 4*3*401))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
aa <- aa[aa$var %in% "cov",]
aa$data <- "H3K4me3"

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GBM.IDHwt/H3K4me3/output_H3K4me3_tss.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("GBM IDHwt", 3*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 3))
bb <- bb[bb$var %in% "cov",]
bb$data <- "H3K4me3"

aa <- rbind(aa,bb)

c <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/H3K4me3/output_H3K4me3_tss.txt")
colnames(c)[1] <- "Position"
cc <-(melt(c,id="Position"))
cc$type <- c(rep("Non-tumor brain", 3*3*401))
cc$type <- as.factor(cc$type)
cc$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 3))
cc <- cc[cc$var %in% "cov",]
cc$data <- "H3K4me3"

aa <- rbind(aa,cc)


ggplot(aa, aes(x = Position, y = value, color = type)) + 
  #geom_density() +
  geom_smooth(span=0.1,method="loess",size=2) + 
  #facet_wrap(~type) + 
  scale_color_manual(values=c("firebrick","darkgreen","grey10","lightslateblue"),                                                                                        labels=c("GCIMP-high","GCIMP-low","GBM IDHwt","Non-tumor brain"), name="Legend") + theme_gdocs(base_size=16) + ggtitle(label = "Genome-wide H3K4me3 across tissues", subtitle = "Centered on TSS") + labs(x = "Position relative to TSS (bp)", y = "ChIP-Fragment Depth (per bp per peak)", colour = "Tissue Type")
#scale_fill_manual(values=c("firebrick","darkgreen","darkgrey","purple","black"))

#### Vista enhancers hg19
#grep ">" vista_enhancer_hg19.txt | awk '{print $1"\t" $6}' > vista_enhancer_hg19_filt.txt
a <- read.table("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/vista_enhancer_hg19_filt.txt", sep="\t")
a$V1 <- gsub(">Human","",a$V1)
a$V1 <- gsub("[|]","",a$V1)
a$V1 <- gsub(":","-",a$V1)
b <- strsplit(a$V1,"-")
b <- data.frame(matrix(unlist(b), nrow=length(b), byrow=T),stringsAsFactors = F)
b$id <- paste0("A",1:nrow(b))
fwrite(b,quote=F,row.names = F,col.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/vista_enhancer_hg19_filt2.txt")

a <- read.table("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/vista_enhancer_hg38.txt")
b <- read.table("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/vista_enhancer_hg19_filt.txt", sep="\t")
a$result.vista <- b$V2
colnames(a) <- c("chrom","start","end","id","result.vista")
aux.GR <- GRanges(seqnames="chr2",IRanges(start=176060241,end=176234660))
aux <- findOverlaps(aux.GR,makeGRangesFromDataFrame(a))

#### GeneEnhancer GeneCards
a <- read.csv("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/genehancer_11032017.csv")
a$score <- gsub(",",".",a$score)
a$score <- as.numeric(a$score)
aux <- strsplit(as.character(a$attributes),";")
#aux <-  data.frame(matrix(unlist(aux), nrow=length(aux), byrow=T),stringsAsFactors = F)
aux <- unlist(aux)
b <- grep("genehancer_id=",aux)
b <- aux[b]
b <- gsub("genehancer_id=","",b)
a$enhancerID <- b
aux <- strsplit(as.character(a$attributes),";")
aux <- lapply(aux, function(x) x[-1])
info <- ldply(aux,function(aux) {
                  aux <- aux[grep("connected_gene=",aux)]
                  aux <- gsub("connected_gene=","",aux)
                  df <- gene.map[gene.map$symbol %in% aux,]
                  df <- data.frame(geneHancer.total.genes=nrow(df),geneHancer.up=nrow(df[df$RNAseq.hg38.summary %in% "Upregulated",]),geneHancer.down=nrow(df[df$RNAseq.hg38.summary %in% "Downregulated",]))
                  return(df)
               })
a <- cbind(a,info)

aux <- rapply(aux,function(aux) paste(aux[grep("connected_gene=",aux)],aux[grep("score=",aux)],sep="#"),how="list")
aux <- lapply(aux, paste, collapse = "%") 
aux <- unlist(aux)
aux <- gsub("connected_gene=","",aux)
aux <- gsub("score=","",aux)
a$connected_genes <- aux
genehancer <- a
colnames(genehancer) <- paste0("geneHancer.",colnames(genehancer))
genehancer <- genehancer[,-9]
save(genehancer,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GeneHancer/genehancer_table.rda")

### map gene symbol to ensembl
aux <- strsplit(as.character(a$attributes),";")
aux <- unlist(aux)
b <- grep("connected_gene=",aux)
b <- aux[b]
b <- gsub("connected_gene=","",b)
gene.map <- data.frame(symbol=unique(b),stringsAsFactors = F)
gene.map2 <- merge(gene.map,genes[,c("gene_name","gene_id")],by.x="symbol",by.y="gene_name",sort=F)
aux <- gene.map2[duplicated(gene.map2$symbol),"symbol"]
gene.map2 <- gene.map2[!(gene.map2$symbol %in% unique(aux)),]
gene.map3 <- merge(gene.map,genes[,c("gene_name","gene_id")],by.x="symbol",by.y="gene_id",sort=F)
gene.map3 <- gene.map3[,c(1,1)]
colnames(gene.map3) <- c("symbol","gene_id")
aux <- gene.map3[duplicated(gene.map3$symbol),"symbol"]
gene.map3 <- gene.map3[!(gene.map3$symbol %in% unique(aux)),]
gene.map2 <- rbind(gene.map2,gene.map3)
gene.map2 <- gene.map2[!duplicated(gene.map2$symbol),]


aux <- unique(gene.map[!(gene.map$symbol %in% gene.map2$symbol),"symbol"])

ensembl <- useMart("ensembl",
                   dataset = "hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id","external_gene_name")
biomart.search <- getBM(attributes = attributes,
                       filters = "external_gene_name",
                       values = aux, 
                       mart = ensembl)

aux <- biomart.search[duplicated(biomart.search$external_gene_name),]
biomart.search <- biomart.search[!(biomart.search$external_gene_name %in% unique(aux$external_gene_name)),]
id <- genes[genes$gene_id %in% biomart.search$ensembl_gene_id,]
biomart.search <- biomart.search[biomart.search$ensembl_gene_id %in% id$gene_id,]
biomart.search <- biomart.search[,c(2,1)]
colnames(biomart.search) <- c("symbol","gene_id")
gene.map2 <- rbind(gene.map2,biomart.search)

gene.map <- gene.map2
gene.map <- merge(gene.map,genes[,c("gene_id","RNAseq.hg38.summary")],by="gene_id")

load("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/DiffBind/enhancer_annot_v22_cnv_wgbs.rda")
aux <- findOverlaps(makeGRangesFromDataFrame(genehancer),makeGRangesFromDataFrame(enhancer))
enhancer.all <- cbind(enhancer[subjectHits(aux),],genehancer[queryHits(aux),c(2,1,4,5,6,9,10,11,12,13)])
aux <- !(enhancer$id %in% unique(enhancer.all$id))
aux <- enhancer[aux,]
a <- matrix(data = NA, nrow = nrow(aux), ncol = 10,
            dimnames = list(1:nrow(aux),colnames(genehancer)[c(2,1,4,5,6,9,10,11,12,13)]))
aux <- cbind(aux,a)
enhancer.all <- rbind(enhancer.all,aux)

enhancer.all <- enhancer.all[enhancer.all$type %in% c("exon","intron","intergenic"),]
enhancer.all <- enhancer.all[with(enhancer.all,order(FDR)),]
save(enhancer.all,gene.map,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/enhancer_all.rda")

aux <- enhancer.all[enhancer.all$status %in% "gain",c(1,2,3,15,5,5)]
aux$strand.1 <- "+"
write.table(aux,quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/Motif_Analysis/HOX_motifs/All_gain.bed",sep="\t")

aux <- enhancer.all[enhancer.all$status %in% "loss",c(1,2,3,15,5,5)]
aux$strand.1 <- "+"
write.table(aux,quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/Motif_Analysis/HOX_motifs/All_loss.bed",sep="\t")

### HOX binding motif output
hox.gain <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/Motif_Analysis/HOX_motifs/output_gain.txt")
colnames(hox.gain)[1] <- "peakID"
colnames(hox.gain)[22] <- "HOXA2.motif"
colnames(hox.gain)[23] <- "HOXA9.motif"
colnames(hox.gain)[24] <- "HOXC9.motif"
colnames(hox.gain)[25] <- "HOXD13.motif"
hox.gain$overlap <- FALSE
hox.gain[hox.gain$HOXA2.motif != "" | hox.gain$HOXA9.motif != "" | hox.gain$HOXC9.motif != "" | hox.gain$HOXD13.motif != "","overlap"] <- TRUE

hox.loss <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/Motif_Analysis/HOX_motifs/output_loss.txt")
colnames(hox.loss)[1] <- "peakID"
colnames(hox.loss)[22] <- "HOXA2.motif"
colnames(hox.loss)[23] <- "HOXA9.motif"
colnames(hox.loss)[24] <- "HOXC9.motif"
colnames(hox.loss)[25] <- "HOXD13.motif"
hox.loss$overlap <- FALSE
hox.loss[hox.loss$HOXA2.motif != "" | hox.loss$HOXA9.motif != "" | hox.loss$HOXC9.motif != "" | hox.loss$HOXD13.motif != "","overlap"] <- TRUE


#Butterfly plot - H3K27ac and H3K4me3
# plot p-value
bt2 <- enhancer.all[,c(15,11,9,30,22,31,17,32,33,34,16,1,2,3,35)]
#bt2 <- enhancer.all[,c(1,2,3,14,15,11,9,30,16,21,22,31,17,32,33,34,23,24,25,27,28,29)]
#bt2 <- na.omit(bt2)
bt2$p.ac <- log10(bt2$FDR)
bt2$p.ac2 <- bt2$p.ac
bt2[bt2$Fold < 0, "p.ac2"] <- -1 * bt2[bt2$Fold < 0, "p.ac"]
bt2$p.score <- bt2$geneHancer.score

bt2$threshold <- "1"
bt2$size <- NA

#f <- subset(bt2, p.ac2 >= -1*log10(0.05)   )
f <- subset(bt2, p.ac2 >= 2 & p.score > 1 )
bt2[rownames(f),"threshold"] <- "2"
#g <- subset(bt2, p.ac2 <= log10(0.05)  )
g <- subset(bt2, p.ac2 <= -2  & p.score > 1)
bt2[rownames(g),"threshold"] <- "3"

size.up <- bt2[bt2$threshold %in% "2","geneHancer.up"]/bt2[bt2$threshold %in% "2","geneHancer.total.genes"]
#size.up <- size.up + 2
size.up <- round(size.up*100)
#min <- min(size.up,na.rm=T)
#max <- max(size.up,na.rm=T)
#size.up <- (((size.up-min)*(10-5))/(max-min)) + 5
#size.up[is.nan(size.up)] <- 1
bt2[bt2$threshold %in% "2","size"] <- size.up

size.down <- bt2[bt2$threshold %in% "3","geneHancer.down"]/bt2[bt2$threshold %in% "3","geneHancer.total.genes"]
#size.down <- size.down + 2
size.down <- round(size.down*100)
#min <- min(size.down,na.rm=T)
#max <- max(size.down,na.rm=T)
#size.down <- (((size.down-min)*(10-5))/(max-min)) + 5
#size.down[is.nan(size.down)] <- 1
bt2[bt2$threshold %in% "3","size"] <- size.down

bt2[is.nan(bt2$size) | is.na(bt2$size),"size"] <- -1
bt2[bt2$size >= 1 & bt2$size <= 50,"size"] <- 1
#bt2[bt2$size >= 11 & bt2$size <= 20,"size"] <- 11
#bt2[bt2$size >= 21 & bt2$size <= 30,"size"] <- 21
#bt2[bt2$size >= 31 & bt2$size <= 40,"size"] <- 31
#bt2[bt2$size >= 41 & bt2$size <= 50,"size"] <- 41
bt2[bt2$size >= 51 & bt2$size <= 75,"size"] <- 51
bt2[bt2$size >= 76 & bt2$size <= 100,"size"] <- 76

bt2$size <- as.factor(bt2$size)
levels(bt2$size)
#levels(bt2$size) <- c("Not Available","0%","1%","2%-10%","11%-20%","21%-30%","31%-40%","41%-50%","51%-75%","76%-100%")
levels(bt2$size) <- c(paste(c("Not Available","0%","1%-50%","51%-75%","76%-100%"),paste0(paste0("(n = ", table(bt2$size)),")")))
#bt2$size <- factor(bt2$size, levels=c("Not Available","0%","1%","2%-10%","11%-20%","21%-30%","31%-40%","41%-50%","51%-75%","76%-100%"))

ggplot(bt2, aes(x = -1*Fold, y =p.score)) + 
  geom_point(data=subset(bt2,threshold %in% c("1")),aes(x = -1*Fold, y =p.score, color = threshold)) +
  geom_point(data=subset(bt2,threshold %in% c("2")),aes(x = -1*Fold, y =p.score, size=size),color = "#800000",shape=21) +
  geom_point(data=subset(bt2,threshold %in% c("3")),aes(x = -1*Fold, y =p.score,  size=size),color = "#006400",shape=21) +
  geom_point(data=subset(bt2,threshold %in% c("2","3")),aes(x = -1*Fold, y =p.score, color = threshold)) + 
  #facet_wrap(~ overlap.me3) + 
  scale_color_manual(values = c("grey48","#F08080","#98FB98"),
                     labels=c("Not significant","Gain of H3K27ac in G-CIMP-low","Loss of H3K27ac in G-CIMP-low")) +
  labs(colour = "Significant (FDR<0.05)", title = "H3K27ac and GeneHancer", y = "GeneHancer score", x = "Fold-Change (H3K27ac) \n G-CIMP low vs. G-CIMP high") +
  theme_gdocs(base_size = 18) + theme(title = element_text(face = "bold"), axis.title = element_text(face = "bold")) + scale_size_discrete(range = c(0,10)) +
  scale_y_continuous(breaks=floor(min(bt2$p.score,na.rm=T)):ceiling(max(bt2$p.score,na.rm=T))) +
  scale_x_continuous(breaks=floor(min(bt2$Fold,na.rm=T)):ceiling(max(bt2$Fold,na.rm=T))) #+
  #facet_wrap(~ CNV.GCIMPlow) #+
  #geom_label_repel(data=bt2[((bt2$p.ac2 < -7 & bt2$p.score > 1) | (bt2$p.ac2 > 7 & bt2$p.score > 1)) & !is.na(bt2$p.ac2) & !is.na(bt2$p.score) & !is.na(bt2$geneHancer.enhancerID),],aes(label=geneHancer.enhancerID),box.padding = unit(1, "lines"),point.padding = unit(1, "lines")) 


#motif analysis homer - enhancer (gain)
a <- bt2[bt2$size %in% c("51%-75% (n = 27)","76%-100% (n = 55)"),]
a <- a[a$status %in% "gain",]
a <- enhancer.all[enhancer.all$id %in% a$id,]
a$strand <- "+"
write.table(a[,c(15,1,2,3,5)],quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/Motif_Analysis/Up/ac_peaks.txt")

#motif analysis homer - TSS (gain)
aux <- strsplit(as.character(a$geneHancer.connected_genes),"%")
aux <- unlist(aux)
aux <- strsplit(aux,"#")
aux <- lapply(aux, function(x) x[1])
aux <- unlist(aux)
aux <- gene.map[gene.map$symbol %in% aux,]
aux <- aux[aux$RNAseq.hg38.summary %in% "Upregulated",]
write.table(aux$gene_id,quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/Motif_Analysis/Up/TSS/gene_id.txt")

#motif analysis homer - enhancer (loss)
a <- bt2[bt2$size %in% c("51%-75% (n = 27)","76%-100% (n = 55)"),]
a <- a[a$status %in% "loss",]
a <- enhancer.all[enhancer.all$id %in% a$id,]
a$strand <- "+"
write.table(a[,c(15,1,2,3,5)],quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/Motif_Analysis/Down/ac_peaks.txt")


#motif analysis homer - TSS (loss)
aux <- strsplit(as.character(a$geneHancer.connected_genes),"%")
aux <- unlist(aux)
aux <- strsplit(aux,"#")
aux <- lapply(aux, function(x) x[1])
aux <- unlist(aux)
aux <- gene.map[gene.map$symbol %in% aux,]
aux <- aux[aux$RNAseq.hg38.summary %in% "Downregulated",]
write.table(aux$gene_id,quote=F,col.names = F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/Motif_Analysis/Down/TSS/gene_id.txt")

aux.GR <- GRanges(seqnames="chr2",IRanges(start=176000000,end=176300000))
aux <- findOverlaps(aux.GR,makeGRangesFromDataFrame(a))
a <- a[subjectHits(aux),]
a <- a[with(a,order(start)),]
aux <- findOverlaps(aux.GR,makeGRangesFromDataFrame(enhancer))
enhancer[subjectHits(aux),]

#remove duplicates
aux <- as.data.frame(aux)
dup <- duplicated(aux$subjectHits)
dup <- aux[dup,]
dup.s <- aux[aux$subjectHits %in% unique(dup$subjectHits),]
enhancer.dup <- enhancer[unique(dup.s$subjectHits),]

b <- a[unique(queryHits(aux)),]
b <- b[with(b,order(-score)),]


bt2 <- genes[,c(8,15,16,18,19,28)]
bt2 <- na.omit(bt2)
bt2$p.ac <- log10(bt2$H3K27ac.FDR)
bt2$p.ac2 <- bt2$p.ac
bt2[bt2$H3K27ac.Fold < 0, "p.ac2"] <- -1 * bt2[bt2$H3K27ac.Fold < 0, "p.ac"]
bt2$p.me3 <- log10(bt2$H3K4me3.FDR)
bt2$p.me <- bt2$p.me3
bt2[bt2$H3K4me3.Fold < 0, "p.me"] <- -1 * bt2[bt2$H3K4me3.Fold < 0, "p.me3"]

bt2$threshold <- "1"

f <- subset(bt2, p.ac2 > -1*log10(0.05)   )
bt2[rownames(f),"threshold"] <- "6"
g <- subset(bt2, p.ac2 < log10(0.05)  )
bt2[rownames(g),"threshold"] <- "6"
h <- subset(bt2, p.me > -1*log10(0.05)  )
bt2[rownames(h),"threshold"] <- "7"
i <- subset(bt2, p.me < log10(0.05)  )
bt2[rownames(i),"threshold"] <- "7"
b <- subset(bt2,p.ac2 < log10(0.05) &  p.me < log10(0.05)  ) 
bt2[rownames(b),"threshold"] <- "2"
c <- subset(bt2,  p.ac2 < log10(0.05) &  p.me > -1*log10(0.05)  )
bt2[rownames(c),"threshold"] <- "3"
d <- subset(bt2, p.ac2 > -1*log10(0.05) &  p.me > -1*log10(0.05)  )
bt2[rownames(d),"threshold"] <- "4"
e <- subset(bt2, p.ac2 > -1*log10(0.05) &  p.me < log10(0.05)  )
bt2[rownames(e),"threshold"] <- "5"
table(bt2$threshold)
bt2$geneID <- rownames(bt2)

t.t <- as.data.frame(table(bt2$threshold))
t.t[,1] <- c("Not Significant", "Depleted in both H3K27ac and H3K4me3","Enriched in both H3K27ac and H3K4me3","Enriched/Depleted only in H3K27ac","Enriched/Depleted only in H3K4me3")
colnames(t.t) <- c("","Count")

star<- ggplot(bt2, aes(x = p.me, y =p.ac2, color = threshold, shape=RNAseq.hg38.summary)) + 
  geom_point(size=2) + 
  #scale_color_manual(values = c("black","darkturquoise", "purple", "darkolivegreen2", "deeppink"),
  #                   labels=c("Not significant","Depleted in both H3K27ac and H3K4me3","Enriched in both H3K27ac and H3K4me3","Enriched/Depleted only in H3K27ac","Enriched/Depleted only in H3K4me3")) +
  labs(colour = "Significant (padj<0.05)", title = "H3K4me3 and H3K27ac - TSS", x = "log10 (FDR) (H3K4me3)", y = "log10 (FDR) (H3K27ac)") +
  annotation_custom(tableGrob(t.t), xmin=10,ymin=10,ymax=20) + theme_gdocs(base_size = 16) +
  geom_label_repel(data=bt2[(bt2$p.ac2 < -2 & bt2$p.me < -1.4) | (bt2$p.ac2 > 4 & bt2$p.me > 3),],aes(label=gene_name),box.padding = unit(1, "lines"),point.padding = unit(1, "lines")) + theme(title = element_text(face = "bold"), axis.title = element_text(face = "bold"))


library(ggplot2); library(grid); library(gridExtra); library(ggthemes)
genes$RNAseq.hg38.summary.reorder <- factor(genes$RNAseq.hg38.summary,levels=c("Downregulated","Upregulated","No.difference","No.info"))

t.rna <- as.data.frame(table(genes$RNAseq.hg38.summary.reorder)); colnames(t.rna) <- c("","Count")
t.k27ac <- as.data.frame(table(genes$H3K27ac.summary)); colnames(t.k27ac) <- c("","Count")
t.k4me3 <- as.data.frame(table(genes$H3K4me3.summary)); colnames(t.k4me3) <- c("","Count")

p.rna <- ggplot(genes, aes(x = RNAseqALL.FC, y = -1*log10(RNAseqALL.pv.adj), color = RNAseq.hg38.summary.reorder)) + geom_point() + annotation_custom(tableGrob(t.rna), xmin=0,ymin=30,ymax =35) + ggtitle(label = "Differential Express Genes (DEG)", subtitle = "G-CIMP low vs G-CIMP high") + labs(x = "Gene Expression Fold Change", y = "Significance (-log10(FDR))", colour = "Gene Exp Summary")
p.rna.theme <- p.rna + theme_gdocs(base_size=18) + scale_colour_tableau("colorblind10")  + theme(title = element_text(face = "bold"), axis.title = element_text(face = "bold"))

p.k4me3 <- ggplot(genes, aes(x = H3K4me3.Fold, y = -1*log10(H3K4me3.FDR), color = H3K4me3.summary)) + geom_point() + annotation_custom(tableGrob(t.k4me3), xmin=0,ymin=25,ymax =30) + ggtitle(label = "Differential Binding H3K4me3", subtitle = "G-CIMP low vs G-CIMP high") + labs(x = "H3K4me3 Fold Change", y = "Significance (-log10(FDR))", colour = "H3K4me3 Summary")
p.k4me3.theme <- p.k4me3 + theme_gdocs(base_size = 18) + scale_colour_tableau("colorblind10")  + theme(title = element_text(face = "bold"), axis.title = element_text(face = "bold"))  + theme(title = element_text(face = "bold"), axis.title = element_text(face = "bold"))

p.k27ac <- ggplot(genes, aes(x = H3K27ac.Fold, y = -1*log10(H3K27ac.FDR), color = H3K27ac.summary)) + geom_point() + annotation_custom(tableGrob(t.k27ac), xmin=0,ymin=15,ymax =20) + ggtitle(label = "Differential Binding H3K27ac", subtitle = "G-CIMP low vs G-CIMP high") + labs(x = "H3K27ac Fold Change", y = "Significance (-log10(FDR))", colour = "H3K27ac Summary")  + theme(title = element_text(face = "bold"), axis.title = element_text(face = "bold"))
p.k27ac.theme <- p.k27ac + theme_gdocs(base_size = 18) + scale_colour_tableau("colorblind10") + theme(title = element_text(face = "bold"), axis.title = element_text(face = "bold"))


grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p.rna.theme, vp = vplayout(1, 1:2))  # key is to define vplayout
print(p.k4me3.theme, vp = vplayout(2, 1))
print(p.k27ac.theme, vp = vplayout(2, 2))



### Profile plot - H3k27ac geneHancer
aux <- bt2[bt2$threshold %in% "2",c(12,13,14,6)]
aux$id <- 2000
aux$strand <- "+"
write.table(aux,quote = F,row.names = F,col.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/ac_gain.bed")

cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer
#bed2pos.pl ac_gain.bed > ac_gain.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/ac_gain.txt hg38 -size 50000 -hist 10 -d GCIMPhigh_H3K27ac_06-0221_woDup/ GCIMPhigh_H3K27ac_DU-6408_woDup/ GCIMPhigh_H3K27ac_DU-6542_woDup/ GCIMPhigh_H3K27ac_DU-7007_woDup/ GCIMPhigh_H3K27ac_DU-7008_woDup/ GCIMPhigh_H3K27ac_DU-7304_woDup/ GCIMPlow_H3K27ac_06-0129_woDup/ GCIMPlow_H3K27ac_06-1805_woDup/ GCIMPlow_H3K27ac_06-2570_woDup/ GCIMPlow_H3K27ac_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/output_H3K27ac_gain.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K4me3/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/ac_gain.txt hg38 -size 50000 -hist 10 -d GCIMPhigh_H3K4me3_06-0221_woDup/ GCIMPhigh_H3K4me3_DU-6408_woDup/ GCIMPhigh_H3K4me3_DU-6542_woDup/ GCIMPhigh_H3K4me3_DU-7007_woDup/ GCIMPhigh_H3K4me3_DU-7008_woDup/ GCIMPhigh_H3K4me3_DU-7304_woDup/ GCIMPlow_H3K4me3_06-0129_woDup/ GCIMPlow_H3K4me3_06-1805_woDup/ GCIMPlow_H3K4me3_06-2570_woDup/ GCIMPlow_H3K4me3_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/output_H3K4me3_gain.txt

a <- bt2[bt2$threshold %in% "2",]
a <- a[a$status %in% "gain",]

#motif analysis homer - TSS (gain)
aux <- strsplit(as.character(a$geneHancer.connected_genes),"%")
aux <- unlist(aux)
aux <- strsplit(aux,"#")
aux <- lapply(aux, function(x) x[1])
aux <- unlist(aux)
aux <- gene.map[gene.map$symbol %in% aux,]
aux <- aux[aux$RNAseq.hg38.summary %in% "Upregulated",]
aux <- genes[genes$gene_id %in% aux$gene_id,]
tss <- promoters(makeGRangesFromDataFrame(aux), upstream = 0, downstream = 1)
tss <- as.data.frame(tss)
tss$id <- rownames(tss)
write.table(tss[,c(1,2,3,6,4,5)],quote=F,row.names = F,col.names = F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/TSS/gene_id_gain.bed",sep="\t")

#bed2pos.pl gene_id_gain.bed > gene_id_gain.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/TSS/gene_id_gain.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K27ac_06-0221_woDup/ GCIMPhigh_H3K27ac_DU-6408_woDup/ GCIMPhigh_H3K27ac_DU-6542_woDup/ GCIMPhigh_H3K27ac_DU-7007_woDup/ GCIMPhigh_H3K27ac_DU-7008_woDup/ GCIMPhigh_H3K27ac_DU-7304_woDup/ GCIMPlow_H3K27ac_06-0129_woDup/ GCIMPlow_H3K27ac_06-1805_woDup/ GCIMPlow_H3K27ac_06-2570_woDup/ GCIMPlow_H3K27ac_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/TSS/output_H3K27ac_TSS_gain.txt

#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K4me3/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/TSS/gene_id_gain.txt hg38 -size 4000 -hist 10 -d GCIMPhigh_H3K4me3_06-0221_woDup/ GCIMPhigh_H3K4me3_DU-6408_woDup/ GCIMPhigh_H3K4me3_DU-6542_woDup/ GCIMPhigh_H3K4me3_DU-7007_woDup/ GCIMPhigh_H3K4me3_DU-7008_woDup/ GCIMPhigh_H3K4me3_DU-7304_woDup/ GCIMPlow_H3K4me3_06-0129_woDup/ GCIMPlow_H3K4me3_06-1805_woDup/ GCIMPlow_H3K4me3_06-2570_woDup/ GCIMPlow_H3K4me3_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/TSS/output_H3K4me3_TSS_gain.txt


library(reshape2)
library(ggplot2)
a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/TSS/output_H3K27ac_TSS_gain.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
aa$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
aa <- aa[aa$var %in% "cov",]
aa$data <- "H3K27ac mark on TSS"

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/TSS/output_H3K4me3_TSS_gain.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*401), rep("low", 4*3*401))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 401), rep("tag5", 401), rep("tag3", 401)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "H3K4me3 mark on TSS"

aa <- rbind(aa,bb)

ggplot(aa, aes(x = Position, y = value, color = type)) + geom_smooth(span=0.01,method="loess",aes(fill=type),alpha=0.1) + facet_wrap(~data,scales="free") + scale_color_manual(values=c("firebrick","darkgreen"),                                                                                        labels=c("GCIMP-high","GCIMP-low"), name="Legend") +
  scale_fill_manual(values=c("firebrick","darkgreen")) + theme_gdocs(base_size = 18) + theme(title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  labs(title = "Histone modifications on TSS of potential target genes (n=150) \n upregulated by nearby active enhancers")


a <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/output_H3K27ac_gain.txt")
colnames(a)[1] <- "Position"
aa <-(melt(a,id="Position"))
aa$type <- c(rep("high", 6*3*5001), rep("low", 4*3*5001))
aa$type <- as.factor(aa$type)
aa$var <- as.factor(rep(c(rep("cov", 5001), rep("tag5", 5001), rep("tag3", 5001)), 10))
aa <- aa[aa$var %in% "cov",]
aa$data <- "H3K27ac mark on potential active enhancer"

b <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/output_H3K4me3_gain.txt")
colnames(b)[1] <- "Position"
bb <-(melt(b,id="Position"))
bb$type <- c(rep("high", 6*3*5001), rep("low", 4*3*5001))
bb$type <- as.factor(bb$type)
bb$var <- as.factor(rep(c(rep("cov", 5001), rep("tag5", 5001), rep("tag3", 5001)), 10))
bb <- bb[bb$var %in% "cov",]
bb$data <- "H3K4me3 mark on potential active enhancer"

aa <- rbind(aa,bb)

ggplot(aa, aes(x = Position, y = value, color = type)) + geom_smooth(span=0.3,method="loess",aes(fill=type),alpha=0.1) + facet_wrap(~data) + scale_color_manual(values=c("firebrick","darkgreen"),                                                                                        labels=c("GCIMP-high","GCIMP-low"), name="Legend") +
  scale_fill_manual(values=c("firebrick","darkgreen")) + theme_gdocs(base_size = 18) + theme(title = element_text(face = "bold"), axis.title = element_text(face = "bold")) +
  labs(title = "Histone modifications on potential active enhancers (n=309) \n based on gain of H3K27ac and GeneHancer score")


aux <- bt2[bt2$threshold %in% "3",c(12,13,14,6)]
aux$id <- 2000
aux$strand <- "+"
write.table(aux,quote = F,row.names = F,col.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/ac_loss.bed")

cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer
#bed2pos.pl ac_loss.bed > ac_loss.txt
#cd /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer
#annotatePeaks.pl /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/ac_loss.txt hg38 -size 2000 -hist 10 -d GCIMPhigh_H3K27ac_06-0221_woDup/ GCIMPhigh_H3K27ac_DU-6408_woDup/ GCIMPhigh_H3K27ac_DU-6542_woDup/ GCIMPhigh_H3K27ac_DU-7007_woDup/ GCIMPhigh_H3K27ac_DU-7008_woDup/ GCIMPhigh_H3K27ac_DU-7304_woDup/ GCIMPlow_H3K27ac_06-0129_woDup/ GCIMPlow_H3K27ac_06-1805_woDup/ GCIMPlow_H3K27ac_06-2570_woDup/ GCIMPlow_H3K27ac_DU-7010_woDup/ > /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/Profile_plot/geneHancer/output_H3K27ac_loss.txt


### GeneHancer com info de transcription factors
load("/media/data2/TCGA_RNAseq/gene_annotation_v22_cnv_wgbs.rda")

#### Aqui eu gerei uma coluna no objeto genes que chama "RNAseq.hg38.summary", que basicamente tem 3 levels: Upregulated, Downregulated, No.difference. Eu classifiquei os genes, usando os dados de RNAseq e fazendo uma analise supervisionada (t-test) nessas 3 classes, baseado em FDR e FC. Imagino que voce tenha isso tambem.. ai soh juntei essa informacao da analise no meu objeto "genes"

a <- read.table("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/enhancer_ids.txt",header=T)

a <- read.csv("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/genehancer_04172018.csv")
a$score <- gsub(",",".",a$score)
a$score <- as.numeric(a$score)
aux <- strsplit(as.character(a$attributes),";")
aux <- unlist(aux)
b <- grep("genehancer_id=",aux)
b <- aux[b]
b <- gsub("genehancer_id=","",b)
a$enhancerID <- b

aux <- strsplit(as.character(a$attributes),";")
aux <- unlist(aux)
b <- grep("connected_gene=",aux)
b <- aux[b]
b <- gsub("connected_gene=","",b)
gene.map <- data.frame(symbol=unique(b),stringsAsFactors = F)
gene.map2 <- merge(gene.map,genes[,c("gene_name","gene_id")],by.x="symbol",by.y="gene_name",sort=F)
aux <- gene.map2[duplicated(gene.map2$symbol),"symbol"]
gene.map2 <- gene.map2[!(gene.map2$symbol %in% unique(aux)),]
gene.map3 <- merge(gene.map,genes[,c("gene_name","gene_id")],by.x="symbol",by.y="gene_id",sort=F)
gene.map3 <- gene.map3[,c(1,1)]
colnames(gene.map3) <- c("symbol","gene_id")
aux <- gene.map3[duplicated(gene.map3$symbol),"symbol"]
gene.map3 <- gene.map3[!(gene.map3$symbol %in% unique(aux)),]
gene.map2 <- rbind(gene.map2,gene.map3)
gene.map2 <- gene.map2[!duplicated(gene.map2$symbol),]


aux <- unique(gene.map[!(gene.map$symbol %in% gene.map2$symbol),"symbol"])

ensembl <- useMart("ensembl",
                   dataset = "hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id","external_gene_name")
biomart.search <- getBM(attributes = attributes,
                        filters = "external_gene_name",
                        values = aux, 
                        mart = ensembl)

aux <- biomart.search[duplicated(biomart.search$external_gene_name),]
biomart.search <- biomart.search[!(biomart.search$external_gene_name %in% unique(aux$external_gene_name)),]
id <- genes[genes$gene_id %in% biomart.search$ensembl_gene_id,]
biomart.search <- biomart.search[biomart.search$ensembl_gene_id %in% id$gene_id,]
biomart.search <- biomart.search[,c(2,1)]
colnames(biomart.search) <- c("symbol","gene_id")
gene.map2 <- rbind(gene.map2,biomart.search)

gene.map <- gene.map2
gene.map <- merge(gene.map,genes[,c("gene_id","RNAseq.hg38.summary")],by="gene_id")


aux <- strsplit(as.character(a$attributes),";")
aux <- lapply(aux, function(x) x[-1])
info <- ldply(aux,function(aux) {
  aux <- aux[grep("connected_gene=",aux)]
  aux <- gsub("connected_gene=","",aux)
  df <- gene.map[gene.map$symbol %in% aux,]
  df <- data.frame(geneHancer.total.genes=nrow(df),geneHancer.up=nrow(df[df$RNAseq.hg38.summary %in% "Upregulated",]),geneHancer.down=nrow(df[df$RNAseq.hg38.summary %in% "Downregulated",]))
  return(df)
})
a <- cbind(a,info)

aux <- rapply(aux,function(aux) paste(aux[grep("connected_gene=",aux)],aux[grep("score=",aux)],sep="#"),how="list")
aux <- lapply(aux, paste, collapse = "%") 
aux <- unlist(aux)
aux <- gsub("connected_gene=","",aux)
aux <- gsub("score=","",aux)
a$connected_genes <- aux
genehancer <- a
colnames(genehancer) <- paste0("geneHancer.",colnames(genehancer))
genehancer <- genehancer[,-9]

enh_ids <- read.table("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/enhancer_ids_041618.txt",header=TRUE,sep="\t",stringsAsFactors = F)
genehancer <- merge(genehancer,enh_ids[,c(4,5,6)],by.x="geneHancer.enhancerID",by.y="GHid")
colnames(genehancer)[(ncol(genehancer)-1):ncol(genehancer)] <- c("geneHancer.cluster.id","geneHancer.is.elite")

save(genehancer,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GeneHancer/genehancer_table_April2018.rda")

load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/cat_5hmc_v2.rda")
aux <- findOverlaps(makeGRangesFromDataFrame(cat.5hmC),makeGRangesFromDataFrame(genehancer))
aux <- as.data.frame(aux)
cat.5hmC$geneHancer <- FALSE
cat.5hmC[unique(aux$queryHits),"geneHancer"] <- TRUE
no.dup <- aux[!(aux$queryHits %in% unique(aux$queryHits[(duplicated(aux$queryHits))])),]
cat.5hmC$geneHancer.id <- NA
cat.5hmC[no.dup$queryHits,"geneHancer.id"] <- as.character(genehancer[no.dup$subjectHits,"geneHancer.enhancerID"])
dup <- aux[aux$queryHits %in% unique(aux$queryHits[(duplicated(aux$queryHits))]),]
dup <- cbind(cat.5hmC[dup$queryHits,c(15)],genehancer[dup$subjectHits,c(1)])
dup <- data.frame(dup,stringsAsFactors = F)
df1 <- aggregate(dup[2], dup[-2], unique)
cat.5hmC[df1$X1,"geneHancer.id"] <- sapply( df1$X2, paste0, collapse="#")

save(cat.5hmC,file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/cat_5hmc_v3.rda")

ggplot(cat.5hmC,aes(status,fill=geneHancer)) + geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent_format()) 

##add the TF to the 5hmC table
enh_tfs <- read.table("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/enhancer_tfs_041618.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE) #20,509 different tissues
enh_tfs[enh_tfs$TF %in% "NSD2","TF"] <- "WHSC1" #ENSG00000109685
enh_tfs[enh_tfs$TF %in% "EMSY","TF"] <- "C11orf30" #ENSG00000158636
enh_tfs[enh_tfs$TF %in% "CAVIN1","TF"] <- "PTRF" #ENSG00000177469
enh_tfs <- merge(enh_tfs,genes[,c("gene_name","RNAseq.hg38.summary")],by.x="TF",by.y="gene_name",all.x=TRUE)
aux <- (mclapply(unique(enh_tfs$enhancer_cluster_id),
                    function(i){
                      aux <- enh_tfs[enh_tfs$enhancer_cluster_id %in% i,]
                      if(nrow(aux[aux$RNAseq.hg38.summary %in% c("Downregulated", "Upregulated"),]) == 0)
                        return(aux)
                      else{
                        aux <- aux[aux$RNAseq.hg38.summary %in% c("Downregulated","Upregulated"),]
                        return(aux)
                      }
                    },mc.cores=12))
a <- melt(aux)
a <- a[,-4]
colnames(a)[4] <- "enhancer_cluster_id"
a <- a[,c(1,4,2,3)]
enh_tfs.s <- a
save(enh_tfs.s,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Enhancer/enhancer_tfs_subset.rda")

# length(unique(enh_tfs.s$enhancer_cluster_id))
# a <- enh_tfs.s[,c(2,1)]
# a1 <- aggregate(a[2], a[-2], unique) #202462
# length(unique(a1$enhancer_cluster_id))
# 
# a <- enh_tfs.s[,c(2,3)]
# a2 <- aggregate(a[2], a[-2], unique) #202462
# length(unique(a2$enhancer_cluster_id))
# 
# a <- enh_tfs.s[,c(2,4)]
# a3 <- aggregate(a[2], a[-2], unique) #202462
# length(unique(a3$enhancer_cluster_id))
# 
# identical(a1$enhancer_cluster_id,a2$enhancer_cluster_id)
# identical(a2$enhancer_cluster_id,a3$enhancer_cluster_id)
# 
# a <- data.frame(enhancer_cluster_id=a1$enhancer_cluster_id,TF=I(a1$TF),Tissue=I(a2$tissues),RNAseq.hg38.summary=I(a3$RNAseq.hg38.summary))
# a <- merge(a,genehancer[,c("geneHancer.enhancerID","geneHancer.is.elite","geneHancer.cluster.id")],by.x="enhancer_cluster_id",by.y="geneHancer.cluster.id")

# cat.5hmC <- merge(cat.5hmC,a,by.x="geneHancer.id",by.y="geneHancer.enhancerID",all.x=T,sort=F)

a1 <- aggregate(TF ~ enhancer_cluster_id, data = enh_tfs.s, paste, collapse = ".")
a2 <- aggregate(tissues ~ enhancer_cluster_id, data = enh_tfs.s, paste, collapse = ".")
a3 <- aggregate(RNAseq.hg38.summary ~ enhancer_cluster_id, data = enh_tfs.s, paste, collapse = ".")
a <- data.frame(enhancer_cluster_id=a1$enhancer_cluster_id,TF=a1$TF,Tissue=a2$tissues,RNAseq.hg38.summary=a3$RNAseq.hg38.summary)

a <- merge(a,genehancer[,c("geneHancer.enhancerID","geneHancer.score","geneHancer.is.elite","geneHancer.cluster.id")],by.x="enhancer_cluster_id",by.y="geneHancer.cluster.id")

cat.5hmC <- merge(cat.5hmC,a,by.x="geneHancer.id",by.y="geneHancer.enhancerID",all.x=T,sort=F)
colnames(cat.5hmC)[31:33] <- c("geneHancer.TF","geneHancer.Tissue","geneHancer.TF.expression")

save(cat.5hmC,file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/cat_5hmc_v4.rda")

require(rtracklayer)
load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/cat_5hmc_v4.rda")
gencode <- readGFF("/media/data2/TCGA_RNAseq/GDC_gencodeV22/gencode.v22.annotation.gtf")
gencode <- gencode[grep("+_PAR_+",gencode$gene_id,invert=T),]
aux <- strsplit(gencode$gene_id, "\\." )
aux <- unlist(lapply(aux, "[[", 1 ))
gencode$gene_id <- aux #Ensembl id list
genes <- gencode[gencode$type %in% "transcript",]

tss <- promoters(makeGRangesFromDataFrame(genes), upstream = 0, downstream = 1)

a <- suppressWarnings(distanceToNearest(makeGRangesFromDataFrame(cat.5hmC),tss,ignore.strand=T))
rownames(cat.5hmC) <- NULL
cat.5hmC$geneID <- NA
cat.5hmC$distance.to.TSS <- NA
cat.5hmC[unique(queryHits(a)),"distance.to.TSS"] <- mcols(a)[,1]
a <- a[mcols(a)[,1] < 2000,]
cat.5hmC$type <- NA
cat.5hmC[unique(queryHits(a)),"type"] <- "TSS" #no duplicates
aux <- as.character(genes[subjectHits(a),"gene_id"])
cat.5hmC[queryHits(a),"geneID"] <- aux

aux <- cat.5hmC[is.na(cat.5hmC$type),]
genes <- gencode[gencode$type %in% "gene",]
a1 <- suppressWarnings(findOverlaps(makeGRangesFromDataFrame(aux),makeGRangesFromDataFrame(genes),ignore.strand=T))
a1 <- as.data.frame(a1)
a1 <- a1[!duplicated(a1$queryHits),]
aux <- aux[a1$queryHits,]
aux$geneID <- as.character(genes[a1$subjectHits,"gene_id"])
genes <- gencode[gencode$type %in% "exon",]
#Remover duplicated baseado em uma condição de outra coluna
# genes.s <- ddply(genes, .(seqid,start), summarize,
#       end=max(end),
#       gene_id=unique(gene_id))
#genes[c(1005424,1005441,1005635),]
a <- suppressWarnings(findOverlaps(makeGRangesFromDataFrame(aux),makeGRangesFromDataFrame(genes),ignore.strand=T))
a <- as.data.frame(a)
dup <- duplicated(a$queryHits)
dup <- a[dup,]
dup.s <- a[a$queryHits %in% unique(dup$queryHits),]
uni <- a[!(a$queryHits %in% unique(dup$queryHits)),]
b <- aux[uni$queryHits,"id"]
cat.5hmC[cat.5hmC$id %in% b,"type"] <- "exon"
cat.5hmC[cat.5hmC$id %in% b,"geneID"] <- as.character(genes[uni$subjectHits,"gene_id"])
for(i in unique(dup.s$queryHits)){
  b <- dup.s[dup.s$queryHits %in% i,"subjectHits"]
  g <- unique(genes[b,"gene_id"])
  if(length(g) > 1)
    g <- paste( g, collapse="#")
  b <- aux[i,"id"]
  cat.5hmC[cat.5hmC$id %in% b,"type"] <- "exon"
  cat.5hmC[cat.5hmC$id %in% b,"geneID"] <- g

}
rownames(aux) <- NULL
b <- aux[!(rownames(aux) %in% a$queryHits),"id"]
cat.5hmC[cat.5hmC$id %in% b,"type"] <- "intron"
cat.5hmC[is.na(cat.5hmC$type),"type"] <- "intergenic"

b <- cat.5hmC[cat.5hmC$id %in% aux$id & !(cat.5hmC$type %in% "exon"),]
a <- aux[aux$id %in% b$id,]
rownames(a) <- as.character(a$id)
identical(rownames(a),as.character(cat.5hmC[cat.5hmC$id %in% aux$id & !(cat.5hmC$type %in% "exon"),"id"]))
cat.5hmC[cat.5hmC$id %in% aux$id & !(cat.5hmC$type %in% "exon"),"geneID"] <- a$geneID

load("/media/data2/TCGA_RNAseq/gene_annotation_v22_cnv_wgbs.rda")
cat.5hmC <- merge(cat.5hmC,genes[,c("gene_id","RNAseq.hg38.summary")],by.x="geneID",by.y="gene_id",all.x=TRUE,sort=FALSE)
colnames(cat.5hmC)[36] <- "gene.expression"
cat.5hmC <- cat.5hmC[,c(16,3:13,28,17,14:15,1,36,18:27,29,2,30:35)]

cat.5hmC <- cat.5hmC[with(cat.5hmC,order(id)),]
save(cat.5hmC,file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/5hmC/cat_5hmc_v5.rda")


##LSD1
load("/media/data2/TCGA_RNAseq/Michele_norm/RNASeq_norm_test.rda")
load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/LGG-GBM.merged.data.FB20150930v3.Rdata")
load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/LGG_GBM_nonTumorBrain_RNAseqNorm.rda")
load("/media/data2/TCGA_RNAseq/gene_annotation_v22_cnv_wgbs.rda")

LSD1 <- RNASeq.NORM.pv["ENSG00000004487",]
aux <- as.character(subset(pd, cartoon %in% c("GCIMP-low"))$case.id)
low.mrna <- LSD1[,substr(colnames(LSD1),1,12) %in% aux]
aux <- as.character(subset(pd, cartoon %in% c("GCIMP-high"))$case.id)
high.mrna <- LSD1[,substr(colnames(LSD1),1,12) %in% aux]
low.mrna <- melt(low.mrna)
low.mrna$Subtype <- "G-CIMP-low"
high.mrna <- melt(high.mrna)
high.mrna$Subtype <- "G-CIMP-high"
all <- rbind(low.mrna,high.mrna)

grob <- grobTree(textGrob("P-value: 0.000003856", x=0.1,  y=0.95, hjust=0))
ggplot(all,aes(Subtype,value,fill=Subtype)) + geom_boxplot(outlier.shape=NA) + scale_fill_manual(values=c("firebrick","darkgreen")) + geom_jitter() + theme_bw() + ggtitle("Expression levels os LSD1") + annotation_custom(grob)

LSD1 <- rna.seq.LGG.GBM.normalBrain[grep("HDAC*|KDM1A",rownames(rna.seq.LGG.GBM.normalBrain)),]
aux <- as.character(subset(pd, cartoon %in% c("GCIMP-low"))$case.id)
low.mrna <- LSD1[,colnames(LSD1) %in% aux]
low.mrna$gene <- rownames(low.mrna)
aux <- as.character(subset(pd, cartoon %in% c("GCIMP-high"))$case.id)
high.mrna <- LSD1[,colnames(LSD1) %in% aux]
high.mrna$gene <- rownames(high.mrna)
aux <- as.character(subset(pd, cartoon %in% c("IDHmut-codel"))$case.id)
codel.mrna <- LSD1[,colnames(LSD1) %in% aux]
codel.mrna$gene <- rownames(codel.mrna)
aux <- as.character(subset(pd, cartoon %in% c("Classic-like"))$case.id)
clas.mrna <- LSD1[,colnames(LSD1) %in% aux]
clas.mrna$gene <- rownames(clas.mrna)
aux <- as.character(subset(pd, cartoon %in% c("Mesenchymal-like"))$case.id)
mes.mrna <- LSD1[,colnames(LSD1) %in% aux]
mes.mrna$gene <- rownames(mes.mrna)
aux <- as.character(subset(pd, cartoon %in% c("LGm6-GBM"))$case.id)
lgm6.mrna <- LSD1[,colnames(LSD1) %in% aux]
lgm6.mrna$gene <- rownames(lgm6.mrna)
aux <- as.character(subset(pd, cartoon %in% c("PA-like"))$case.id)
pa.mrna <- LSD1[,colnames(LSD1) %in% aux]
pa.mrna$gene <- rownames(pa.mrna)
normal.mrna <- LSD1[,1:5]
normal.mrna$gene <- rownames(normal.mrna)
low.mrna <- melt(low.mrna,id.vars="gene")
low.mrna$Subtype <- "G-CIMP-low"
high.mrna <- melt(high.mrna,id.vars="gene")
high.mrna$Subtype <- "G-CIMP-high"
codel.mrna <- melt(codel.mrna,id.vars="gene")
codel.mrna$Subtype <- "IDHmut-codel"
clas.mrna <- melt(clas.mrna,id.vars="gene")
clas.mrna$Subtype <- "Classic-like"
mes.mrna <- melt(mes.mrna,id.vars="gene")
mes.mrna$Subtype <- "Mesenchymal-like"
lgm6.mrna <- melt(lgm6.mrna,id.vars="gene")
lgm6.mrna$Subtype <- "LGm6-GBM"
pa.mrna <- melt(pa.mrna,id.vars="gene")
pa.mrna$Subtype <- "PA-like"
normal.mrna <- melt(normal.mrna,id.vars="gene")
normal.mrna$Subtype <- "Non-tumor brain"

all <- rbind(low.mrna,high.mrna,codel.mrna,clas.mrna,mes.mrna,lgm6.mrna,pa.mrna,normal.mrna)
all$Subtype <- factor(all$Subtype)
all$Subtype <- factor(all$Subtype, levels = c("Non-tumor brain","G-CIMP-low","G-CIMP-high","IDHmut-codel","Classic-like","Mesenchymal-like","LGm6-GBM","PA-like"))

a <- read.table("/media/data1/thais.projects/Glioma.Rec/Risk_p_by_samples.txt",header=TRUE,sep="\t")
rownames(a) <- as.character(a$Patient.ID)
a <- a[intersect(as.character(a$Patient.ID),as.character(all$variable)),]
all$Risk <- "No.info"
rownames(all) <- as.character(all$variable)
all[intersect(as.character(a$Patient.ID),as.character(all$variable)),"Risk"] <- as.character(a[intersect(as.character(a$Patient.ID),as.character(all$variable)),"Risk.Factor"])

all$gene <- factor(all$gene)
a <- genes[genes$gene_name %in% as.character(all$gene),]
rownames(a) <- as.character(a$gene_name)
a <- a[levels(all$gene),]
a$RNAseqALL.pv.adj <- round(a$RNAseqALL.pv.adj,digits=3)


grob <- grobTree(textGrob("P-value: 0.000003856", x=0.1,  y=0.95, hjust=0))

ggplot(all,aes(Subtype,value,fill=Subtype)) + 
  geom_boxplot(outlier.shape=NA) + 
  scale_fill_manual(values=c("grey","darkgreen","firebrick","purple","orange","yellow","blue","cyan")) + 
  #geom_jitter(aes(color=Risk)) + 
  geom_jitter() +
  scale_color_manual(values=c("black","grey","blue")) + 
  theme_bw() + 
  ggtitle("Expression levels os LSD1 and HDACs") + 
  facet_wrap(~ gene,scales="free") +
  geom_text(label=as.character(a$RNAseqALL.pv.adj), 
            colour="black", vjust=1,hjust=1)
  #annotation_custom(grob)