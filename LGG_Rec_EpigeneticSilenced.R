#### LGG Recurrent - Epigenetically Silenced Method
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGGRecurrent_DNAmeth.Rda") # LGGrecurrent - DNA methylation obj
lggRec_mRNA_IDs <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2_IDS.txt")
lggRec_mRNA_IDs <- subset(lggRec_mRNA_IDs, Comment..TCGA.Data.Type..1 %in% "RSEM_genes")

lggRec_mRNA_IDs$TCGAid <- substr(as.character(lggRec_mRNA_IDs$Comment..TCGA.Barcode.),1,16)
lggRec_mRNA_IDs.sel <- lggRec_mRNA_IDs[lggRec_mRNA_IDs$TCGAid %in% substr(colnames(LGGrecurrent)[5:36],1,16),]

files <- list.files("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/Level_3/")
files <- grep("*rsem.genes.normalized_results",files)
#setdiff(files,lggRec_mRNA_IDs.sel$Derived.Data.File) OK

LGGrec.mrna = NULL
rownames(lggRec_mRNA_IDs.sel) <- as.character(lggRec_mRNA_IDs.sel$Derived.Data.File)
setwd("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/Level_3/")
z = 1 #first file
while (!is.na(files[z])) {
  aux <- lggRec_mRNA_IDs.sel[as.character(files[z]),"Comment..TCGA.Barcode."]
  if(!is.na(aux)){ #do not read metadata files and non-tumor samples
    if(is.null(LGGrec.mrna)){ #read the first sample
      LGGrec.mrna = read.table(files[z],sep="\t",header=TRUE) 
      colnames(LGGrec.mrna)[2] <- as.character(aux)
      LGGrec.mrna <- LGGrec.mrna[,c(1,2)]
    } 
    else{
      aux2 = read.table(files[z], header=TRUE,sep="\t")
      colnames(aux2)[2] = as.character(aux)
      LGGrec.mrna = merge(LGGrec.mrna, aux2[1:2], by="gene_id")
    }
  }
  z = z + 1 #next file
}#end while
# lggRec 32 samples

LGGrec.mrna$entrezID <- str_extract(LGGrec.mrna$gene_id,"[^|]*$")
LGGrec.mrna$gene_symbol <- str_extract(LGGrec.mrna$gene_id,"^*[^|]*")
rownames(LGGrec.mrna) <- as.character(LGGrec.mrna$entrezID)


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#### how to retrieve data from TxDb objects: http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf
#### All TxDb objects are backed by a SQLite database that manages genomic locations and the relationships between pre-processed mRNA transcripts, exons, protein coding sequences, and their related gene identiers.
txdb.gene <- select(TxDb.Hsapiens.UCSC.hg19.knownGene,keys=as.character(LGGrec.mrna$entrezID),columns=c("TXNAME","TXCHROM","TXSTART","TXEND","TXSTRAND"),keytype="GENEID")
txdb.gene <- na.omit(txdb.gene)
txdb.gene <- subset(txdb.gene,TXCHROM %in% paste0("chr",c(1:22,"X","Y")))
rownames(txdb.gene) <- NULL

library(GenomicRanges)
probe.info <- subset(LGGrecurrent[,1:4],Chromosome %in% c(1:22,"X","Y"))
gene.GR <- GRanges(seqnames = as.character(txdb.gene$TXCHROM),ranges = IRanges(start = txdb.gene$TXSTART, end=txdb.gene$TXEND), strand = txdb.gene$TXSTRAND, idGene = txdb.gene$GENEID, txName = txdb.gene$TXNAME)

probe.info.gr <- GRanges(seqnames = paste0("chr",probe.info$Chromosome), ranges = IRanges(start = probe.info$Genomic_Coordinate, end = probe.info$Genomic_Coordinate), probeID = probe.info$Composite.Element.REF)
distance <- as.data.frame(distanceToNearest(probe.info.gr,gene.GR)) #closest gene to each 450k probe

gene.order.by.distance <- txdb.gene[distance$subjectHits,]
gene.order.by.distance$distance <- as.matrix(distance$distance)
LGGrec.nearest.gene <- subset(LGGrecurrent,Chromosome %in% c(1:22,"X","Y"))
LGGrec.nearest.gene <- cbind(LGGrec.nearest.gene, gene.order.by.distance[,c("GENEID","distance")])
probe.info <- LGGrec.nearest.gene[,c(1:4,37:38)]
colnames(LGGrec.nearest.gene) <- substr(colnames(LGGrec.nearest.gene),1,16)
LGGrec.nearest.gene <- LGGrec.nearest.gene[,substr(colnames(LGGrec.mrna)[2:33],1,16)]
LGGrec.nearest.gene <- cbind(probe.info,LGGrec.nearest.gene)
colnames(LGGrec.mrna)[2:33] <- substr(colnames(LGGrec.mrna)[2:33],1,16)
cpg.gene <- merge(LGGrec.nearest.gene,LGGrec.mrna[,2:35],by.x="GENEID",by.y="entrezID")
dimnames(cpg.gene)[[1]] <- paste(cpg.gene[,"Composite.Element.REF"], cpg.gene[,"gene_symbol"], sep=".")
#(d[,7:38],d[,39:70]) #mDNA methylation x Gene expression
cpg.gene <- na.omit(cpg.gene)

#aux <- subset(codel, hypomethylated_enhancer_enrichment %in% "Enhancer depleted" & Tumor.status %in% c("Primary","First Recurrence"))
aux <- subset(codel, Tumor.status %in% c("Primary","First Recurrence"))
meth.g <- cpg.gene[,7:38] #only tumor
names(meth.g) <- substr(names(meth.g), 1, 16)
expr.g <- cpg.gene[,39:70] #only tumor
names(expr.g) <- substr(names(expr.g), 1, 16)
#identical(colnames(meth.g),colnames(expr.g))
meth.g <- meth.g[,rownames(aux)]
expr.g <- expr.g[,rownames(aux)]

epiGene.LGGrec <- list()
j=1
for(i in 1:dim(cpg.gene)[1]) {
  #for(i in 138:250) {
  
  meth.gu <- names(meth.g)[meth.g[i, ] < 0.3]
  meth.gm <- names(meth.g)[meth.g[i, ] >= 0.3]
  # if(fisher.test(table(yy))$p.value < 0.05){
  if(length(meth.gu)>1 & length(meth.gm)>1){
    expr.mean.gm <- mean(as.numeric(expr.g[i, meth.gm]), na.rm=T)
    if(expr.mean.gm < 1.28*sd(as.numeric(expr.g[i, meth.gu]), na.rm=T)){
      expr.mean.gu <- mean(as.numeric(expr.g[i, meth.gu]), na.rm = T)
      yy <- cbind(t(meth.g[i, ] >= 0.3), t(expr.g[i, ] < expr.mean.gu))
      colnames(yy) <- c("meth.gt0.3", "exp.lt.meanUnmeth")
      yy <- as.data.frame(yy)
      gene <- rownames(meth.g[i,])
      if(as.matrix(table(yy$meth.gt0.3, yy$exp.lt.meanUnmeth))[2,2]/as.numeric(table(yy$meth.gt0.3)[2])>0.8){
        yy$geneName <- gene
        yy$id <- rownames(yy)
        yy$methValues <- as.numeric(meth.g[i, ])
        yy$exprValues <- as.numeric(expr.g[i, ])  
        yy$group <- as.character(codel[rownames(yy),"hypomethylated_enhancer_enrichment"])
        #info <- metadata[rownames(yy),c("cluster.meth","IDH.status","codel1p19q","cluster.meth.gbm","cluster.meth.lgg","tumor.type")]
        #yy$tumor.type <- as.character(info$tumor.type)
        #yy$codel1p19q <- as.character(info$codel1p19q)
        #yy <- merge(metadata, yy, by.x = "TCGAIDnew", by.y = "id")
        yy$status <- ifelse(substr(rownames(yy),14,15)=="01","primary","recurrent")
        epiGene.LGGrec[[j]] <- as.data.frame(yy)
        names(epiGene.LGGrec)[[j]] <- as.character(gene)
        j <- j+1
        cat("found a ES", j, "\n")
      }
      else{
        cat("no step3", "\n")
      }
      
    }else{
      cat("no ES group", "\n")
    }
  }else{
    cat("no methylated group", i, "\n")
  }
}  


uu.nl <- epiGene.LGGrec[[1]]
clusters <- c("primary","recurrent")
#clusters <- "LGm1"


enrichment <- lapply(clusters, function(cluster, epiGene.ss=epiGene.LGGrec) {
  enrichment.group <- unlist(mclapply(epiGene.ss, function(uu.nl, cluster.group=cluster) {
#enrichment <- unlist(mclapply(epiGene.LGGrec,function(uu.nl, cluster.group=clusters){ 
   # uu.nl <- subset(uu.nl, group %in% "Enriched")
    uu.nl$es <- uu.nl$meth.gt0.3 & uu.nl$exp.lt.meanUnmeth
    uu.nl$Kclass <- NA
    uu.nl$Kclass <- uu.nl$status == cluster.group
    es.table <- table(uu.nl$Kclass, uu.nl$es)#, dnn=c("Class", "ES"))
    #print()
    if((es.table[2,2]/rowSums(es.table)[2] > 0.5) & (es.table[2,2]/colSums(es.table)[2] > 0.5)) {
      return(fisher.test(table(uu.nl$Kclass, uu.nl$es))$p.value)
    } else {
      return(NULL)
    }
    #}
    #}  
  }, mc.cores=1))
  return(enrichment.group)
})


names(enrichment) <- c("primary","recurrent")
#enrichment.adj <- lapply(enrichment, p.adjust)
#enrichmentList.sig <- lapply(enrichment.adj, function(x) {
#  x <- x[x < 0.05]
#})
enrichmentList.sig <- lapply(enrichment, function(x) {
  x <- x[x < 0.05]
})

a <- names(c(enrichmentList.sig$primary,enrichmentList.sig$recurrent))
#aux <- subset(codel, hypomethylated_enhancer_enrichment %in% "Enhancer depleted")
meth.g <- cpg.gene[,7:38] #only tumor
names(meth.g) <- substr(names(meth.g), 1, 16)
expr.g <- cpg.gene[,39:70] #only tumor
names(expr.g) <- substr(names(expr.g), 1, 16)
#meth.g <- meth.g[,rownames(aux)]
#expr.g <- expr.g[,rownames(aux)]

LGGrec.meth.pri <- meth.g[as.character(a),rownames(subset(aux,Tumor.status=="Primary"))]
order.pri.meth <- hclust(dist(t(LGGrec.meth.pri)))
LGGrec.meth.rec <- meth.g[as.character(a),rownames(subset(aux,Tumor.status!="Primary"))]
order.rec.meth <- hclust(dist(t(LGGrec.meth.rec)))

LGGrec.meth <- cbind(LGGrec.meth.pri[,order.pri.meth$order],LGGrec.meth.rec[,order.rec.meth$order])

#clab <- c(rep("deeppink",5),rep("darkturquoise",7))
rlab <- c(rep("red",107),rep("green",8)) #summary(enrichmentList.dep)

#order.rec.row <- hclust(dist(LGGrec.meth.rec))
clab <- codel
rownames(clab) <- as.character(clab$bcr.tcga)
clab$Patient.Id <- as.factor(clab$Patient.Id)
clab$Tumor.status <- as.factor(clab$Tumor.status)
clab$IDH1_1p19qCodeletion <- as.factor(clab$IDH1_1p19qCodeletion)
clab$hypomethylated_enhancer_enrichment <- as.factor(clab$hypomethylated_enhancer_enrichment)
clab$Histology <- as.factor(clab$Histology)
clab$Grade <- as.factor(clab$Grade)
levels(clab$Patient.Id) <- c("slateblue4", "violet", "turquoise3", "yellow2", "snow3", "purple", "tan1", "slategray1", "yellowgreen", "pink4", "red3", "paleturquoise4", "darkgreen", "black")
clab$Primary_to_recurrence <- as.numeric(clab$Primary_to_recurrence)
clab$Primary_to_recurrence <- jet.colors(100)[ceiling((as.numeric(clab$Primary_to_recurrence)-min(as.numeric(clab$Primary_to_recurrence),na.rm=T))/(max(as.numeric(clab$Primary_to_recurrence),na.rm=T)-min(as.numeric(clab$Primary_to_recurrence),na.rm=T))*(length(jet.colors(100))-1))+1]
clab$Primary_to_recurrence[is.na(clab$Primary_to_recurrence)] <- "white"
levels(clab$Tumor.status) <- c("deeppink","darkturquoise","brown")
levels(clab$Grade) <- c("gray66","gray40","black")
levels(clab$Histology) <- c("blue","green","red","blue","red","yellow","green","red")
levels(clab$IDH1_1p19qCodeletion) <- c("indianred1","hotpink4","lightblue4")
levels(clab$hypomethylated_enhancer_enrichment) <- c("orange3","royalblue2","white")


z <- seq(min(clab$Primary_to_recurrence,na.rm=TRUE), max(clab$Primary_to_recurrence,na.rm=TRUE), length = 200) #min(lgg.gbm.27.450.order) and #max(lgg.gbm.27.450.order)
n <- max(clab$Primary_to_recurrence,na.rm=TRUE)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/key_TimeToRecurrence2.pdf",width=10, height=05,bg = "transparent")
image(matrix(z, ncol = 1), col = jet.colors(75),
      xaxt = "n", yaxt = "n", main = "Overlap count Key")
box()
par(usr = c(0, n, 0, n))
axis(1, at = c(0,max(z)), c(min(z),max(z)))
dev.off()

clab <- clab[colnames(LGGrec.meth),]

aux <- hclust(dist(LGGrec.meth[names(enrichmentList.sig$primary),]))
aux2 <- hclust(dist(LGGrec.meth[names(enrichmentList.sig$recurrent),]))
aux <- c(aux$order,(aux2$order+max(aux$order)))

png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/Heatmap_ES_meth_ALL2.png",res=400,width=3000,height=3000)
heatmap.plus.sm(as.matrix(LGGrec.meth[aux,]),
                col = jet.colors(75),
                scale = "none",
                labRow = str_extract(rownames(LGGrec.meth[aux,]),"^*[^.]*"),
                cexRow = 0.8,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = as.matrix(clab[,c(2,3,4,6,11,12,13)]),
               # margins = c(1, 6)
                RowSideColors = cbind(rlab,rlab)
                
)
dev.off()

aux <- subset(codel, hypomethylated_enhancer_enrichment %in% "Enhancer depleted")
#aux <- codel
a <- names(c(enrichmentList.sig$primary,enrichmentList.sig$recurrent))
LGGrec.mrna <- expr.g[a,colnames(LGGrec.meth)]
LGGrec.mrna <- expr.g[a,as.character(aux$bcr.tcga)]
LGGrec.mrna$meanPri <- apply(LGGrec.mrna[,rownames(subset(aux,Tumor.status=="Primary"))],1,mean,na.rm=TRUE)
LGGrec.mrna$meanRec <- apply(LGGrec.mrna[,rownames(subset(aux,Tumor.status!="Primary"))],1,mean,na.rm=TRUE)
LGGrec.mrna <- LGGrec.mrna[,c("meanPri","meanRec")]
LGGrec.mrna <- LGGrec.mrna + 1
LGGrec.mrna <- log2(LGGrec.mrna) 

aux <- hclust(dist(LGGrec.meth[names(enrichmentList.sig$primary),]))
aux2 <- hclust(dist(LGGrec.meth[names(enrichmentList.sig$recurrent),]))
aux <- c(aux$order,(aux2$order+max(aux$order)))

png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/Heatmap_ES_mRNA_ALL2.png",res=400,width=3000,height=3000)
heatmap.plus.sm(as.matrix(LGGrec.mrna[rownames(LGGrec.meth),]),
                col=greenred(1000),
                scale = "none",
                #trace = "none",
                labRow = str_extract(rownames(LGGrec.mrna[rownames(LGGrec.meth),]),"[^.]*$"),
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = cbind(c("darkturquoise","deeppink"),c("darkturquoise","deeppink")),
                RowSideColors = cbind(rlab,rlab)
                
)
dev.off()

### Epiplot
mean.meth.ES.pri <- apply(meth.g[names(c(enrichmentList.sig$primary)),], 2, mean, na.rm=TRUE)
mean.meth.ES.rec <- apply(meth.g[names(c(enrichmentList.sig$recurrent)),], 2, mean, na.rm=TRUE)

mean.expr.ES.pri <- apply(expr.g[names(c(enrichmentList.sig$primary)),], 2, mean, na.rm=TRUE)
mean.expr.ES.rec <- apply(expr.g[names(c(enrichmentList.sig$recurrent)),], 2, mean, na.rm=TRUE)

epiplot.ES.pri <- cbind(mean.meth.ES.pri,mean.expr.ES.pri)
epiplot.ES.rec <- cbind(mean.meth.ES.rec,mean.expr.ES.rec)

epiplot.ES.pri <- as.data.frame(epiplot.ES.pri)
epiplot.ES.rec <- as.data.frame(epiplot.ES.rec)

colnames(epiplot.ES.pri) <- c("methValues","exprValues")
colnames(epiplot.ES.rec) <- c("methValues","exprValues")


epiplot.ES.pri <- merge(epiplot.ES.pri,yy[,c("id","status")],by=0,all.x=T)
epiplot.ES.rec <- merge(epiplot.ES.rec,yy[,c("id","status")],by=0,all.x=T)

epiplot.ES.pri$group <- "ES.pri"
epiplot.ES.rec$group <- "ES.rec"

epiplot <- rbind(epiplot.ES.pri,epiplot.ES.rec)
g1 <- ggplot(epiplot, aes(methValues, exprValues)) +
  stat_smooth(method = "lm", fullrange=FALSE, na.rm=TRUE, se=FALSE, size=1) +
  geom_vline(xintercept = 0.3, linetype = "longdash", colour="black") +
  #geom_hline(yintercept = expr.mean.gu, linetype = "longdash", colour="black") +
  geom_point(size=3, aes(colour = factor(status))) +
  scale_colour_manual(values = c("deeppink","darkturquoise")) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_log10() +
  #scale_size("cluster.gbm") +
  ylab("mean RNAseq") +
  xlab(paste("DNA Methylation Values\n")) + 
  labs(title = "", colour = "Status") +
  #theme(legend.position = "none") +
  facet_grid(. ~ group,scales="free", space="free")
ggsave(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/epiplot.pdf")


### Enhancer Groups (Enriched, depleted, unknown)
grp <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/LGGrec_EnhancerGroups.txt")
yy$ID <- substr(rownames(yy),9,12)
yy$barcode <- rownames(yy)
grp <- merge(grp,yy[,c("barcode","ID")],by="ID")
rownames(grp) <- as.character(grp$barcode)
grp$Status <- substr(grp$barcode,14,16)

### Clinical Data - 12-08-2015
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/LGGrec_ClinicalData_12082015.Rda")
codel <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/LGG_codel.txt")
rownames(codel) <- as.character(codel$bcr.tcga)
codel$id <- substr(codel$bcr.tcga,9,12)
codel <- subset(codel, id %in% as.character(TCGA_LGG$Patient.Id))
TCGA_LGG$chr.del <- NA
aux <- codel[codel$X19q < -0.3 & codel$X1p > -0.3,"id"]
TCGA_LGG[TCGA_LGG$Patient.Id %in% aux,"chr.del"] <- "Del.19q"
#aux <- codel[codel$X19q > -0.3 & codel$X1p < -0.3,"id"]
#TCGA_LGG[TCGA_LGG$Patient.Id %in% aux,"chr.del"] <- "Del.1p"
aux <- codel[codel$subtype %in% "IDHmut-codel","id"] 
TCGA_LGG[TCGA_LGG$Patient.Id %in% aux,"chr.del"] <- "Codel"
write.table(TCGA_LGG,file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/LGG_codel.txt",quote=F,row.names=F,sep="\t")


### KEY
z <- seq(min(LGGrec.mrna,na.rm=TRUE), max(LGGrec.mrna,na.rm=TRUE), length = 200) #min(lgg.gbm.27.450.order) and #max(lgg.gbm.27.450.order)
n <- max(LGGrec.mrna,na.rm=TRUE)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/key_RNAseqMean.pdf",width=10, height=05,bg = "transparent")
image(matrix(z, ncol = 1), col = greenred(75),
      xaxt = "n", yaxt = "n", main = "Overlap count Key")
box()
par(usr = c(0, n, 0, n))
axis(1, at = c(0,max(z)), c(min(z),max(z)))
dev.off()

z <- seq(min(LGGrec.meth,na.rm=TRUE), max(LGGrec.meth,na.rm=TRUE), length = 200) #min(lgg.gbm.27.450.order) and #max(lgg.gbm.27.450.order)
n <- max(LGGrec.meth,na.rm=TRUE)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/key_DNAMeth.pdf",width=10, height=05,bg = "transparent")
image(matrix(z, ncol = 1), col = jet.colors(75),
      xaxt = "n", yaxt = "n", main = "Overlap count Key")
box()
par(usr = c(0, n, 0, n))
axis(1, at = c(0,max(z)), c(min(z),max(z)))
dev.off()


######## Radiation - RNAseq x DNA methylation

###Radiation cases: 6397,7304,5963,5870,5965,6404,A7RV,6407,5872
###Radiation cases: 6397,7304,5963,5870,5965,6404,A7RV,6407,5872


load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGGRecurrent_DNAmeth.Rda") # LGGrecurrent - DNA methylation obj
lggRec_mRNA_IDs <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2_IDS.txt")
lggRec_mRNA_IDs <- subset(lggRec_mRNA_IDs, Comment..TCGA.Data.Type..1 %in% "RSEM_genes_normalized")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/LGG_rad_pvalues.rda")


lggRec_mRNA_IDs$TCGAid <- substr(as.character(lggRec_mRNA_IDs$Comment..TCGA.Barcode.),1,16)
lggRec_mRNA_IDs.sel <- lggRec_mRNA_IDs[lggRec_mRNA_IDs$TCGAid %in% substr(colnames(LGGrecurrent)[5:36],1,16),]


files <- list.files("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/RNASeqV2_new/UNC__IlluminaHiSeq_RNASeqV2/Level_3")
files <- files[grep("*rsem.genes.normalized_results",files)]



LGGrec.mrna = NULL
rownames(lggRec_mRNA_IDs.sel) <- as.character(lggRec_mRNA_IDs.sel$Derived.Data.File)
setwd("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/RNASeqV2_new/UNC__IlluminaHiSeq_RNASeqV2/Level_3")
z = 1 #first file
while (!is.na(files[z])) {
  aux <- lggRec_mRNA_IDs.sel[as.character(files[z]),"Comment..TCGA.Barcode."]
  if(!is.na(aux)){ #do not read metadata files and non-tumor samples
    if(is.null(LGGrec.mrna)){ #read the first sample
      LGGrec.mrna = read.table(files[z],sep="\t",header=TRUE) 
      colnames(LGGrec.mrna)[2] <- as.character(aux)
      LGGrec.mrna <- LGGrec.mrna[,c(1,2)]
    } 
    else{
      aux2 = read.table(files[z], header=TRUE,sep="\t")
      colnames(aux2)[2] = as.character(aux)
      LGGrec.mrna = merge(LGGrec.mrna, aux2[1:2], by="gene_id")
    }
  }
  z = z + 1 #next file
}#end while
# lggRec 32 samples

### map gene ID to genomic coordinates
gene.location <- read.delim("/dados/ResearchProjects/thais/TCGA/Normals/mRNA_Brain/Galaxy1-[Homo_sapiens_genes_(GRCh37.p13)].tabular")
gene.location <- gene.location[gene.location$Chromosome.Name %in% c(1:22,"X","Y"),]
gene.location <- gene.location[!is.na(gene.location$EntrezGene.ID),]
gene.location <- gene.location[!duplicated(gene.location$EntrezGene.ID),] #the duplicates are different transcripts, not different coordinates
library(GenomicRanges)
gene.GR <- GRanges(seqnames = paste0("chr",gene.location$Chromosome.Name),ranges = IRanges(start = gene.location$Gene.Start..bp., end=gene.location$Gene.End..bp.), strand = gene.location$Strand, symbol = gene.location$Associated.Gene.Name, EntrezID = gene.location$EntrezGene.ID)
probe.info <- GRanges(seqnames = paste0("chr",LGGrecurrent$Chromosome), ranges = IRanges(start = LGGrecurrent$Genomic_Coor, end = LGGrecurrent$Genomic_Coor), probeID = LGGrecurrent$Composite.Element.REF)
distance <- as.data.frame(distanceToNearest(probe.info,gene.GR)) #closest gene to each 450k probe
rownames(gene.location) <- NULL
gene.order.by.distance <- gene.location[distance$subjectHits,]
gene.order.by.distance$distance <- as.matrix(distance$distance)
Recurrent.nearest.gene <- cbind(LGGrecurrent[,1:4], gene.order.by.distance[,c("Associated.Gene.Name","EntrezGene.ID","distance")])

LGG.rad.s <- subset(LGG.rad.volcano, p.value < 0.04 & (diffmean < -0.22 | diffmean > 0.34))
LGG.rad.s <- cbind(LGG.rad.s,Recurrent.nearest.gene[rownames(LGG.rad.s),c("Associated.Gene.Name","EntrezGene.ID")])

a <- LGGrecurrent[rownames(LGG.rad.s),]
LGG.rad.GR <- GRanges(seqnames = paste0("chr",a$Chromosome), ranges = IRanges(start = a$Genomic_Coordinate-3000, end = a$Genomic_Coordinate+3000), probeID = a$Composite.Element.REF)

#### find the top 10 genes that precede a CpG site (450k)
LGGrec.mrna$entrezID <- str_extract(LGGrec.mrna$gene_id,"[^|]*$")
LGGrec.mrna$gene_symbol <- str_extract(LGGrec.mrna$gene_id,"^*[^|]*")
rownames(LGGrec.mrna) <- as.character(LGGrec.mrna$entrezID)
gene.location.RNAseq <- subset(gene.location, EntrezGene.ID %in% rownames(LGGrec.mrna))
genes.RNAseq.GR <- GRanges(seqnames = paste0("chr",gene.location.RNAseq$Chromosome.Name),ranges = IRanges(start = gene.location.RNAseq$Gene.Start..bp., end=gene.location.RNAseq$Gene.End..bp.), strand = 
                             gene.location.RNAseq$Strand, symbol = gene.location.RNAseq$Associated.Gene.Name, EntrezID = gene.location.RNAseq$EntrezGene.ID)
rownames(gene.location.RNAseq) <- as.character(gene.location.RNAseq$EntrezGene.ID)
genes.precede <- matrix(data=NA,nrow=133,ncol=10,dimnames=list(LGG.rad.GR$probeID,paste("Precede.round",1:10,sep=".")))

for(i in 1:10){
  genes.RNAseq.GR <- GRanges(seqnames = paste0("chr",gene.location.RNAseq$Chromosome.Name),ranges = IRanges(start = gene.location.RNAseq$Gene.Start..bp., end=gene.location.RNAseq$Gene.End..bp.), strand = gene.location.RNAseq$Strand, symbol = gene.location.RNAseq$Associated.Gene.Name, EntrezID = gene.location.RNAseq$EntrezGene.ID)
  aux <- precede(LGG.rad.GR,genes.RNAseq.GR)
  genes.precede[,i] <- as.character(gene.location.RNAseq[aux,"EntrezGene.ID"])
  gene.location.RNAseq <- gene.location.RNAseq[-aux,]
  
}


gene.location.RNAseq <- subset(gene.location, EntrezGene.ID %in% rownames(LGGrec.mrna))
genes.RNAseq.GR <- GRanges(seqnames = paste0("chr",gene.location.RNAseq$Chromosome.Name),ranges = IRanges(start = gene.location.RNAseq$Gene.Start..bp., end=gene.location.RNAseq$Gene.End..bp.), strand = 
                             gene.location.RNAseq$Strand, symbol = gene.location.RNAseq$Associated.Gene.Name, EntrezID = gene.location.RNAseq$EntrezGene.ID)
rownames(gene.location.RNAseq) <- as.character(gene.location.RNAseq$EntrezGene.ID)
genes.follow <- matrix(data=NA,nrow=133,ncol=10,dimnames=list(LGG.rad.GR$probeID,paste("Follow.round",1:10,sep=".")))

for(i in 1:10){
  genes.RNAseq.GR <- GRanges(seqnames = paste0("chr",gene.location.RNAseq$Chromosome.Name),ranges = IRanges(start = gene.location.RNAseq$Gene.Start..bp., end=gene.location.RNAseq$Gene.End..bp.), strand = gene.location.RNAseq$Strand, symbol = gene.location.RNAseq$Associated.Gene.Name, EntrezID = gene.location.RNAseq$EntrezGene.ID)
  aux <- follow(LGG.rad.GR,genes.RNAseq.GR)
  genes.follow[,i] <- as.character(gene.location.RNAseq[aux,"EntrezGene.ID"])
  gene.location.RNAseq <- gene.location.RNAseq[-aux,]
  
}

genes.all <- cbind(genes.precede,genes.follow)

library(reshape)

#color -> radiation
#shape -> status primary and recurrent
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/gbm.lgg.recurrent.meta.rda")
rownames(gbm.lgg.recurrent.meta) <- as.character(gbm.lgg.recurrent.meta$Patient.ID)
gbm.lgg.recurrent.meta$radiation <- "N"
gbm.lgg.recurrent.meta[c("TCGA-DU-6397","TCGA-DU-7304","TCGA-FG-5963","TCGA-DU-5870","TCGA-FG-5965","TCGA-DU-6404","TCGA-TQ-A7RV","TCGA-DU-6407","TCGA-DU-5872"),"radiation"] <- "Y"
meta <- subset(gbm.lgg.recurrent.meta,type %in% "LGG" & Age.Diagnosis < 70)

setwd("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/RNASeqV2_new/plots")
for(i in 1:nrow(genes.all)){
  exp <- LGGrec.mrna[as.character(genes.all[i,]),substr(colnames(LGGrec.mrna),1,12) %in% rownames(meta) & substr(colnames(LGGrec.mrna),14,16) != "02B"]
  exp$entrezID <- as.character(genes.all[i,])
  dna.m <- LGGrecurrent[rownames(genes.all)[i],substr(colnames(LGGrecurrent),1,12) %in% rownames(meta) & substr(colnames(LGGrecurrent),14,16) != "02B"]
  dna.m <- t(dna.m)
  rownames(dna.m) <- substr(rownames(dna.m),1,16)
  df <- melt(exp, id = "entrezID")
  df$status <- substr(df$variable,14,15)
  df$variable <- substr(df$variable,1,16)
  df <- merge(df,dna.m,by.x="variable",by.y=0)
  df$variable <- substr(df$variable,1,12)
  df <- merge(df,meta[,c("Patient.ID","radiation")],by.x="variable",by.y="Patient.ID")
  #df$entrezID <- as.character(gene.location.RNAseq[df$entrezID,"Associated.Gene.Name"])
  filename <- paste0(rownames(genes.all)[i],".pdf")
  df$value <- df$value + 1
  df$value <- log2(df$value)
  colnames(df)[5] <- "probe"
  ggplot(df, aes(x= probe, y=value, color=factor(radiation))) + 
    geom_point(size=3,aes(shape=factor(status)))+ 
    stat_smooth(method = "lm", fullrange=FALSE, na.rm=TRUE, se=FALSE, size=0.5,colour="black") +
    geom_vline(xintercept = 0.5, linetype = "longdash", colour="black") +
    facet_wrap( ~ entrezID,ncol=5)+ #scales="free_y"
    scale_colour_manual(values = c("forestgreen","darkorchid1"), 
                        labels=c("No Radiation","Radiation"),
                        name="Biological Group") + 
    scale_shape_discrete(name  ="Sample Type",
                         breaks=c("01", "02"),
                         labels=c("Primary", "First Recurrence")) +
    ylab("RNASeq Expression (log2)") +
    xlab("DNA Methylation Values") + 
    scale_x_continuous(limits=c(0,1),breaks=c(0,0.25,0.5,0.75,1)) + 
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  ggsave(filename = filename)
}


probes <- c("cg20852945","cg01690062","cg01715901","cg02162404","cg03950246","cg04237608","cg08439109")
genes <- c("4192","317761","11227","2888","4016","4237","30009")
aux <- LGGrecurrent[probes,substr(colnames(LGGrecurrent),14,15) %in% "01"]
aux.Y <- aux[,substr(colnames(aux),1,12) %in% rownames(subset(meta, radiation == "Y"))]
aux.N <- aux[,substr(colnames(aux),1,12) %in% rownames(subset(meta, radiation == "N"))]
p <- cbind(aux.Y,aux.N)
clab <- c(rep("plum1",ncol(aux.Y)),rep("antiquewhite4",ncol(aux.N)))
aux <- LGGrecurrent[probes,substr(colnames(LGGrecurrent),14,16) %in% "02A"]
aux <- aux[,substr(colnames(aux),1,12) %in% substr(colnames(p),1,12)]
clab <- cbind(clab,c(rep("slateblue4",ncol(p)),rep("darkslategray",ncol(aux))))
meth <- cbind(p,aux)

ordem <- c("6397-01A","7304-01A","5963-01A","5870-01A","5965-01B","6404-01A","A7RV-01A","6407-01A","5872-01A","A7CF-01A","A8XE-01A","A7RK-01A","A4MT-01A","6397-02A","7304-02A","5963-02A","5870-02A","5965-02A","6404-02A","A7RV-02A","6407-02A","5872-02A","A7CF-02A","A8XE-02A","A7RK-02A","A4MT-02A")

colnames(meth) <- substr(colnames(meth),9,16)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/RNASeqV2_new/plots/Heatmap_DNAm_7p.pdf")
heatmap.plus.sm(as.matrix(meth[,ordem]),
                col = jet.colors(75),
                scale = "none",
                #labRow = str_extract(rownames(LGGrec.meth[aux,]),"^*[^.]*"),
                cexRow = 0.8,
               # labCol = substr(colnames(p),6,16),
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab
                # margins = c(1, 6)
                #RowSideColors = cbind(rlab,rlab)
)
dev.off()

aux <- LGGrec.mrna[genes,substr(colnames(LGGrec.mrna),14,15) %in% "01"]
aux.Y <- aux[,substr(colnames(aux),1,12) %in% rownames(subset(meta, radiation == "Y"))]
aux.N <- aux[,substr(colnames(aux),1,12) %in% rownames(subset(meta, radiation == "N"))]
mean.Y <- apply(aux.Y,1,mean,na.rm=TRUE)
mean.N <- apply(aux.N,1,mean,na.rm=TRUE)
p <- cbind(mean.Y,mean.N)
aux <- LGGrec.mrna[genes,substr(colnames(LGGrec.mrna),14,16) %in% "02A"]
aux.Y <- aux[,substr(colnames(aux),1,12) %in% rownames(subset(meta, radiation == "Y"))]
aux.N <- aux[,substr(colnames(aux),1,12) %in% rownames(subset(meta, radiation == "N"))]
mean.Y <- apply(aux.Y,1,mean,na.rm=TRUE)
mean.N <- apply(aux.N,1,mean,na.rm=TRUE)
p <- cbind(p,mean.Y,mean.N)
#clab <- c(rep("plum1",ncol(mean.Y)),rep("antiquewhite4",ncol(mean.N)))
clab <- c("plum1","antiquewhite4","plum1","antiquewhite4")


p <- p + 1
p <- log2(p)


pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/RNASeqV2_new/plots/Heatmap_mRNA_7p_mean.pdf")
heatmap.plus.sm(as.matrix(p),
                col = greenred(1000),
                scale = "row",
                #labRow = str_extract(rownames(LGGrec.meth[aux,]),"^*[^.]*"),
                cexRow = 0.8,
                # labCol = substr(colnames(p),6,16),
                Colv = NA,
                Rowv = NA,
                ColSideColors = cbind(clab,clab)
                # margins = c(1, 6)
                #RowSideColors = cbind(rlab,rlab)
)
dev.off()

## Key
z <- seq(min(p,na.rm=TRUE), max(p,na.rm=TRUE), length = 200) #min(lgg.gbm.27.450.order) and #max(lgg.gbm.27.450.order)
n <- max(p,na.rm=TRUE)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/RNASeqV2_new/plots/key_RNAseq_7p_mean.pdf",width=10, height=05,bg = "transparent")
image(matrix(z, ncol = 1), col = greenred(75),
      xaxt = "n", yaxt = "n", main = "Overlap count Key")
box()
par(usr = c(0, n, 0, n))
axis(1, at = c(0,max(z)), c(min(z),max(z)))
dev.off()


#### Gviz
library(Gviz)
aux <- LGGrec.mrna[genes,substr(colnames(LGGrec.mrna),14,15) %in% "01"]
aux.Y <- aux[,substr(colnames(aux),1,12) %in% rownames(subset(meta, radiation == "Y"))]
aux.N <- aux[,substr(colnames(aux),1,12) %in% rownames(subset(meta, radiation == "N"))]
id.p.Y <- substr(colnames(aux.Y),1,16)
id.p.N <- substr(colnames(aux.N),1,16)
aux <- LGGrec.mrna[genes,substr(colnames(LGGrec.mrna),14,16) %in% "02A"]
aux.Y <- aux[,substr(colnames(aux),1,12) %in% rownames(subset(meta, radiation == "Y"))]
aux.N <- aux[,substr(colnames(aux),1,12) %in% rownames(subset(meta, radiation == "N"))]
id.r.Y <- substr(colnames(aux.Y),1,16)
id.r.N <- substr(colnames(aux.N),1,16)

probes <- c("cg20852945","cg01690062","cg01715901","cg02162404","cg03950246","cg04237608","cg08439109")


probes.df <- LGGrecurrent[probes,1:4]
probes.df <- probes.df[3,]
probes.GR <- GRanges(seqnames = paste0("chr",probes.df$Chromosome),ranges = IRanges(start = probes.df$Genomic_Coordinate, end=probes.df$Genomic_Coordinate), probeID = probes.df$Composite.Element.REF, EntrezID = genes[3])
genome(probes.GR) <- "hg19"
atrack <- AnnotationTrack(probes.GR, name = "CpG")
plotTracks(atrack)
gtrack <- GenomeAxisTrack()
plotTracks(list(gtrack, atrack))
gen <- genome(probes.GR)
chr <- as.character(unique(seqnames(probes.GR)))
itrack <- IdeogramTrack(genome = gen, chromosome = chr)
plotTracks(list(itrack, gtrack, atrack))
library("biomaRt")
#mart = useMart("ensembl")
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",  path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
ensembl = useDataset("hsapiens_gene_ensembl", mart = mart)
gene.location.df <- getBM(attributes=c("ensembl_exon_id","exon_chrom_start","exon_chrom_end","ensembl_transcript_id","ensembl_gene_id","chromosome_name","strand","start_position","end_position"), filters = "entrezgene", values = genes[3], mart = ensembl)
gene.location.df <- gene.location.df[,c(6,2,3,7,5,1,4)]
colnames(gene.location.df) <- c("chromosome","start","end","strand","gene","exon","transcript")
gene.location.df$chromosome <- paste0("chr",gene.location.df$chromosome)
gene.location.df$strand <- ifelse(gene.location.df$chromosome == "1","+","-")

grtrack <- GeneRegionTrack(gene.location.df, genome = gen, chromosome = chr, name = "Gene")
plotTracks(list(itrack, gtrack, atrack, grtrack))


mean.Y.p <- mean(as.numeric(LGGrecurrent[probes[3],substr(colnames(LGGrecurrent),1,16) %in% id.p.Y]),na.rm=TRUE)
mean.N.p <- mean(as.numeric(LGGrecurrent[probes[3],substr(colnames(LGGrecurrent),1,16) %in% id.p.N]),na.rm=TRUE)
mean.Y.r <- mean(as.numeric(LGGrecurrent[probes[3],substr(colnames(LGGrecurrent),1,16) %in% id.r.Y]),na.rm=TRUE)
mean.N.r <- mean(as.numeric(LGGrecurrent[probes[3],substr(colnames(LGGrecurrent),1,16) %in% id.r.N]),na.rm=TRUE)
probes.df <- LGGrecurrent[probes,1:4]
probes.df <- probes.df[3,]
probes.GR <- GRanges(seqnames = paste0("chr",probes.df$Chromosome),ranges = IRanges(start = probes.df$Genomic_Coordinate, end=probes.df$Genomic_Coordinate),mean.Y.p=mean.Y.p,mean.N.p=mean.N.p,mean.Y.r=mean.Y.r,mean.N.r=mean.N.r)

genome(probes.GR) <- "hg19"


dTrack <- DataTrack(probes.GR, name = "DNA methylation")
#plotTracks(list(itrack, gtrack, atrack, grtrack,dTrack), groups = c("Rad.p","NonRad.p","Rad.r","NonRad,r"), type = c("a", "p", "confint"))

plotTracks(list(itrack, gtrack, atrack, grtrack,dTrack), groups = c("Rad.p","NonRad.p","Rad.r","NonRad.r"), type = c("a", "p", "confint"),legend=TRUE,ylim=c(0,1))




######### GVIZ
library(Gviz)
aux <- LGGrec.mrna[genes,substr(colnames(LGGrec.mrna),14,15) %in% "01"]
aux.Y <- aux[,substr(colnames(aux),1,12) %in% rownames(subset(meta, radiation == "Y"))]
aux.N <- aux[,substr(colnames(aux),1,12) %in% rownames(subset(meta, radiation == "N"))]
id.p.Y <- substr(colnames(aux.Y),1,16)
id.p.N <- substr(colnames(aux.N),1,16)
aux <- LGGrec.mrna[genes,substr(colnames(LGGrec.mrna),14,16) %in% "02A"]
aux.Y <- aux[,substr(colnames(aux),1,12) %in% rownames(subset(meta, radiation == "Y"))]
aux.N <- aux[,substr(colnames(aux),1,12) %in% rownames(subset(meta, radiation == "N"))]
id.r.Y <- substr(colnames(aux.Y),1,16)
id.r.N <- substr(colnames(aux.N),1,16)


probes <- c("cg20852945","cg01690062","cg01715901","cg02162404","cg03950246","cg04237608","cg08439109")


probes.df <- LGGrecurrent[probes,1:4]
probes.df <- probes.df[3,]
probes.GR <- GRanges(seqnames = paste0("chr",probes.df$Chromosome),ranges = IRanges(start = probes.df$Genomic_Coordinate, end=probes.df$Genomic_Coordinate), probeID = probes.df$Composite.Element.REF, EntrezID = genes[3])

library(biomaRt)

mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",  path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
ensembl = useDataset("hsapiens_gene_ensembl", mart = mart)
gene.location.df <- getBM(attributes=c("ensembl_exon_id","exon_chrom_start","exon_chrom_end","ensembl_transcript_id","ensembl_gene_id","chromosome_name","strand","start_position","end_position"), filters = "entrezgene", values = genes[3], mart = ensembl)
gene.location.df <- gene.location.df[,c(6,2,3,7,5,1,4)]
colnames(gene.location.df) <- c("chromosome","start","end","strand","gene","exon","transcript")
gene.location.df$chromosome <- paste0("chr",gene.location.df$chromosome)
gene.location.df$strand <- ifelse(gene.location.df$chromosome == "1","+","-")

range.cpg <- GRanges(seqnames = paste0("chr",probes.df$Chromosome),ranges = IRanges(start = min(gene.location.df$start), end=probes.df$Genomic_Coordinate))
all.cpg <- GRanges(seqnames = paste0("chr",LGGrecurrent$Chromosome),ranges = IRanges(start = LGGrecurrent$Genomic_Coordinate, end=LGGrecurrent$Genomic_Coordinate))
overlap.cpg <- as.data.frame(findOverlaps(all.cpg,range.cpg))
overlap.cpg <- LGGrecurrent[overlap.cpg$queryHits,]
overlap.cpg <- na.omit(overlap.cpg)




#grtrack <- GeneRegionTrack(gene.location.df, genome = gen, chromosome = chr, name = "Gene")
#plotTracks(list(itrack, gtrack, atrack, grtrack))


mean.Y.p <- apply(overlap.cpg[,substr(colnames(overlap.cpg),1,16) %in% id.p.Y],1,mean, na.rm=TRUE)
mean.N.p <- apply(overlap.cpg[,substr(colnames(overlap.cpg),1,16) %in% id.p.N],1,mean, na.rm=TRUE)
mean.Y.r <- apply(overlap.cpg[,substr(colnames(overlap.cpg),1,16) %in% id.r.Y],1,mean, na.rm=TRUE)
mean.N.r <- apply(overlap.cpg[,substr(colnames(overlap.cpg),1,16) %in% id.r.N],1,mean, na.rm=TRUE)


#probes.df <- cbind(overlap.cpg[,3:4],overlap.cpg[,substr(colnames(overlap.cpg),1,16) %in% id.p.Y],overlap.cpg[,substr(colnames(overlap.cpg),1,16) %in% id.p.N],overlap.cpg[,substr(colnames(overlap.cpg),1,16) %in% id.r.Y],overlap.cpg[,substr(colnames(overlap.cpg),1,16) %in% id.r.N])

probes.df <- cbind(overlap.cpg[,3:4],overlap.cpg[,substr(colnames(overlap.cpg),1,16) %in% id.p.Y],overlap.cpg[,substr(colnames(overlap.cpg),1,16) %in% id.r.Y])

probes.df$Chromosome <- paste0("chr",probes.df$Chromosome)
#probes.GR <- GRanges(seqnames = paste0("chr",probes.df$Chromosome),ranges = IRanges(start = probes.df$Genomic_Coordinate, end=probes.df$Genomic_Coordinate),mean.Y.p=mean.Y.p,mean.N.p=mean.N.p,mean.Y.r=mean.Y.r,mean.N.r=mean.N.r)
probes.GR <- makeGRangesFromDataFrame(probes.df,keep.extra.columns=T,seqnames.field="Chromosome",start.field="Genomic_Coordinate",end.field="Genomic_Coordinate",ignore.strand=T)
genome(probes.GR) <- "hg19"


atrack <- AnnotationTrack(probes.GR, name = "CpG",stacking="dense",featureAnnotation = "feature", fontcolor.feature = 1,cex.feature = 0.5)
feature(atrack) <- rownames(probes.df)


pdf("teste.pdf")
#plotTracks(list(itrack,atrack,dTrack,knownGenes), groups = c(rep("Rad.p",length(id.p.Y)),rep("NonRad.p",length(id.p.N)),rep("Rad.r",length(id.r.Y)),rep("NonRad.r",length(id.r.N))), type = c("a", "p", "confint"),showSampleNames = TRUE,ylim=c(0,1),legend=T)

plotTracks(atrack )
dev.off()

#plotTracks(atrack)
gtrack <- GenomeAxisTrack()
#plotTracks(list(gtrack, atrack))
gen <- genome(probes.GR)
chr <- as.character(unique(seqnames(probes.GR)))
itrack <- IdeogramTrack(genome = gen, chromosome = chr)
#plotTracks(list(itrack, gtrack, atrack))
#probes.GR.mean <- GRanges(seqnames = paste0("chr",probes.df$Chromosome),ranges = IRanges(start = probes.df$Genomic_Coordinate, end=probes.df$Genomic_Coordinate),mean.Y.p=mean.Y.p,mean.N.p=mean.N.p,mean.Y.r=mean.Y.r,mean.N.r=mean.N.r)

#genome(probes.GR.mean) <- "hg19"


dTrack <- DataTrack(probes.GR, name = "DNA methylation")
#plotTracks(list(itrack, gtrack, atrack, grtrack,dTrack), groups = c("Rad.p","NonRad.p","Rad.r","NonRad,r"), type = c("a", "p", "confint"))

                                      
knownGenes <- UcscTrack(genome = "hg19", chromosome = "chr2",track = "knownGene", from = min(probes.df$Genomic_Coordinate), to = max(probes.df$Genomic_Coordinate), trackType = "GeneRegionTrack",rstarts = "exonStarts", rends = "exonEnds", strand = "strand",collapseTranscripts=TRUE, name = "UCSC Genes",gene="name",symbol="name", transcript="name",id="name",transcriptAnnotation = "gene")

#ucsf <- UcscTrack(genome = "hg19", chromosome = "chr2",track = "ucsfBrainMethyl", from = min(probes.df$Genomic_Coordinate), to = max(probes.df$Genomic_Coordinate),trackType = "AnnotationTrack", start = "chromStart",end = "chromEnd", id = "name", shape = "box",fill = "#006400", name = "UCSF Brain Methyl")

#miRNA <- UcscTrack(genome = "hg19", chromosome = "chr2",track = "wgRna", from = min(probes.df$Genomic_Coordinate), to = max(probes.df$Genomic_Coordinate),trackType = "GeneRegionTrack", start = "chromStart",end = "chromEnd", id = "name", shape = "box",fill = "#006400", name = "miRNA")

ht <- HighlightTrack(trackList = list(atrack,dTrack), start = 158297683, end=158297685,chromosome = 2)
ht.g <- HighlightTrack(trackList = knownGenes, start = 158114110, end=158170723,chromosome = 2)
#ht.g <- HighlightTrack(trackList = knownGenes, start = 158297683, end=158297685,chromosome = 2)

setwd("/dados/ResearchProjects/thais/TCGA/LGG.GBM")

pdf("teste2.pdf")
#plotTracks(list(itrack,atrack,dTrack,knownGenes), groups = c(rep("Rad.p",length(id.p.Y)),rep("NonRad.p",length(id.p.N)),rep("Rad.r",length(id.r.Y)),rep("NonRad.r",length(id.r.N))), type = c("a", "p", "confint"),showSampleNames = TRUE,ylim=c(0,1),legend=T)

plotTracks(list(itrack,gtrack,atrack,dTrack,knownGenes), groups = c(rep("Rad.p",length(id.p.Y)),rep("Rad.r",length(id.r.Y))), type = c("a","confint"),ylim=c(0,1),legend=T,labelPos = "below")
dev.off()

