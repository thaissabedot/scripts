

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


LGGrec.mrna$entrezID <- str_extract(LGGrec.mrna$gene_id,"[^|]*$")
LGGrec.mrna$gene_symbol <- str_extract(LGGrec.mrna$gene_id,"^*[^|]*")
rownames(LGGrec.mrna) <- as.character(LGGrec.mrna$entrezID)
gene.location.RNAseq <- subset(gene.location, EntrezGene.ID %in% rownames(LGGrec.mrna))
genes.RNAseq.GR <- GRanges(seqnames = paste0("chr",gene.location.RNAseq$Chromosome.Name),ranges = IRanges(start = gene.location.RNAseq$Gene.Start..bp., end=gene.location.RNAseq$Gene.End..bp.), strand = 
                             gene.location.RNAseq$Strand, symbol = gene.location.RNAseq$Associated.Gene.Name, EntrezID = gene.location.RNAseq$EntrezGene.ID)
rownames(gene.location.RNAseq) <- as.character(gene.location.RNAseq$EntrezGene.ID)
genes.precede <- matrix(data=NA,nrow=133,ncol=10,dimnames=list(LGG.rad.GR$probeID,paste("Precede.round",1:10,sep=".")))

#### find the top 10 genes upstream a CpG site (450k)
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

#### find the top 10 genes downstream a CpG site (450k)
for(i in 1:10){
  genes.RNAseq.GR <- GRanges(seqnames = paste0("chr",gene.location.RNAseq$Chromosome.Name),ranges = IRanges(start = gene.location.RNAseq$Gene.Start..bp., end=gene.location.RNAseq$Gene.End..bp.), strand = gene.location.RNAseq$Strand, symbol = gene.location.RNAseq$Associated.Gene.Name, EntrezID = gene.location.RNAseq$EntrezGene.ID)
  aux <- follow(LGG.rad.GR,genes.RNAseq.GR)
  genes.follow[,i] <- as.character(gene.location.RNAseq[aux,"EntrezGene.ID"])
  gene.location.RNAseq <- gene.location.RNAseq[-aux,]
  
}

genes.all <- cbind(genes.precede,genes.follow)

library(reshape)
library(scales)
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
  colnames(df)[5] <- "probe"
  ggplot(df, aes(x= probe, y=log2(value), color=factor(radiation))) + 
    geom_point(size=0.9,aes(shape=factor(status)))+ 
    facet_wrap( ~ entrezID,ncol=5)+ #scales="free_y"
    scale_colour_manual(values = c("black","red")) + 
    scale_x_continuous(limits=c(0,1),breaks=c(0,0.25,0.5,0.75,1)) + 
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  ggsave(filename = filename)
}


