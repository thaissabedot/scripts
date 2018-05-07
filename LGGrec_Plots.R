#### find the top -10/+10 genes that are around a CpG site (450k)
gene.location <- read.delim("/dados/ResearchProjects/thais/TCGA/Normals/mRNA_Brain/Galaxy1-[Homo_sapiens_genes_(GRCh37.p13)].tabular")
gene.location <- gene.location[gene.location$Chromosome.Name %in% c(1:22,"X","Y"),]
gene.location <- gene.location[!is.na(gene.location$EntrezGene.ID),]
gene.location <- gene.location[!duplicated(gene.location$EntrezGene.ID),] #the duplicates are different transcripts


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#### how to retrieve data from TxDb objects: http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf
#### All TxDb objects are backed by a SQLite database that manages genomic locations and the relationships between pre-processed mRNA transcripts, exons, protein coding sequences, and their related gene identiers.
txdb.gene <- select(TxDb.Hsapiens.UCSC.hg19.knownGene,keys=as.character(LGGrec.mrna$entrezID),columns=c("TXNAME","TXCHROM","TXSTART","TXEND","TXSTRAND"),keytype="GENEID")
txdb.gene <- na.omit(txdb.gene)
txdb.gene <- subset(txdb.gene,TXCHROM %in% paste0("chr",c(1:22,"X","Y")))
rownames(txdb.gene) <- NULL

a <- LGGrecurrent[rownames(LGG.rad.s),] #LGGrecurrent -> DNA methylation matrix. LGG.rad.s -> 133 CpG sites related to radiation
LGG.rad.GR <- GRanges(seqnames = paste0("chr",a$Chromosome), ranges = IRanges(start = a$Genomic_Coordinate, end = a$Genomic_Coordinate), probeID = a$Composite.Element.REF)
LGGrec.mrna$entrezID <- str_extract(LGGrec.mrna$gene_id,"[^|]*$") #LGGrec.mrna -> Gene expression matrix
LGGrec.mrna$gene_symbol <- str_extract(LGGrec.mrna$gene_id,"^*[^|]*")
rownames(LGGrec.mrna) <- as.character(LGGrec.mrna$entrezID)
gene.location.RNAseq <- subset(gene.location, EntrezGene.ID %in% rownames(LGGrec.mrna))
genes.RNAseq.GR <- GRanges(seqnames = paste0("chr",gene.location.RNAseq$Chromosome.Name),ranges = IRanges(start = gene.location.RNAseq$Gene.Start..bp., end=gene.location.RNAseq$Gene.End..bp.), strand = gene.location.RNAseq$Strand, symbol = gene.location.RNAseq$Associated.Gene.Name, EntrezID = gene.location.RNAseq$EntrezGene.ID)
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
#shape -> status primary/recurrent
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/gbm.lgg.recurrent.meta.rda")
rownames(gbm.lgg.recurrent.meta) <- as.character(gbm.lgg.recurrent.meta$Patient.ID) #metadata
gbm.lgg.recurrent.meta$radiation <- "N" #adding a column to radiation status
gbm.lgg.recurrent.meta[c("TCGA-DU-6397","TCGA-DU-7304","TCGA-FG-5963","TCGA-DU-5870","TCGA-FG-5965","TCGA-DU-6404","TCGA-TQ-A7RV","TCGA-DU-6407","TCGA-DU-5872"),"radiation"] <- "Y" #patients with radiation
meta <- subset(gbm.lgg.recurrent.meta,type %in% "LGG")

aux <- subset(meta, radiation == "Y" & Age.Diagnosis < 70)
setwd("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/RNASeqV2_new/plots")
ES.rec.Rad <- NULL


for(i in 1:nrow(genes.all)){
  exp <- LGGrec.mrna[as.character(genes.all[i,]),substr(colnames(LGGrec.mrna),1,12) %in% rownames(meta)]
  exp <- log2(exp+1)
  exp$entrezID <- as.character(genes.all[i,])
  mean.Rad.primary <- data.frame(mean=apply(exp[,substr(colnames(exp),1,15) %in% paste(as.character(aux$Patient.ID),"-01",sep="")],1,mean,na.rm=TRUE),entrezID=as.character(exp$entrezID))
  if(is.null(ES.rec.Rad)){
    ES.rec.Rad <- matrix((rowSums((exp[,substr(colnames(exp),1,16) %in% paste(as.character(aux$Patient.ID),"-02A",sep="")]) > mean.Rad.primary$mean) > 6 & rowSums(LGGrecurrent[rownames(genes.all)[i],substr(colnames(LGGrecurrent),1,16) %in% paste(as.character(aux$Patient.ID),"-02A",sep="")] < 0.5) > 6 & rowSums(LGGrecurrent[rownames(genes.all)[i],substr(colnames(LGGrecurrent),1,15) %in% paste(as.character(aux$Patient.ID),"-01",sep="")] > 0.5) > 6),dimnames=list(NULL,rownames(genes.all)[i]))
    
  }
  else {
    ES.rec.Rad <- cbind(ES.rec.Rad,matrix((rowSums((exp[,substr(colnames(exp),1,16) %in% paste(as.character(aux$Patient.ID),"-02A",sep="")]) > mean.Rad.primary$mean) > 6 & rowSums(LGGrecurrent[rownames(genes.all)[i],substr(colnames(LGGrecurrent),1,16) %in% paste(as.character(aux$Patient.ID),"-02A",sep="")] < 0.5) > 6 & rowSums(LGGrecurrent[rownames(genes.all)[i],substr(colnames(LGGrecurrent),1,15) %in% paste(as.character(aux$Patient.ID),"-01",sep="")] > 0.5) > 6),dimnames=list(NULL,rownames(genes.all)[i])))
  }

  dna.m <- LGGrecurrent[rownames(genes.all)[i],substr(colnames(LGGrecurrent),1,12) %in% rownames(meta)]
  dna.m <- t(dna.m)
  rownames(dna.m) <- substr(rownames(dna.m),1,16)
  df <- melt(exp, id = "entrezID")
  df$status <- substr(df$variable,14,16) #primary, first recurrence or second recurrence
  df[df$status %in% "01B","status"] <- "01A" #change 01B to 01A
  df$variable <- substr(df$variable,1,16)
  df <- merge(df,dna.m,by.x="variable",by.y=0) #merge DNA methylation with gene expression
  df$variable <- substr(df$variable,1,12)
  df <- merge(df,meta[,c("Patient.ID","radiation")],by.x="variable",by.y="Patient.ID")
  filename <- paste0(rownames(genes.all)[i],".pdf")
  colnames(df)[5] <- "probe"
  df$entrezID <- as.factor(df$entrezID)
  ggplot(df, aes(x= probe, y=value, color=factor(radiation))) + 
    geom_point(size=4,aes(shape=factor(status)))+ 
    stat_smooth(method = "lm", fullrange=FALSE, na.rm=TRUE, se=FALSE, size=0.5,colour="black") +
    geom_vline(xintercept = 0.5, linetype = "longdash", colour="black") +
    facet_wrap( ~ entrezID,ncol=5)+
    geom_hline(aes(yintercept = mean), mean.Rad.primary, linetype = "longdash", colour="black") +
    scale_colour_manual(values = c("forestgreen","darkorchid1"), #check if the order is right
                        labels=c("No Radiation","Radiation"), #check if the order is right
                        name="Biological Group") + 
    scale_shape_discrete(name  ="Sample Type",
                         breaks=c("01A", "02A", "02B"),
                         labels=c("Primary", "First Recurrence", "Second Recurrence")) +
    ylab("RNASeq Expression (log2)") +
    xlab("DNA Methylation Values") + 
    scale_x_continuous(limits=c(0,1),breaks=c(0,0.25,0.5,0.75,1)) + 
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  ggsave(filename = filename, height = 10, width = 18)
}

