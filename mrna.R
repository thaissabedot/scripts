id.mrna <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RNAseq_05-11-14/gbm_file_manifest.txt")
rownames(id.mrna) <- id.mrna$File.Name
gbm.mrna = NULL
pattern.primary = 'TCGA-([0-9A-Z]{2})-([0-9A-Z]{1,})-([0-1A-Z]{3})-([0-9A-Z]{3})-([0-9A-Z]{4}-([0-9]{2}))' #barcode
files = list.files("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RNAseq_05-11-14/gbm_RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/GBM_Genes_normalized") 
setwd("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RNAseq_05-11-14/gbm_RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/GBM_Genes_normalized")
z = 1 #first file
while (!is.na(files[z])) {
  aux <- id.mrna[as.character(files[z]),"Barcode"]
  c <- str_extract(as.character(aux), pattern.primary)
    if(!is.na(c)){ #do not read metadata files and non-tumor samples
        if(is.null(gbm.mrna)){ #read the first sample
          gbm.mrna = read.table(files[z],sep="\t",header=TRUE) 
          colnames(gbm.mrna)[2] <- as.character(c)
          gbm.mrna <- gbm.mrna[,c(1,2)]
       } 
        else{
          aux = read.table(files[z], header=TRUE,sep="\t")
          colnames(aux)[2] = c
          gbm.mrna = merge(gbm.mrna, aux[1:2], by="gene_id")
        }
      }
    z = z + 1 #next file
}#end while
# gbm -> 154 tumor samples

#duplicated(substr(colnames(gbm.mrna),1,12)) #TCGA-06-0211 AND TCGA-06-0156
gbm.mrna <- gbm.mrna[,-c(77,92)] #remove duas duplicatas
colnames(gbm.mrna)[52] <- "TCGA-06-0211-01A-01R-1849-01"
colnames(gbm.mrna)[71] <- "TCGA-06-0156-01A-03R-1849-01"

lgg.mrna = NULL
id.mrna <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RNAseq_05-11-14/lgg_file_manifest.txt")
rownames(id.mrna) <- id.mrna$File.Name
files = list.files("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RNAseq_05-11-14/lgg_RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/LGG_Genes_normalized") 
setwd("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RNAseq_05-11-14/lgg_RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/LGG_Genes_normalized")
z = 1 #first file
while (!is.na(files[z])) {
  aux <- id.mrna[as.character(files[z]),"Barcode"]
  c <- str_extract(as.character(aux), pattern.primary)
  if(!is.na(c)){ #do not read metadata files and non-tumor samples
    if(is.null(lgg.mrna)){ #read the first sample
      lgg.mrna = read.table(files[z],sep="\t",header=TRUE) 
      colnames(lgg.mrna)[2] <- as.character(c)
      lgg.mrna <- lgg.mrna[,c(1,2)]
    }
    else{
      aux = read.table(files[z], header=TRUE,sep="\t")
      colnames(aux)[2] = c
      lgg.mrna = merge(lgg.mrna, aux[1:2], by="gene_id")
    }
  }
  z = z + 1 #next file
}#end while
# lgg -> 514 tumor samples

#save(lgg.mrna,gbm.mrna,file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/RNAseq_05-11-14/mRNA.Rda")
mrna <- merge(lgg.mrna,gbm.mrna,by="gene_id")
rownames(mrna) <- as.character(mrna$gene_id)
colnames(mrna) <- substr(colnames(mrna),1,12)
mrna <- mrna[,-which(names(mrna) %in% c("TCGA-HT-A61A"))] #remove this sample
entrezID <- str_extract(mrna$gene_id,"[^|]*$")
rownames(mrna) <- entrezID
#entrezID <- str_extract(mrna$gene_id,"*[^|]*")


### Control TCGA DNA methylation
ids <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Control_EScalls/TCGAnormalsID.unique2.txt",sep=",")
ids$TCGAnew <- substr(ids$TCGAID,1,12)
aux <- subset(ids, Disease %in% c("BLCA","BRCA","COAD","HNSC","KIRC","KIRP","LIHC","LUAD","PRAD","THCA","UCEC"))
aux$Disease <-  factor(aux$Disease)
write.table(unique(as.character(subset(aux, Disease %in% levels(aux$Disease)[9])$TCGAnew))[1:10],file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Control_EScalls/idParaDownload.txt",quote=F,row.names=F,sep="\t") #de 1 a 11
normal.dnameth = NULL
folders = list.files(path="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Control_EScalls/",pattern = "DNA_Methylation")
z = 1 #first folder
pattern = 'TCGA-([0-9A-Z]{2})-([0-9A-Z]{1,})-([0-9A-Z]{3})-([0-9A-Z]{3})-([0-9A-Z]{4}-([0-9]{2}))' #barcode
beginning <- "/dados/ResearchProjects/thais/TCGA/LGG.GBM/Control_EScalls/"
end <- "/JHU_USC__HumanMethylation450/Level_3"
while (!is.na(folders[z])) {
  path <- (paste(c(paste(c(paste(c(beginning,folders[z]), collapse=''),""),collapse=''),end),collapse=''))
  files = list.files(path)
  i <- 1
  while (!is.na(files[i])) {
    patientID <- str_extract(files[i], pattern)
    if(is.null(normal.dnameth)){ #read the first sample
      normal.dnameth = read.table(paste(path,files[i],sep="/"),header=TRUE,quote="",comment="",sep="\t", skip=1) 
      normal.dnameth = subset(normal.dnameth,select=c(Composite.Element.REF,Gene_Symbol,Chromosome,Genomic_Coordinate,Beta_value))
      colnames(normal.dnameth)[5] = patientID
    }
    else{
      aux = read.table(paste(path,files[i],sep="/"),header=TRUE,quote="",comment="",sep="\t", skip=1) 
      colnames(aux)[2] = patientID
      normal.dnameth = merge(normal.dnameth, aux[,1:2], by="Composite.Element.REF")
    }
    i = i + 1 #next file
  }#end while files
  z = z + 1 #next folder
}#end while folders
#27k sample (COAD)
#file <- "/dados/ResearchProjects/thais/TCGA/LGG.GBM/Control_EScalls/COAD_DNA_Methylation/JHU_USC__HumanMethylation27/Level_3/jhu-usc.edu_COAD.HumanMethylation27.2.lvl-3.TCGA-AA-3516-11A-01D-1551-05.txt"
file <- "/dados/ResearchProjects/thais/TCGA/LGG.GBM/Control_EScalls/COAD_DNA_Methylation/jhu-usc.edu_COAD.HumanMethylation450.1.lvl-3.TCGA-A6-2681-11A-01D-1551-05.txt"
patientID <- str_extract(file, pattern)
aux <- read.table(file,header=TRUE,quote="",comment="",sep="\t", skip=1)
colnames(aux)[2] = patientID
normal.dnameth = merge(normal.dnameth, aux[,1:2], by="Composite.Element.REF")
normal.dnameth = normal.dnameth[,c(1:24,114,25:113)] #deixar na ordem


#RNAseq Control

normal.mrna = NULL
folders = list.files(path="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Control_EScalls/",pattern = "RNASeqV2")
beginning <- "/dados/ResearchProjects/thais/TCGA/LGG.GBM/Control_EScalls/"
end <- "/UNC__IlluminaHiSeq_RNASeqV2/Level_3"
metadata <- "/UNC__IlluminaHiSeq_RNASeqV2/file_manifest.txt"
z = 1 #first folder
while (!is.na(folders[z])) {
  path <- (paste(c(paste(c(paste(c(beginning,folders[z]), collapse=''),""),collapse=''),end),collapse=''))
  files = list.files(path)
  id.mrna <- read.delim((paste(c(paste(c(paste(c(beginning,folders[z]), collapse=''),""),collapse=''),metadata),collapse='')))
  rownames(id.mrna) <- id.mrna$File.Name
  i <- 1
  while (!is.na(files[i])) {
    patientID <- id.mrna[as.character(files[i]),"Barcode"]
    if(is.null(normal.mrna)){ #read the first sample
      normal.mrna = read.table(paste(path,files[i],sep="/"),header=TRUE,quote="",comment="",sep="\t") 
      colnames(normal.mrna)[2] = as.character(patientID)
    }
    else{
      aux = read.table(paste(path,files[i],sep="/"),header=TRUE,quote="",comment="",sep="\t")
      colnames(aux)[2] = as.character(patientID)
      normal.mrna = merge(normal.mrna, aux[,1:2], by="gene_id")
    }
    i = i + 1 #next file
  }#end while files
  z = z + 1 #next folder
}#end while folders

normal.mrna <- normal.mrna[,-c(3,4)]
#save(normal.dnameth,normal.mrna,file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Control_EScalls/mRNA.Rda")
#rownames(normal.mrna) <- as.character(normal.mrna$gene_id)
colnames(normal.mrna) <- substr(colnames(normal.mrna),1,12)
colnames(normal.dnameth) <- substr(colnames(normal.dnameth),1,12)
#identical(normal.mrna$gene_id,mrna$gene_id)
normal.mrna$entrezID <- as.matrix(entrezID)
rownames(normal.mrna) <- entrezID
### FIM Control TCGA

mrna <- merge(lgg.mrna,gbm.mrna,by="gene_id")
rownames(mrna) <- as.character(mrna$gene_id)
colnames(mrna) <- substr(colnames(mrna),1,12)
mrna <- mrna[,-which(names(mrna) %in% c("TCGA-HT-A61A"))] #remove this sample
entrezID <- str_extract(mrna$gene_id,"[^|]*$")
#entrezID[entrezID == "?"] <- str_extract(mrna$gene_id[1:29],"[^|]*$")
#entrezID[entrezID == "SLC35E2"] <- c("SLC35E2.728661","SLC35E2.9906")
mrna$entrezID <- as.matrix(entrezID)
symbol <- str_extract(mrna$gene_id,"^*[^|]*")
symbol[symbol == "?"] <- str_extract(mrna$gene_id[1:29],"[^|]*$")
symbol[symbol == "SLC35E2"] <- c("SLC35E2.728661","SLC35E2.9906")
rownames(mrna) <- symbol
mrna$symbol <- as.matrix(symbol)

mrna.all <- merge(mrna[,-1], normal.mrna[,-1],by="entrezID")
rownames(mrna.all) <- as.character(mrna.all$symbol)
mrna.all <- mrna.all[,-c(669)] #777 samples: 667 tumor + 110 normal

#######Normalize RNAseq data - Stefano M. Pagnotta 
#Tumor and Normal samples
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RNAseq_05-11-14/geneInfo.RData") 
library("EDASeq")
geneInfo <- na.omit(geneInfo)
wwhich <- which(rownames(geneInfo) %in% rownames(mrna.all))
geneInfo <- geneInfo[wwhich,]
mrna <- mrna.all[rownames(geneInfo),]
id <- mrna$entrezID
mrna <- mrna[,-1] #id column
#geneInfo <- geneInfo[rownames(rna.seq.lgg.gbm),]
ffData  <- as.data.frame(geneInfo)
#rna.seq.lgg.gbm.norm <- floor(as.matrix(rna.seq.lgg.gbm[,1:636]))
mrna <- floor(as.matrix(mrna))
tmp <- newSeqExpressionSet(mrna, featureData = ffData)
fData(tmp)[, "gcContent"] <- as.numeric(geneInfo[, "gcContent"])
tmp <- withinLaneNormalization(tmp, "gcContent", which = "upper", offset = TRUE)
tmp <- betweenLaneNormalization(tmp, which = "upper", offset = TRUE)
########################
normCounts <-  log(mrna + .1) + offst(tmp)
normCounts <-  floor(exp(normCounts) - .1)
library(steFunctions) # R package attached that contains the quantileNormalization function and other stuff
tmp <- t(quantileNormalization(t(normCounts)))
normCounts <- floor(tmp)
rna.seq.norm <- as.data.frame(normCounts)
#rna.seq.lgg.gbm.norm$entrezID <- as.character(rna.seq.lgg.gbm$entrezID)

#rna.seq.lgg.gbm <- rna.seq.norm[,colnames(rna.seq.norm) %in% colnames(LGG.GBM)]
rna.seq.norm$entrezID <- as.matrix(id)
#identical(rownames(mrna),rownames(rna.seq.norm)) OK

#save(normal.mrna,rna.seq.lgg.gbm,mrna,rna.seq.norm,file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/RNAseq_05-11-14/rnaSeqNorm.Rda")

### map gene ID to genomic coordinates
gene.location <- read.delim("/dados/ResearchProjects/thais/TCGA/Normals/mRNA_Brain/Galaxy1-[Homo_sapiens_genes_(GRCh37.p13)].tabular")
gene.location <- gene.location[gene.location$Chromosome.Name %in% c(1:22,"X","Y"),]
gene.location <- gene.location[!is.na(gene.location$EntrezGene.ID),]
gene.location <- gene.location[!duplicated(gene.location$EntrezGene.ID),] #the duplicates are different transcripts, not different coordinates
library(GenomicRanges)
gene.GR <- GRanges(seqnames = paste0("chr",gene.location$Chromosome.Name),ranges = IRanges(start = gene.location$Gene.Start..bp., end=gene.location$Gene.End..bp.), strand = gene.location$Strand, symbol = gene.location$Associated.Gene.Name, EntrezID = gene.location$EntrezGene.ID)
probe.info <- GRanges(seqnames = paste0("chr",LGG.GBM.new$Chromosome), ranges = IRanges(start = LGG.GBM.new$Genomic_Coordinate, end = LGG.GBM.new$Genomic_Coordinate), probeID = LGG.GBM.new$Composite.Element.REF)
distance <- as.data.frame(distanceToNearest(probe.info,gene.GR)) #closest gene to each 450k probe
rownames(gene.location) <- NULL
gene.order.by.distance <- gene.location[distance$subjectHits,]
gene.order.by.distance$distance <- as.matrix(distance$distance)
GBM.LGG.27.450k.nearest.gene <- cbind(LGG.GBM.new[,1:4],LGG.GBM)
GBM.LGG.27.450k.nearest.gene <- cbind(GBM.LGG.27.450k.nearest.gene, gene.order.by.distance[,c("Associated.Gene.Name","EntrezGene.ID","distance")])
info <- GBM.LGG.27.450k.nearest.gene[,c(1:4,937:939)]
colnames(GBM.LGG.27.450k.nearest.gene) <- substr(colnames(GBM.LGG.27.450k.nearest.gene),1,12)
GBM.LGG.27.450k.nearest.gene <- GBM.LGG.27.450k.nearest.gene[,colnames(GBM.LGG.27.450k.nearest.gene) %in% substr(colnames(rna.seq.norm)[1:667],1,12)]
GBM.LGG.27.450k.nearest.gene <- cbind(info,GBM.LGG.27.450k.nearest.gene)

#cpg.gene <- merge(GBM.LGG.27.450k.nearest.gene,rna.seq.lgg.gbm,by.x="EntrezGene.ID",by.y="entrezID")
cpg.gene <- merge(GBM.LGG.27.450k.nearest.gene,rna.seq.norm[,colnames(GBM.LGG.27.450k.nearest.gene)[8:643]],by.x="Associated.Gene.Name",by.y=0) 
#25978 probes - 1387 samples
dimnames(cpg.gene)[[1]] <- paste(cpg.gene[,"Composite.El"], cpg.gene[,"Associated.Gene.Name"], sep=".")
#(d[,8:584],d[,586:1162]) #mDNA methylation x Gene expression

meth.g <- cpg.gene[,8:643] #only tumor
names(meth.g) <- substr(names(meth.g), 1, 12)
expr.g <- cpg.gene[,644:1279] #only tumor
names(expr.g) <- substr(names(expr.g), 1, 12)
#identical(colnames(meth.g),colnames(expr.g))



normalTCGA.mrna <- rna.seq.norm[,colnames(rna.seq.norm) %in% colnames(normal.mrna)]
normalTCGA.mrna$entrezID <- as.matrix(id)
#map dna methylation probes to gene (normal)
gene.GR <- GRanges(seqnames = paste0("chr",gene.location$Chromosome.Name),ranges = IRanges(start = gene.location$Gene.Start..bp., end=gene.location$Gene.End..bp.), strand = gene.location$Strand, symbol = gene.location$Associated.Gene.Name, EntrezID = gene.location$EntrezGene.ID)
probe.info <- GRanges(seqnames = paste0("chr",normal.dnameth$Chromosome), ranges = IRanges(start = normal.dnameth$Genomic_Coor, end = normal.dnameth$Genomic_Coor), probeID = normal.dnameth$Composite.El)
distance <- as.data.frame(distanceToNearest(probe.info,gene.GR)) #closest gene to each 450k probe
rownames(gene.location) <- NULL
gene.order.by.distance <- gene.location[distance$subjectHits,]
gene.order.by.distance$distance <- as.matrix(distance$distance)
normalTCGA.dnaMeth.nearest.gene <- cbind(normal.dnameth, gene.order.by.distance[,c("Associated.Gene.Name","EntrezGene.ID","distance")])
info <- normalTCGA.dnaMeth.nearest.gene[,c(1:4,115:117)]
colnames(normalTCGA.dnaMeth.nearest.gene) <- substr(colnames(normalTCGA.dnaMeth.nearest.gene),1,12)
normalTCGA.dnaMeth.nearest.gene <- normalTCGA.dnaMeth.nearest.gene[,colnames(normalTCGA.mrna)[1:110]]
normalTCGA.dnaMeth.nearest.gene <- cbind(info,normalTCGA.dnaMeth.nearest.gene)
normalTCGA.dnaMeth.nearest.gene <- subset(normalTCGA.dnaMeth.nearest.gene, Composite.El %in% as.character(GBM.LGG.27.450k.nearest.gene$Composite.El))


cpg.gene.normal <- merge(normalTCGA.dnaMeth.nearest.gene,normalTCGA.mrna,by.x="EntrezGene.ID",by.y="entrezID")
dimnames(cpg.gene.normal)[[1]] <- paste(cpg.gene.normal[,"Composite.El"], cpg.gene.normal[,"Associated.Gene.Name"], sep=".") 

#identical(rownames(cpg.gene),rownames(cpg.gene.normal)) OK

meth.g.norm <- cpg.gene.normal[,8:117]
names(meth.g.norm) <- substr(names(meth.g.norm), 1, 12)
expr.g.norm <- cpg.gene.normal[,118:227]
names(expr.g.norm) <- substr(names(expr.g.norm), 1, 12)
epiGene.RNAseq <- list()
j=1
#metadata <- subset(metadata.1129.samples.20141030,Study == "TCGA")
a <- subset(pd,case.id %in% colnames(meth.g))
rownames(a) <- as.character(a$case.id)

for(i in 1:dim(cpg.gene.normal)[1]) {
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
      aux <- cbind(t(meth.g.norm[rownames(meth.g[i,]), ] >= 0.3), t(expr.g.norm[rownames(meth.g[i,]), ] < expr.mean.gu))
      colnames(aux) <- c("meth.gt0.3", "exp.lt.meanUnmeth")
      yy <- rbind(yy, aux)
      yy <- as.data.frame(yy)
      gene <- rownames(meth.g.norm[i,])
      if(as.matrix(table(yy$meth.gt0.3[1:636], yy$exp.lt.meanUnmeth[1:636]))[2,2]/as.numeric(table(yy$meth.gt0.3[1:636])[2])>0.8){
        yy$geneName <- gene
        yy$id <- rownames(yy)
        yy$methValues <- NA
        yy$exprValues <- NA
        yy$methValues[1:636] <- as.numeric(meth.g[i, ])
        yy$exprValues[1:636] <- as.numeric(expr.g[i, ])  
        #ic <-248+(491*(k-1))
        #fim <- 247+(491*k)
        yy$methValues[637:nrow(yy)] <- as.vector(as.matrix(t(meth.g.norm[rownames(meth.g[i,]), ])))
        yy$exprValues[637:nrow(yy)] <- as.vector(as.matrix(t(expr.g.norm[rownames(meth.g[i,]), ])))
        #yy$methValues[248:nrow(yy)] <- as.numeric(meth.g.norm[index[k], ])
        #yy$exprValues[248:nrow(yy)] <- as.numeric(expr.g.norm[index[k], ])
        info <- a[rownames(yy)[1:636],c("cluster.meth","IDH.status","Codel.1p.19q","cluster.meth.gbm","cluster.meth.lgg","tumor.type")]
        yy$IDH.status <- "WT"
        yy$codel1p19q <- "non-codel"
        yy$cluster <- NA
        yy$cluster.gbm <- NA
        yy$cluster.lgg <- NA
        yy$tumor.type <- NA
        yy$IDH.status[1:636] <- as.character(info$IDH.status)
        yy$cluster.gbm[1:636] <- as.character(info$cluster.meth.gbm)
        yy$cluster.lgg[1:636] <- as.character(info$cluster.meth.lgg)
        yy$tumor.type[1:636] <- as.character(info$tumor.type)
        yy$codel1p19q[1:636] <- as.character(info$Codel.1p.19q)
        yy$cluster[1:636] <- as.character(info$cluster.meth)
        #yy <- merge(metadata, yy, by.x = "TCGAIDnew", by.y = "id")
        epiGene.RNAseq[[j]] <- as.data.frame(yy)
        names(epiGene.RNAseq)[[j]] <- as.character(gene)
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

#save(epiGene.RNAseq,file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/RNAseq_05-11-14/epiGene_RNAseq.Rda")

uu.nl <- epiGene.RNAseq[[1]]
clusters <- c(paste0("LGm", 1:6))
#clusters <- "LGm1"


enrichment.normal <- lapply(clusters, function(cluster, epiGene.ss=epiGene.RNAseq) {
  enrichment.group <- unlist(mclapply(epiGene.ss, function(uu.nl, cluster.group=cluster) {
    control <- (!uu.nl$meth.gt0.3[637:length(uu.nl$meth.gt0.3)] & !uu.nl$exp.lt.meanUnmeth[637:length(uu.nl$exp.lt.meanUnmeth)])
    uu.nl$es <- NA
    uu.nl$es[1:636] <- uu.nl$meth.gt0.3[1:636] & uu.nl$exp.lt.meanUnmeth[1:636]
    uu.nl$Kclass <- NA
    uu.nl$Kclass <- uu.nl$cluster == cluster.group
    es.table <- table(uu.nl$Kclass, uu.nl$es)#, dnn=c("Class", "ES"))
    #print()
    if((es.table[2,2]/rowSums(es.table)[2] > 0.5) & (es.table[2,2]/colSums(es.table)[2] > 0.5)  & (table(control)[2]/(table(control)[2] + table(control)[1]) > 0.5 | is.na(table(control)[2]/(table(control)[2] + table(control)[1])))) {
      return(fisher.test(table(uu.nl$Kclass, uu.nl$es))$p.value)
    } else {
      return(NULL)
    }
    #}
    #}  
  }, mc.cores=1))
  return(enrichment.group)
})

names(enrichment.normal) <- c("LGm1", "LGm2", "LGm3", "LGm4", "LGm5", "LGm6")
enrichment.normal.adj <- lapply(enrichment.normal, p.adjust, method="BH")
enrichmentList.sig <- lapply(enrichment.normal.adj, function(x) {
  x <- x[x < 0.05]
})

a <- c(enrichmentList.sig$LGm2)
names(a) <- str_extract(names(a),"[^.]*$")
a <- as.matrix(a)
write.table(unique(rownames(a)),file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/RNAseq_05-11-14/gene.LGm2.txt",quote=F,row.names=F,sep="\t")

gene.names.LGm2 <- unlist(strsplit(names(enrichmentList.sig$LGm2),"[.]"))
gene.names.LGm3 <- unlist(strsplit(names(enrichmentList.sig$LGm3),"[.]"))
even_indexes<-seq(2,length(gene.names.LGm2),2)
even_indexes2<-seq(2,length(gene.names.LGm3),2)
gene.names.LGm2 <- gene.names.LGm2[even_indexes]
gene.names.LGm3 <- gene.names.LGm3[even_indexes2]
genes.names.LGm2 <- c(gene.names.LGm2,gene.names.LGm3)
genes.names.LGm2 <- unique(genes.names.LGm2)
#gene.names.LGm2[988] <- "AC091878.1"
#gene.names.LGm2.a <- unlist(strsplit(names(enrichmentList.sig$LGm6)[1220:1261],"[.]"))
#even_indexes<-seq(2,length(gene.names.LGm2.a),2)
#gene.names.LGm2.a <- gene.names.LGm2.a[even_indexes]
#gene.names.LGm2.a <- unique(gene.names.LGm2.a)
write.table(gene.names.LGm2,file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/RNAseq_05-11-14/gene.names.LGm5.txt",quote=F,row.names=F,sep="\t")

##Heatmap ES calls

a <- c(names(enrichmentList.sig$LGm2),names(enrichmentList.sig$LGm3),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm5))
a <- c(names(enrichmentList.sig$LGm3),names(enrichmentList.sig$LGm4)) #Soh ES call
#a <- c(names(enrichmentList.sig$LGm3),names(enrichmentList.sig$LGm4),names(enrichmentList.rev$LGm2),names(enrichmentList.rev$LGm5)) #es calls and awake
b <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/probesLGm2.txt",header=F)
a <- c(names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm3),names(sort(enrichmentList.sig$LGm2)[1:15]),as.character(b$V1))
probe.meth <- str_extract(a,"*[^.]*")
probe.mrna <- str_extract(a,"[^.]*$")


a <- names(enrichmentList.sig$LGm2)
probe.meth <- str_extract(a,"*[^.]*")
aux <- subset(metadata.1129.samples.20150112, cluster.meth %in% c("LGm1","LGm2","LGm3"))
lgm <- LGG.GBM.250914[as.character(probe.meth),as.character(aux$id)]
lgm$mean <- apply(lgm,1,mean,na.rm=TRUE)
aux <- subset(metadata.1129.samples.20150112, cluster.meth %in% c("LGm4","LGm5","LGm6"))
not.lgm <- LGG.GBM.250914[as.character(probe.meth),as.character(aux$id)]
not.lgm$mean <- apply(not.lgm,1,mean,na.rm=TRUE)
ESg3 <- data.frame(lgm1.to.3.mean=lgm$mean,lgm4.to.6.mean=not.lgm$mean)
rownames(ESg3) <- rownames(lgm)
save(ESg3,file="ESg3_mean.Rda")


###LGG.GBM.new.noXY.s <- LGG.GBM.250914[ as.character(probe.meth),5:936]
###LGG.GBM.new.order <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
###colnames(LGG.GBM.new.order) <- substr(colnames(LGG.GBM.new.order),1,12)
###LGG.GBM.new.order <- LGG.GBM.new.order[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
###a <- colnames(LGG.GBM.new.order) %in% colnames(rna.seq.lgg.gbm)
###teste <- LGG.GBM.new.order[,a]
#LGG.GBM.new.order <- rna.seq.lgg.gbm[as.character(unique(probe.mrna)),colnames(teste)]
######LGG.GBM.new.order <- rna.seq.lgg.gbm[as.character(probe.mrna),colnames(teste)]
#LGG.GBM.new.order <- na.omit(LGG.GBM.new.order)
#LGG.GBM.new.order <- subset(LGG.GBM.new.order,select=colnames(GBM.LGG.27.450k.nearest.gene)[8:643])

c <- cpg.gene
d <- cpg.gene.normal
rownames(c) <- as.character(c$Composite.El)
rownames(d) <- as.character(d$Composite.El)
LGG.GBM.new.order <- c[as.character(probe.meth),644:1279]
aux <- as.character(c[as.character(probe.meth),"Associated.Gene.Name"])
#rownames(LGG.GBM.new.order) <- as.character(aux)
normals.sel <- d[as.character(probe.meth),118:227]
#rownames(normals.sel) <- as.character(aux)
a <- metadata.1129.samples.20150312[substr(colnames(LGG.GBM.new.order),1,12),]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$flag.categories),as.matrix(a$cluster.cnv),as.matrix(a$mut.TP53),as.matrix(a$cnv.EGFR),as.matrix(a$mut.BRAF),as.matrix(a$histology),as.matrix(a$grade),as.matrix(a$subtype))
clab.6 <- clab.6[,c(1,20,19,3,17,16,4,21,14,15,10,11)]
rownames(normalTCGA.dnaMeth.nearest.gene) <- as.character(normalTCGA.dnaMeth.nearest.gene$Composite.El)
####normals.sel <- normalTCGA.mrna[rownames(LGG.GBM.new.order),1:110]
norm.label <- matrix("white", nrow = 110, ncol = 12)
norm.label[,10] <- "darkblue"
clab.6 <- rbind(norm.label,clab.6)
LGG.GBM.new.order <- cbind(normals.sel,LGG.GBM.new.order)
LGG.GBM.new.order <- LGG.GBM.new.order + 1
LGG.GBM.new.order <- log2(LGG.GBM.new.order)
LGG.GBM.new.order <- LGG.GBM.new.order + 1
colnames(LGG.GBM.new.order) <- substr(colnames(LGG.GBM.new.order),1,12)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_RNAseqBoth2.pdf")
heatmap.plus.sm(as.matrix(LGG.GBM.new.order),
                col=rev(redgreen(1000)),
                scale = "none",
                #trace = "none",
                labRow = as.character(aux),
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab.6
                #RowSideColors = rlab
                
)
dev.off()

c <- cpg.gene
rownames(c) <- as.character(c$Composite.El)
LGG.GBM.new.order <- c[rownames(order),644:1279]
aux <- as.character(c[rownames(order),"Associated.Gene.Name"])
#b <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/probesLGm2.txt",header=F)
#a <- c(rownames(aux),names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm3),names(sort(enrichmentList.sig$LGm2)[1:15]),as.character(b$V1))
#probe.meth <- str_extract(a,"*[^.]*")
colnames(LGG.GBM.new.order) <- substr(colnames(LGG.GBM.new.order),1,12)
a <- metadata.1129.samples.20150312[metadata.1129.samples.20150312$id %in% colnames(LGG.GBM.new.order),]
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGm_subgroups.Rda")
a <- merge(a,new.pdata,by="id")
lgm1.hypo <- subset(a, new.subgroup == "LGm1_hypo")
lgm1.hypo.s <- LGG.GBM.new.order[,as.character(lgm1.hypo$id)]
mean <- apply(lgm1.hypo.s, 1, mean, na.rm=TRUE)
lgm1 <- subset(a, new.subgroup == "LGm1")
lgm1.s <- LGG.GBM.new.order[,as.character(lgm1$id)]
mean <- cbind(mean,apply(lgm1.s, 1, mean, na.rm=TRUE))
lgm2 <- subset(a, new.subgroup == "LGm2")
lgm2.s <- LGG.GBM.new.order[,as.character(lgm2$id)]
mean <- cbind(mean,apply(lgm2.s, 1, mean, na.rm=TRUE))
lgm3 <- subset(a, new.subgroup == "LGm3")
lgm3.s <- LGG.GBM.new.order[,as.character(lgm3$id)]
mean <- cbind(mean,apply(lgm3.s, 1, mean, na.rm=TRUE))
lgm4 <- subset(a, new.subgroup == "LGm4")
lgm4.s <- LGG.GBM.new.order[,as.character(lgm4$id)]
mean <- cbind(mean,apply(lgm4.s, 1, mean, na.rm=TRUE))
lgm5 <- subset(a, new.subgroup == "LGm5")
lgm5.s <- LGG.GBM.new.order[,as.character(lgm5$id)]
mean <- cbind(mean,apply(lgm5.s, 1, mean, na.rm=TRUE))
lgm6 <- subset(a, new.subgroup == "LGm6")
lgm6.s <- LGG.GBM.new.order[,as.character(lgm6$id)]
mean <- cbind(mean,apply(lgm6.s, 1, mean, na.rm=TRUE))
mean <- as.data.frame(mean)
colnames(mean) <- c("LGm1-hypo","LGm1","LGm2","LGm3","LGm4","LGm5","LGm6")
a <- metadata.1129.samples.20150312[1:6,]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$flag.categories),as.matrix(a$cluster.cnv),as.matrix(a$mut.TP53),as.matrix(a$cnv.EGFR),as.matrix(a$mut.BRAF),as.matrix(a$histology),as.matrix(a$grade),as.matrix(a$subtype))
clab.6 <- clab.6[,c(1,20,19,18,3,17,16,4,21,14,15,10,11)]
clab.6[,13] <-c("green","red","purple","orange","yellow","blue")
mean <- mean + 1
teste <- c("green3","green","red","purple","orange","yellow","blue")
#png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_MeanRNA3.png",res=1200,width=8000,height=8000)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_MeanRNA3.pdf")
heatmap.plus.sm(as.matrix(mean),
                col=rev(redgreen(1000)),
                #scale = "none",
                #trace = "none",
                labRow = as.character(aux),
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                #ColSideColors = clab.6
                ColSideColors = cbind(teste,teste)
                #RowSideColors = rlab
                
)
dev.off()

tert.expr.g <- read.delim("TERTexpr_grps.txt")
rownames(tert.expr.g) <- as.character(tert.expr.g$id)

b <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/probesLGm2.txt",header=F)
aux <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/EReg5.txt")
a <- c(as.character(aux$id),names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm3),names(sort(enrichmentList.sig$LGm2)[1:15]),as.character(b$V1))
probe.meth <- str_extract(a,"*[^.]*")

LGG.GBM.new.noXY.s <- LGG.GBM.250914[ as.character(probe.meth),5:936]
LGG.GBM.new.order2 <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
colnames(LGG.GBM.new.order2) <- substr(colnames(LGG.GBM.new.order2),1,12)
LGG.GBM.new.order2 <- LGG.GBM.new.order2[,c(1:63,64:327,405:531,532:682,683:932,328:404)]

ESg1 <- LGG.GBM.new.order2[as.character(b$V1),]
ESg1.p <- hclust(dist(ESg1))
probe.meth <- str_extract(names(sort(enrichmentList.sig$LGm2)[1:15]),"*[^.]*")
ESg2 <- LGG.GBM.new.order2[as.character(probe.meth),]
ESg2.p <- hclust(dist(ESg2))
probe.meth <- str_extract(names(enrichmentList.sig$LGm3),"*[^.]*")
ESg3 <- LGG.GBM.new.order2[as.character(probe.meth),]
ESg3.p <- hclust(dist(ESg3))
probe.meth <- str_extract(c(names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4)),"*[^.]*")
ESg4 <- LGG.GBM.new.order2[as.character(probe.meth),]
ESg4.p <- hclust(dist(ESg4))
probe.meth <- str_extract(as.character(aux$id),"*[^.]*")
ESg5 <- LGG.GBM.new.order2[as.character(probe.meth),]
ESg5.p <- hclust(dist(ESg5))


#order <- rbind(ESg1[ESg1.p$order,],ESg2[ESg2.p$order,],ESg3[ESg3.p$order,],ESg4[ESg4.p$order,],ESg5[ESg5.p$order,])
#a <- metadata.1129.samples.20150312[colnames(order),]
order <- rbind(ESg5[ESg5.p$order,],ESg4[ESg4.p$order,],ESg3[ESg3.p$order,],ESg2[ESg2.p$order,],ESg1[ESg1.p$order,])
a <- metadata.1129.samples.20150312[colnames(order),]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$flag.categories),as.matrix(a$cluster.cnv),as.matrix(a$mut.TP53),as.matrix(a$cnv.EGFR),as.matrix(a$mut.BRAF),as.matrix(a$histology),as.matrix(a$grade),as.matrix(a$subtype))
#clab.6 <- metadata.1129.samples.20150312[colnames(order),c("cluster.meth","histology")]
d <- normal.dnameth
rownames(d) <- as.character(d$Composite.El)
d <- d[rownames(order),5:114]
norm.label <- matrix("white", nrow = 110, ncol = 21)
norm.label[,14] <- "darkblue"
order <- cbind(d,order)
clab.6 <- rbind(norm.label,clab.6)
rownames(clab.6) <- colnames(order)
clab <- merge(clab.6, tert.expr.g[,c(1,3)],all.x=T,by.x=0,by.y="id")
levels(clab$grp) <- c("black","purple")
clab$grp <- as.character(clab$grp)
clab[is.na(clab$grp),"grp"] <- "white"
rownames(clab) <- as.character(clab$Row.names)
clab <- clab[colnames(order),]
levels(clab$cluster.meth) <- c("green","red","purple","orange","yellow","blue")
clab$cluster.meth <- as.character(clab$cluster.meth)
clab[is.na(clab$cluster.meth),"cluster.meth"] <- "white"
clab$histology <- as.factor(clab$histology)
levels(clab$histology) <- c("red","purple","cyan","green")
clab$histology <- as.character(clab$histology)
clab[is.na(clab$histology),"histology"] <- "white"
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_DNAMeth_new2.png",res=1200,width=8000,height=8000)
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_DNAMeth_new.pdf")
heatmap.plus.sm(as.matrix(order[1:2,]),
                col = jet.colors(75),
                scale = "none",
                #labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = as.matrix(clab[,c(3,4,2)]),
                margins = c(1, 6)
                #RowSideColors = rlab
                
)
dev.off()

LGG.GBM.new.noXY.s <- LGG.GBM.250914[ as.character(probe.meth),5:936]
LGG.GBM.new.order2 <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
colnames(LGG.GBM.new.order2) <- substr(colnames(LGG.GBM.new.order2),1,12)
#a <- colnames(LGG.GBM.new.order) %in% colnames(GBM.LGG.27.450k.nearest.gene)
#teste <- LGG.GBM.new.order[,a]
#LGG.GBM.new.order <- subset(LGG.GBM.new.order,select=colnames(GBM.LGG.27.450k.nearest.gene)[8:643])
a <- metadata.1129.samples.20150312[substr(colnames(LGG.GBM.new.order2),1,12),]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$flag.categories),as.matrix(a$cluster.cnv),as.matrix(a$mut.TP53),as.matrix(a$cnv.EGFR),as.matrix(a$mut.BRAF),as.matrix(a$histology),as.matrix(a$grade),as.matrix(a$subtype))
clab.6 <- clab.6[,c(1,20,19,18,3,17,16,4,21,14,15,10,11)]
LGG.GBM.new.order2 <- LGG.GBM.new.order2[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
clab.6 <- clab.6[c(1:63,64:327,405:531,532:682,683:932,328:404),]
rownames(normalTCGA.dnaMeth.nearest.gene) <- as.character(normalTCGA.dnaMeth.nearest.gene$Composite.El)
normals.sel <- normalTCGA.dnaMeth.nearest.gene[rownames(LGG.GBM.new.order2),8:117]
norm.label <- matrix("white", nrow = 110, ncol = 13)
norm.label[,10] <- "darkblue"
clab.6 <- rbind(norm.label,clab.6)
rownames(clab.6) <- colnames(LGG.GBM.new.order2)
clab <- merge(clab.6, tert.expr.g[,c(1,3)],all.x=T,by.x=0,by.y="id")
levels(clab$grp) <- c("black","purple")
clab$grp <- as.character(clab$grp)
clab[is.na(clab$grp),"grp"] <- "white"
LGG.GBM.new.order2 <- cbind(normals.sel,LGG.GBM.new.order2)
rownames(clab) <- as.character(clab$Row.names)
clab <- clab[colnames(LGG.GBM.new.order2),]
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_DNAMeth_tert.png",res=1200,width=8000,height=8000)
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_DNAMeth_tert.pdf")
heatmap.plus.sm(as.matrix(LGG.GBM.new.order2),
                col = jet.colors(75),
                scale = "none",
                #labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = as.matrix(clab[,15:14]),
                margins = c(1, 6)
                #RowSideColors = rlab
                
)
dev.off()


oGE.2 <- oGE[,michele.row.order$rowInd]
oGE.2 <- t(oGE.2)
cc.col2 <- cc.col
cc.col2[,4] <- str_replace(cc.col2[,4], "black","darkblue")
cc.col2[,4] <- str_replace(cc.col2[,4], "gray","cyan")
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_LabelsMichele.pdf", width=30, height=20,bg = "transparent")
#layout(matrix(c(1,2,3,4),2,2, byrow=T), widths=c(2,2), heights=c(3,1), respect=T)
heatmap.plus.sm(oGE.2[1:2,],
                col = greenred(75),
                scale = "row",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = cc.col2,
                margins = c(1, 6)
                #RowSideColors = rlab
                
)
heatmap.plus.sm(as.matrix(LGG.GBM.new.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab.6
                #RowSideColors = rlab
                
)
dev.off()

## Key
z <- seq(min(lgg.gbm[,455:932],na.rm=TRUE), max(lgg.gbm[,455:932],na.rm=TRUE), length = 200) #min(lgg.gbm.27.450.order) and #max(lgg.gbm.27.450.order)
n <- max(lgg.gbm[,455:932],na.rm=TRUE)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/key.pdf",width=10, height=05,bg = "transparent")
image(matrix(z, ncol = 1), col = jet.colors(75),
      xaxt = "n", yaxt = "n", main = "Overlap count Key")
box()
par(usr = c(0, n, 0, n))
axis(1, at = c(0,max(z)), c(min(z),max(z)))
dev.off()

z <- seq(min(order,na.rm=TRUE), max(order,na.rm=TRUE), length = 200) #min(lgg.gbm.27.450.order) and #max(lgg.gbm.27.450.order)
n <- max(order,na.rm=TRUE)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/key_DNAMeth.pdf",width=10, height=05,bg = "transparent")
image(matrix(z, ncol = 1), col = jet.colors(75),
      xaxt = "n", yaxt = "n", main = "Overlap count Key")
box()
par(usr = c(0, n, 0, n))
axis(1, at = c(0,max(z)), c(min(z),max(z)))
dev.off()

png("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_2b.png", bg = "transparent", width = 8000, height = 8000, res = 800)
layout(matrix(c(1,2),1,2, byrow=T), widths = c(1.51, 1), heights = c(1, 1.75))
#par(mfrow=c(1,2)+0.1)
image(t(as.matrix(LGG.GBM.new.order)),
      col = jet.colors(75),axes=FALSE)
image(oGE.2,
      col = greenred(75),axes=FALSE)
dev.off()



pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_Michele.pdf")
michele.row.order <- heatmap.2(t(oGE),
                col = greenred(75),
                scale = "row",
                #labRow = NA,
                labCol = NA,
                Colv = F,
                #Rowv = NA,
                #RowSideColors = rlab
                
)
dev.off()

## ES reverse

meth.g.norm <- cpg.gene.normal[,8:117]
names(meth.g.norm) <- substr(names(meth.g.norm), 1, 12)
expr.g.norm <- cpg.gene.normal[,118:227]
names(expr.g.norm) <- substr(names(expr.g.norm), 1, 12)

meth.g <- cpg.gene[,8:643] #only tumor
names(meth.g) <- substr(names(meth.g), 1, 12)
expr.g <- cpg.gene[,644:1279] #only tumor
names(expr.g) <- substr(names(expr.g), 1, 12)
epiGene.rev <- list()

j=1
metadata <- subset(metadata.1129.samples.20141030,Study == "TCGA")
a <- subset(metadata.1129.samples.20141030,id %in% colnames(meth.g))

for(i in 1:dim(cpg.gene.normal)[1]) {
  #for(i in 138:250) {
  
  meth.gu <- names(meth.g)[meth.g[i, ] < 0.3]
  meth.gm <- names(meth.g)[meth.g[i, ] >= 0.3]
  # if(fisher.test(table(yy))$p.value < 0.05){
  if(length(meth.gu)>1 & length(meth.gm)>1){
    expr.mean.gu <- mean(as.numeric(expr.g[i, meth.gu]), na.rm=T)
    if(expr.mean.gu > 1.28*sd(as.numeric(expr.g[i, meth.gm]), na.rm=T)){
      expr.mean.gm <- mean(as.numeric(expr.g[i, meth.gm]), na.rm = T)
      yy <- cbind(t(meth.g[i, ] < 0.3), t(expr.g[i, ] > expr.mean.gm))
      colnames(yy) <- c("meth.gt0.3", "exp.lt.meanUnmeth")
      aux <- cbind(t(meth.g.norm[rownames(meth.g[i,]), ] > 0.3), t(expr.g.norm[rownames(meth.g[i,]), ] < expr.mean.gm))
      colnames(aux) <- c("meth.gt0.3", "exp.lt.meanUnmeth")
      yy <- rbind(yy, aux)
      yy <- as.data.frame(yy)
      gene <- rownames(meth.g.norm[i,])
      if(as.matrix(table(yy$meth.gt0.3[1:636], yy$exp.lt.meanUnmeth[1:636]))[2,2]/as.numeric(table(yy$meth.gt0.3[1:636])[2])>0.8){
        yy$geneName <- gene
        yy$id <- rownames(yy)
        yy$methValues <- NA
        yy$exprValues <- NA
        yy$methValues[1:636] <- as.numeric(meth.g[i, ])
        yy$exprValues[1:636] <- as.numeric(expr.g[i, ])  
        #ic <-248+(491*(k-1))
        #fim <- 247+(491*k)
        yy$methValues[637:nrow(yy)] <- as.vector(as.matrix(t(meth.g.norm[rownames(meth.g[i,]), ])))
        yy$exprValues[637:nrow(yy)] <- as.vector(as.matrix(t(expr.g.norm[rownames(meth.g[i,]), ])))
        #yy$methValues[248:nrow(yy)] <- as.numeric(meth.g.norm[index[k], ])
        #yy$exprValues[248:nrow(yy)] <- as.numeric(expr.g.norm[index[k], ])
        info <- metadata[rownames(yy)[1:636],c("cluster.meth","IDH.status","codel1p19q","cluster.meth.gbm","cluster.meth.lgg","tumor.type")]
        yy$IDH.status <- "WT"
        yy$codel1p19q <- "non-codel"
        yy$cluster <- NA
        yy$cluster.gbm <- NA
        yy$cluster.lgg <- NA
        yy$tumor.type <- NA
        yy$cluster.gbm[1:636] <- as.character(info$cluster.meth.gbm)
        yy$cluster.lgg[1:636] <- as.character(info$cluster.meth.lgg)
        yy$tumor.type[1:636] <- as.character(info$tumor.type)
        yy$codel1p19q[1:636] <- as.character(info$codel1p19q)
        yy$cluster[1:636] <- as.character(info$cluster.meth)
        #yy <- merge(metadata, yy, by.x = "TCGAIDnew", by.y = "id")
        epiGene.rev[[j]] <- as.data.frame(yy)
        names(epiGene.rev)[[j]] <- as.character(gene)
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

uu.nl <- epiGene.rev[[1]]
#uu.nl <- epiGene.rev[["cg24881834.ME1"]]
clusters <- c(paste0("LGm", 1:6))


enrichment.normal <- lapply(clusters, function(cluster, epiGene.ss=epiGene.rev) {
  enrichment.group <- unlist(mclapply(epiGene.ss, function(uu.nl, cluster.group=cluster) {
    #control <- (!uu.nl$meth.gt0.3[637:length(uu.nl$meth.gt0.3)] & !uu.nl$exp.lt.meanUnmeth[637:length(uu.nl$exp.lt.meanUnmeth)])
    control <- sum(uu.nl$methValues[637:length(uu.nl$meth.gt0.3)] > 0.3 & uu.nl$exp.lt.meanUnmeth[637:length(uu.nl$exp.lt.meanUnmeth)])
    uu.nl$es <- NA
    uu.nl$es[1:636] <- uu.nl$meth.gt0.3[1:636] & uu.nl$exp.lt.meanUnmeth[1:636]
    uu.nl$Kclass <- NA
    uu.nl$Kclass <- uu.nl$cluster == cluster.group
    es.table <- table(uu.nl$Kclass, uu.nl$es)#, dnn=c("Class", "ES"))
    #print()
    if((es.table[2,2]/rowSums(es.table)[2] > 0.5) & (es.table[2,2]/colSums(es.table)[2] > 0.5)  & control > 55) {
      return(fisher.test(table(uu.nl$Kclass, uu.nl$es))$p.value)
    } else {
      return(NULL)
    }
    #}
    #}  
  }, mc.cores=1))
  return(enrichment.group)
})



names(enrichment.normal) <- c("LGm1", "LGm2", "LGm3", "LGm4", "LGm5", "LGm6")
enrichment.normal.adj <- lapply(enrichment.normal, p.adjust, method="BH")
enrichmentList.rev <- lapply(enrichment.normal.adj, function(x) {
  x <- x[x < 0.05]
})

##Heatmap ES calls
a <- c(names(enrichmentList.rev$LGm2),names(enrichmentList.rev$LGm5))
#a <- c(names(enrichmentList.sig$LGm3),names(enrichmentList.sig$LGm4))
probe.meth <- str_extract(a,"*[^.]*")
probe.mrna <- str_extract(a,"[^.]*$")



LGG.GBM.new.noXY.s <- LGG.GBM.250914[ as.character(probe.meth),5:936]
LGG.GBM.new.order <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
colnames(LGG.GBM.new.order) <- substr(colnames(LGG.GBM.new.order),1,12)
LGG.GBM.new.order <- LGG.GBM.new.order[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
a <- colnames(LGG.GBM.new.order) %in% colnames(rna.seq.lgg.gbm)
teste <- LGG.GBM.new.order[,a]
LGG.GBM.new.order <- rna.seq.lgg.gbm[as.character(unique(probe.mrna)),colnames(teste)]
LGG.GBM.new.order <- na.omit(LGG.GBM.new.order)
#LGG.GBM.new.order <- subset(LGG.GBM.new.order,select=colnames(GBM.LGG.27.450k.nearest.gene)[8:643])
a <- metadata.1129.samples.20150112[substr(colnames(LGG.GBM.new.order),1,12),]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$flag.categories),as.matrix(a$cluster.cnv),as.matrix(a$mut.TP53),as.matrix(a$cnv.EGFR),as.matrix(a$mut.BRAF),as.matrix(a$histology),as.matrix(a$grade),as.matrix(a$subtype))
clab.6 <- clab.6[,c(1,20,19,18,3,17,16,4,21,14,15,10,11)]
rownames(normalTCGA.dnaMeth.nearest.gene) <- as.character(normalTCGA.dnaMeth.nearest.gene$Composite.El)
normals.sel <- normalTCGA.mrna[rownames(LGG.GBM.new.order),1:110]
norm.label <- matrix("white", nrow = 110, ncol = 12)
norm.label[,10] <- "darkblue"
clab.6 <- rbind(norm.label,clab.6)
LGG.GBM.new.order <- cbind(normals.sel,LGG.GBM.new.order)
LGG.GBM.new.order <- LGG.GBM.new.order + 1
LGG.GBM.new.order <- log2(LGG.GBM.new.order)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_RNAseqRev.pdf")
heatmap.plus.sm(as.matrix(LGG.GBM.new.order),
                col=rev(redgreen(75)),
                scale = "none",
                #trace = "none",
                #labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab.6
                #RowSideColors = rlab
                
)
dev.off()

LGG.GBM.new.noXY.s <- LGG.GBM.250914[ as.character(probe.meth),5:936]
LGG.GBM.new.order2 <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
colnames(LGG.GBM.new.order2) <- substr(colnames(LGG.GBM.new.order2),1,12)
#a <- colnames(LGG.GBM.new.order) %in% colnames(GBM.LGG.27.450k.nearest.gene)
#teste <- LGG.GBM.new.order[,a]
#LGG.GBM.new.order <- subset(LGG.GBM.new.order,select=colnames(GBM.LGG.27.450k.nearest.gene)[8:643])
a <- metadata.1129.samples.20150112[substr(colnames(LGG.GBM.new.order2),1,12),]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$flag.categories),as.matrix(a$cluster.cnv),as.matrix(a$mut.TP53),as.matrix(a$cnv.EGFR),as.matrix(a$mut.BRAF),as.matrix(a$histology),as.matrix(a$grade),as.matrix(a$subtype))
clab.6 <- clab.6[,c(1,20,19,18,3,17,16,4,21,14,15,10,11)]
LGG.GBM.new.order2 <- LGG.GBM.new.order2[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
clab.6 <- clab.6[c(1:63,64:327,405:531,532:682,683:932,328:404),]
rownames(normalTCGA.dnaMeth.nearest.gene) <- as.character(normalTCGA.dnaMeth.nearest.gene$Composite.El)
normals.sel <- normalTCGA.dnaMeth.nearest.gene[rownames(LGG.GBM.new.order2),8:117]
norm.label <- matrix("white", nrow = 110, ncol = 13)
norm.label[,10] <- "darkblue"
clab.6 <- rbind(norm.label,clab.6)
LGG.GBM.new.order2 <- cbind(normals.sel,LGG.GBM.new.order2)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_DNAMethRev.pdf")
heatmap.plus.sm(as.matrix(LGG.GBM.new.order2),
                col = jet.colors(75),
                scale = "none",
                #labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab.6,
                margins = c(1,6)
                #RowSideColors = rlab
                
)
dev.off()



###### OLD

cpg.gene.norm <- t(apply(cpg.gene[,255:501], 1, function(x) { (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) }))
cpg.gene.norm[is.nan(cpg.gene.norm)] <- 0
cpg.gene.norm <- cbind(cpg.gene[,1:254],cpg.gene.norm)
cpg.gene.normal.norm <- t(apply(cpg.gene.normal[,496:986], 1, function(x) {  (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) }))
cpg.gene.normal.norm[is.nan(cpg.gene.normal.norm)] <- 0
cpg.gene.normal.norm <- cbind(cpg.gene.normal[,1:495],cpg.gene.normal.norm)
#cpg.gene.norm <- cpg.gene.norm[c(names(tent.enrichmentList.sig$LGm3),names(tent.enrichmentList.sig$LGm4),names(tent.enrichmentList.sig$LGm5)),]
#index <- unlist(lapply(c(names(tent.enrichmentList.sig$LGm3),names(tent.enrichmentList.sig$LGm4),names(tent.enrichmentList.sig$LGm5)), function (x) { grep(x,rownames(cpg.gene.normal)) }))
#cpg.gene.normal.norm <- cpg.gene.normal[index,]

meth.g.norm <- cpg.gene.normal.norm[,3:493]
colnames(meth.g.norm) <- substr(colnames(meth.g.norm), 1, (nchar(names(meth.g.norm))-2))
expr.g.norm <- cpg.gene.normal.norm[,496:986]
colnames(expr.g.norm) <- substr(colnames(expr.g.norm), 1, (nchar(names(expr.g.norm))-2))

meth.g <- cpg.gene.norm[,8:254]
names(meth.g) <- substr(names(meth.g), 1, 12)
expr.g <- cpg.gene.norm[,255:501]
names(expr.g) <- substr(names(expr.g), 1, 12)
epiGene.2 <- list()
j=1
rownames(metadata) <- metadata$TCGAIDnew
for(i in 1:dim(cpg.gene.norm)[1]) {
  #for(i in 138:142) {
  
  meth.gu <- names(meth.g)[meth.g[i, ] < 0.3]
  meth.gm <- names(meth.g)[meth.g[i, ] >= 0.3]
  index <- grep(rownames(meth.g[i,]),rownames(cpg.gene.normal.norm))
  # if(fisher.test(table(yy))$p.value < 0.05){
  if(length(meth.gu)>1 & length(meth.gm)>1 & length(index) > 0){
    expr.mean.gm <- mean(as.numeric(expr.g[i, meth.gm]), na.rm=T)
    if(expr.mean.gm < 1.28*sd(as.numeric(expr.g[i, meth.gu]), na.rm=T)){
      expr.mean.gu <- mean(as.numeric(expr.g[i, meth.gu]), na.rm = T)
      yy <- cbind(t(meth.g[i, ] >= 0.3), t(expr.g[i, ] < expr.mean.gu))
      colnames(yy) <- c("meth.gt0.3", "exp.lt.meanUnmeth")
      for(k in 1:length(index)){
        aux <- cbind(t(meth.g.norm[index[k], ] >= 0.3), t(expr.g.norm[index[k], ] < expr.mean.gu))
        colnames(aux) <- c("meth.gt0.3", "exp.lt.meanUnmeth")
        yy <- rbind(yy, aux)
        yy <- as.data.frame(yy)
        gene <- rownames(meth.g.norm[index[k],])
      }
      if(as.matrix(table(yy$meth.gt0.3[1:247], yy$exp.lt.meanUnmeth[1:247]))[2,2]/as.numeric(table(yy$meth.gt0.3[1:247])[2])>0.8){
        yy$geneName <- gene
        yy$id <- rownames(yy)
        yy$methValues <- NA
        yy$exprValues <- NA
        yy$methValues[1:247] <- as.numeric(meth.g[i, ])
        yy$exprValues[1:247] <- as.numeric(expr.g[i, ])  
        #ic <-248+(491*(k-1))
        #fim <- 247+(491*k)
        yy$methValues[248:nrow(yy)] <- as.vector(as.matrix(t(meth.g.norm[index, ])))
        yy$exprValues[248:nrow(yy)] <- as.vector(as.matrix(t(expr.g.norm[index, ])))
        #yy$methValues[248:nrow(yy)] <- as.numeric(meth.g.norm[index[k], ])
        #yy$exprValues[248:nrow(yy)] <- as.numeric(expr.g.norm[index[k], ])
        info <- metadata[rownames(yy)[1:247],c("cluster.meth","IDH.status","codel1p19q")]
        yy$IDH.status <- "WT"
        yy$codel1p19q <- "non-codel"
        yy$cluster <- NA
        ix <- grep("CRBLM",rownames(yy))
        yy$cluster[ix] <- "Cerebellum"
        ix <- grep("FCTX",rownames(yy))
        yy$cluster[ix] <- "Frontal.cortex"
        ix <- grep("PONS",rownames(yy))
        yy$cluster[ix] <- "Pons"
        ix <- grep("TCTX",rownames(yy))
        yy$cluster[ix] <- "Temporal.cortex"
        yy$IDH.status[1:247] <- as.character(info$IDH.status)
        yy$codel1p19q[1:247] <- as.character(info$codel1p19q)
        yy$cluster[1:247] <- as.character(info$cluster.meth)
        #yy <- merge(metadata, yy, by.x = "TCGAIDnew", by.y = "id")
        epiGene.2[[j]] <- as.data.frame(yy)
        names(epiGene.2)[[j]] <- as.character(gene)
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


uu.nl <- epiGene.2[[1]]
clusters <- c(paste0("LGm", 1:6))
#normal <- epiGene.2
#i<-1

enrichment.normal <- lapply(clusters, function(cluster, epiGene.ss=epiGene.2) {
  enrichment.group <- unlist(mclapply(epiGene.ss, function(uu.nl, cluster.group=cluster) {
    control <- (!uu.nl$meth.gt0.3[248:length(uu.nl$meth.gt0.3)] & !uu.nl$exp.lt.meanUnmeth[248:length(uu.nl$exp.lt.meanUnmeth)])
    uu.nl$es <- NA
    uu.nl$es[1:247] <- uu.nl$meth.gt0.3[1:247] & uu.nl$exp.lt.meanUnmeth[1:247]
    uu.nl$Kclass <- NA
    uu.nl$Kclass <- uu.nl$cluster == cluster.group
    es.table <- table(uu.nl$Kclass, uu.nl$es)#, dnn=c("Class", "ES"))
    if((es.table[2,2]/rowSums(es.table)[2] > 0.5) & (es.table[2,2]/colSums(es.table)[2] > 0.5)  & (table(control)[2]/(table(control)[2] + table(control)[1]) > 0.5)) {
      return(fisher.test(table(uu.nl$Kclass, uu.nl$es))$p.value)
    } else {
      return(NULL)
    }
    #}
    #}  
  }, mc.cores=1))
  return(enrichment.group)
})



names(enrichment.normal) <- c("LGm1", "LGm2", "LGm3", "LGm4", "LGm5", "LGm6")
enrichment.normal.adj <- lapply(enrichment.normal, p.adjust, method="BH")
enrichmentList.sig <- lapply(enrichment.normal.adj, function(x) {
  x <- x[x < 0.05]
})

v <- lapply(enrichmentList.sig, function(x) {
  x <- x[x < 1e-20]
})

es.genes <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/ES_genes.txt")
es.genes$gene.name <- paste(paste(".",es.genes$gene.name,sep=""),".",sep="")
for(i in 1:360){
  teste <- grep(es.genes$gene.name[i],c(names(enrichmentList.sig$LGm2),names(enrichmentList.sig$LGm3),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm5)))
  print(teste)
}
metadata <- metadata[colnames(heatmap.es.order),]
clab.6 <- colors.label(as.matrix(metadata$K6),as.matrix(metadata$DNA.Methylation.Clusters),as.matrix(metadata$MethylationCluster),as.matrix(metadata$tumor.type),as.matrix(metadata$age),as.matrix(metadata$IDH.status),as.matrix(metadata$Cluster),as.matrix(metadata$cluster.mRNA),as.matrix(metadata$COCCluster),as.numeric(metadata$leukocyte),as.matrix(metadata$rna.estimates))
clab.6 <- clab.6[c(856:905,1:333,794:855,334:489,490:708,709:793),]
heatmap.es.order<- heatmap.es.order[,c(856:905,1:333,794:855,334:489,490:708,709:793)]
probe.ex <- rownames(heatmap.es.order)[rownames(heatmap.es.order) %nin% rownames(norms)]
norms <- rbind(norms,NA)
rownames(norms)[392868] <- probe.ex
norm.sel <- norms[rownames(heatmap.es.order),]
e <- merge(norm.sel,heatmap.es.order,by=0)
norm.label <- matrix("white", nrow = 77, ncol = 11)
heatmap.es.order <- cbind(norm.sel,heatmap.es.order)
clab.6 <- rbind(norm.label,clab.6)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_905_s/heatmap_ES3_normal.png",res=500,width=2000,height=2000)
heatmap.plus.sm(as.matrix(heatmap.es.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab.6,
                RowSideColors = cbind(rlab,rlab)
)
dev.off()

epiplot <- epiGene.RNAseq[["cg24719575.REST"]]
epiplot[is.na(epiplot$cluster),"cluster"] <- "Non-tumor"
epiplot[is.na(epiplot$tumor.type),"tumor.type"] <- "Non-tumor"
epiplot$exprValues <- epiplot$exprValues + 1
aux <- subset(epiplot[1:636,], meth.gt0.3==FALSE)
expr.mean.gu <- mean(as.numeric(aux[, "exprValues"]), na.rm = T)
png("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RNAseq_05-11-14/x.png",res=200,width=4000,height=1000)
p <- ggplot(epiplot, aes(methValues, exprValues))
p + 
  stat_smooth(method = "lm", fullrange=FALSE, na.rm=TRUE, se=FALSE, size=1) +
  geom_vline(xintercept = 0.3, linetype = "longdash", colour="black") +
  geom_hline(yintercept = expr.mean.gu, linetype = "longdash", colour="black") +
  geom_point(size=3, aes(colour = factor(cluster), shape=tumor.type,size=factor(cluster.gbm))) +
  scale_colour_manual(values = c("LGm1" = "green","LGm2" = "red","LGm3" = "purple", "LGm4" = "orange", "LGm5" = "yellow", "LGm6" = "blue", "Non-tumor" = "gray")) +
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.3, 0.5, 0.7, 1)) +
  scale_y_log10() +
  scale_size("cluster.gbm") +
  ylab("RNA seq expression (log2)") +
  xlab(paste("DNA Methylation Values\n")) + 
  labs(title = "DNA Methylation vs Expression", colour = "DNA Meth. Cluster") 
dev.off()

pdata <- read.delim("LGG-GBM.merged.data.txt")
aux <- epiGene.RNAseq[["cg06190732.RP11-986E7.7"]]
aux[is.na(aux$cluster),"cluster"] <- "Non-tumor"
aux[is.na(aux$tumor.type),"tumor.type"] <- "Non-tumor"
#aux[is.na(aux$IDH.status) & aux$tumor.type %in% "Non-tumor","IDH.status"] <- "WT"
#aux[1:636,"IDH.status"] <- as.character(metadata.1129.samples.20150312[rownames(metadata.1129.samples.20150312) %in% as.character(rownames(aux)),"IDH.status"])
aux[1:636,"IDH.status"] <- as.character(pdata[as.character(pdata$case.id) %in% as.character(rownames(aux)),"IDH.status"])
aux[637:746,"IDH.status"] <- "WT"
b <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/probesLGm2.txt",header=F)
#c <- subset(volcano, threshold == 3) #Hypo no LGm6
c <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/EReg5.txt")
c <- as.character(c$id)
a <- c(c,names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm3),names(sort(enrichmentList.sig$LGm2)[1:15]),as.character(b$V1))
d <- cpg.gene
rownames(d) <- as.character(d$Composite.El)
probe.meth <- str_extract(a,"[a-z0-9]*[^.]*")
d <- d[probe.meth,]
d <- na.omit(d)
probes <- paste(as.character(d$Composite.El),as.character(d$Associated.Gene.Name),sep=".")
#expr.v <- sapply(epiGene.RNAseq, `[[`, 'exprValues')
#expr.v <- as.data.frame(expr.v)
#rownames(expr.v) <- rownames(aux)
expr.lgm <- d[as.character(probe.meth),644:1279]
colnames(expr.lgm) <- substr(colnames(expr.lgm),1,12)
expr.lgm <- expr.lgm[,rownames(aux)[1:636]]
expr.lgm <- na.omit(expr.lgm)
##falta arrumar aqui tudo probeid, nao gene..no d
mean.ESg1.expr <- apply(expr.lgm[57:65,], 2, mean, na.rm=TRUE)
mean.ESg2.expr <- apply(expr.lgm[42:56,], 2, mean, na.rm=TRUE)
mean.ESg3.expr <- apply(expr.lgm[27:41,], 2, mean, na.rm=TRUE)
mean.ESg4.expr <- apply(expr.lgm[12:26,], 2, mean, na.rm=TRUE)
mean.ESg5.expr <- apply(expr.lgm[1:11,], 2, mean, na.rm=TRUE)
#mean.EReg <- c(mean.ESg1,mean.ESg2,mean.ESg3,mean.ESg4,mean.ESg5)
#mean.EReg.expr <- as.matrix(mean.EReg)

#meth.v <- sapply(epiGene.RNAseq, `[[`, 'methValues')
#meth.v <- as.data.frame(meth.v)
#rownames(meth.v) <- rownames(aux)
#meth.lgm <- meth.v[,as.character(probes)]
meth.lgm <- d[rownames(expr.lgm),8:643]
colnames(meth.lgm) <- substr(colnames(meth.lgm),1,12)
meth.lgm <- meth.lgm[,rownames(aux)[1:636]]
mean.ESg1.meth <- apply(meth.lgm[57:65,], 2, mean, na.rm=TRUE)
mean.ESg2.meth <- apply(meth.lgm[42:56,], 2, mean, na.rm=TRUE)
mean.ESg3.meth <- apply(meth.lgm[27:41,], 2, mean, na.rm=TRUE)
mean.ESg4.meth <- apply(meth.lgm[12:26,], 2, mean, na.rm=TRUE)
mean.ESg5.meth <- apply(meth.lgm[1:11,], 2, mean, na.rm=TRUE)
#mean.EReg.meth <- c(mean.ESg1,mean.ESg2,mean.ESg3,mean.ESg4,mean.ESg5)
#mean.EReg.meth <- as.matrix(mean.EReg)
#meth.lgm$mean.group1 <- apply(meth.lgm[,names(enrichmentList.sig$LGm2)], 1, mean, na.rm=TRUE)
#meth.lgm$mean.group2 <- apply(meth.lgm[,names(enrichmentList.sig$LGm3)], 1, mean, na.rm=TRUE)
#meth.lgm$mean.group3 <- apply(meth.lgm[,c(names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm5))], 1, mean, na.rm=TRUE)

epiplot.group1 <- cbind(mean.ESg1.meth,mean.ESg1.expr)
epiplot.group2 <- cbind(mean.ESg2.meth,mean.ESg2.expr)
epiplot.group3 <- cbind(mean.ESg3.meth,mean.ESg3.expr)
epiplot.group4 <- cbind(mean.ESg4.meth,mean.ESg4.expr)
epiplot.group5 <- cbind(mean.ESg5.meth,mean.ESg5.expr)
epiplot.group1 <- as.data.frame(epiplot.group1)
epiplot.group2 <- as.data.frame(epiplot.group2)
epiplot.group3 <- as.data.frame(epiplot.group3)
epiplot.group4 <- as.data.frame(epiplot.group4)
epiplot.group5 <- as.data.frame(epiplot.group5)
colnames(epiplot.group1) <- c("methValues","exprValues")
colnames(epiplot.group2) <- c("methValues","exprValues")
colnames(epiplot.group3) <- c("methValues","exprValues")
colnames(epiplot.group4) <- c("methValues","exprValues")
colnames(epiplot.group5) <- c("methValues","exprValues")
rownames(pdata) <- as.character(pdata$case.id)
epiplot.group1 <- merge(epiplot.group1,pdata[,c("cluster.meth","IDH.status")],by=0,all.x=T)
epiplot.group2 <- merge(epiplot.group2,pdata[,c("cluster.meth","IDH.status")],by=0,all.x=T)
epiplot.group3 <- merge(epiplot.group3,pdata[,c("cluster.meth","IDH.status")],by=0,all.x=T)
epiplot.group4 <- merge(epiplot.group4,pdata[,c("cluster.meth","IDH.status")],by=0,all.x=T)
epiplot.group5 <- merge(epiplot.group5,pdata[,c("cluster.meth","IDH.status")],by=0,all.x=T)
epiplot.group1$group <- "group1"
epiplot.group2$group <- "group2"
epiplot.group3$group <- "group3"
epiplot.group4$group <- "group4"
epiplot.group5$group <- "group5"
epiplot <- rbind(epiplot.group1,epiplot.group2,epiplot.group3,epiplot.group4,epiplot.group5)
g1 <- ggplot(epiplot, aes(methValues, exprValues)) +
  stat_smooth(method = "lm", fullrange=FALSE, na.rm=TRUE, se=FALSE, size=1) +
  geom_vline(xintercept = 0.3, linetype = "longdash", colour="black") +
  #geom_hline(yintercept = expr.mean.gu, linetype = "longdash", colour="black") +
  geom_point(size=3, aes(colour = factor(cluster.meth), shape=IDH.status)) +
  scale_colour_manual(values = c("LGm1" = "green","LGm2" = "red","LGm3" = "purple", "LGm4" = "orange", "LGm5" = "yellow", "LGm6" = "blue", "Non-tumor" = "gray")) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_log10() +
  #scale_size("cluster.gbm") +
  ylab("mean RNAseq") +
  xlab(paste("DNA Methylation Values\n")) + 
  labs(title = "", colour = "DNA Meth. Cluster") +
  #theme(legend.position = "none") +
  facet_grid(. ~ group,scales="free", space="free")
ggsave(file="epiplot2.pdf")


epiplot.group1 <- cbind(meth.lgm$mean.group1,expr.lgm$mean.group1)
epiplot.group2 <- cbind(meth.lgm$mean.group2,expr.lgm$mean.group2)
epiplot.group3 <- cbind(meth.lgm$mean.group3,expr.lgm$mean.group3)
epiplot.group1 <- as.data.frame(epiplot.group1)
epiplot.group2 <- as.data.frame(epiplot.group2)
epiplot.group3 <- as.data.frame(epiplot.group3)
rownames(epiplot.group1) <- rownames(meth.lgm)
rownames(epiplot.group2) <- rownames(meth.lgm)
rownames(epiplot.group3) <- rownames(meth.lgm)
colnames(epiplot.group1) <- c("methValues","exprValues")
colnames(epiplot.group2) <- c("methValues","exprValues")
colnames(epiplot.group3) <- c("methValues","exprValues")
epiplot.group1$cluster <- as.character(aux$cluster)
epiplot.group2$cluster <- as.character(aux$cluster)
epiplot.group3$cluster <- as.character(aux$cluster)
epiplot.group1$tumor.type <- as.character(aux$tumor.type)
epiplot.group2$tumor.type <- as.character(aux$tumor.type)
epiplot.group3$tumor.type <- as.character(aux$tumor.type)
epiplot.group1$IDH.status <- as.character(aux$IDH.status)
epiplot.group2$IDH.status <- as.character(aux$IDH.status)
epiplot.group3$IDH.status <- as.character(aux$IDH.status)

epiplot.group1$group <- "group1"
epiplot.group2$group <- "group2"
epiplot.group3$group <- "group3"
epiplot <- rbind(epiplot.group1,epiplot.group2,epiplot.group3)
#epiplot$exprValues <- epiplot$exprValues + 1
#aux <- subset(epiplot[1:636,], meth.gt0.3==FALSE)
#expr.mean.gu <- mean(as.numeric(aux[, "exprValues"]), na.rm = T)
g1 <- ggplot(epiplot, aes(methValues, exprValues)) +
  stat_smooth(method = "lm", fullrange=FALSE, na.rm=TRUE, se=FALSE, size=1) +
  geom_vline(xintercept = 0.3, linetype = "longdash", colour="black") +
  #geom_hline(yintercept = expr.mean.gu, linetype = "longdash", colour="black") +
  geom_point(size=3, aes(colour = factor(cluster), shape=IDH.status)) +
  scale_colour_manual(values = c("LGm1" = "green","LGm2" = "red","LGm3" = "purple", "LGm4" = "orange", "LGm5" = "yellow", "LGm6" = "blue", "Non-tumor" = "gray")) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_log10() +
  #scale_size("cluster.gbm") +
  ylab("mean RNAseq") +
  xlab(paste("DNA Methylation Values\n")) + 
  labs(title = "LGm1/LGm2/LGm3", colour = "DNA Meth. Cluster") +
  #theme(legend.position = "none") +
  facet_grid(. ~ group,scales="free", space="free")
ggsave(file="epiplot2.pdf")

g2 <- ggplot(epiplot.group2, aes(methValues, exprValues)) +
  stat_smooth(method = "lm", fullrange=FALSE, na.rm=TRUE, se=FALSE, size=1) +
  geom_vline(xintercept = 0.3, linetype = "longdash", colour="black") +
  geom_hline(yintercept = expr.mean.gu, linetype = "longdash", colour="black") +
  geom_point(size=3, aes(colour = factor(cluster), shape=tumor.type)) +
  scale_colour_manual(values = c("LGm1" = "green","LGm2" = "red","LGm3" = "purple", "LGm4" = "orange", "LGm5" = "yellow", "LGm6" = "blue", "Non-tumor" = "gray")) +
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.3, 0.5, 0.7, 1)) +
  #scale_y_log10() +
  #scale_size("cluster.gbm") +
  #ylab("mean RNAseq") +
  #xlab(paste("DNA Methylation Values\n")) + 
  labs(title = "LGm3", colour = "DNA Meth. Cluster") +
  theme(legend.position = "none")

g3 <- ggplot(epiplot.group3, aes(methValues, exprValues)) +
  stat_smooth(method = "lm", fullrange=FALSE, na.rm=TRUE, se=FALSE, size=1) +
  geom_vline(xintercept = 0.3, linetype = "longdash", colour="black") +
  geom_hline(yintercept = expr.mean.gu, linetype = "longdash", colour="black") +
  geom_point(size=3, aes(colour = factor(cluster), shape=tumor.type)) +
  scale_colour_manual(values = c("LGm1" = "green","LGm2" = "red","LGm3" = "purple", "LGm4" = "orange", "LGm5" = "yellow", "LGm6" = "blue", "Non-tumor" = "gray")) +
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.3, 0.5, 0.7, 1)) +
  #scale_y_log10() +
  #scale_size("cluster.gbm") +
  #ylab("RNA seq expression (log2)") +
 # xlab(paste("DNA Methylation Values\n")) + 
  labs(title = "LGm4/LGm5", colour = "DNA Meth. Cluster") +
  theme(legend.position = "none")

pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RNAseq_05-11-14/es_group.pdf", width=25, height=20,bg = "transparent")
multiplot(g1,g2,g3)
dev.off()

### Heatmap Sturm DNA Methylation
a <- c(names(enrichmentList.sig$LGm2),names(enrichmentList.sig$LGm3),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm5))
#a <- c(names(enrichmentList.rev$LGm2),names(enrichmentList.rev$LGm5))
probe.meth <- str_extract(a,"*[^.]*")
sturm.samples <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Sturm_GBM_GSE36278.txt",skip=139)
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/gbm450.Rda")
#sturm.samples <- na.omit(sturm.samples)
sturm.samples.nocontrol <- sturm.samples[,-c(77:82)]
rownames(sturm.samples.nocontrol) <- sturm.samples.nocontrol$ID_REF
colnames(sturm.samples.nocontrol)[2:137] <- as.character(rownames(RF.Sturm.on.GBM)[c(1:75,82:142)])
sturm.tcga <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Sturm_TCGA.txt")
colnames(gbm.450k) <- substr(colnames(gbm.450k),1,12)
sturm.tcga.samples <- gbm.450k[,as.character(sturm.tcga$id)]
sturm.tcga.samples$probe.id <- as.character(gbm.450k[,1])
sturm.210 <- merge(sturm.samples.nocontrol,sturm.tcga.samples,by.x="ID_REF",by.y="probe.id")
rownames(sturm.210) <- as.character(sturm.210$ID_REF)
b <- metadata.1129.samples.20150112[as.character(colnames(sturm.210)[2:211]),]
idh <- subset(b, GBM.Cluster.Sturm.et.al == "IDH")
k27 <- subset(b, GBM.Cluster.Sturm.et.al == "K27")
G34 <- subset(b, GBM.Cluster.Sturm.et.al == "G34")
pdgfra <- subset(b, GBM.Cluster.Sturm.et.al == "RTK I 'PDGFRA'")
mesen <- subset(b, GBM.Cluster.Sturm.et.al == "Mesenchymal")
classic <- subset(b, GBM.Cluster.Sturm.et.al == "RTK II 'Classic'")
order <- c(as.character(idh$id),as.character(k27$id),as.character(G34$id),as.character(pdgfra$id),as.character(mesen$id),as.character(classic$id))
b$GBM.Cluster.Sturm.et.al <- factor(b$GBM.Cluster.Sturm.et.al)
levels(b$GBM.Cluster.Sturm.et.al) <- c("blue","red","green","yellow","orange","brown")
b$cluster.meth.all <- factor(b$cluster.meth.all)
levels(b$cluster.meth.all) <- c("green","red","purple","orange","yellow","blue")
sturm.hmp <- sturm.210[as.character(probe.meth),2:211]
#png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmap_EA_Sturm_Order.png",res=800,width=4000,height=4000)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_Sturm_Order.png",res=800,width=4000,height=4000)
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmap_ES_Sturm_Labels.pdf")
heatmap.plus.sm(as.matrix(sturm.hmp[,order]),
                col = jet.colors(75),
                scale = "none",
               #labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = as.matrix(b[order,c("cluster.meth.all","GBM.Cluster.Sturm.et.al")])
                #RowSideColors = rlab
                
)
dev.off()

##Survival Sturm
b <- metadata.1129.samples.20141209[as.character(colnames(sturm.210)[2:211]),]
b$vital <- factor(b$vital)
b$gender <- factor(b$gender)
levels(b$gender) <- c("female","male","female","male")
pd = subset(b, GBM.Cluster.Sturm.et.al %in% c("IDH", "RTK II 'Classic'"))
table(pd$GBM.Cluster.Sturm.et.al) 

library(survival)
pd$s = pd$vital == "dead"
pd$os = pd$os.months
pd$GBM.Cluster.Sturm.et.al <- factor(pd$GBM.Cluster.Sturm.et.al)
f = formula(Surv(os,s) ~ GBM.Cluster.Sturm.et.al)
fit = survfit(f, data=pd)
plot(fit, col=c("red","brown"))
legend("bottomleft", levels(pd$GBM.Cluster.Sturm.et.al), fill=c("red","brown"))
diff = survdiff(f, data=pd)
p = round(1 - pchisq(diff$chisq, length(diff$n) - 1), 4)
text(locator(1), sprintf("P=%s", p))

coxph(f, data=pd) ## Should be similar to survdiff (logrank test)

f2 = formula(Surv(os,s) ~ age)
f2.5 = formula(Surv(os,s) ~ GBM.Cluster.Sturm.et.al + age)
model1 = coxph(f2, data=pd) ## survival by age
model2 = coxph(f2.5, data=pd) ## By age and CIMPclass
anova(model2,model1) #does adding cimpclass to age alone add a significant predictor??

##### Heatmap paper
colors.label <- function(lgg.cluster.new,gbm.cluster.new,p19q,age,cluster,tumor,IDH,rna.dna.c,rna.cluster,coc.c,leukocyte,estimates,sturm,anno,cnv,tp53,egfr,braf,histology,grade,subtype){
  
  for(i in 1:length(cluster)){
    if(is.na(cluster[i]))
      cluster[i] <- "white"
    else{ 
      if(cluster[i] == "LGm1")
        cluster[i] <- "green"
      else if(cluster[i] == "LGm2")
        cluster[i] <- "red"
      else if(cluster[i] == "LGm3")
        cluster[i] <- "purple"
      else if(cluster[i] == "LGm4")
        cluster[i] <- "orange"
      else if(cluster[i] == "LGm5")
        cluster[i] <- "yellow"
      else if(cluster[i] == "LGm6")
        cluster[i] <- "blue"
      else 
        cluster[i] <- "white"
    }
  }
  
  for(i in 1:length(tp53)){
    if(is.na(tp53[i]))
      tp53[i] <- "white"
    else{
      if(tp53[i] == "Mut")
        tp53[i] <- "black"
      else if(tp53[i] == "WT")
        tp53[i] <- "grey"
      else
        tp53[i] <- "white"
      
    }
  }
  
  
  for(i in 1:length(egfr)){
    if(is.na(egfr[i]))
      egfr[i] <- "white"
    else{
      if(egfr[i] == "Hemi-del")
        egfr[i] <- "blue"
      # else if(egfr[i] == "Homo-del")
      #  egfr[i] <- "cadetblue4"
      else if(egfr[i] == "Low-amp")
        egfr[i] <- "pink"
      else if(egfr[i] == "High-amp")
        egfr[i] <- "red"
      else if(egfr[i] == "Neutral")
        egfr[i] <- "gray"
      else
        egfr[i] <- "white"
      
    }
  }
  
  for(i in 1:length(braf)){
    if(is.na(braf[i]))
      braf[i] <- "white"
    else{
      if(braf[i] == "Mut")
        braf[i] <- "black"
      else if(braf[i] == "WT")
        braf[i] <- "grey"
      else
        braf[i] <- "white"
      
    }
  }
  
  
  for(i in 1:length(histology)){
    if(is.na(histology[i]))
      histology[i] <- "white"
    else{
      if(histology[i] == "astrocytoma")
        histology[i] <- "red"
      else if(histology[i] == "glioblastoma")
        histology[i] <- "purple"
      else if(histology[i] == "oligoastrocytoma")
        histology[i] <- "cyan"
      else if(histology[i] == "oligodendroglioma")
        histology[i] <- "green3"
      #else if(histology[i] == "pilocytic astrocytoma")
      # histology[i] <- "deeppink3"
      else
        histology[i] <- "white"
      
    }
  }
  
  for(i in 1:length(grade)){
    if(is.na(grade[i]))
      grade[i] <- "white"
    else{
      if(grade[i] == "G2")
        grade[i] <- "green3"
      else if(grade[i] == "G3")
        grade[i] <- "red"
      else if(grade[i] == "G4")
        grade[i] <- "blue"
      else
        grade[i] <- "white"
      
    }
  }
  
  # for(i in 1:length(gbm.cluster)){
  #  if(is.na(gbm.cluster[i]))
  #   gbm.cluster[i] <- "white"
  #else{
  #      if(gbm.cluster[i] == "G-CIMP")
  #       gbm.cluster[i] <- "yellow"
  #    else if(gbm.cluster[i] == "M1")
  #     gbm.cluster[i] <- "green"
  #  else if(gbm.cluster[i] == "M2")
  #        gbm.cluster[i] <- "black"
  #     else if(gbm.cluster[i] == "M3")
  #      gbm.cluster[i] <- "blue"
  #   else if(gbm.cluster[i] == "M4")
  #    gbm.cluster[i] <- "red"
  #      else if(gbm.cluster[i] == "M6")
  #       gbm.cluster[i] <- "darksalmon"
  #    else 
  #     gbm.cluster[i] <- "white"
  #}
  #}
  
  for(i in 1:length(cnv)){
    if(is.na(cnv[i]))
      cnv[i] <- "white"
    else{ 
      if(cnv[i] == "LGc1")
        cnv[i] <- "red"
      else if(cnv[i] == "LGc2")
        cnv[i] <- "purple"
      else if(cnv[i] == "LGc3")
        cnv[i] <- "cyan"
      else if(cnv[i] == "LGc4")
        cnv[i] <- "magenta"
      else if(cnv[i] == "LGc5")
        cnv[i] <- "orange"
      else if(cluster[i] == "LGc6")
        cnv[i] <- "green3"
      else 
        cnv[i] <- "white"
    }
  }
  
  
  
  for(i in 1:length(gbm.cluster.new)){
    if(is.na(gbm.cluster.new[i]))
      gbm.cluster.new[i] <- "white"
    else{
      if(gbm.cluster.new[i] == "Classical")
        gbm.cluster.new[i] <- "red"
      else if(gbm.cluster.new[i] == "G-CIMP")
        gbm.cluster.new[i] <- "yellow"
      else if(gbm.cluster.new[i] == "Mesenchymal")
        gbm.cluster.new[i] <- "green"
      else if(gbm.cluster.new[i] == "Neural")
        gbm.cluster.new[i] <- "blue"
      else if(gbm.cluster.new[i] == "Proneural")
        gbm.cluster.new[i] <- "orange"
      else 
        gbm.cluster.new[i] <- "white"
    }
  }
  
  #  for(i in 1:length(lgg.cluster)){
  #   if(is.na(lgg.cluster[i]))
  #    lgg.cluster[i] <- "white"
  #  else{
  #   if(lgg.cluster[i] == "M5")
  #    lgg.cluster[i] <- "darkorchid4"
  #  else if(lgg.cluster[i] == "M1")
  #   lgg.cluster[i] <- "black"
  #  else if(lgg.cluster[i] == "M2")
  #      lgg.cluster[i] <- "red"
  #   else if(lgg.cluster[i] == "M3")
  #   lgg.cluster[i] <- "blue"
  #else if(lgg.cluster[i] == "M4")
  #        lgg.cluster[i] <- "forestgreen"
  #     else 
  #      lgg.cluster[i] <- "white"
  # }
  #}
  
  for(i in 1:length(subtype)){
    if(is.na(subtype[i]))
      subtype[i] <- "white"
    else{ 
      if(subtype[i] == "Classical")
        subtype[i] <- "red"
      else if(subtype[i] == "G-CIMP")
        subtype[i] <- "black"
      else if(subtype[i] == "IDHmut-codel")
        subtype[i] <- "cyan"
      else if(subtype[i] == "IDHmut-non-codel")
        subtype[i] <- "tomato"
      else if(subtype[i] == "IDHwt")
        subtype[i] <- "gold"
      else if(subtype[i] == "Mesenchymal")
        subtype[i] <- "green"
      else if(subtype[i] == "Neural")
        subtype[i] <- "blue"
      else if(subtype[i] == "Proneural")
        subtype[i] <- "orange"
      else 
        subtype[i] <- "white"
    }
  }
  
  
  for(i in 1:length(lgg.cluster.new)){
    if(is.na(lgg.cluster.new[i]))
      lgg.cluster.new[i] <- "white"
    else{
      if(lgg.cluster.new[i] == "M5")
        lgg.cluster.new[i] <- "darkorchid4"
      else if(lgg.cluster.new[i] == "M1")
        lgg.cluster.new[i] <- "black"
      else if(lgg.cluster.new[i] == "M2")
        lgg.cluster.new[i] <- "red"
      else if(lgg.cluster.new[i] == "M3")
        lgg.cluster.new[i] <- "blue"
      else if(lgg.cluster.new[i] == "M4")
        lgg.cluster.new[i] <- "forestgreen"
      else 
        lgg.cluster.new[i] <- "white"
    }
  }
  
  
  for(i in 1:length(sturm)){
    if(is.na(sturm[i]))
      sturm[i] <- "white"
    else{
      if(sturm[i] == "IDH")
        sturm[i] <- "red3"
      else if(sturm[i] == "Mesenchymal")
        sturm[i] <- "deeppink"
      else if(sturm[i] == "RTK I 'PDGFRA'")
        sturm[i] <- "tan3"
      else if(sturm[i] == "RTK II 'Classic'")
        sturm[i] <- "darkslateblue"
      else if(sturm[i] == "G34")
        sturm[i] <- "green"
      else if(sturm[i] == "K27")
        sturm[i] <- "cyan"
      else 
        sturm[i] <- "white"
    }
  }
  
  for(i in 1:length(tumor)){
    if(is.na(tumor[i]))
      tumor[i] <- "white"
    else{ 
      if(tumor[i] == "LGG")
        tumor[i] <- "black"
      else if(tumor[i] == "GBM")
        tumor[i] <- "gray"
      else
        tumor[i] <- "khaki3"
    }
    
  }
  
  for(i in 1:length(IDH)){
    if(is.na(IDH[i]))
      IDH[i] <- "white"
    else{
      if(IDH[i] == "Mutant")
        IDH[i] <- "black"
      else if(IDH[i] == "WT")
        IDH[i] <- "grey"
      else
        IDH[i] <- "white"
      
    }
  }
  
  for(i in 1:length(p19q)){
    if(is.na(p19q[i]))
      p19q[i] <- "white"
    else{
      if(p19q[i] == "codel")
        p19q[i] <- "darkblue"
      else if(p19q[i] == "non-codel")
        p19q[i] <- "cyan"
      else
        p19q[i] <- "white"
      
    }
  }
  
  #  for(i in 1:length(idh2)){
  #   if(is.na(idh2[i]))
  #    idh2[i] <- "white"
  # else{ 
  #      if(idh2[i] == "R172G")
  #       idh2[i] <- "black"
  #    else if(idh2[i] == "R172K")
  #     idh2[i] <- "black"
  #      else if(idh2[i] == "R172M")
  #       idh2[i] <- "black"
  #    else if(idh2[i] == "R172S")
  #     idh2[i] <- "black"
  #  else if(idh2[i] == "R172W")
  #   idh2[i] <- "black"
  #      else if(idh2[i] == "WT")
  #       idh2[i] <- "gray"
  #    else 
  #     idh2[i] <- "white"
  #  }
  #}
  
  # for(i in 1:length(idh1)){
  #  if(is.na(idh1[i]))
  #   idh1[i] <- "white"
  # else{ 
  #    if(idh1[i] == "R132C")
  #     idh1[i] <- "black"
  #  else if(idh1[i] == "R132G")
  #        idh1[i] <- "black"
  #     else if(idh1[i] == "R132H")
  #      idh1[i] <- "black"
  #   else if(idh1[i] == "R132S")
  #        idh1[i] <- "black"
  #     else if(idh1[i] == "WT")
  #      idh1[i] <- "gray"
  #      else 
  #       idh1[i] <- "white"
  #  }
  #}
  
  estimates <- jet.colors(100)[ceiling((estimates-min(estimates,na.rm=T))/(max(estimates,na.rm=T)-min(estimates,na.rm=T))*(length(jet.colors(100))-1))+1]
  estimates[is.na(estimates)] <- "white"
  leukocyte <- jet.colors(100)[ceiling((leukocyte-min(leukocyte,na.rm=T))/(max(leukocyte,na.rm=T)-min(leukocyte,na.rm=T))*(length(jet.colors(100))-1))+1]
  leukocyte[is.na(leukocyte)] <- "white"
  age <- jet.colors(100)[ceiling((age-min(age,na.rm=T))/(max(age,na.rm=T)-min(age,na.rm=T))*(length(jet.colors(100))-1))+1]
  age[is.na(age)] <- "white"
  
  for(i in 1:length(rna.dna.c)){
    if(is.na(rna.dna.c[i]))
      rna.dna.c[i] <- "white"
    else{
      if(rna.dna.c[i] == "Classical")
        rna.dna.c[i] <- "red"
      else if(rna.dna.c[i] == "G-CIMP")
        rna.dna.c[i] <- "black"
      else if(rna.dna.c[i] == "Mesenchymal")
        rna.dna.c[i] <- "green"
      else if(rna.dna.c[i] == "Neural")
        rna.dna.c[i] <- "blue"
      else if(rna.dna.c[i] == "Proneural")
        rna.dna.c[i] <- "orange"
      else 
        rna.dna.c[i] <- "white"
    }
  }
  
  
  for(i in 1:length(rna.cluster)){
    if(is.na(rna.cluster[i]))
      rna.cluster[i] <- "white"
    else{
      if(rna.cluster[i] == "LGr1")
        rna.cluster[i] <- "yellow"
      else if(rna.cluster[i] == "LGr2")
        rna.cluster[i] <- "cyan"
      else if(rna.cluster[i] == "LGr3")
        rna.cluster[i] <- "blue"
      else if(rna.cluster[i] == "LGr4")
        rna.cluster[i] <- "red"
      else if(rna.cluster[i] == "unclassified")
        rna.cluster[i] <- "brown"
      else 
        rna.cluster[i] <- "white"
    }
  }
  
  
  for(i in 1:length(coc.c)){
    if(is.na(coc.c[i]))
      coc.c[i] <- "white"
    else{
      if(coc.c[i] == "coc1")
        coc.c[i] <- "darkblue"
      else if(coc.c[i] == "coc2")
        coc.c[i] <- "darkcyan"
      else if(coc.c[i] == "coc3")
        coc.c[i] <- "cyan1"
      else
        coc.c[i] <- "white"
      
    }
  }
  
  
  
  for(i in 1:length(anno)){
    if(is.na(anno[i]))
      anno[i] <- "white"
    else{
      if(anno[i] == "Case submitted is found to be a recurrence after submission")
        anno[i] <- "yellow"
      else if(anno[i] == "Case submitted is found to be a recurrence after submission/Item does not meet study protocol")
        anno[i] <- "yellow"
      else if(anno[i] == "History of unacceptable prior treatment related to a prior/other malignancy")
        anno[i] <- "yellow"
      else if(anno[i] == "Item in special subset")
        anno[i] <- "gray"
      else if(anno[i] == "Molecular analysis outside specification")
        anno[i] <- "gray"
      else if(anno[i] == "Neoadjuvant therapy")
        anno[i] <- "cyan"
      else if(anno[i] == "Neoadjuvant therapy/Case submitted is found to be a recurrence after submission")
        anno[i] <- "yellow"
      else if(anno[i] == "Normal tissue origin incorrect")
        anno[i] <- "gray"
      else if(anno[i] == "Item may not meet study protocol")
        anno[i] <- "gray"
      else if(anno[i] == "Pathology outside specification")
        anno[i] <- "gray"
      else if(anno[i] == "Prior malignancy")
        anno[i] <- "yellow"
      else if(anno[i] == "Prior malignancy/History of acceptable prior treatment related to a prior/other malignancy")
        anno[i] <- "yellow"
      else if(anno[i] == "Qualified in error")
        anno[i] <- "yellow"
      else if(anno[i] == "Qualified in error/Item in special subset")
        anno[i] <- "gray"
      else if(anno[i] == "Subject withdrew consent")
        anno[i] <- "magenta"
      
      else 
        anno[i] <- "white"
    }
  }
  
  
  
  #cluster <- cbind(age,anno,p19q,idh1,idh2,IDH,sturm,rna.dna.c,gbm.cluster.new,gbm.cluster,coc.c,lgg.cluster.new,lgg.cluster,rna.cluster,cluster,new.cl,estimates,leukocyte,tumor)
  cluster <- cbind(age,anno,p19q,IDH,sturm,rna.dna.c,gbm.cluster.new,coc.c,lgg.cluster.new,rna.cluster,cluster,estimates,leukocyte,tumor,cnv,tp53,egfr,braf,histology,grade,subtype)
  colnames(cluster) <- c("Age","Annotation","1p19q codel","IDH","Sturm et al","GBM Clusters [Brennan, 2013]","GBM DNA Meth Predicted Cluster","LGG COC [NEJM, unpublished]","LGG Predicted Clusters", "RNA Clusters","Cluster meth","RNA Leukocyte", "DNA Meth Leukocyte", "Sample Type","Cluster CNV","TP53","EGFR CNV","BRAF","Histology","Grade","Subtype")
  #colnames(cluster) <- c("Age","Annotation","1p19q codel","IDH1","IDH2", "IDH","Sturm et al","GBM Clusters [Brennan, 2013]","GBM DNA Meth Predicted Cluster","GBM DNA Meth Cl [Brennan, 2013]","COC Clusters [NEJM, unpublished]","LGG Predicted Clusters", "LGG Clusters [NEJM, unpublished]", "RNA Clusters","Cluster meth","New Master Cluster","RNA Leukocyte", "DNA Meth Leukocyte", "Sample Type")
  
  return(cluster)
}

source("/dados/ResearchProjects/thais/heatmap.plus.R")



LGG.GBM.new.noXY.s <- LGG.GBM.new.noXY[rownames(full.heatmap.orderd),5:936]
LGG.GBM.new.order <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
a <- metadata.1129.samples.20150312[substr(colnames(LGG.GBM.new.order),1,12),]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$flag.categories),as.matrix(a$cluster.cnv),as.matrix(a$mut.TP53),as.matrix(a$cnv.EGFR),as.matrix(a$mut.BRAF),as.matrix(a$histology),as.matrix(a$grade),as.matrix(a$subtype))
clab.6 <- clab.6[,c(8,1,20,19,3,17,16,4,21,14,15,10,11)]
LGG.GBM.new.order <- LGG.GBM.new.order[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
clab.6 <- clab.6[c(1:63,64:327,405:531,532:682,683:932,328:404),]
normals.sel <- norms[rownames(LGG.GBM.new.order),]
norm.label <- matrix("white", nrow = 77, ncol = 12)
norm.label[,10] <- "darkblue"
clab.6 <- rbind(norm.label,clab.6)
LGG.GBM.new.order <- cbind(normals.sel,LGG.GBM.new.order)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/normals77.png",res=800,width=2000,height=2000)
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmap_Labels.pdf", width=30, height=20,bg = "transparent")
heatmap.plus.sm(as.matrix(LGG.GBM.new.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA
                #ColSideColors = clab.6[,c(9,11,12)]
                #RowSideColors = rlab
                
)
dev.off()


### Multiplot
info <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/pdata.516lgg.606gbm.txt")
pd <- info
pd.m = subset(pd, !is.na(pd$cluster.meth))
pd.r <- subset(pd, !is.na(pd$cluster.expr) & cluster.expr != "unclassified")
pd.m$cluster.meth <- factor(pd.m$cluster.meth) 
pd.r$cluster.expr <- factor(pd.r$cluster.expr) 

library(survival)

pd.m$s = pd.m$vital == "dead"
pd.m$os = pd.m$os.months

pd.r$s = pd.r$vital == "dead"
pd.r$os = pd.r$os.months

f.m = formula(Surv(os,s) ~ cluster.meth)
f.r = formula(Surv(os,s) ~ cluster.expr)

fit.m = survfit(f.m, data=pd.m)
fit.r = survfit(f.r, data=pd.r)

surv.m <- ggsurv(fit.m, surv.col = c("green","red","purple", "orange", "yellow","blue"),ylab = "Percent Probability\nOf Survival", xlab = "Time After Diagnosis (Months)", main="DNA Methylation Clusters",back.white=T,cens.col="lightgray")


surv.r <- ggsurv(fit.r, surv.col = c("yellow","cyan","blue", "red"), ylab = "Percent Probability\nOf Survival", xlab = "Time After Diagnosis (Months)", main="RNA Clusters",back.white=T,cens.col="lightgray")


library(ggplot2)
#ages <- na.omit(ages)
a <- pd[,c("cluster.meth","case.id","tumor.type")]
a <- na.omit(a)
meanLGm1 <- apply(LGG.GBM[,substr(colnames(LGG.GBM),1,12) %in% as.character(subset(a,cluster.meth %in% "LGm1")$case.id)],2,mean,na.rm=T)
meanLGm2 <- apply(LGG.GBM[,substr(colnames(LGG.GBM),1,12) %in% as.character(subset(a,cluster.meth %in% "LGm2")$case.id)],2,mean,na.rm=T)
meanLGm3 <- apply(LGG.GBM[,substr(colnames(LGG.GBM),1,12) %in% as.character(subset(a,cluster.meth %in% "LGm3")$case.id)],2,mean,na.rm=T)
meanLGm4 <- apply(LGG.GBM[,substr(colnames(LGG.GBM),1,12) %in% as.character(subset(a,cluster.meth %in% "LGm4")$case.id)],2,mean,na.rm=T)
meanLGm5 <- apply(LGG.GBM[,substr(colnames(LGG.GBM),1,12) %in% as.character(subset(a,cluster.meth %in% "LGm5")$case.id)],2,mean,na.rm=T)
meanLGm6g <- apply(LGG.GBM[,substr(colnames(LGG.GBM),1,12) %in% as.character(subset(a,cluster.meth %in% "LGm6" & tumor.type == "GBM")$case.id)],2,mean,na.rm=T)
meanLGml <- apply(LGG.GBM[,substr(colnames(LGG.GBM),1,12) %in% as.character(subset(a,cluster.meth %in% "LGm6" & tumor.type == "LGG")$case.id)],2,mean,na.rm=T)
meanNormal <- apply(norms[rownames(norms) %in% rownames(LGG.GBM),],2,mean,na.rm=T)
mean <- data.frame(mean=c(meanLGm1,meanLGm2,meanLGm3,meanLGm4,meanLGm5,meanLGm6g,meanLGml,meanNormal),cluster=c(rep("LGm1",length(meanLGm1)),rep("LGm2",length(meanLGm2)),rep("LGm3",length(meanLGm3)),rep("LGm4",length(meanLGm4)),rep("LGm5",length(meanLGm5)),rep("LGm6-GBM",length(meanLGm6g)),rep("LGm6-LGG",length(meanLGml)),rep("normal",length(meanNormal))))
mean <- na.omit(mean)

age.m <- ggplot(mean, aes(factor(cluster), mean)) + 
  geom_boxplot(aes(fill = factor(cluster)), notchwidth=0.25) + 
  geom_jitter(height = 0, position = position_jitter(width = .1), size=2)  + 
  scale_fill_manual(values=c("green","red","purple", "orange", "yellow","blue","lightblue","gray"
  )) +
  #scale_fill_manual(values = rev(c("darkorchid","darkgreen","blue","red","DARKGREY"))) +
  ylab("Mean DNA Methylation") +
  xlab("DNA Methylation Clusters") + 
  #  facet_grid(. ~ tumor.type) +
  labs(title = "DNA Methylation Clusters by DNA Methylation", fill = "DNA Meth. Cluster")# +
#theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white"), panel.border = element_rect(colour = "white"), legend.position="none")
ggsave("age.pdf")
#ages <- na.omit(ages)
age.r <- ggplot(pd.r, aes(factor(cluster.expr), age)) + 
  geom_boxplot(aes(fill = factor(cluster.expr)), notchwidth=0.25) + 
  geom_jitter(height = 0, position = position_jitter(width = .1), size=2)  + 
  scale_fill_manual(values=c("yellow","cyan","blue", "red"
  ), labels = c("LGr1", "LGr2", "LGr3","LGr4"
  )) +
  #scale_fill_manual(values = rev(c("darkorchid","darkgreen","blue","red","DARKGREY"))) +
  ylab("Age At Diagnosis") +
  xlab("RNA Clusters") + 
  #  facet_grid(. ~ tumor.type) +
  labs(title = "RNA Clusters by Age At Diagnosis", fill = "RNA Cluster")# +



pdf("survival_age.pdf")
multiplot(surv.m,surv.r,age.m,age.r,cols=2)
dev.off()


##### p value
diff.m = survdiff(f.m, data=pd.m)
p = round(1 - pchisq(diff.m$chisq, length(diff.m$n) - 1), 4)
coxph(f.m, data=pd.m) 

diff.r = survdiff(f.r, data=pd.r)
p = round(1 - pchisq(diff.r$chisq, length(diff.r$n) - 1), 4)
coxph(f.r, data=pd.r) 




f2 = formula(Surv(os,s) ~ age)

f2.5 = formula(Surv(os,s) ~ CIMPclass + age)

model1 = coxph(f2, data=pd) ## survival by age

model2 = coxph(f2.5, data=pd) ## By age and CIMPclass

anova(model2,model1) #does adding cimpclass to age alone add a significant predictor??



f2 = formula(Surv(os,s) ~ age + codel1p19q)

f2.5 = formula(Surv(os,s) ~ CIMPclass + age + codel1p19q)

model1 = coxph(f2, data=pd) ## survival by age

model2 = coxph(f2.5, data=pd) ## By age and CIMPclass

anova(model2,model1) #does adding cimpclass to age + codel add a significant predictor??

####### Volcano LGm4+LGm5 x LGm6
#Beta Value Diff M1 x M2
M1 <- metadata.1129.samples.20150112[metadata.1129.samples.20150112$cluster.meth %in% "LGm4",]
M2 <- metadata.1129.samples.20150112[metadata.1129.samples.20150112$cluster.meth %in% "LGm6",]
#Selecting the Samples
data <- LGG.GBM.250914[,5:936]
M1 <- data[,as.character(M1$id)]
M2 <- data[,as.character(M2$id)]
volcano <- cbind(M1,M2)

values <- t(na.omit(volcano))
values <- data.frame(values)
#values.names <- colnames(values)

require(exactRankTests)
require(parallel)
system.time(w.p.values <- unlist(mclapply(values,
                                          function(probe) {
                                            zz <- wilcox.exact(probe[1:151],probe[152:228], exact=TRUE)
                                            z <- zz$p.value
                                            return(z)
                                          }, mc.cores=8)))

#lggs.p.value <- as.matrix(w.p.values)
#gbms.p.value <- as.matrix(w.p.values)
w.p.values.adj <- p.adjust(w.p.values, method = "BH")

volcano$meanM1 <- apply(volcano[,1:151],1,mean,na.rm=T)
volcano$meanM2 <- apply(volcano[,152:228],1,mean,na.rm=T)
volcano$DiffMean <- volcano$meanM1 - volcano$meanM2
volcano <- na.omit(volcano)
volcano$p.value.adj <- w.p.values.adj
volcano$p.value <- w.p.values

volcano$threshold <- "1"
a <- subset(volcano, p.value.adj < 0.05)
b <- subset(a, DiffMean < -0.5) #hyper no da direita
c <- subset(a, DiffMean < 0 & DiffMean > -0.5)
volcano[rownames(b),"threshold"] <- "2"
volcano[rownames(c),"threshold"] <- "4"
b <- subset(a, DiffMean > 0.5) #hyper no da esquerda
c <- subset(a, DiffMean > 0 & DiffMean < 0.5)
volcano[rownames(b),"threshold"] <- "3"
volcano[rownames(c),"threshold"] <- "4"
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/volcano_LGm4xLGm6_0.5.png",res=200,width=2000,height=2000)
ggplot(data=volcano, aes(x=DiffMean, y=-1*log10(p.value.adj), colour=threshold)) +
  geom_point() +
  xlim(c(-1,1)) + ylim(c(0,48)) +
  xlab("DNA Methylation Diff") + ylab("-1 * log10 of the Significance") +
  labs(title = "Volcano Plot") +
  scale_color_manual(breaks=c("1","3","4"), # color scale (for points)
                     values=c("black",  "red","grey"),
                     labels=c("Not Significant","Hypermethylated in LGm4","Not significant"),
                     name="Legend") 
dev.off()

aux <- cpg.gene
rownames(aux) <- as.character(cpg.gene$Composite.El)
a <- aux[rownames(b),]
a <- a[,1:7]
a$p.value.adj <- as.character(b$p.value.adj)
write.table(a, file="Hyper_LGm4_CutOff_at_0.5.txt",quote=F,row.names=F,sep="\t")


#### Mur validation
load("mur.DNA_Meth.Rda")
valid.clinical <- gset.mur.p[intersect(colnames(gset.mur.m),rownames(gset.mur.p)),]
valid.mur <- gset.mur.m[as.character(probe.meth),intersect(colnames(gset.mur.m),rownames(gset.mur.p))] #54 probes #46 samples

aux <- subset(valid.clinical, CIMP == "CIMP+")
cimp.p <- valid.mur[,as.character(aux$geo_accession)]
hc.cimp.p <- hclust(dist(t(cimp.p)))
aux <- subset(valid.clinical, CIMP == "CIMP-")
cimp.n <- valid.mur[,as.character(aux$geo_accession)]
hc.cimp.n <- hclust(dist(t(cimp.n)))
aux <- subset(valid.clinical, CIMP == "CD-CIMP+")
cimp.cd <- valid.mur[,as.character(aux$geo_accession)]
hc.cimp.cd <- hclust(dist(t(cimp.cd)))

mur.order <- cbind(cimp.p[,hc.cimp.p$order],cimp.n[,hc.cimp.n$order],cimp.cd[,hc.cimp.cd$order])
valid.clinical <- valid.clinical[colnames(mur.order),]
levels(valid.clinical$CIMP) <- c("pink","khaki1","lightgreen") #blue = codel #green cimp+ #red cimp-
levels(valid.clinical$IDHmutated) <- c("grey","black")

aux <- LGG.GBM.250914[,5:936]
aux <- colnames(aux[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]])
aux <- aux[c(1:63,64:327,405:531,532:682,683:932,328:404)]

rownames(all.450k) <- as.character(all.450k$Composite.Element.REF)
meanLGml <- apply(all.450k[,substr(colnames(all.450k),1,12) %in% as.character(subset(pd,cluster.meth %in% "LGm6" & tumor.type == "LGG")$case.id)],1,mean,na.rm=T)
a <- as.data.frame(meanLGml)
a$id <- rownames(a)
a<-a[with(a,order(meanLGml)),]
a <- na.omit(a)
lgg.gbm <- all.450k[rownames(a),5:649]
colnames(lgg.gbm) <- substr(colnames(lgg.gbm),1,12)
b <- aux %in% colnames(lgg.gbm)
b <- aux[b]
lgg.gbm <- lgg.gbm[,b]
clab <- pd[b,c("cluster.meth","tumor.type")]
c <- subset(clab, cluster.meth %in% "LGm6" & tumor.type == "LGG")
d <- subset(clab, cluster.meth %in% "LGm6" & tumor.type == "GBM")
clab <- subset(clab, cluster.meth != "LGm6")
clab <- rbind(clab,c)
clab <- rbind(clab,d)
levels(clab$cluster.meth) <- c("green","red","purple","orange","yellow","blue")
clab$tumor.type <- as.factor(clab$tumor.type)
levels(clab$tumor.type) <- c("gray","black")
lgg.gbm <- lgg.gbm[,430:645]
e <- lgg.gbm[,rownames(clab)[430:645]]
rlab <- data.frame(CpG=rep("white",nrow(e)),row.names=rownames(e),stringsAsFactors=F)
rlab[rownames(dat.lgg.gbm.new.noXY.dic.oecg),"CpG"] <- "black"
#clab <- c(rep("green",49),rep("red",256),rep("purple",124),rep("orange",72),rep("yellow",103),rep("blue",41))
norm <- as.data.frame(norms)
norm <- norm[rownames(norm) %in% rownames(e),]
e <- merge(e,norm,by=0,all.x=T,sort=FALSE)

rownames(e) <- as.character(e$Row.names)
e <- e[,-1]
e <- e[rownames(lgg.gbm),]

aux <- data.frame(cluster.meth=rep("white",77),tumor.type=rep("white",77))
clab <- rbind(clab,aux)

rlab <- data.frame(CpG=rep("white",nrow(e)),row.names=rownames(e),stringsAsFactors=F)
rlab[rownames(dat.lgg.gbm.new.noXY.dic.oecg),"CpG"] <- "black"
rlab <- cbind(rlab,rlab)

png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/t450.png",res=800,width=4000,height=4000)
heatmap.plus.sm(as.matrix(e),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
               ColSideColors = as.matrix(clab[430:722,]),
               RowSideColors = as.matrix(rlab),
                margins = c(1, 6)
                #RowSideColors = rlab
                
)
dev.off()

library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO
gset <- getGEO("GSE36245", GSEMatrix =TRUE)[[1]]
gset.expr <- exprs(gset)
gset.pdata <- pData(gset)

gene.location <- read.delim("/dados/ResearchProjects/thais/TCGA/Normals/mRNA_Brain/Galaxy1-[Homo_sapiens_genes_(GRCh37.p13)].tabular")
gene.location <- gene.location[gene.location$Chromosome.Name %in% c(1:22,"X","Y"),]
gene.location <- gene.location[!is.na(gene.location$EntrezGene.ID),]
gene.location <- gene.location[!duplicated(gene.location$EntrezGene.ID),] #the duplicates are different transcripts, not different
library(biomaRt)
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
affy.map <- getBM(attributes=c("affy_hg_u133_plus_2","chromosome_name","start_position","end_position","strand"),filters="affy_hg_u133_plus_2",values=rownames(gset.expr),mart=ensembl)
a <- regmatches(affy.map$strand, regexpr("^[^[:digit:]]*", affy.map$strand))
for (i in 1:length(a)){
  if(a[i] == "")
    a[i] <- "+"
}

affy.map$strandNew <- as.character(a)
Affy.genes <- GRanges(seqnames = paste("chr", as.character(affy.map$chromosome_name), sep=""), ranges = IRanges(start = affy.map$start_position, end = affy.map$end_position), strand = affy.map$strandNew,probe.id=as.character(affy.map$affy_hg_u133_plus_2)) #Transformando meu objeto em um obj do tipo GRanges
gene.GR <- GRanges(seqnames = paste0("chr",gene.location$Chromosome.Name),ranges = IRanges(start = gene.location$Gene.Start..bp., end=gene.location$Gene.End..bp.), strand = gene.location$Strand, symbol = gene.location$Associated.Gene.Name, EntrezID = gene.location$EntrezGene.ID)
###########parei aqui
distance <- as.data.frame(distanceToNearest(Affy.genes,gene.GR)) #closest gene to each 450k probe
rownames(gene.location) <- NULL
gene.order.by.distance <- gene.location[distance$subjectHits,]
gene.order.by.distance$distance <- as.matrix(distance$distance)
affy.map.nearest.gene <- cbind(affy.map, gene.order.by.distance[,c("Associated.Gene.Name","EntrezGene.ID","distance")])
sturm.expr <- merge(affy.map.nearest.gene[,c(1,7,8)],gset.expr,by.x=1,by.y=0)
c <- cpg.gene
rownames(c) <- as.character(c$Composite.El)
aux <- as.character(c[as.character(probe.meth),"Associated.Gene.Name"])
sturm.expr$Associated.Gene.Name <- as.character(sturm.expr$Associated.Gene.Name)
a <- intersect(aux,sturm.expr$Associated.Gene.Name)
sturm.hmp <- subset(sturm.expr, Associated.Gene.Name %in% as.character(a))

sturm.hmp <- sturm.hmp[with(sturm.hmp,order(Associated.Gene.Name)),]
idh <- subset(gset.pdata, characteristics_ch1.5 %in% "subgroup: IDH")
idh.s <- sturm.hmp[,as.character(idh$geo_accession)]
mean <- apply(idh.s, 1, mean, na.rm=TRUE)
g34 <- subset(gset.pdata, characteristics_ch1.5 %in% "subgroup: G34")
g34.s <- sturm.hmp[,as.character(g34$geo_accession)]
mean <- cbind(mean,apply(g34.s, 1, mean, na.rm=TRUE))
k27 <- subset(gset.pdata, characteristics_ch1.5 %in% "subgroup: K27")
k27.s <- sturm.hmp[,as.character(k27$geo_accession)]
mean <- cbind(mean,apply(k27.s, 1, mean, na.rm=TRUE))
mes <- subset(gset.pdata, characteristics_ch1.5 %in% "subgroup: Mesenchymal")
mes.s <- sturm.hmp[,as.character(mes$geo_accession)]
mean <- cbind(mean,apply(mes.s, 1, mean, na.rm=TRUE))
rtk1 <- subset(gset.pdata, characteristics_ch1.5 %in% "subgroup: RTK I")
rtk1.s <- sturm.hmp[,as.character(rtk1$geo_accession)]
mean <- cbind(mean,apply(rtk1.s, 1, mean, na.rm=TRUE))
rtk2 <- subset(gset.pdata, characteristics_ch1.5 %in% "subgroup: RTK II")
rtk2.s <- sturm.hmp[,as.character(rtk2$geo_accession)]
mean <- cbind(mean,apply(rtk2.s, 1, mean, na.rm=TRUE))
colnames(mean) <- c("IDH","G34","K27","Mesenchymal","RTKI","RTKII")
mean <- mean[,c(1,3,2,5,4,6)]
a <- metadata.1129.samples.20150112[1:6,]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$flag.categories),as.matrix(a$cluster.cnv),as.matrix(a$mut.TP53),as.matrix(a$cnv.EGFR),as.matrix(a$mut.BRAF),as.matrix(a$histology),as.matrix(a$grade),as.matrix(a$subtype))
clab.6 <- clab.6[,c(1,20,19,18,3,17,16,4,21,14,15,10,11)]
clab.6[,13] <-c("red","green3","blue","orange","yellow","brown")
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_sturmexpr2.png",res=800,width=4000,height=4000)
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_MeanRNA2.pdf")
heatmap.plus.sm(as.matrix(log2(mean)),
                col=rev(redgreen(1000)),
                #scale = "none",
                #trace = "none",
                labRow = as.character(sturm.hmp$Associated.Gene.Name),
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab.6
                #RowSideColors = rlab
                
)
dev.off()

###Survival Validation (Mur + DKFZ + Turcan)
#identical(rownames(LGG.GBM.new.order2),rownames(order))
#b <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/probesLGm2.txt",header=F)
#a <- c(names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm3),names(sort(enrichmentList.sig$LGm2)[1:15]),as.character(b$V1))
turcan.sv <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/turcan_clinical.txt")
turcan.sv$title <- paste("Tumor ID", as.character(turcan.sv$ID.number),sep=" ") #GAMBIARRA
turcan.sv <- merge(turcan.sv,turcan.clinical[,c("title","geo_accession","cluster.meth.all")],by="title")
turcan.sv$Study <- "Turcan"
rownames(turcan.sv) <- as.character(turcan.sv$geo_accession)

mur_clinical <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Mur_clinical.txt")
mur_id <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/mur_id.txt")
mur_clinical_46 <- merge(mur_id,mur_clinical,by="Case")
mur_clinical_46$Study <- "Mur"
sturm.samples <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Sturm_GBM_GSE36278.txt",skip=139)
sturm.samples <- na.omit(sturm.samples)
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/mur.DNA_Meth.Rda")
mur.46 <- gset.mur.m[,as.character(gset.mur.p$geo_accession)]

gset <- getGEO("GSE44684", GSEMatrix =TRUE)[[1]] #PA meth
pa.meth <- exprs(gset)
pa.pdata <- pData(gset)

gset <- getGEO("GSE36278", GSEMatrix =TRUE)[[1]] #Sturm meth
sturm.pdata.meth <- pData(gset)
sturm.ids <- subset(sturm.pdata.meth, characteristics_ch1 == "tissue: primary glioblastoma")
sturm.samples.nocontrol <- sturm.samples[,-c(77:82)]
rownames(sturm.samples.nocontrol) <- sturm.samples.nocontrol$ID_REF
sturm.samples.nocontrol <- sturm.samples.nocontrol[,-1]
colnames(sturm.samples.nocontrol) <- as.character(sturm.ids$geo_accession)

pa.meth.nocontrol <- pa.meth[,1:61]
sturm.pa.mur <- merge(sturm.samples.nocontrol,pa.meth.nocontrol,by=0)
rownames(sturm.pa.mur) <- as.character(sturm.pa.mur$Row.names)
sturm.pa.mur <- sturm.pa.mur[,-1]
sturm.pa.mur <- merge(sturm.pa.mur,mur.46,by=0)
rownames(sturm.pa.mur) <- as.character(sturm.pa.mur$Row.names)
sturm.pa.mur <- sturm.pa.mur[,-1]

ESg <- function(lgg.gbm, sturm.pa.mur){
  teste <- NULL
  for (i in 1:length(probes)){	
    mean.g <- mean(as.numeric(lgg.gbm[probes[i],as.character(id$id)]))
    sd.g <- sd(as.numeric(lgg.gbm[probes[i],as.character(id$id)]))
    if(is.null(teste))
      teste <- sturm.pa.mur[probes[i],] >= (mean.g - sd.g)
    else{
      aux <- sturm.pa.mur[probes[i],] >= (mean.g - sd.g)
      teste <- rbind(teste,aux)
    }
  }
  return(teste)
}

LGG.GBM.new.order2 <- LGG.GBM.250914[rownames(dat.lgg.gbm.new.noXY.dic.oecg),5:936]
#probes <- c(names(sort(enrichmentList.sig$LGm2)[1:15]),names(enrichmentList.sig$LGm3),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm5))
validation.set.noNA <- na.omit(validation.set) #136 Sturm, 61 PA, 46 Mur, 81 Turcan

b <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/probesLGm2.txt",header=F)
probes <- as.character(b$V1)
id <- subset(metadata.1129.samples.20150112,cluster.meth %in% c("LGm2","LGm3"))
#probes <- probes[probes %in% rownames(sturm.pa.mur.noNA)]
ESg1 <- ESg(LGG.GBM.250914,validation.set)

probes <- names(sort(enrichmentList.sig$LGm2)[1:15])
probes <- str_extract(probes,"*[^.]*")
id <- subset(metadata.1129.samples.20150112,cluster.meth %in% c("LGm1","LGm2","LGm3"))
#probes <- probes[probes %in% rownames(validation.set.noNA)]
ESg2 <- ESg(LGG.GBM.250914,validation.set)


probes <- names(enrichmentList.sig$LGm3)
probes <- str_extract(probes,"*[^.]*")
id <- subset(metadata.1129.samples.20150112,cluster.meth %in% c("LGm3"))
#probes <- probes[probes %in% rownames(sturm.pa.mur.noNA)]
ESg3 <- ESg(LGG.GBM.250914,validation.set)


probes <- c(names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4))
probes <- str_extract(probes,"*[^.]*")
id <- subset(metadata.1129.samples.20150112,cluster.meth %in% c("LGm4","LGm5"))
#probes <- probes[probes %in% rownames(sturm.pa.mur.noNA)]
ESg4 <- ESg(LGG.GBM.250914,validation.set)

#colSums(ESg2) > 8
#a <- as.matrix(colSums(ESg2) > 8)
#colnames(a) <- "Group"
#b <- merge(a,pa.sturm.RF.s1,by.x=0,by.y="nonTCGA.ID")

a <- as.data.frame(colSums(ESg4,na.rm=T) > 8 & colSums(ESg3,na.rm=T) < 8 & colSums(ESg2,na.rm=T) < 8)
a$type <- "ESg4"
a$id <- rownames(a)

b <- as.data.frame(colSums(ESg3,na.rm=T) > 8 & colSums(ESg2,na.rm=T) > 8 & colSums(ESg4,na.rm=T) < 8 & colSums(ESg1,na.rm=T) > 4)
b$type <- "ESg3"
b$id <- rownames(b)

c <- as.data.frame(colSums(ESg2,na.rm=T) > 8 & colSums(ESg3,na.rm=T) < 8 & colSums(ESg4,na.rm=T) < 8 & colSums(ESg1,na.rm=T) > 4)
c$type <- "ESg2"
c$id <- rownames(c)

d <- as.data.frame(colSums(ESg2,na.rm=T) > 8 & colSums(ESg1,na.rm=T) < 4)
d$type <- "ESg1"
d$id <- rownames(d)

colnames(a)[1] <- "Group"
colnames(b)[1] <- "Group"
colnames(c)[1] <- "Group"
colnames(d)[1] <- "Group"

t <- intersect(c[c$Group == "FALSE","id"],b[b$Group == "FALSE","id"])
t <- intersect(t,a[a$Group == "FALSE","id"])
t <- intersect(t,d[d$Group == "FALSE","id"])

#load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/pa.sturm.RF.s1.rda")
e <- rbind(c[c$Group == "TRUE",],b[b$Group == "TRUE",],a[a$Group == "TRUE",],d[d$Group == "TRUE",])
aux <- a[as.character(t),]
aux$type <- "Unk"
e <- rbind(e,aux)
e$type <- as.factor(e$type)
#d$type2 <- paste(d$Group,d$type,sep=".")
#d[d$id %in% as.character(t),"type2"] <- "Unk"
#d$type2 <- as.factor(d$type2)
#levels(d$type2)[c(2,4:6,8)] <- "Unk"


#### Suvival com ESg
clinical <- pa.sturm.RF.s1[,c("nonTCGA.ID","Study","Sturm.Death","Sturm.OS.months","Study.Class.original","Age.at.diagnosis.years","Gender")]
mur_cl <- mur_clinical_46[,c("geo_accession","Study","Status","OS.months","CIMP","Age","Gender")]
colnames(clinical) <- c("geo_accession","Study","Status","OS.months","CIMP","Age","Gender")
turcan.sv1 <- turcan.sv[,c("geo_accession","Study","Deceased","Follow.Up..Months.","CIMP.phenotype","Age.at.Diagnosis","Gender")]
colnames(turcan.sv1) <- c("geo_accession","Study","Status","OS.months","CIMP","Age","Gender")
clinical.St_PA_Mur_Turcan <- rbind(clinical,mur_cl)
clinical.St_PA_Mur_Turcan <- rbind(clinical.St_PA_Mur_Turcan,turcan.sv1)
clinical.St_PA_Mur_Turcan$Status <- as.factor(clinical.St_PA_Mur_Turcan$Status)
levels(clinical.St_PA_Mur_Turcan$Status) <- c("0","1","0","1",NA,"0","1") #1 dead, 0 alive, YES dead, NO alive
levels(clinical.St_PA_Mur_Turcan$Gender) <- c("female","male","female","male","female","female","male")
clinical.St_PA_Mur_Turcan$CIMP <- as.character(clinical.St_PA_Mur_Turcan$CIMP)
all <- merge(e,clinical.St_PA_Mur_Turcan,by.x="id",by.y="geo_accession")
all$Status <- as.character(all$Status)
all$Status <- as.numeric(all$Status)
all$OS.months <- as.numeric(all$OS.months)
rownames(all) <- as.character(all$id)
aux <- subset(all, OS.months > 60)
all[rownames(aux),"OS.months"] <- 60 #colocar a coluna com a data do ultimo contato com o paciente
all[rownames(aux), "Status"] <- "0" #colocar a coluna com informacao sobre o estado do paciente 

#### Survival com LGm
b <- metadata.1129.samples.20150312[,c("cluster.meth.all","id","Sturm.PA.merged.clusters","Study","gender")]
mur.cl <- mur.RF[,c("RFtype","geo_accession","CIMP","Gender")]
colnames(mur.cl) <- c("cluster.meth.all","id","Sturm.PA.merged.clusters","gender")
mur.cl$Study <- "Mur"
levels(turcan.sv$CIMP.phenotype) <- c("Turcan.CIMP+","Turcan.CIMP-")
turcan.cl <- turcan.sv[,c("cluster.meth.all","geo_accession","CIMP.phenotype","Gender")]
colnames(turcan.cl) <- c("cluster.meth.all","id","Sturm.PA.merged.clusters","gender")
turcan.cl$Study <- "Turcan"
c <- rbind(b,mur.cl)
c <- rbind(c,turcan.cl)
rownames(c) <- as.character(c$id)
b <- subset(c, Study != "TCGA")
b$gender <- as.factor(b$gender)
levels(b$gender) <- c("female","female","male","male")

b$cluster.meth.all <- as.character(b$cluster.meth.all)
b <- merge(b, all, by="id")
b$Status <- as.numeric(b$Status)
b$OS.months <- as.numeric(b$OS.months)
rownames(b) <- as.character(b$id)
aux <- subset(b, OS.months > 60)
b[rownames(aux),"OS.months"] <- 60 #colocar a coluna com a data do ultimo contato com o paciente
b[rownames(aux), "Status"] <- "0" #colocar a coluna com informacao sobre o estado do paciente 


#ggsurv(fit.m,ylab = "Percent Probability\nOf Survival", xlab = "Time After Diagnosis (Months)", main="DNA Methylation Clusters",back.white=T,cens.col="lightgray")
b$Status <- as.numeric(b$Status)
c <- subset(b, Study.x != "PA")
c$type <- as.factor(c$cluster.meth.all)
levels(c$type) <- c("LGm1","LGm2.3","LGm2.3","LGm4.5","LGm4.5","LGm6")
f.m = formula(Surv(OS.months,Status) ~ type)
fit.m = survfit(f.m, data=c)
c$Age <- as.numeric(c$Age)


age.quint <- quantile(c$Age, probs=seq(0,1,.25), na.rm=T) #colocar a coluna com a idade no momento do diagnostico
#check values and make whole number version
age.quint <- c(2,23,42,52,84) #verificar quais valores estao em age.quint
agecat <- cut(c$Age, age.quint)
fdata <- c
fdata$age2 <- agecat
fdata$sex <- fdata$Gender #colocar a coluna com o sexo do paciente
fdata$group <- as.factor(fdata$type) #colocar a coluna com o resultado do Consensus Cluster
#library(survival)
data(uspop2)
refpop <- uspop2[as.character(floor(min(age.quint)):ceiling(max(age.quint))), , "2000"]
temp <- as.numeric(cut(floor(min(age.quint)):ceiling(max(age.quint)), age.quint, include.lowest=T))
pi.us<- tapply(refpop, list(temp[row(refpop)], col(refpop)), sum)/sum(refpop)
tab2 <- with(fdata, table(age2, sex, group))/ nrow(fdata)
pi.us2 <- pi.us[,2:1] # this puts male/female in the right order (check data)
us.wt <- rep(pi.us2, length(levels(fdata$group)))/ tab2
fdata$sex <- ifelse(fdata$sex == "female", 1, 2) #convert female to 1 and male to 2
index <- with(fdata, cbind(as.numeric(age2), as.numeric(sex), as.numeric(group)))
fdata$uswt <- us.wt[index]
a <- subset(fdata, type %in% c("LGm4.5","LGm6"))
a$type <- factor(a$type)
gg.s<-coxph(Surv(as.numeric(OS.months), as.numeric(Status)) ~ type, data = a, weights = uswt)

#png(filename = "/dados/ResearchProjects/thais/TCGA/LGG.GBM/survival_validation_5y_all_LGmred.png", bg="white", res=300, width=2000, height=2000)
pdf(file = "/dados/ResearchProjects/thais/TCGA/LGG.GBM/survival_validation_5y_all_LGmred.pdf")
plot(fit.m, 
     lwd=4,
     #lty = c(2,1,2,1,1),
     col=c("green",
       "firebrick4", #colocar a cor de acordo com a ordem de 1 a 6
           "orange3",
           "blue"
           
           
     ),
     main="Kaplan-Meier Overall Survival Curves\nValidation Set\nDKFZ GBM, Mur and Turcan Samples ", 
     xlab="TIME SINCE DIAGNOSIS (MONTHS)", 
     ylab="PERCENT PROBABILITY OF SURVIVAL", 
     yscale=100,
     # xscale=365.25, 
     bg="black"
)
box(col="black", lwd=3);
#abline(h=0, lwd=4, lty=1)
legend("bottomleft", 
       legend=c(
         paste("LGm1 (n=",fit.m$n[1],")",sep=""),
         paste("LGm2.3 (n=",fit.m$n[2],")",sep=""),
         paste("LGm4.5 (n=",fit.m$n[3],")",sep=""),
         paste("LGm6 (n=",fit.m$n[4],")",sep="")
         
       ),
       #lty = c(1,3,1,3,1),
       col=c("green",
             "firebrick4",
             "orange3",
             "blue"
             
       ),
       lwd=3,
       title="Legend",
       box.lwd=3,
       bg="white")
dev.off()

ggplot(c, aes(x = factor(cluster.meth.all), fill = Sturm.PA.merged.clusters)) + geom_bar() 

load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGm_subgroups.Rda")
x <- merge(metadata.1129.samples.20150312,new.pdata,by="id")
aux <- subset(x, os.months > 60)
x[rownames(aux),"os.months"] <- 60 #colocar a coluna com a data do ultimo contato com o paciente
x[rownames(aux), "vital"] <- "0" #colocar a coluna com informacao sobre o estado do paciente 
x$s <- x$vital == "1"
x$type <- as.factor(x$new.subgroup)
x$cluster.meth <- as.factor(x$cluster.meth)
levels(x$type) <- c("LGm1.2.3","LGm1.hypo","LGm1.2.3","LGm1.2.3","LGm4.5","LGm4.5","LGm6")
levels(x$cluster.meth) <- c("LGm1","LGm2-3","LGm2-3","LGm4-5","LGm4-5","LGm6")

age.quint <- quantile(x$age, probs=seq(0,1,.25), na.rm=T) #colocar a coluna com a idade no momento do diagnostico
#check values and make whole number version
age.quint <- c(10,38,51,62,89) #verificar quais valores estao em age.quint
agecat <- cut(x$age, age.quint)
fdata <- x
fdata$age2 <- agecat
fdata$sex <- as.factor(fdata$gender) #colocar a coluna com o sexo do paciente
fdata$group <- as.factor(fdata$type) #colocar a coluna com o resultado do Consensus Cluster
#library(survival)
data(uspop2)
refpop <- uspop2[as.character(floor(min(age.quint)):ceiling(max(age.quint))), , "2000"]
temp <- as.numeric(cut(floor(min(age.quint)):ceiling(max(age.quint)), age.quint, include.lowest=T))
pi.us<- tapply(refpop, list(temp[row(refpop)], col(refpop)), sum)/sum(refpop)
tab2 <- with(fdata, table(age2, sex, group))/ nrow(fdata)
pi.us2 <- pi.us[,2:1] # this puts male/female in the right order (check data)
us.wt <- rep(pi.us2, length(levels(fdata$group)))/ tab2
fdata$sex <- ifelse(fdata$sex == "female", 1, 2) #convert female to 1 and male to 2
index <- with(fdata, cbind(as.numeric(age2), as.numeric(sex), as.numeric(group)))
fdata$uswt <- us.wt[index]
a <- subset(fdata, cluster.meth %in% c("LGm4-5","LGm6"))
a$cluster.meth <- factor(a$cluster.meth)
gg.s<-coxph(Surv(as.numeric(os.months), as.numeric(vital)) ~ cluster.meth, data = fdata, weights = uswt)

f.m = formula(Surv(os.months,s) ~ cluster.meth)
fit.lgm = survfit(f.m, data=x)



#png(filename = "/dados/ResearchProjects/thais/TCGA/LGG.GBM/survival_ES_5y2.png", bg="white", res=300, width=2000, height=2000)
pdf(file = "/dados/ResearchProjects/thais/TCGA/LGG.GBM/survival_ES_5y.pdf")
plot(fit.lgm, 
     lwd=4,
     #lty = c(2,1,2,1,1),
     col=c("green",
           "firebrick4",
            "orange3",
           "blue"
           
           
     ),
     main="Kaplan-Meier Overall Survival Curves\nTCGA Samples\nLGG and GBM ", 
     xlab="TIME SINCE DIAGNOSIS (MONTHS)", 
     ylab="PERCENT PROBABILITY OF SURVIVAL", 
     yscale=100,
     # xscale=365.25, 
     bg="black"
)
box(col="black", lwd=3);
#abline(h=0, lwd=4, lty=1)
legend("bottomleft", 
       legend=c(
         paste("LGm1 (n=",fit.lgm$n[1],")",sep=""),
         paste("LGm2-3 (n=",fit.lgm$n[2],")",sep=""),
         paste("LGm4-5 (n=",fit.lgm$n[3],")",sep=""),
         paste("LGm6 (n=",fit.lgm$n[4],")",sep="")
         
       ),
       #lty = c(1,3,1,3,1),
       col=c("green",
             "firebrick4",
             "orange3",
             "blue"
             
             
             
       ),
       lwd=3,
       title="Legend",
       box.lwd=3,
       bg="white")
dev.off()


####Age Validation
library(ggplot2)
x <- subset(metadata.1129.samples.20150112,Study=="TCGA")
x$type <- as.factor(x$cluster.meth)
levels(x$type) <- c("LGm1","LGm2-3","LGm2-3","LGm4-5","LGm4-5","LGm6")
p <- ggplot(x, aes(factor(type), age))
png("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/ES_age_LGm2.png",res=200,width=4000,height=2000)
p + geom_boxplot(aes(fill = factor(type)), notchwidth=0.25) + 
  geom_jitter(height = 0, position = position_jitter(width = .1), size=2)  + 
  scale_fill_manual(values=c("green","firebrick4","orange3", "blue"
  ), labels = c("LGm1", "LGm2-3", "LGm4-5","LGm6"
  )) +
  ylim(0,90) +
  ylab("Age At Diagnosis") +
  xlab("DNA Methylation Clusters") + 
  labs(title = "TCGA Samples\nLGG and GBM", fill = "Legend")
dev.off()

t.test(x[x$type %in% "LGm1","age"],x[x$type %in% "LGm6","age"])


b <- metadata.1129.samples.20150112[,c("cluster.meth.all","id","Sturm.PA.merged.clusters","Study")]
mur.cl <- mur.RF[,c("RFtype","geo_accession","CIMP")]
mur.cl$Study <- "Mur"
colnames(mur.cl) <- c("cluster.meth.all","id","Sturm.PA.merged.clusters","Study")
turcan.cl <- turcan.clinical[,c("cluster.meth.all","geo_accession","characteristics_ch1.2")]
turcan.cl$Study <- "Turcan"
colnames(turcan.cl) <- c("cluster.meth.all","id","Sturm.PA.merged.clusters","Study")

c <- rbind(b,mur.cl)
c <- rbind(c,turcan.cl)
rownames(c) <- as.character(c$id)
b <- subset(c, Study != "TCGA")

b$Sturm.PA.merged.clusters <- factor(b$Sturm.PA.merged.clusters)
ggplot(b, aes(x = factor(cluster.meth.all), fill = factor(Sturm.PA.merged.clusters))) + geom_bar() +  scale_fill_manual(values=c("cornflowerblue","chartreuse4","yellow","black","lightgreen","blue4","orchid","green3","gold3","darkorange2","firebrick4","red","lightblue4")) + ylab("Counts") +
  xlab("DNA Methylation Clusters") + 
  labs(title = "Distribution of Previous Labels across DNA Methylation Clusters", fill = "Legend") 
ggsave("distr_cl.png")

b$cluster.meth.all <- as.character(b$cluster.meth.all)
b <- merge(b, all, by="id")
b$Status <- as.numeric(b$Status)
b$OS.months <- as.numeric(b$OS.months)
rownames(b) <- as.character(b$id)
b$Status <- as.numeric(b$Status)
c <- subset(b, Study.x != "PA")
c$type <- as.factor(c$cluster.meth.all)
levels(c$type) <- c("LGm1","LGm2.3","LGm2.3","LGm4.5","LGm4.5","LGm6")
c$Age <- as.numeric(c$Age)
c <- subset(c, Age > 10)
p <- ggplot(c, aes(factor(type), Age))
png("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/ES_age_SPM2.png",res=200,width=4000,height=2000)
p + geom_boxplot(aes(fill = factor(type)), notchwidth=0.25) + 
  geom_jitter(height = 0, position = position_jitter(width = .1), size=2)  + 
  scale_fill_manual(values=c("green","firebrick4","orange3", "blue"
  ), labels = c("LGm1", "LGm2-3", "LGm4-5","LGm6"
  )) +
  
  ylim(0,90) +
  ylab("Age At Diagnosis") +
  xlab("Random Forest Groups") + 
  labs(title = "Validation Set\nDKFZ GBM, Mur and Turcan Samples", fill = "Legend")
dev.off()

t.test(c[c$type %in% "LGm2.3","Age"],c[c$type %in% "LGm1","Age"])


###Heatmap ES Validation
a <- c(names(enrichmentList.sig$LGm2),names(enrichmentList.sig$LGm3),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm5))
#a <- c(names(enrichmentList.rev$LGm2),names(enrichmentList.rev$LGm5))
probe.meth <- str_extract(a,"*[^.]*")
#sturm.samples <- na.omit(sturm.samples)
sturm.samples.nocontrol <- sturm.samples[,-c(77:82)]
rownames(sturm.samples.nocontrol) <- sturm.samples.nocontrol$ID_REF
colnames(sturm.samples.nocontrol)[2:137] <- as.character(rownames(RF.Sturm.on.GBM)[c(1:75,82:142)])
sturm.tcga <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Sturm_TCGA.txt")
colnames(gbm.450k) <- substr(colnames(gbm.450k),1,12)
sturm.tcga.samples <- gbm.450k[,as.character(sturm.tcga$id)]
sturm.tcga.samples$probe.id <- as.character(gbm.450k[,1])
sturm.210 <- merge(sturm.samples.nocontrol,sturm.tcga.samples,by.x="ID_REF",by.y="probe.id")
rownames(sturm.210) <- as.character(sturm.210$ID_REF)
b <- metadata.1129.samples.20150112[as.character(colnames(sturm.210)[2:211]),]
idh <- subset(b, GBM.Cluster.Sturm.et.al == "IDH")
k27 <- subset(b, GBM.Cluster.Sturm.et.al == "K27")
G34 <- subset(b, GBM.Cluster.Sturm.et.al == "G34")
pdgfra <- subset(b, GBM.Cluster.Sturm.et.al == "RTK I 'PDGFRA'")
mesen <- subset(b, GBM.Cluster.Sturm.et.al == "Mesenchymal")
classic <- subset(b, GBM.Cluster.Sturm.et.al == "RTK II 'Classic'")
order <- c(as.character(idh$id),as.character(k27$id),as.character(G34$id),as.character(pdgfra$id),as.character(mesen$id),as.character(classic$id))
b$GBM.Cluster.Sturm.et.al <- factor(b$GBM.Cluster.Sturm.et.al)
levels(b$GBM.Cluster.Sturm.et.al) <- c("blue","red","green","yellow","orange","brown")
b$cluster.meth.all <- factor(b$cluster.meth.all)
levels(b$cluster.meth.all) <- c("green","red","purple","orange","yellow","blue")
sturm.hmp <- sturm.210[as.character(probe.meth),2:211]
#png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmap_EA_Sturm_Order.png",res=800,width=4000,height=4000)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_Sturm_Order.png",res=800,width=4000,height=4000)
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmap_ES_Sturm_Labels.pdf")
heatmap.plus.sm(as.matrix(sturm.hmp[,order]),
                col = jet.colors(75),
                scale = "none",
                #labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = as.matrix(b[order,c("cluster.meth.all","GBM.Cluster.Sturm.et.al")])
                #RowSideColors = rlab
                
)
dev.off()



gset <- getGEO("GSE44684", GSEMatrix =TRUE)[[1]] #PA
pa.meth <- exprs(gset)
pa.pdata <- pData(gset)


pa.meth.nocontrol <- pa.meth[,1:61]
validation.set <- merge(sturm.samples.nocontrol[,-1],pa.meth.nocontrol,by=0)
rownames(validation.set) <- as.character(validation.set$Row.names)
b <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/probesLGm2.txt",header=F)
a <- c(names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm3),names(sort(enrichmentList.sig$LGm2)[1:15]),as.character(b$V1))
probe.meth <- str_extract(a,"*[^.]*")
validation.set <- validation.set[as.character(probe.meth),-1] #136 DKFZ, 61 PA
b <- subset(metadata.1129.samples.20150312, Study %in% "PA" | Study %in% "Sturm")

aux <- subset(b, cluster.meth.all == "LGm1")
lgm1.p <- validation.set[,as.character(aux$id)]
hc.lgm1.p <- hclust(dist(t(lgm1.p)))
aux <- subset(b, cluster.meth.all == "LGm2") 
lgm2.p <- validation.set[,as.character(aux$id)]
aux <- subset(b, cluster.meth.all == "LGm3") #No LGm3 so tem 1 amostra
lgm3.p <- as.matrix(validation.set[,as.character(aux$id)])
colnames(lgm3.p) <- as.character(aux$id)
hc.lgm2.p <- hclust(dist(t(lgm2.p)))
aux <- subset(b, cluster.meth.all == "LGm4")
lgm4.p <- validation.set[,as.character(aux$id)]
hc.lgm4.p <- hclust(dist(t(lgm4.p)))
aux <- subset(b, cluster.meth.all == "LGm5")
lgm5.p <- validation.set[,as.character(aux$id)]
hc.lgm5.p <- hclust(dist(t(lgm5.p)))
aux <- subset(b, cluster.meth.all == "LGm6")
lgm6.p <- validation.set[,as.character(aux$id)]
hc.lgm6.p <- hclust(dist(t(lgm6.p)))

order <- cbind(lgm1.p[,hc.lgm1.p$order],lgm2.p[,hc.lgm2.p$order],as.matrix(lgm3.p),lgm4.p[,hc.lgm4.p$order],lgm5.p[,hc.lgm5.p$order],lgm6.p[,hc.lgm6.p$order])
validation.info <- b[colnames(order),]
validation.info$cluster.meth.all <- as.factor(validation.info$cluster.meth.all)
validation.info$Sturm.PA.merged.clusters <- as.factor(validation.info$Sturm.PA.merged.clusters)
levels(validation.info$cluster.meth.all) <- c("green","red","purple","orange","yellow","blue")
levels(validation.info$Sturm.PA.merged.clusters) <- c("cornflowerblue","chartreuse4","blue4","brown1","green3","gold3","darkorange2","firebrick4") #PA1 PA2 G34 IDH K27 Mesenc RTKI RTKII 

png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_SturmPA.png",res=800,width=4000,height=4000)
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmap_ES_Sturm_Labels.pdf")
heatmap.plus.sm(as.matrix(order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = as.matrix(validation.info[,c("cluster.meth.all","Sturm.PA.merged.clusters")])
                #RowSideColors = rlab
                
)
dev.off()


gset <- getGEO("GSE44971", GSEMatrix =TRUE)[[1]] #PA exprs
pa.exprs <- exprs(gset)
pa.pdata.exprs <- pData(gset)

id <- subset(pa.pdata.exprs, characteristics_ch1.4 %in% c("dna methylation subgroup: 1","dna methylation subgroup: 2")) #40 de 58
pa.exprs <- pa.exprs[,as.character(id$geo_accession)]


gset <- getGEO("GSE36245", GSEMatrix =TRUE)[[1]]
gset.sturm.exprs <- exprs(gset)
gset.sturm.pdata <- pData(gset)


id <- subset(gset.sturm.pdata, characteristics_ch1.5 != "subgroup: NA") #30 de 46
sturm.exprs <- gset.sturm.exprs[,as.character(id$geo_accession)]
validation.exprs <- merge(sturm.exprs,pa.exprs,by=0)
rownames(validation.exprs) <- as.character(validation.exprs$Row.names)
validation.exprs <- validation.exprs[,-1]
#validation.set <- validation.set[as.character(probe.meth),-1] #30 DKFZ, 40 PA


gene.location <- read.delim("/dados/ResearchProjects/thais/TCGA/Normals/mRNA_Brain/Galaxy1-[Homo_sapiens_genes_(GRCh37.p13)].tabular")
gene.location <- gene.location[gene.location$Chromosome.Name %in% c(1:22,"X","Y"),]
gene.location <- gene.location[!is.na(gene.location$EntrezGene.ID),]
gene.location <- gene.location[!duplicated(gene.location$EntrezGene.ID),] #the duplicates are different transcripts, not different
library(biomaRt)
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
affy.map <- getBM(attributes=c("affy_hg_u133_plus_2","chromosome_name","start_position","end_position","strand"),filters="affy_hg_u133_plus_2",values=rownames(validation.exprs),mart=ensembl)
a <- regmatches(affy.map$strand, regexpr("^[^[:digit:]]*", affy.map$strand))
for (i in 1:length(a)){
  if(a[i] == "")
    a[i] <- "+"
}

affy.map$strandNew <- as.character(a)
Affy.genes <- GRanges(seqnames = paste("chr", as.character(affy.map$chromosome_name), sep=""), ranges = IRanges(start = affy.map$start_position, end = affy.map$end_position), strand = affy.map$strandNew,probe.id=as.character(affy.map$affy_hg_u133_plus_2)) #Transformando meu objeto em um obj do tipo GRanges
gene.GR <- GRanges(seqnames = paste0("chr",gene.location$Chromosome.Name),ranges = IRanges(start = gene.location$Gene.Start..bp., end=gene.location$Gene.End..bp.), strand = gene.location$Strand, symbol = gene.location$Associated.Gene.Name, EntrezID = gene.location$EntrezGene.ID)
distance <- as.data.frame(distanceToNearest(Affy.genes,gene.GR)) #closest gene to each 450k probe
rownames(gene.location) <- NULL
gene.order.by.distance <- gene.location[distance$subjectHits,]
gene.order.by.distance$distance <- as.matrix(distance$distance)
affy.map.nearest.gene <- cbind(affy.map, gene.order.by.distance[,c("Associated.Gene.Name","EntrezGene.ID","distance")])
exprs.data <- merge(affy.map.nearest.gene[,c(1,7,8)],validation.exprs,by.x=1,by.y=0)
c <- cpg.gene
rownames(c) <- as.character(c$Composite.El)
aux <- as.character(c[as.character(probe.meth),"Associated.Gene.Name"])
exprs.data$Associated.Gene.Name <- as.character(exprs.data$Associated.Gene.Name)
a <- intersect(aux,exprs.data$Associated.Gene.Name)
exprs.hmp <- subset(exprs.data, Associated.Gene.Name %in% as.character(a))

gset <- getGEO("GSE36278", GSEMatrix =TRUE)[[1]] #PA exprs
sturm.pdata.meth <- pData(gset)


c <- merge(pa.pdata[,1:2],pa.pdata.exprs[,1:2],by="title")
d <- merge(sturm.pdata.meth[,1:2],gset.sturm.pdata[,1:2],by="title")
c <- rbind(c,d)
b <- metadata.1129.samples.20150112[as.character(c$geo_accession.x),]

aux <- subset(b, cluster.meth.all == "LGm1")
aux <- subset(c, geo_accession.x %in% as.character(aux$id))$geo_accession.y
lgm1.p <- exprs.hmp[,as.character(aux)]
mean <- apply(lgm1.p, 1, mean, na.rm=TRUE)
aux <- subset(b, cluster.meth.all == "LGm2") 
aux <- subset(c, geo_accession.x %in% as.character(aux$id))$geo_accession.y
lgm2.p <- exprs.hmp[,as.character(aux)]
mean <- cbind(mean,apply(lgm2.p, 1, mean, na.rm=TRUE))
aux <- subset(b, cluster.meth.all == "LGm4")
aux <- subset(c, geo_accession.x %in% as.character(aux$id))$geo_accession.y
lgm4.p <- exprs.hmp[,as.character(aux)]
mean <- cbind(mean,apply(lgm4.p, 1, mean, na.rm=TRUE))
aux <- subset(b, cluster.meth.all == "LGm5")
aux <- subset(c, geo_accession.x %in% as.character(aux$id))$geo_accession.y
lgm5.p <- exprs.hmp[,as.character(aux)]
mean <- cbind(mean,apply(lgm5.p, 1, mean, na.rm=TRUE))
aux <- subset(b, cluster.meth.all == "LGm6")
aux <- subset(c, geo_accession.x %in% as.character(aux$id))$geo_accession.y
lgm6.p <- exprs.hmp[,as.character(aux)]
mean <- cbind(mean,apply(lgm6.p, 1, mean, na.rm=TRUE))
colnames(mean) <- c("LGm1","LGm2","LGm4","LGm5","LGm6")

order <- c(colnames(lgm1.p),colnames(lgm2.p),colnames(lgm4.p),colnames(lgm5.p),colnames(lgm6.p))
rownames(c) <- as.character(c$geo_accession.y)
c <- c[as.character(order),]
b <- b[as.character(c$geo_accession.x),]
b$cluster.meth.all <- as.factor(b$cluster.meth.all)
b$Sturm.PA.merged.clusters <- as.factor(b$Sturm.PA.merged.clusters)
levels(b$cluster.meth.all) <- c("green","red","orange","yellow","blue")
levels(b$Sturm.PA.merged.clusters) <- c("cornflowerblue","chartreuse4","blue4","brown1","green3","gold3","darkorange2","firebrick4") #PA1 PA2 G34 IDH K27 Mesenc RTKI RTKII 

clab <- c("green","red","orange","yellow","blue")
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_sturmPa_exprs.png",res=800,width=4000,height=4000)
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_MeanRNA2.pdf")
heatmap.plus.sm(as.matrix(log2(mean)),
                col=rev(redgreen(1000)),
                #scale = "none",
                #trace = "none",
                labRow = as.character(exprs.hmp$Associated.Gene.Name),
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = cbind(clab,clab)
                #RowSideColors = rlab
                
)
dev.off()

###Random Forest (mur,Sturm,PA,turcan)
load("Sturm_DNA_meth.rda") #450k
mur.46s.t$id <- as.character(rownames(mur.46s.t))
mur.RF <- merge(mur_clinical_46, mur.46s.t[,c("RFtype","id")],by.x="geo_accession",by.y="id")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/mur.DNA_Meth.Rda")
mur.meth <- gset.mur.m[,as.character(gset.mur.p$geo_accession)]

validation.set <- merge(sturm.meth.tumor,pa.meth.nocontrol,by=0)
rownames(validation.set) <- as.character(validation.set$Row.names)
validation.set <- validation.set[,-1]
validation.set <- merge(validation.set,mur.meth,by=0)
rownames(validation.set) <- as.character(validation.set$Row.names)
validation.set <- validation.set[,-1]
validation.set <- merge(validation.set,turcan.tumor.bv,by=0)
rownames(validation.set) <- as.character(validation.set$Row.names)
validation.set <- validation.set[,-1]
b <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/probesLGm2.txt",header=F)
aux <- subset(volcano, threshold == 3) #Hypo no LGm6
a <- c(rownames(aux),names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm3),names(sort(enrichmentList.sig$LGm2)[1:15]),as.character(b$V1))
probe.meth <- str_extract(a,"*[^.]*")
#validation.set.p <- validation.set[as.character(probe.meth),] #136 DKFZ, 61 PA, 46 Mur, 81 Turcan
validation.set.p <- validation.set[rownames(order),] #136 DKFZ, 61 PA, 46 Mur, 81 Turcan
b <- metadata.1129.samples.20150312[,c("cluster.meth.all","id","Sturm.PA.merged.clusters","histology","age","Study")]
mur.cl <- mur.RF[,c("RFtype","geo_accession","CIMP","Diagnosis","Age")]
mur.cl$Study <- "Mur"
colnames(mur.cl) <- c("cluster.meth.all","id","Sturm.PA.merged.clusters","histology","age","Study")
turcan.cl <- turcan.clinical[,c("cluster.meth.all","CIMP.phenotype","Pathology","Age.at.Diagnosis")]
turcan.cl$Study <- "Turcan"
turcan.cl$geo_accession <- rownames(turcan.cl)
turcan.cl <- turcan.cl[,c(1,6,2,3,4,5)]
colnames(turcan.cl) <- c("cluster.meth.all","id","Sturm.PA.merged.clusters","histology","age","Study")

c <- rbind(b,mur.cl)
c <- rbind(c,turcan.cl)
rownames(c) <- as.character(c$id)
b <- subset(c, Study != "TCGA")

aux <- subset(b, cluster.meth.all == "LGm1")
lgm1.p <- validation.set.p[,as.character(aux$id)]
hc.lgm1.p <- hclust(dist(t(lgm1.p)))
aux <- subset(b, cluster.meth.all == "LGm2") 
lgm2.p <- validation.set.p[,as.character(aux$id)]
hc.lgm2.p <- hclust(dist(t(lgm2.p)))
aux <- subset(b, cluster.meth.all == "LGm3") 
lgm3.p <- validation.set.p[,as.character(aux$id)]
hc.lgm3.p <- hclust(dist(t(lgm3.p)))
aux <- subset(b, cluster.meth.all == "LGm4")
lgm4.p <- validation.set.p[,as.character(aux$id)]
hc.lgm4.p <- hclust(dist(t(lgm4.p)))
aux <- subset(b, cluster.meth.all == "LGm5")
lgm5.p <- validation.set.p[,as.character(aux$id)]
hc.lgm5.p <- hclust(dist(t(lgm5.p)))
aux <- subset(b, cluster.meth.all == "LGm6")
lgm6.p <- validation.set.p[,as.character(aux$id)]
hc.lgm6.p <- hclust(dist(t(lgm6.p)))

order <- cbind(lgm1.p[,hc.lgm1.p$order],lgm2.p[,hc.lgm2.p$order],lgm3.p[,hc.lgm3.p$order],lgm4.p[,hc.lgm4.p$order],lgm5.p[,hc.lgm5.p$order],lgm6.p[,hc.lgm6.p$order])
validation.info <- b[colnames(order),]
#validation.info$RF.cluster.meth <- as.factor(validation.info$RF.cluster.meth)
#validation.info$Original.Clusters <- as.factor(validation.info$Original.Clusters)
#validation.info$Age <- as.numeric(validation.info$Age)
#levels(validation.info$RF.cluster.meth) <- c("green","red","purple","orange","yellow","blue")
#levels(validation.info$Original.Clusters) <- c("cornflowerblue","chartreuse4","yellow","black","lightgreen","blue4","orchid","green3","gold3","darkorange2","firebrick4","red","lightblue4")
validation.info[validation.info$age == "ND ","age"] <- NA
validation.info$age <- as.numeric(validation.info$age)
validation.info[validation.info$age <= 21 & !is.na(validation.info$age) & validation.info$Study != "PA","histology"] <- "Ped.glioma"
#validation.info <- b[colnames(order),]
validation.info$cluster.meth.all <- as.factor(validation.info$cluster.meth.all)
validation.info$histology <- as.factor(validation.info$histology)
levels(validation.info$histology) <- c("GBM","OligoA","OligoA","OligoD","other","Ped.glioma","Astro","other","other","GBM","OligoA","OligoD","other","PA")
validation.info$Sturm.PA.merged.clusters <- as.factor(validation.info$Sturm.PA.merged.clusters)
levels(validation.info$cluster.meth.all) <- c("green","red","purple","orange","yellow","blue")
levels(validation.info$Sturm.PA.merged.clusters) <- c("cornflowerblue","chartreuse4","yellow","black","lightgreen","blue4","orchid","green3","gold3","darkorange2","firebrick4","red","lightblue4") #PA1 PA2 CD-CIMP+ CIMP+ CIMP- G34 IDH K27 Mesenc RTKI RTKII CIMP+ CIMP-
levels(validation.info$histology) <- c("purple","cyan","green","grey","darkblue","red","orange")


####outra coisa
a <- merge(validation.info,all,by="id")
a$type <- as.factor(a$type) 
levels(a$type) <- c("green","aquamarine3","purple","orange3","blue")
rownames(a) <- as.character(a$id)
a <- a[colnames(order),]
c <- cpg.gene
genes <- as.character(c[as.character(probe.meth),"Associated.Gene.Name"])
names <- paste(rownames(LGG.GBM.new.order2),genes,sep=".")
#png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_SturmPAMurTurcan5.png",res=1600,width=8000,height=8000)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_SturmPAMurTurcan5.pdf")
heatmap.plus.sm(as.matrix(order[1:2,]),
                col = jet.colors(75),
                scale = "none",
                #labRow = NA,
                #cexRow = 0.4,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = as.matrix(validation.info[,c("histology","Sturm.PA.merged.clusters","cluster.meth.all")])
                #RowSideColors = rlab
                
)
dev.off()


gset.turcan <- getGEO("GSE30339", GSEMatrix =TRUE)[[1]]
gset.turcan.m <- exprs(gset.turcan)
gset.turcan.p <- pData(gset.turcan)

id <- subset(gset.turcan.p,characteristics_ch1 == "tissue type: frozen tissue for primary tumor")
turcan.tumor <- gset.turcan.m[,as.character(id$geo_accession)]

turcan.tumor.bv <- apply(turcan.tumor, c(1,2), function(x){y <- exp(x)/(1+exp(x)); return(y)})

id <- turcan.sv[rownames(turcan.1281.t.bv),]
turcan.clinical <- cbind(id,cluster.meth.all=as.matrix(lgg_gbm.RF.to.turcan.pred))

b <- metadata.1129.samples.20150112[,c("cluster.meth.all","id","Sturm.PA.merged.clusters","Study")]
mur.cl <- mur.RF[,c("RFtype","geo_accession","CIMP")]
mur.cl$Study <- "Mur"
colnames(mur.cl) <- c("cluster.meth.all","id","Sturm.PA.merged.clusters","Study")
turcan.cl <- turcan.clinical[,c("cluster.meth.all","geo_accession","characteristics_ch1.2")]
turcan.cl$Study <- "Turcan"
colnames(turcan.cl) <- c("cluster.meth.all","id","Sturm.PA.merged.clusters","Study")

c <- rbind(b,mur.cl)
c <- rbind(c,turcan.cl)
#b[b$cluster.meth.all]
rownames(c) <- as.character(c$id)
b <- subset(c, Study != "TCGA")

aux <- subset(b, cluster.meth.all == "LGm1")
lgm1.p <- validation.set.p[,as.character(aux$id)]
hc.lgm1.p <- hclust(dist(t(lgm1.p)))
aux <- subset(b, cluster.meth.all == "LGm2") 
lgm2.p <- validation.set.p[,as.character(aux$id)]
hc.lgm2.p <- hclust(dist(t(lgm2.p)))
aux <- subset(b, cluster.meth.all == "LGm3") 
lgm3.p <- validation.set.p[,as.character(aux$id)]
hc.lgm3.p <- hclust(dist(t(lgm3.p)))
aux <- subset(b, cluster.meth.all == "LGm4")
lgm4.p <- validation.set.p[,as.character(aux$id)]
hc.lgm4.p <- hclust(dist(t(lgm4.p)))
aux <- subset(b, cluster.meth.all == "LGm5")
lgm5.p <- validation.set.p[,as.character(aux$id)]
hc.lgm5.p <- hclust(dist(t(lgm5.p)))
aux <- subset(b, cluster.meth.all == "LGm6")
lgm6.p <- validation.set.p[,as.character(aux$id)]
hc.lgm6.p <- hclust(dist(t(lgm6.p)))

order <- cbind(lgm1.p[,hc.lgm1.p$order],lgm2.p[,hc.lgm2.p$order],lgm3.p[,hc.lgm3.p$order],lgm4.p[,hc.lgm4.p$order],lgm5.p[,hc.lgm5.p$order],lgm6.p[,hc.lgm6.p$order])
validation.info <- b[colnames(order),]
validation.info$cluster.meth.all <- as.factor(validation.info$cluster.meth.all)
validation.info$Sturm.PA.merged.clusters <- as.factor(validation.info$Sturm.PA.merged.clusters)
levels(validation.info$cluster.meth.all) <- c("green","red","purple","orange","yellow","blue")
levels(validation.info$Sturm.PA.merged.clusters) <- c("cornflowerblue","chartreuse4","lightcoral","khaki1","lightgreen","blue4","brown1","green3","gold3","darkorange2","firebrick4") #PA1 PA2 G34 IDH K27 Mesenc RTKI RTKII 

png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_SturmPAMur.png",res=800,width=4000,height=4000)
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmap_ES_Sturm_Labels.pdf")
heatmap.plus.sm(as.matrix(order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = as.matrix(validation.info[,c("cluster.meth.all","Sturm.PA.merged.clusters")])
                #RowSideColors = rlab
                
)
dev.off()


########### LGm6
aux <- metadata.1129.samples.20150312[metadata.1129.samples.20150312$cluster.meth %in% c("LGm1"),]
M1 <- sample(aux$id,10)
aux <- metadata.1129.samples.20150312[metadata.1129.samples.20150312$cluster.meth %in% c("LGm2"),]
M2 <- sample(aux$id,43)
aux <- metadata.1129.samples.20150312[metadata.1129.samples.20150312$cluster.meth %in% c("LGm3"),]
M3 <- sample(aux$id,21)
aux <- metadata.1129.samples.20150312[metadata.1129.samples.20150312$cluster.meth %in% c("LGm4"),]
M4 <- sample(aux$id,25)
aux <- metadata.1129.samples.20150312[metadata.1129.samples.20150312$cluster.meth %in% c("LGm5"),]
M5 <- sample(aux$id,41)
id <- c(M1,M2,M3,M4,M5)

M6 <- metadata.1129.samples.20150312[metadata.1129.samples.20150312$cluster.meth %in% "LGm6",]

data <- LGG.GBM.250914[,5:936]
M1 <- data[,as.character(id)]
M2 <- data[,as.character(M6$id)]
volcano <- cbind(M1,M2)

values <- t(na.omit(volcano))
values <- data.frame(values)
require(exactRankTests)
require(parallel)
w.p.values <- unlist(mclapply(values,
                              function(probe) {
                                zz <- wilcox.exact(as.matrix(probe[1:140]),as.matrix(probe[141:217]),exact=TRUE,alternative = "two.sided") 
                                z <- zz$p.value
                                return(z)
                              }, mc.cores=8))

w.p.values.adj <- p.adjust(w.p.values, method = "BH")


volcano$meanM1 <- apply(volcano[,1:140],1,mean,na.rm=T)
volcano$meanM2 <- apply(volcano[,141:217],1,mean,na.rm=T)
volcano$DiffMean <- volcano$meanM1 - volcano$meanM2
volcano <- na.omit(volcano)
volcano$p.value.adj <- w.p.values.adj
volcano$p.value <- w.p.values

volcano$threshold <- "1"
a <- subset(volcano, p.value.adj < 10^-21)
b <- subset(a, DiffMean < -0.35) #hyper no da direita
c <- subset(a, DiffMean < 0 & DiffMean > -0.30)
volcano[rownames(b),"threshold"] <- "2"
volcano[rownames(c),"threshold"] <- "4"
b <- subset(a, DiffMean > 0.33) #hyper no da esquerda
c <- subset(a, DiffMean > 0 & DiffMean < 0.33)
volcano[rownames(b),"threshold"] <- "3"
volcano[rownames(c),"threshold"] <- "4"
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/volcano_LGm6xAll_0.5.png",res=200,width=2000,height=2000)
ggplot(data=volcano, aes(x=DiffMean, y=-1*log10(p.value.adj), colour=threshold)) +
  geom_point() +
  xlim(c(-1,1)) + ylim(c(0,48)) +
  xlab("DNA Methylation Diff") + ylab("-1 * log10 of the Significance") +
  labs(title = "Volcano Plot") +
  scale_color_manual(breaks=c("1","3","4"), # color scale (for points)
                     values=c("black",  "red","grey"),
                     labels=c("Not Significant","Hypomethylated in LGm6","Not significant"),
                     name="Legend") 
dev.off()

aux <- subset(volcano, threshold == 3) #Hypo no LGm6
c <- cpg.gene
rownames(c) <- as.character(c$Composite.El)
LGG.GBM.new.order <- c[as.character(probe.meth),644:1279]
colnames(LGG.GBM.new.order) <- substr(colnames(LGG.GBM.new.order),1,12)
a <- metadata.1129.samples.20150312[substr(colnames(LGG.GBM.new.order),1,12),]
id6 <- a[a$cluster.meth %in% c("LGm6"),"id"]
id.others <- a[a$cluster.meth %in% c("LGm1","LGm2","LGm3","LGm4","LGm5"),"id"]
sd6 <- apply(LGG.GBM.new.order[,as.character(id6)],1,sd,na.rm=T)
sd.others <- apply(LGG.GBM.new.order[,as.character(id.others)],1,sd,na.rm=T)

lgg.gbm <- LGG.GBM.250914[rownames(aux),5:936]
lgg.gbm <- lgg.gbm[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
colnames(lgg.gbm) <- substr(colnames(lgg.gbm),1,12)


a <- metadata.1129.samples.20150312[substr(colnames(lgg.gbm),1,12),]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$flag.categories),as.matrix(a$cluster.cnv),as.matrix(a$mut.TP53),as.matrix(a$cnv.EGFR),as.matrix(a$mut.BRAF),as.matrix(a$histology),as.matrix(a$grade),as.matrix(a$subtype))
clab.6 <- clab.6[,c(1,20,19,18,3,17,16,4,21,14,15,10,11)]
lgg.gbm <- lgg.gbm[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
lgg.gbm <- lgg.gbm[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
clab.6 <- clab.6[c(1:63,64:327,405:531,532:682,683:932,328:404),]

png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap.png",res=800,width=4000,height=4000)
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmap_ES_Sturm_Labels.pdf")
heatmap.plus.sm(as.matrix(lgg.gbm),
                col = jet.colors(75),
                scale = "none",
                #labRow = NA,
                labCol = NA,
                Colv = NA,
               # Rowv = NA
                ColSideColors = as.matrix(clab.6[,12:13])
                #RowSideColors = rlab
                
)
dev.off()

## Supplemental - Heatmap ALL
b <- metadata.1129.samples.20150312[,c("cluster.meth.all","id","Sturm.PA.merged.clusters","Study")]
mur.cl <- mur.RF[,c("RFtype","geo_accession","CIMP")]
mur.cl$Study <- "Mur"
colnames(mur.cl) <- c("cluster.meth.all","id","Sturm.PA.merged.clusters","Study")
turcan.cl <- turcan.clinical[,c("cluster.meth.all","geo_accession","characteristics_ch1.2")]
turcan.cl$Study <- "Turcan"
colnames(turcan.cl) <- c("cluster.meth.all","id","Sturm.PA.merged.clusters","Study")
c <- rbind(b,mur.cl)
c <- rbind(c,turcan.cl)
rownames(c) <- as.character(c$id)

LGG.GBM.new.noXY.s <- LGG.GBM.new.noXY[rownames(full.heatmap.orderd) ,5:936]
all.1256 <- merge(LGG.GBM.new.noXY.s,validation.set,by=0,all.x=T)
rownames(all.1256) <- as.character(all.1256$Row.names)
all.1256 <- all.1256[,-1]

aux <- subset(c, cluster.meth.all == "LGm1")
lgm1 <- all.1256[,as.character(aux$id)]
hc.lgm1 <- hclust(dist(t(lgm1)))
aux <- subset(c, cluster.meth.all == "LGm2")
lgm2 <- all.1256[,as.character(aux$id)]
hc.lgm2 <- hclust(dist(t(lgm2)))
aux <- subset(c, cluster.meth.all == "LGm3")
lgm3 <- all.1256[,as.character(aux$id)]
hc.lgm3 <- hclust(dist(t(lgm3)))
aux <- subset(c, cluster.meth.all == "LGm4")
lgm4 <- all.1256[,as.character(aux$id)]
hc.lgm4 <- hclust(dist(t(lgm4)))
aux <- subset(c, cluster.meth.all == "LGm5")
lgm5 <- all.1256[,as.character(aux$id)]
hc.lgm5 <- hclust(dist(t(lgm5)))
aux <- subset(c, cluster.meth.all == "LGm6")
lgm6 <- all.1256[,as.character(aux$id)]
hc.lgm6 <- hclust(dist(t(lgm6)))

all.order <- cbind(lgm1[,hc.lgm1$order],lgm2[,hc.lgm2$order],lgm3[,hc.lgm3$order],lgm4[,hc.lgm4$order],lgm5[,hc.lgm5$order],lgm6[,hc.lgm6$order])
info <- c[colnames(all.order),]
info$cluster.meth.all <- as.factor(info$cluster.meth.all)
info$Sturm.PA.merged.clusters <- as.factor(info$Sturm.PA.merged.clusters)
info$Study <- as.factor(info$Study)
levels(info$cluster.meth.all) <- c("green","red","purple","orange","yellow","blue")
#levels(info$Sturm.PA.merged.clusters) <- c("cornflowerblue","chartreuse4","lightgoldenrod","black","lightgreen","blue4","orchid","green3","green","red","purple","orange","yellow","blue","gold3","darkorange2","firebrick4","red","lightblue4") #PA1 PA2 CD-CIMP+ CIMP+ CIMP- G34 IDH K27 Mesenc RTKI RTKII C
levels(info$Sturm.PA.merged.clusters) <- c("cornflowerblue","chartreuse4","lightgoldenrod","black","lightgreen","blue4","orchid","green3","lightcyan4","lightcyan4","lightcyan4","lightcyan4","lightcyan4","lightcyan4","gold3","darkorange2","firebrick4","red","lightblue4")
levels(info$Study) <- c("magenta","lightpink4","salmon","turquoise","tomato4")

all.order <- all.order[rownames(full.heatmap.orderd),]
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_All.png",res=1600,width=8000,height=8000)
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_SturmPAMurTurcan4.pdf")
heatmap.plus.sm(as.matrix(all.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                #cexRow = 0.4,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = as.matrix(info[,c("Study","Sturm.PA.merged.clusters","cluster.meth.all")])
                #RowSideColors = rlab
                
)
dev.off()


##### LGm1 Hypogroup
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/ids.newdiscovery.450.rda")
ids.newdiscovery.450 <- as.data.frame(ids.newdiscovery.450)
ids.newdiscovery.450$newID <- substr(ids.newdiscovery.450$ids.newdiscovery.450,1,12)
aux <- metadata.1129.samples.20150312[as.character(ids.newdiscovery.450$newID),c("tumor.type","Methylation.Platform")] #5 GBMs and 245 LGGs
info <- merge(ids.newdiscovery.450,aux,by.x="newID",by.y=0,all.x=T)
aux <- subset(info, tumor.type == "GBM")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/gbm450.Rda")
gbms <- gbm.450k[,as.character(aux$ids.newdiscovery.450)]
aux <- subset(info, tumor.type == "LGG")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/lgg.Rda")
lggs <- lgg.450k[,as.character(aux$ids.newdiscovery.450)]
LGG.GBM.newdiscovey.450 <- cbind(gbms,lggs)
LGG.GBM.newdiscovey.450 <- cbind(lgg.450k[,1:4],LGG.GBM.newdiscovey.450)
rownames(LGG.GBM.newdiscovey.450) <- as.character(LGG.GBM.newdiscovey.450$Composite.Element.REF)
save(LGG.GBM.newdiscovey.450,file="LGG.GBM.newdiscovey.450.rda")

write.table(d[,c(2,6)],file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Todos_ES.txt",quote=F,row.names=F,sep="\t") #de 1 a 11

aux <- as.data.frame(sort(c(enrichmentList.sig$LGm5,enrichmentList.sig$LGm4)))
aux$ProbeID <- str_extract(as.character(rownames(aux)),"*[^.]*")
aux$Gene <- str_extract(as.character(rownames(aux)),"[^.]*$")

write.table(aux[,c(2,1)],file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/EReg52.txt",quote=F,row.names=F,sep="\t") #de 1 a 11


######### ROC Curves (RF)
library(ROCR)
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/PA.onLGG.GBM.Rda") #pa.probes.pa.prob
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Sturm.onLGG.GBM.Rda") #sturm.probes.sturm.prob
#lgg_gbm.RF.to.turcan.prob
#lgg_gbm.RF.to.mur.prob

pr <- predict(lgg_gbm.RF.to.mur, myTest, type="prob")[,6]
type <- factor(ifelse(myTest$RFtype == "LGm6", "LGm6", "notLGm6"))
pred = prediction(pr, type)
t1 = performance(pred,"tpr","fpr")
t2 = performance(pred,"tpr","fpr")
t3 = performance(pred,"tpr","fpr")
t4 = performance(pred,"tpr","fpr")
t5 = performance(pred,"tpr","fpr")
t6 = performance(pred,"tpr","fpr")
plot(t6,main="ROC Curve for Random Forest",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")


#aux <- validation.info[rownames(sturm.probes.sturm.prob),"cluster.meth.all"]
#aux <- validation.info[rownames(pa.probes.pa.prob),"cluster.meth.all"]
#aux <- validation.info[rownames(lgg_gbm.RF.to.mur.prob),"cluster.meth.all"]
aux <- metadata.1129.samples.20150312[substr(rownames(myTest),1,12),"cluster.meth.all"]
type1 <- factor(ifelse(aux == "LGm1", "LGm1", "notLGm1"))
type2 <- factor(ifelse(aux == "LGm2", "LGm2", "notLGm2"))
type3 <- factor(ifelse(aux == "LGm3", "LGm3", "notLGm3"))
type4 <- factor(ifelse(aux == "LGm4", "LGm4", "notLGm4"))
type5 <- factor(ifelse(aux == "LGm5", "LGm5", "notLGm5"))
type6 <- factor(ifelse(aux == "LGm6", "LGm6", "notLGm6"))

t1 <- roc(response = type1,
          #predictor = sturm.probes.sturm.prob[, 1],
          #predictor = pa.probes.pa.prob[, 1],
          predictor = sturm.probes.lgg_gbm.pred[, 1],
          levels = rev(levels(type1)))
t2 <- roc(response = type2,
          #predictor = sturm.probes.sturm.prob[, 2],
          #predictor = pa.probes.pa.prob[, 2],
          predictor = sturm.probes.lgg_gbm.pred[, 2],
          levels = rev(levels(type2)))
t3 <- roc(response = type3,
          #predictor = sturm.probes.sturm.prob[, 3],
          #predictor = pa.probes.pa.prob[, 3],
          predictor = sturm.probes.lgg_gbm.pred[, 3],
          levels = rev(levels(type3)))
t4 <- roc(response = type4,
          #predictor = sturm.probes.sturm.prob[, 4],
          #predictor = pa.probes.pa.prob[, 4],
          predictor = sturm.probes.lgg_gbm.pred[, 4],
          levels = rev(levels(type4)))
t5 <- roc(response = type5,
          #predictor = sturm.probes.sturm.prob[, 5],
          #predictor = pa.probes.pa.prob[, 5],
          predictor = sturm.probes.lgg_gbm.pred[, 5],
          levels = rev(levels(type5)))
t6 <- roc(response = type6,
          #predictor = sturm.probes.sturm.prob[, 6],
          #predictor = pa.probes.pa.prob[, 6],
          predictor = sturm.probes.lgg_gbm.pred[, 6],
          levels = rev(levels(type6)))

data.roc <- data.frame(sens=t1$sensitivities,spec=(1-t1$specificities),group="LGm1")
aux <- data.frame(sens=t2$sensitivities,spec=(1-t2$specificities),group="LGm2")
data.roc <- rbind(data.roc,aux)
aux <- data.frame(sens=t3$sensitivities,spec=(1-t3$specificities),group="LGm3")
data.roc <- rbind(data.roc,aux)
aux <- data.frame(sens=t4$sensitivities,spec=(1-t4$specificities),group="LGm4")
data.roc <- rbind(data.roc,aux)
aux <- data.frame(sens=t5$sensitivities,spec=(1-t5$specificities),group="LGm5")
data.roc <- rbind(data.roc,aux)
aux <- data.frame(sens=t6$sensitivities,spec=(1-t6$specificities),group="LGm6")
data.roc <- rbind(data.roc,aux)

ggplot(data.roc,aes(spec,sens,color=group)) +
  geom_line(size = 1, alpha = 0.7)+
  scale_colour_manual(values = c("LGm1" = "green","LGm2" = "red","LGm3" = "purple", "LGm4" = "orange", "LGm5" = "yellow", "LGm6" = "blue")) +
  #scale_colour_manual("green","LGm2" = "red","LGm3" = "purple", "LGm4" = "orange", "LGm5" = "yellow", "LGm6" = "blue")) +
  labs(title= "ROC curve\nSturm et al., 2012", 
       x = "False Positive Rate (1-Specificity)", 
       y = "True Positive Rate (Sensitivity)",
       color = "Legend")


plot(t1, legacy.axes = T)
lines(t2$sensitivities ~ t2$specificities, col = 2)
lines(t3$sensitivities ~ t3$specificities, col = 3)
lines(t4$sensitivities ~ t4$specificities, col = 4)
lines(t5$sensitivities ~ t5$specificities, col = 5)
lines(t6$sensitivities ~ t6$specificities, col = 6)


####### LGGs 450k NEJM probes
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/lgg.Rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_NEJM_11k_Probes.Rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/metadata_LGG_GBM_PA_Sturm-20150312.Rda")
rownames(lgg.450k) <- as.character(lgg.450k$Composite.Element.REF)
meta.lgg <- metadata.1129.samples.20150312[rownames(metadata.1129.samples.20150312) %in% substr(colnames(lgg.450k),1,12),]
colnames(lgg.450k)[5:520] <- substr(colnames(lgg.450k)[5:520],1,12)
lgg.450k <- lgg.450k[as.character(LGG.puritycorrection.HeatmapOrder.rows.ID),]
lgg.450k <- na.omit(lgg.450k) #37 probes were excluded

aux <- subset(meta.lgg, cluster.meth == "LGm1")
lgm1 <- lgg.450k[,as.character(aux$id)]
hc.lgm1 <- hclust(dist(t(lgm1)))
aux <- subset(meta.lgg, cluster.meth == "LGm2")
lgm2 <- lgg.450k[,as.character(aux$id)]
hc.lgm2 <- hclust(dist(t(lgm2)))
aux <- subset(meta.lgg, cluster.meth == "LGm3")
lgm3 <- lgg.450k[,as.character(aux$id)]
hc.lgm3 <- hclust(dist(t(lgm3)))
aux <- subset(meta.lgg, cluster.meth == "LGm4")
lgm4 <- lgg.450k[,as.character(aux$id)]
hc.lgm4 <- hclust(dist(t(lgm4)))
aux <- subset(meta.lgg, cluster.meth == "LGm5")
lgm5 <- lgg.450k[,as.character(aux$id)]
hc.lgm5 <- hclust(dist(t(lgm5)))
aux <- subset(meta.lgg, cluster.meth == "LGm6")
lgm6 <- lgg.450k[,as.character(aux$id)]
hc.lgm6 <- hclust(dist(t(lgm6)))

lgg.450k <- lgg.450k[,5:520]

lgg.order <- cbind(lgm1[,hc.lgm1$order],lgm2[,hc.lgm2$order],lgm3[,hc.lgm3$order],lgm4[,hc.lgm4$order],lgm5[,hc.lgm5$order],lgm6[,hc.lgm6$order])
meta.order <- meta.lgg[colnames(lgg.order),]
levels(meta.order$cluster.meth) <- c("green","red","purple","orange","yellow","blue") 
levels(meta.order$cluster.meth.lgg) <- c("black","red","blue","forestgreen","darkorchid4") 
meta.order$subtype <- factor(meta.order$subtype)
levels(meta.order$subtype) <- c("cyan","tomato","gold")

png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_LGG_NEJM.png",res=800,width=4000,height=4000)
heatmap.plus.sm(as.matrix(lgg.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = as.matrix(meta.order[,c("cluster.meth","cluster.meth.lgg","subtype")]),
                margins = c(1, 6)
                #RowSideColors = rlab
                
)
dev.off()

aux <- subset(metadata.1129.samples.20150312, cluster.meth == "LGm6" & IDH.status == "WT")
a <- read.delim("TERTexpr_grps.txt")
a <- merge(aux,a,by.x=0,by.y="id",all.x=TRUE)

load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/pdata.panglioma.Rda")
aux <- merge(metadata.1129.samples.20150312,pdata.panglioma,by="id")
aux <- subset(aux, new.anno %in% c("IDH wt","IDH wt good"))

### Mutation - NADA (Quase - NF1)
j <- NULL  
for(i in 64:200){
j <- c(j,fisher.test(table(aux$new.anno, aux[,i]))$p.value)
}

### Amp/Del - Neutral x Non-neutral
j <- NULL  
for(i in 286:328){
  levels(aux[,i])[levels(aux[,i]) != "Neutral"] <- "Del"
  j <- c(j,fisher.test(table(aux$new.anno, aux[,i]))$p.value)
}

### Colunas -> 315,316,317.. 326 e 327 (cromossomo 13q e 19q)

### CNV - Neutral x Non-neutral
j <- NULL  
for(i in 201:256){
  levels(aux[,i])[levels(aux[,i]) != "Neutral"] <- "Non-neutral"
  j <- c(j,fisher.test(table(aux$cluster.meth.lgg, aux[,i]))$p.value)
}

### Colunas -> 205,232,233, 242 e 243 (cromossomo 13q e 19q)

lm1 <- subset(aux, cluster.meth.lgg == "M1")
lm1 <- lgg.450k[,as.character(lm1$id)]
lm1.m <- apply(lm1,2,mean,na.rm=T)

lm4 <- subset(aux, cluster.meth.lgg == "M4")
lm4 <- lgg.450k[,as.character(lm4$id)]
lm4.m <- apply(lm4,2,mean,na.rm=T)


library(ggplot2)
x <- data.frame(mean = c(lm1.m,lm4.m))
x$type <- NA
x$type[1:13] <- "Lm1"
x$type[14:26] <- "Lm4"

#colnames(x) <- c("Cerebellum", "Lm1","PA","LGm1+LGm2+LGm3+LGm4+LGm5")
png("/dados/ResearchProjects/thais/TCGA/LGG.GBM/meanM1xM4_GW.png")
p <- ggplot(x, aes(factor(type), mean))
p + geom_boxplot(aes(fill = factor(type)), notchwidth=0.25) +
  geom_jitter(height = 0, position = position_jitter(width = .1), size=2)  +
  scale_fill_manual(values=c("black","darkgreen"
  ), labels = c("Lm1", "Lm4"
  )) +
  ylab("Mean Methylation Values") +
  xlab("Lm Groups in LGm6") +  ylim(c(0,0.6)) +
  labs(title = "Mean Methylation Values (Genome-Wide) By Groups", fill = "Groups") #+
  #theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white"), panel.border = element_rect(colour = "white"), legend.position="none")
dev.off()

y <- x #GW


### camila

CodingVariants()

load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Annotation.rda")
a <- as.data.frame(loc_hg19,row.names = NULL)
b <- as.data.frame(enhancers.cg_GR)
c <- merge(a,b,by=c("seqnames","start","end"),all.x=T)
c$cgID <- as.factor(c$cgID)
ts <- function(c){
  for(i in 1:nlevels(c$cgID)){
    aux <- subset(c,cgID %in% levels(c$cgID)[i])
    
  }
}


c <- read.delim("LGG-GBM.merged.data.txt")
a <- subset(metadata.1129.samples.20150312, cluster.meth %in% c("LGm4","LGm5","LGm6") & IDH.status == "WT")
a$cluster.meth <- factor(a$cluster.meth)
a$IDH.status <- factor(a$IDH.status)
p.m <- ggplot(a, aes(factor(cluster.meth), ESTIMATEScore)) +
geom_boxplot(aes(fill = factor(cluster.meth)), notchwidth=0.25) +
geom_jitter(height = 0, position = position_jitter(width = .1), size=2)  +
scale_fill_manual(values=c("orange", "yellow","blue"
), labels = c("LGm4", "LGm5", "LGm6"
)) +
#scale_fill_manual(values = rev(c("darkorchid","darkgreen","blue","red","DARKGREY"))) +
ylab("Purity Estimates Scores") +
xlab("DNA Methylation Clusters") +
facet_grid(. ~ tumor.type) +
labs(title = "DNA Methylation Clusters by Purity Scores", fill = "DNA Meth. Cluster")
p.m
ggsave(filename = "PurityEst-IDHwt.pdf")

b <- subset(a,cluster.meth %in% "LGm6")
c <- subset(b, tumor.type %in% "LGG")
d <- subset(b, tumor.type %in% "GBM")

clab.6 <- pdata[rownames(oGEz),c("tumor.type","histology","subtype","cluster.meth","cluster.expr")]
levels(clab.6$cluster.meth) <- c("green","red","purple","orange","yellow","blue")
levels(clab.6$subtype) <- c("red","black","cyan","tomato","gold","green","blue","orange")
levels(clab.6$histology) <- c("red","purple","cyan","green3")
levels(clab.6$tumor.type) <- c("grey","black")
levels(clab.6$cluster.expr) <- c("yellow","cyan","blue","red","white")
#png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmap_Michele.png",res=1600,width=8000,height=8000)
pdf(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmap_MicheleLABELS.pdf")
heatmap.plus.sm(as.matrix(t(oGEz[,1:5])),
                col=greenred(1000),
                scale = "none",
                #trace = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                #Rowv = NA,
                ColSideColors = as.matrix(clab.6)
                #RowSideColors = rlab
                
)
dev.off()

################# REGEX NOVO
#"[^.]*$" - Gene -> Depois do ponto, sem inclui-lo
#"^*[^.]*" - Probe -> Antes do ponto, sem inclui-lo



### novo pdata - usar sempre este
load("LGG-GBM.merged.data_20150716.Rdata")

### all 450k
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/lgg.Rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/gbm450.Rda")
all.450k <- merge(gbm.450k,lgg.450k[,-c(2,3,4)],by="Composite.Element.REF")

lgg.mirna = NULL
pattern.primary = 'TCGA-([0-9A-Z]{2})-([0-9A-Z]{1,})-([0-1A-Z]{3})-([0-9A-Z]{3})-([0-9A-Z]{4}-([0-9]{2}))' #barcode
files = list.files("/home/thais/LGG-GBM/miRNA/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3") 
setwd("/home/thais/LGG-GBM/miRNA/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3")
z = 1 #first file
while (!is.na(files[z])) {
  c <- str_extract(as.character(files[z]), pattern.primary)
  if(!is.na(c)){ #do not read metadata files and non-tumor samples
    if(is.null(lgg.mirna)){ #read the first sample
      lgg.mirna = read.table(files[z],sep="\t",header=TRUE) 
      colnames(lgg.mirna)[2] <- as.character(c)
      lgg.mirna <- lgg.mirna[,c(1,2)]
    } 
    else{
      aux = read.table(files[z], header=TRUE,sep="\t")
      colnames(aux)[2] = c
      lgg.mirna = merge(lgg.mirna, aux[,c(1,2)], by="miRNA_ID")
    }
  }
  z = z + 1 #next file
}#end while

#### microRNA - miRNA
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/lgg.mirna_RPM.rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/lgg.mirna_RC.rda")
LGm6.LGG.ids <- read.table("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGm6_LGG.txt")

#LGm6.LGG.mi <- lgg.mirna.RC[,substr(colnames(lgg.mirna.RC),1,12) %in% LGm6.LGG.ids$V1]
#id.others <- lgg.mirna.RC[,!(substr(colnames(lgg.mirna.RC),1,12) %in% LGm6.LGG.ids$V1)]

aux <- subset(ids, cluster.meth %in% "LGm6-PAlike")
LGm6.LGG.mi <- lgg.mirna.RC[,substr(colnames(lgg.mirna.RC),1,12) %in% aux$case.id]
aux <- subset(ids, cluster.meth != "LGm6-PAlike")
id.others <- lgg.mirna.RC[,substr(colnames(lgg.mirna.RC),1,12) %in% aux$case.id]

a <- cbind(id.others,LGm6.LGG.mi)
rownames(a) <- as.character(lgg.mirna.RC$miRNA_ID)
library(samr)
y <- c(rep(1,ncol(id.others)),rep(2,ncol(LGm6.LGG.mi)))

#### colocar na.rm = TRUE em todas as funcoes hihihi (fonte: http://stackoverflow.com/questions/17418640/is-it-possible-to-set-na-rm-to-true-globally)
Funs <- Filter(is.function,sapply(ls(baseenv()),get,baseenv()))
na.rm.f <- names(Filter(function(x) any(names(formals(args(x)))%in% 'na.rm'),Funs))
ll <- lapply(na.rm.f,function(x)
{
  tt <- get(x)
  ss = body(tt)
  if (class(ss)!="{") ss = as.call(c(as.name("{"), ss))
  if(length(ss) < 2) print(x)
  else{
    if (!length(grep("na.rm = TRUE",ss[[2]],fixed=TRUE))) {
      ss = ss[c(1,NA,2:length(ss))]
      ss[[2]] = parse(text="na.rm = TRUE")[[1]]
      body(tt)=ss
      (unlockBinding)(x,baseenv())
      assign(x,tt,envir=asNamespace("base"),inherits=FALSE)
      lockBinding(x,baseenv())
    }
  }
})


diff.mirna <- SAMseq(a,y,resp.type = "Two class unpaired",fdr.output = 0.05, geneid = as.character(lgg.mirna.RC$miRNA_ID))
genes.low <- as.data.frame(diff.mirna$siggenes.table$genes.lo)
genes.up <- as.data.frame(diff.mirna$siggenes.table$genes.up)
volcano.miRNA <- rbind(genes.low,genes.up)
rownames(volcano.miRNA) <- as.character(volcano.miRNA$`Gene Name`)

####### Volcano microRNA
volcano.miRNA$`Fold Change` <- as.numeric(as.character(volcano.miRNA$`Fold Change`))
volcano.miRNA$`q-value(%)` <- as.numeric(as.character(volcano.miRNA$`q-value(%)`))
volcano.miRNA$threshold <- "1"
volcano.miRNA$log2FC <- volcano.miRNA$`Fold Change` + 0.01
volcano.miRNA$log2FC <- log2(volcano.miRNA$`Fold Change`)
volcano.miRNA[is.infinite(volcano.miRNA$log2FC),"log2FC"] <- 0
a <- subset(volcano.miRNA, `q-value(%)` < 5.0)
b <- subset(a, log2FC < -5) #hyper no da direita
c <- subset(a, log2FC < 0 & log2FC > -5)
volcano.miRNA[rownames(b),"threshold"] <- "2"
volcano.miRNA[rownames(c),"threshold"] <- "4"
b <- subset(a, log2FC > 5) #hyper no da esquerda
c <- subset(a, log2FC > 0 & log2FC < 5)
volcano.miRNA[rownames(b),"threshold"] <- "3"
volcano.miRNA[rownames(c),"threshold"] <- "4"
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_905_s/volcano_LGG.png",res=200,width=2000,height=2000)
ggplot(data=volcano.miRNA, aes(x=log2FC, y=`q-value(%)`,label=`Gene Name`,colour=threshold)) +
  geom_point() + geom_text() +
  #xlim(c(-1,1)) + ylim(c(0,35)) +
  xlab("log2 Fold-change") + ylab("q-value(%)") +
  labs(title = "Volcano Plot") +
  scale_color_manual(values=c("black", "green", "red","grey"),
                     labels=c("Not Significant","Hypermethylated","Hypermethylated","Not significant"),
                     name="Legend") 
dev.off()



####### Wilcoxon test - LGm6-LGG-IDHwt x todo mundo
aux <- subset(ids, cluster.meth %in% "LGm6-PAlike")
LGm6.LGG.mi <- lgg.mirna.RPM[,substr(colnames(lgg.mirna.RPM),1,12) %in% aux$case.id]
aux <- subset(ids, cluster.meth != "LGm6-PAlike")
id.others <- lgg.mirna.RPM[,substr(colnames(lgg.mirna.RPM),1,12) %in% aux$case.id]

volcano <- cbind(id.others,LGm6.LGG.mi)
rownames(volcano) <- as.character(lgg.mirna.RPM$miRNA_ID)
values <- t(na.omit(volcano))
values <- data.frame(values)


require(exactRankTests)
require(parallel)
system.time(w.p.values <- unlist(mclapply(values,
                                          function(probe) {
                                            zz <- wilcox.exact(probe[1:68],probe[69:94], exact=TRUE)
                                            z <- zz$p.value
                                            return(z)
                                          }, mc.cores=8)))

miRNA.pvalue.adj <- p.adjust(w.p.values, method = "BH")
save(miRNA.pvalue.adj, file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/miRNA_pvalues.rda")

load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/miRNA_pvalues.rda")

LGm6.LGG.mi_median <- apply(LGm6.LGG.mi,1,median,na.rm=TRUE)
others.mi_median <- apply(id.others,1,median,na.rm=TRUE)
subset.mirna <- data.frame(median_LGm6.LGG = LGm6.LGG.mi_median, median_others = others.mi_median, row.names = lgg.mirna.RPM$miRNA_ID, p.adj = miRNA.pvalue.adj)

a <- subset(subset.mirna, median_LGm6.LGG > 50 & median_others > 50 & p.adj < 0.05)
b <- as.data.frame(diff.mirna$siggenes.table$genes.up)
c <- as.data.frame(diff.mirna$siggenes.table$genes.lo)
colnames(b)[2] <- "teste"
colnames(c)[2] <- "teste"
d <- subset(b, teste  %in% rownames(a))
e <- subset(c, teste  %in% rownames(a))

#heatmap.miRNA <- volcano[c(as.character(d$teste),as.character(e$teste)),]
heatmap.miRNA <- a[as.character(subset(volcano.miRNA, threshold %in% c("3"))$`Gene Name`),]
#heatmap.miRNA <- volcano

clab <- as.matrix(c(rep("darkseagreen",68),rep("blue",26)))

heatmap.miRNA <- heatmap.miRNA + 1
#png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/testeSturm.png",res=300,width=2000,height=2000)
pdf("heatmap_miRNA2.pdf")
heatmap.plus.sm(as.matrix(log2(heatmap.miRNA)),
                col = colorpanel(75, "yellow", "black", "darkturquoise"),
                scale = "row",
                #labRow = NA,
                labCol = NA,
                Colv = NA,
                #Rowv = NA,
                ColSideColors = cbind(clab,clab),
                cexRow = 0.7
                #RowSideColors = rlab
                
)
dev.off()

##### CNV by Chromosome
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG-GBM.CNVbyChr.FB20150205.Rdata")
others.id <- subset(pd, cluster.meth %in% c("LGm4","LGm5","LGm6") & IDH.status %in% "WT" & !(case.id %in% LGm6.LGG.ids$V1))
a <- data.frame(case.id=as.character(LGm6.LGG.ids$V1),cluster.meth=rep("LGm6-PAlike",length(LGm6.LGG.ids$V1)))
ids <- rbind(others.id[,c("case.id","cluster.meth")],a)
cna.LGm <- merge(cnadf,ids,by.x="case",by.y="case.id")

ggplot(cna.LGm, aes(chr.arm, fill=factor(cna)))  +  
  geom_bar() + 
  scale_y_continuous(labels = percent) +  
  facet_wrap( ~ cluster.meth, scales = "free") +
  scale_fill_manual(values=c("cadetblue4","bisque","firebrick"))

#heatmap - CNV by Chr
png(filename = "CNVbyChr_Heatmap.png", bg="white", res=200, width=4000, height=2000)
ggplot(cna.LGm, aes(y=factor(cluster.meth),x=case)) + 
  geom_tile(aes(fill=factor(cna))) +
  #scale_fill_manual(values=c("cadetblue1","brown4","cadetblue4","red","purple","cyan","magenta","orange","green3","green3","red","black","green","red","purple","orange","yellow","blue","yellow","gray","cyan","magenta","blue","green3","red","brown1","black","cornsilk2","gray","white"),guide = guide_legend(title = "Legend")) + 
  scale_fill_manual(values=c("darkslategray4","bisque1","brown4")) +
  theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white"),axis.text.x = element_blank(),panel.margin = unit(0, "lines")) + #, strip.background=element_rect(fill=c("green","red","purple","orange","yellow","blue")),strip.text=element_text(size=12,color=c("white"))) + 
  facet_grid(cluster.meth ~ chr.arm,space="free",scales="free") +
  xlab("TCGA LGG and GBM") + 
  ylab("Gene Name") +
  guides(col = guide_legend(nrow = 8))
dev.off()

##### CNV by Gene
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG-GBM.CNVbyGene.FB20150205.Rdata")
aux <- subset(ids, cluster.meth %in% "LGm6-PAlike")
aux <- cn[rownames(cn) %in% as.character(aux$case.id),]
hc.LGm6.PA <- as.hclust(agnes(daisy(data.frame(aux), metric="gower"), method="ward")) #Fazer a clusterização
cnv.s <- aux[hc.LGm6.PA$order,] #Ordena as colunas com o resultado da clusterização

aux <- subset(ids, cluster.meth != "LGm6-PAlike")
aux <- cn[rownames(cn) %in% as.character(aux$case.id),]
hc.others <- as.hclust(agnes(daisy(data.frame(aux), metric="gower"), method="ward")) #Fazer a clusterização
cnv.s2 <- aux[hc.others$order,] #Ordena as colunas com o resultado da clusterização

all <- rbind(cnv.s,cnv.s2) #label + dados
all <- as.data.frame(all)
all$case.id <- rownames(all)
all.m <- melt(all,id="case.id")

all.lb <- merge(all.m,pd[,c("cluster.meth.lgg","cluster.meth.gbm","case.id","cluster.meth")],by="case.id")
all.lb$case.id <- factor(all.lb$case.id)
all.lb$case.id <- factor(all.lb$case.id, levels = rownames(all))
a <- subset(pd,IDH.status=="WT" & age <=40)
all.lb$case.id <- as.factor(all.lb$case.id)
aux <- matrix()
png(filename = "CNVbyGene.png", bg="white", res=200, width=4000, height=2000)
ggplot(all.lb, aes(y=factor(variable),x=case.id)) + 
  geom_tile(aes(fill=factor(value))) +
  #scale_fill_manual(values=c("cadetblue1","brown4","cadetblue4","red","purple","cyan","magenta","orange","green3","green3","red","black","green","red","purple","orange","yellow","blue","yellow","gray","cyan","magenta","blue","green3","red","brown1","black","cornsilk2","gray","white"),guide = guide_legend(title = "Legend")) + 
  scale_fill_manual(values=c("darkslategray4","darkslategray2","bisque1","brown2","brown4")) +
  theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white"),axis.text.x = element_blank(),axis.text.y = element_blank(),panel.margin = unit(0, "lines")) + #, strip.background=element_rect(fill=c("green","red","purple","orange","yellow","blue")),strip.text=element_text(size=12,color=c("white"))) + 
  facet_grid(age.status ~ cluster.meth,space="free",scales="free") +
  xlab("TCGA LGG and GBM") + 
  ylab("Gene Name") +
  guides(col = guide_legend(nrow = 8))
dev.off()


##### Mutation - LGm6
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG-GBM.Mutation.FB20150708.Rdata")
aux <- subset(ids, cluster.meth %in% "LGm6-PAlike")
aux <- mut[rownames(mut) %in% as.character(aux$case.id),]
hc.LGm6.PA <- as.hclust(agnes(daisy(data.frame(aux), metric="gower"), method="ward")) #Fazer a clusterização
mut.s <- aux[hc.LGm6.PA$order,] #Ordena as colunas com o resultado da clusterização

aux <- subset(ids, cluster.meth != "LGm6-PAlike")
aux <- mut[rownames(mut) %in% as.character(aux$case.id),]
hc.others <- as.hclust(agnes(daisy(data.frame(aux), metric="gower"), method="ward")) #Fazer a clusterização
mut.s2 <- aux[hc.others$order,] #Ordena as colunas com o resultado da clusterização

all.mut <- rbind(mut.s,mut.s2) #label + dados
all.mut <- as.data.frame(all.mut)
all.mut$case.id <- rownames(all.mut)
all.2 <- melt(all.mut,id="case.id")

all.melt <- merge(all.2,ids,by="case.id")
all.melt$case.id <- factor(all.melt$case.id)
all.melt$case.id <- factor(all.melt$case.id, levels = rownames(all.mut))
png(filename = "Mut_ByGene.png", bg="white", res=200, width=4000, height=2000)
ggplot(all.melt, aes(y=factor(variable),x=case.id)) + 
  geom_tile(aes(fill=factor(value))) +
  #scale_fill_manual(values=c("cadetblue1","brown4","cadetblue4","red","purple","cyan","magenta","orange","green3","green3","red","black","green","red","purple","orange","yellow","blue","yellow","gray","cyan","magenta","blue","green3","red","brown1","black","cornsilk2","gray","white"),guide = guide_legend(title = "Legend")) + 
  scale_fill_manual(values=c("gray","black")) +
  theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white"),axis.text.x = element_blank(),axis.text.y = element_blank(),panel.margin = unit(0, "lines")) + #, strip.background=element_rect(fill=c("green","red","purple","orange","yellow","blue")),strip.text=element_text(size=12,color=c("white"))) + 
  facet_grid(. ~ cluster.meth,space="free",scales="free") +
  xlab("TCGA LGG and GBM") + 
  ylab("Gene Name") +
  guides(col = guide_legend(nrow = 8))
dev.off()


##### ES only for LGm6-PA like

uu.nl <- epiGene.RNAseq[[1]]
clusters <- c("LGm6","others")
#clusters <- "LGm1"


enrichment.normal <- lapply(clusters, function(cluster, epiGene.ss=epiGene.RNAseq) {
  enrichment.group <- unlist(mclapply(epiGene.ss, function(uu.nl, cluster.group=cluster) {
    control <- (!uu.nl$meth.gt0.3[637:length(uu.nl$meth.gt0.3)] & !uu.nl$exp.lt.meanUnmeth[637:length(uu.nl$exp.lt.meanUnmeth)])
    uu.nl <- subset(uu.nl, IDH.status %in% "WT")
    uu.nl[uu.nl$cluster %in% c("LGm1","LGm2","LGm3","LGm4","LGm5") | (uu.nl$cluster %in% "LGm6" & uu.nl$tumor.type %in% "GBM"),"cluster"] <- "others"
    uu.nl$es <- NA
    uu.nl$es[1:length(uu.nl$es)] <- uu.nl$meth.gt0.3[1:length(uu.nl$es)] & uu.nl$exp.lt.meanUnmeth[1:length(uu.nl$es)]
    uu.nl$Kclass <- NA
    uu.nl$Kclass <- uu.nl$cluster == cluster.group
    es.table <- table(uu.nl$Kclass, uu.nl$es)#, dnn=c("Class", "ES"))
    #print()
    if(dim(es.table)[1] > 1 & dim(es.table)[2] > 1){
      # print("teste")
      if((es.table[2,2]/rowSums(es.table)[2] > 0.5) & (es.table[2,2]/colSums(es.table)[2] > 0.5)  & (table(control)[2]/(table(control)[2] + table(control)[1]) > 0.2 | is.na(table(control)[2]/(table(control)[2] + table(control)[1])))) {
        return(fisher.test(table(uu.nl$Kclass, uu.nl$es))$p.value)
      } else {
        return(NULL)
      }
    }
    #}  
  }, mc.cores=8))
  return(enrichment.group)
})

### nothing. Notice that the cutoff for the normals were reduced from 50% to 20%, still nothing


# Wilcoxon-test -> IDHwt from LGm4-5-6
####### Volcano LGm4+LGm5 x LGm6
#Beta Value Diff M1 x M2
rownames(LGG.GBM.new) <- as.character(LGG.GBM.new$Composite.Element.REF)
data <- LGG.GBM.new[, colnames(LGG.GBM.new) %in% colnames(dat.lgg.gbm.new.noXY.dic.oecg)] #[1] 25978   936

M1 <- pd[pd$cluster.meth %in% "LGm4" & !is.na(pd$cluster.meth) & !is.na(pd$IDH.status) & pd$IDH.status %in% "WT",]
M2 <- pd[pd$cluster.meth %in% "LGm5" & !is.na(pd$cluster.meth) & !is.na(pd$IDH.status) & pd$IDH.status %in% "WT",]
#Selecting the Samples
M1 <- data[,substr(colnames(data),1,12) %in% as.character(M1$case.id)]
M2 <- data[,substr(colnames(data),1,12) %in% as.character(M2$case.id)]
volcano <- cbind(M1,M2)

values <- t(na.omit(volcano))
values <- data.frame(values)
#values.names <- colnames(values)

LGm4xLGm5 <- values
save(LGm4xLGm5,file='LGm4xLGm5.rda')




require(exactRankTests)
require(parallel)
dim(LGm6.LGGxLGm6.GBM)
system.time(w.p.values <- unlist(mclapply(LGm6.LGGxLGm6.GBM,
                                          function(probe) {
                                            zz <- wilcox.exact(probe[1:41],probe[42:67], exact=TRUE)
                                            z <- zz$p.value
                                            return(z)
                                          }, mc.cores=2)))

#lggs.p.value <- as.matrix(w.p.values)



#LGm4 x LGm5 - 1:136 137:362
#LGm4 x LGm6-LGG - 1:136 137:162
#LGm4 x LGm6-GBM - 1:136 137:177
#LGm5 x LGm6-GBM - 1:226 227:267
#LGm5 x LGm6-LGG - 1:226 227:252
#LGm6-GBM x LGm6-LGG - 1:41 42:67


LGm5xLGm6.LGG <- t(LGm5xLGm6.LGG)
LGm5xLGm6.LGG <- as.data.frame(LGm5xLGm6.LGG)
LGm5xLGm6.LGG$p.value <- w.p.values
LGm5xLGm6.LGG$p.value.adj <- p.adjust(LGm5xLGm6.LGG$p.value, method="BH")
LGm5xLGm6.LGG$meanM1 <- apply(LGm5xLGm6.LGG[,1:226],1,mean,na.rm=T)
LGm5xLGm6.LGG$meanM2 <- apply(LGm5xLGm6.LGG[,227:252],1,mean,na.rm=T)
LGm5xLGm6.LGG$DiffMean <- LGm5xLGm6.LGG$meanM1 - LGm5xLGm6.LGG$meanM2
dim(subset(LGm5xLGm6.LGG, p.value.adj < 0.05 & abs(DiffMean) > 0.3))  
save(LGm5xLGm6.LGG,file="LGm5xLGm6-LGG_pvalue.rda")

my.wilcox.test.p.value <- function(...) {
  obj<-try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error"))
    return(NA)
  else
    return(obj$p.value)
}

wilcoxt <- apply(volcano, 1, function(i) my.wilcox.test.p.value(i[1:204],i[205:244]) )
volcano$p.value = wilcoxt
volcano$p.value.adj = p.adjust(volcano$p.value, method = "fdr")

p4x5 <- rownames(subset(LGm4xLGm5, p.value.adj < 0.05 & abs(DiffMean) > 0.3))
p4x6l <- rownames(subset(LGm4xLGm6.LGG, p.value.adj < 0.05 & abs(DiffMean) > 0.3))
p4x6g <- rownames(subset(LGm4xLGm6.GBM, p.value.adj < 0.05 & abs(DiffMean) > 0.3))
p5x6 <- rownames(subset(LGm5xLGm6.GBM, p.value.adj < 0.05 & abs(DiffMean) > 0.3))
p6lx6g <- rownames(subset(LGm6.LGGxLGm6.GBM, p.value.adj < 0.05 & abs(DiffMean) > 0.3))



p <- c(p4x5,p4x6g,p4x6l,p5x6,p6lx6g,p5x6l)
rlab <- data.frame(cor=c(rep("pink",length(p4x5)),rep("pink2",length(p4x6g)),rep("pink4",length(p4x6l)),rep("purple",length(p5x6)),rep("palevioletred3",length(p6lx6g)),rep("palevioletred1",length(p5x6l))),id=p)
p <- unique(p)
p <- p6lx6g
rlab <- rlab[!duplicated(rlab$id),]
LGG.GBM.new.noXY.s <- LGG.GBM.new
rownames(LGG.GBM.new.noXY.s) <- as.character(LGG.GBM.new.noXY.s$Composite.Element.REF)
LGG.GBM.new.noXY.s <- LGG.GBM.new.noXY.s[,colnames(dat.lgg.gbm.new.noXY.dic.oecg)]
LGG.GBM.new.noXY.s <- LGG.GBM.new.noXY.s[ rownames(dat.lgg.gbm.new.noXY.dic.oecg),]
LGG.GBM.new.order2 <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
LGG.GBM.new.order2 <- LGG.GBM.new.order2[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
colnames(LGG.GBM.new.order2) <- substr(colnames(LGG.GBM.new.order2),1,12)
a <- pd
rownames(a) <- as.character(a$case.id)
a <- a[colnames(LGG.GBM.new.order2),]
levels(a$cluster.meth) <- c("green","red","purple","orange","yellow","blue")
a$tumor.type <- as.factor(a$tumor.type)
levels(a$tumor.type) <- c("gray","black") #LGG black
a$has.HM450 <- as.factor(a$has.HM450)
levels(a$has.HM450) <- c("gray","black") #black yes
a$treat.adjuvant.radiotherapy <- as.factor(a$treat.adjuvant.radiotherapy)
levels(a$treat.adjuvant.radiotherapy) <- c("gray","black") #black yes
#a <- colnames(LGG.GBM.new.order) %in% colnames(GBM.LGG.27.450k.nearest.gene)
#teste <- LGG.GBM.new.order[,a]
a <- a[856:932,c("cluster.meth","tumor.type","age","has.HM450","treat.adjuvant.radiotherapy")]
a$age <- jet.colors(100)[ceiling((a$age-min(a$age,na.rm=T))/(max(a$age,na.rm=T)-min(a$age,na.rm=T))*(length(jet.colors(100))-1))+1]
LGG.GBM.new.order2 <- as.matrix(LGG.GBM.new.order2[,856:932])
aux <- subset(a, tumor.type == "gray")
aux2 <- subset(a, tumor.type == "black")
a <- a[c(rownames(aux),rownames(aux2)),]

  
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/heatmap_IDHwtProbes2.png",res=300,width=2000,height=2000)
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_DNAMeth_tert.pdf")
heatmap.plus.sm(as.matrix(LGm6.450k[rownames(subset(LGm6.450k, p.value.adj < 0.05 & abs(DiffMean) > 0.3)),1:38]),
                col = jet.colors(75),
                scale = "none",
                #labRow = NA,
                labCol = NA,
                Colv = NA,
                # Rowv = NA,
               # ColSideColors = as.matrix(a),
               ColSideColors = cbind(c(rep("lightblue",26),rep("blue",12)),c(rep("lightblue",26),rep("blue",12))),
                #RowSideColors = cbind(as.character(rlab$cor),as.character(rlab$cor)),
                margins = c(1, 6)
                #RowSideColors = rlab
                
)
dev.off()

aux <- subset(pd, cluster.meth %in% "LGm6" & tumor.type %in% "LGG" & IDH.status %in% "WT")
aux2 <- subset(pd, cluster.meth %in% "LGm6" & tumor.type %in% "GBM" & IDH.status %in% "WT")

rownames(all.450k) <- as.character(all.450k$Composite.Element.REF)
LGm6.450k.LGG <- all.450k[,substr(colnames(all.450k),1,12) %in% as.character(aux$case.id)]
LGm6.450k.GBM <- all.450k[,substr(colnames(all.450k),1,12) %in% as.character(aux2$case.id)]
LGm6.450k <- cbind(LGm6.450k.LGG,LGm6.450k.)

aux <- subset(pd, cluster.meth %in% "LGm6" & IDH.status %in% "WT")

