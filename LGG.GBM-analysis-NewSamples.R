
#LGG New Samples. Batches: 
x = 1 #batch
lgg.patients = NULL
normal = NULL
control = NULL
i=0
end = ".12.0"
beginning <- "/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/jhu-usc.edu_LGG.HumanMethylation450.Level_3." #folder
pattern = 'TCGA-([0-9A-Z]{2})-([0-9A-Z]{1,})-([0-9A-Z]{3})-([0-9A-Z]{3})-([0-9A-Z]{4}-([0-9]{2}))' #barcode
pattern.primary = 'TCGA-([0-9A-Z]{2})-([0-9A-Z]{1,})-([0-1A-Z]{3})-([0-9A-Z]{3})-([0-9A-Z]{4}-([0-9]{2}))' #barcode
while (x<17){ #from 3.12.10.0 ate 3.16.10.0
  dir = (paste(c(paste(c(paste(c(beginning,x), collapse=''),""),collapse=''),end),collapse=''))
  files = list.files(dir,pattern="txt",full.names=T) 
  z = 1 #first file
  while (!is.na(files[z])) {
    a = strsplit(files[z],"/") 
    b = a[[1]][[9]] #file name
    c <- str_extract(b, pattern.primary)
    patientID <- c
    if(!is.na(c)){ #do not read metadata files
      if (substr(c,14,15) > "19"){ #control samples
        if(is.null(control)){ #read the first file
          control = read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1) 
          control = subset(control,select=c(Composite.Element.REF,Gene_Symbol,Chromosome,Genomic_Coordinate,Beta_value))
          colnames(control)[5] = patientID
        }
        else{
          aux = read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1)
          colnames(aux)[2] = patientID
          control = merge(control, aux[,1:2], by.x ="Composite.Element.REF", by.y = "Composite.Element.REF") 
        }
      } #end if control
      else if(substr(c,14,15) > "09"){ #normal samples
        if(is.null(normal)){ #read the first file
          normal = read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1) 
          normal = subset(normal,select=c(Composite.Element.REF,Gene_Symbol,Chromosome,Genomic_Coordinate,Beta_value)) 
          colnames(normal)[5] = patientID
        }
        else{
          aux = read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1)
          colnames(aux)[2] = patientID
          normal = merge(normal, aux[,1:2], by.x ="Composite.Element.REF", by.y = "Composite.Element.REF")
        }
      }
      else { #tumor samples
        if(is.null(lgg.patients)){ #read the first sample
          lgg.patients = read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1) 
          lgg.patients = subset(lgg.patients,select=c(Composite.Element.REF,Gene_Symbol,Chromosome,Genomic_Coordinate,Beta_value)) 
          colnames(lgg.patients)[5] = patientID
        }
        else{
          aux = read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1)
          colnames(aux)[2] = patientID
          lgg.patients = merge(lgg.patients, aux[,1:2], by.x ="Composite.Element.REF", by.y = "Composite.Element.REF")
        }
      }
    } #end !is.na(c)
    z = z + 1 #next file
  }#end while
  print(x)
  x = x + 1 #next folder
}
#save(lgg.patients,control,file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/lgg_new_samples.Rda")
library(stringr)
pattern.primary <- 'TCGA-([0-9A-Z]{2})-([0-9A-Z]{1,})-([0-1A-Z]{3})-([0-9A-Z]{3})-([0-9A-Z]{4}-([0-9]{2}))'
metadata <- lgg.patients[,1:4]
primary <- str_extract(colnames(lgg.patients)[5:ncol(lgg.patients)], pattern.primary) #remove the recurrent samples
primary <- na.omit(primary) #5 recurrents
lgg.patients <-lgg.patients[,primary]
lgg.patients <- cbind(metadata,lgg.patients) #220 samples
load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG.GBM.losangeles/lgg.gbm.normals.plusAnno.rda")
load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Brennan2013_GBM27K-450K//gbm.27.450.Brennan2013.rda") ##SNPs have already been removed
load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG.GBM.losangeles/lgg_289_samples.Rda")
#primary <- str_extract(colnames(dat)[5:ncol(dat)], pattern.primary) #no recurrent
LGG.merge <- merge(dat,lgg.patients[,c(1,5:224)], by.x="Composite.Element.REF", by.y="Composite.Element.REF")
GBM.LGG.27.450k <- merge(LGG.merge, GBM27.450.brennan2013, by.x = "Composite.Element.REF", by.y = "TargetID", all.x = F)
GBM.LGG.27.450k <- GBM.LGG.27.450k[,-c(514:516)]
dimnames(GBM.LGG.27.450k)[[1]] <- GBM.LGG.27.450k$Composite.Element.REF #dim 909
#GBM.LGG.27.450k.noX <- subset(GBM.LGG.27.450k, Chromosome.x != "X")
#GBM.LGG.27.450k.noXY <- subset(GBM.LGG.27.450k.noX, Chromosome.x != "Y")
GBM.LGG.27.450k.noXY <- GBM.LGG.27.450k[GBM.LGG.27.450k$Chromosome.x != "X" & GBM.LGG.27.450k$Chromosome.x != "Y",]


## dichotomizing the data
GBM.LGG.27.450k.noXY.dic <- GBM.LGG.27.450k.noXY[ ,5:909]
GBM.LGG.27.450k.noXY.dic <- GBM.LGG.27.450k.noXY.dic>0.3
storage.mode(GBM.LGG.27.450k.noXY.dic) <- "numeric"
## Use probes methylated more than 10% of the tumor samples
#GBM.LGG.27.450k.noXY.dic <- GBM.LGG.27.450k.noXY.dic[(rowSums(GBM.LGG.27.450k.noXY[,5:909])/ncol(GBM.LGG.27.450k.noXY[,5:909]))*100 >10,]
GBM.LGG.27.450k.noXY.dic <- GBM.LGG.27.450k.noXY.dic[(rowSums(GBM.LGG.27.450k.noXY.dic[,5:909])/ncol(GBM.LGG.27.450k.noXY.dic[,5:909]))*100 >10,]

#head(FDb.HM450.oecg.3kb)
FDb.HM450.oecg.3kb.s <- subset(FDb.HM450.oecg.3kb, oecg.3kb>1.25)
GBM.LGG.27.450k.noXY.dic.oecg <- GBM.LGG.27.450k.noXY.dic[rownames(GBM.LGG.27.450k.noXY.dic) %in% (FDb.HM450.oecg.3kb.s$probeid),]

#TODO:
## need to select tissue specific probes (>0.3)
avg.lgg.gbm.27.450.noXY <- apply(GBM.LGG.27.450k.noXY[,5:909], 1, mean, na.rm=TRUE)
avg.norm.lgg.gbm.27.450.noXY <- merge(as.data.frame(avg.lgg.gbm.27.450.noXY), avg.normals, by = 0, all.x = T)
#library(LSD)
#comparisonplot(avg.norm.lgg.gbm.27.450.noXY$avg.normals, avg.norm.lgg.gbm.27.450.noXY$avg.lgg.gbm.27.450.noXY, pimp=TRUE)
lgg.gbm.27.450.tumor.specific <- (subset(avg.norm.lgg.gbm.27.450.noXY, avg.normals < 0.3))

##selecting tumor specific probes for lgg.gbm
dat.lgg.gbm.27.450.noXY.dic.oecg <- GBM.LGG.27.450k.noXY.dic.oecg[(rownames(GBM.LGG.27.450k.noXY.dic.oecg) %in% lgg.gbm.27.450.tumor.specific$Row.names), ]

library(ConsensusClusterPlus)
title="/dados/ResearchProjects/thais/TCGA/LGG.GBM/ConsensusCluster.LGG.GBM.27.450_905_samples"
setwd(title)
cc.out.purity.lgg.gbm.27.450 <- ConsensusClusterPlus(d=dat.lgg.gbm.27.450.noXY.dic.oecg, 
                                                     maxK=10, 
                                                     reps=1000, 
                                                     pItem=0.8, 
                                                     pFeature=1,
                                                     title=title, 
                                                     clusterAlg="hc",
                                                     distance="binary",
                                                     innerLinkage = "ward", 
                                                     finalLinkage = "ward",
                                                     seed=1022,
                                                     plot="pdf",
                                                     #writeTable=TRUE,
                                                     verbose = TRUE)
ccICL = calcICL(cc.out.purity.lgg.gbm.27.450, plot='pdf')
cc.out.purity.lgg.gbm.27.450.2 <- cbind(as.matrix(cc.out.purity.lgg.gbm.27.450[[2]][[3]]), 
                                        as.matrix(cc.out.purity.lgg.gbm.27.450[[3]][[3]]), 
                                        as.matrix(cc.out.purity.lgg.gbm.27.450[[4]][[3]]), 
                                        as.matrix(cc.out.purity.lgg.gbm.27.450[[5]][[3]]), 
                                        as.matrix(cc.out.purity.lgg.gbm.27.450[[6]][[3]]),
                                        as.matrix(cc.out.purity.lgg.gbm.27.450[[7]][[3]]),
                                        as.matrix(cc.out.purity.lgg.gbm.27.450[[8]][[3]]),
                                        as.matrix(cc.out.purity.lgg.gbm.27.450[[9]][[3]]),
                                        as.matrix(cc.out.purity.lgg.gbm.27.450[[10]][[3]]))

cc.out.purity.lgg.gbm.27.450.2 <- as.data.frame(cc.out.purity.lgg.gbm.27.450.2)
names(cc.out.purity.lgg.gbm.27.450.2) <- c("K2",  "K3",  "K4" , "K5" , "K6" , "K7",  "K8" , "K9" , "K10")

#cc.out.purity.lgg.gbm.27.450.2$TumorType <- "UNK"
#cc.out.purity.lgg.gbm.27.450.2[1:509,"TumorType"] <- "LGG"
#cc.out.purity.lgg.gbm.27.450.2[510:905,"TumorType"] <- "GBM"
cc.out.purity.lgg.gbm.27.450.2$TCGAIDnew <- substr(rownames(cc.out.purity.lgg.gbm.27.450.2),1,12)
cc.out.purity.lgg.gbm.27.450.2$TCGAID <- rownames(cc.out.purity.lgg.gbm.27.450.2)
#509 LGG 396 GBM

##### New Clinical + Mutation Status data (Synapse) 22/07/2014
#sem fazer o consensus
cc.out.purity.lgg.gbm.27.450.2 <- data.frame(TCGAID = colnames(GBM.LGG.27.450k)[5:909])
cc.out.purity.lgg.gbm.27.450.2$TCGAIDnew <- substr(cc.out.purity.lgg.gbm.27.450.2$TCGAID,1,12)
#fim sem fazer o cc
mutation.clinical.data.synapse <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/new_mut_info.txt")
mutation.clinical.cc <- merge(cc.out.purity.lgg.gbm.27.450.2,mutation.clinical.data.synapse,by.x="TCGAIDnew",by.y="id",all.x=T)
rownames(mutation.clinical.cc) <- mutation.clinical.cc$TCGAID
#mutation.clinical.cc$IDH1.status <- str_replace(mutation.clinical.cc$IDH1.status, "R132C|R132G|R132H|R132S|Atypical","Mut")
agecat <- cut(mutation.clinical.cc$age, c(10.90,38.80,52,63,89.30))
mutation.clinical.cc$age.strat <- as.character(agecat)
#table(mutation.clinical.cc$IDH1.status,mutation.clinical.cc$cluster.meth)
#table(mutation.clinical.cc$mut.IDH1,mutation.clinical.cc$cluster.meth)

#parei aqui, depois dar crtl z
lgg.gbm.anno.27.450 <- merge(mutation.clinical.cc, brennan.class.s[,c("DNA.Methylation.Clusters","TCGAIDnew")], by = "TCGAIDnew", all.x =T)
rownames(lgg.gbm.anno.27.450) <- lgg.gbm.anno.27.450$TCGAID
xlgg.gbm.anno.27.450 <- merge(lgg.gbm.anno.27.450, clinical.cluster.genetics[,c("Tumor" ,"MethylationCluster","COCCluster")], by.x = "TCGAIDnew", by.y = "Tumor", all.x = T)
rownames(xlgg.gbm.anno.27.450) <- xlgg.gbm.anno.27.450$TCGAID
xlgg.gbm.anno.27.450 <- merge(xlgg.gbm.anno.27.450,gbm.molecular.subtype[,c("sample.id","Cluster")], by.x = "TCGAIDnew", by.y = "sample.id", all.x = T)
rownames(xlgg.gbm.anno.27.450) <- xlgg.gbm.anno.27.450$TCGAID
#colnames(xlgg.gbm.anno.27.450)[c(61,62)] <- c("GBM.methylation.cl","LGG.methylation.c")
xxMODr <- merge(xlgg.gbm.anno.27.450,leukocyte, by.x=0, by.y= "x",all.x=T)
rownames(xxMODr) <- xxMODr$Row.names
xxMODr <- xxMODr[,-1]
leukocytecat <- cut(xxMODr$leukocyte, c(0.0136,0.10340341,0.15807936,0.24441133,0.77))
xxMODr$leukocyte.strat <- leukocytecat
metadata <- merge(xxMODr,rna.estimates,by.x="TCGAIDnew",by.y=0,all.x=T)
colnames(metadata)[266] <- "rna.estimates"
rownames(metadata) <- metadata$TCGAID
#metadata$cluster.mrna.new <- NA
#metadata[metadata$cluster.mRNA %in% "yellow","cluster.mrna.new"] <- "LGr1"
#metadata[metadata$cluster.mRNA %in% "gray","cluster.mrna.new"] <- "LGr2"
#metadata[metadata$cluster.mRNA %in% "cyan","cluster.mrna.new"] <- "LGr3"
#metadata[metadata$cluster.mRNA %in% "magenta","cluster.mrna.new"] <- "LGr4"
#metadata[metadata$cluster.mRNA %in% "blue","cluster.mrna.new"] <- "LGr5"
#metadata[metadata$cluster.mRNA %in% "green3","cluster.mrna.new"] <- "LGr6"
#metadata[metadata$cluster.mRNA %in% "red","cluster.mrna.new"] <- "LGr7"
#metadata[metadata$cluster.mRNA %in% "brown","cluster.mrna.new"] <- "Undefined"
sturm.et.al <- read.xlsx("/dados/ResearchProjects/thais/TCGA/LGG.GBM/mmc2.xlsx",1)
metadata <- merge(metadata, sturm.et.al[,c(1,2)],by.x="TCGAIDnew",by.y="Sample.ID",all.x=T)
rownames(metadata) <- metadata$TCGAID

#RNA sequencing data for LGG + GBM (from synapse, Aug 4th)
rna <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/ge.all617.unc.rsem.matrix.txt")
colnames(rna) <- str_replace_all(colnames(rna),"[.]","-")
rna.seq.lgg.gbm <- rna[,colnames(rna) %in% as.character(metadata$TCGAIDnew)]
entrezID <- str_extract(rownames(rna.seq.lgg.gbm),"[^|]*$")
rna.seq.lgg.gbm$entrezID <- as.matrix(entrezID)
#library(biomaRt)
#ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
#gene.location <- getBM(attributes=c("chromosome_name","start_position","end_position","entrezgene"),filters="entrezgene",values=rna.seq.lgg.gbm$entrezID,mart=ensembl)
#b <- gene.location[gene.location$chromosome_name %in% c(1:22,"X","Y"),]
### map gene ID to genomic coordinates
gene.location <- read.delim("/dados/ResearchProjects/thais/TCGA/Normals/mRNA_Brain/Galaxy1-[Homo_sapiens_genes_(GRCh37.p13)].tabular")
gene.location <- gene.location[gene.location$Chromosome.Name %in% c(1:22,"X","Y"),]
gene.location <- gene.location[!is.na(gene.location$EntrezGene.ID),]
gene.location <- gene.location[!duplicated(gene.location$EntrezGene.ID),] #the duplicates are different transcripts, not different coordinates
library(GenomicRanges)
gene.GR <- GRanges(seqnames = paste0("chr",gene.location$Chromosome.Name),ranges = IRanges(start = gene.location$Gene.Start..bp., end=gene.location$Gene.End..bp.), strand = gene.location$Strand, symbol = gene.location$Associated.Gene.Name, EntrezID = gene.location$EntrezGene.ID)
probe.info <- GRanges(seqnames = paste0("chr",GBM.LGG.27.450k$Chromosome.x), ranges = IRanges(start = GBM.LGG.27.450k$Genomic_Coordinate.x, end = GBM.LGG.27.450k$Genomic_Coordinate.x), probeID = GBM.LGG.27.450k$Composite.Element.REF)
distance <- as.data.frame(distanceToNearest(probe.info,gene.GR)) #closest gene to each 450k probe
rownames(gene.location) <- NULL
gene.order.by.distance <- gene.location[distance$subjectHits,]
gene.order.by.distance$distance <- as.matrix(distance$distance)
GBM.LGG.27.450k.nearest.gene <- cbind(GBM.LGG.27.450k, gene.order.by.distance[,c("Associated.Gene.Name","EntrezGene.ID","distance")])
info <- GBM.LGG.27.450k.nearest.gene[,c(1:4,910:912)]
colnames(GBM.LGG.27.450k.nearest.gene) <- substr(colnames(GBM.LGG.27.450k.nearest.gene),1,12)
GBM.LGG.27.450k.nearest.gene <- GBM.LGG.27.450k.nearest.gene[,colnames(rna.seq.lgg.gbm)[1:577]]
GBM.LGG.27.450k.nearest.gene <- cbind(info,GBM.LGG.27.450k.nearest.gene)
b <- metadata[metadata$TCGAIDnew %in% colnames(GBM.LGG.27.450k.nearest.gene),] #select the metadata information for the samples that have both DNA methylation and RNA seq data
#to adjust the proportion of LGGs:
#calculate the GBM proportion in the original data (metadata)
#match the proportion with the 'b' object (randomly select)
#LGm1: orig: 68% LGG (34) now: 15
#LGm2: orig: 96% LGG (34) now: 24
#LGm3: OK
#LGm4: orig: 14% LGG (34) now: 8
#LGm5: orig: 19% LGG (34) now: 12
#LGm6: orig: 36% LGG (34) now: 8
aux <- b[b$cluster.meth == "LGm1" & b$tumor.type == "LGG",]
b <- b[!(b$cluster.meth == "LGm1" & b$tumor.type == "LGG"),]
aux <- aux[1:15,]
b <- rbind(b,aux)
aux <- b[b$cluster.meth == "LGm2" & b$tumor.type == "LGG",]
b <- b[!(b$cluster.meth == "LGm2" & b$tumor.type == "LGG"),]
aux <- aux[100:123,]
b <- rbind(b,aux)
aux <- b[b$cluster.meth == "LGm4" & b$tumor.type == "LGG",]
b <- b[!(b$cluster.meth == "LGm4" & b$tumor.type == "LGG"),]
aux <- aux[10:17,]
b <- rbind(b,aux)
aux <- b[b$cluster.meth == "LGm5" & b$tumor.type == "LGG",]
b <- b[!(b$cluster.meth == "LGm5" & b$tumor.type == "LGG"),]
aux <- aux[15:26,]
b <- rbind(b,aux)
aux <- b[b$cluster.meth == "LGm6" & b$tumor.type == "LGG",]
b <- b[!(b$cluster.meth == "LGm6" & b$tumor.type == "LGG"),]
aux <- aux[1:8,]
b <- rbind(b,aux)
GBM.LGG.27.450k.nearest.gene.sel <- GBM.LGG.27.450k.nearest.gene[,b$TCGAIDnew]
GBM.LGG.27.450k.nearest.gene.sel <- cbind(GBM.LGG.27.450k.nearest.gene[,1:7],GBM.LGG.27.450k.nearest.gene.sel)
rna.seq.lgg.gbm.sel <- rna.seq.lgg.gbm[,b$TCGAIDnew]
rna.seq.lgg.gbm.sel <- cbind(rna.seq.lgg.gbm$entrezID,rna.seq.lgg.gbm.sel)

cpg.gene <- merge(na.omit(GBM.LGG.27.450k.nearest.gene.sel),rna.seq.lgg.gbm.sel,by.x="EntrezGene.ID",by.y="rna.seq.lgg.gbm$entrezID")
dimnames(cpg.gene)[[1]] <- paste(cpg.gene[,"Composite.Element.REF"], cpg.gene[,"Associated.Gene.Name"], sep=".")
#(d[,8:584],d[,586:1162]) #mDNA methylation x Gene expression

meth.g <- cpg.gene[,8:254]
names(meth.g) <- substr(names(meth.g), 1, 12)
expr.g <- cpg.gene[,255:501]
names(expr.g) <- substr(names(expr.g), 1, 12)
epiGene <- list()
j=1
for(i in 1:dim(cpg.gene)[1]) {
  
  meth.gu <- names(meth.g)[meth.g[i, ] < 0.3]
  meth.gm <- names(meth.g)[meth.g[i, ] >= 0.3]
  # if(fisher.test(table(yy))$p.value < 0.05){
  if(length(meth.gu)>1 & length(meth.gm)>1){
    expr.mean.gm <- mean(as.numeric(expr.g[i, meth.gm]), na.rm=T)
    if(expr.mean.gm < 1.28*sd(as.numeric(expr.g[i, meth.gu]), na.rm=T)){
      expr.mean.gu <- mean(as.numeric(expr.g[i, meth.gu]), na.rm = T)
      yy <- cbind(t(meth.g[i, ] >= 0.3), t(expr.g[i, ] < expr.mean.gu))
      yy <- as.data.frame(yy)
      gene <- colnames(yy)[1]
      colnames(yy) <- c("meth.gt0.3", "exp.lt.meanUnmeth")
      if(as.matrix(table(yy$meth.gt0.3, yy$exp.lt.meanUnmeth))[2,2]/as.numeric(table(yy$meth.gt0.3)[2])>0.8){
        yy$geneName <- gene
        yy$id <- rownames(yy)
        yy$methValues <- as.numeric(meth.g[i, ])
        yy$exprValues <- as.numeric(expr.g[i, ])    
        yy <- merge(metadata, yy, by.x = "TCGAIDnew", by.y = "id")
        epiGene[[j]] <- as.data.frame(yy)
        names(epiGene)[[j]] <- as.character(gene)
        j <- j+1
        cat("found a ES", j, "\n")
      }else{
        cat("no step3", "\n")
      }
    }else{
      cat("no ES group", "\n")
    }
  }else{
    cat("no methylated group", i, "\n")
  }
}


### PULEI ESSA PARTE
promoter.probes <- subset(cpg.gene, distance <= 5000)
promoter.probes <- promoter.probes[,c("Composite.Element.REF","Gene_Symbol.x","EntrezGene.ID","Associated.Gene.Name")]
gene.id.prom <- intersect(names(epiGene), promoter.probes$Associated.Gene.Name)
gene.id.prom <- promoter.probes[promoter.probes$Associated.Gene.Name %in% gene.id.prom,]

all.probes <- subset(cpg.gene, distance <= 5000)
all.probes <- all.probes[,c("Composite.Element.REF","Gene_Symbol.x","EntrezGene.ID","Associated.Gene.Name")]
all.probes.id <- all.probes$Associated.Gene.Name
all.probes.id <- all.probes[all.probes$Associated.Gene.Name %in% all.probes.id,]  

gene.id.prom.count <- data.frame(gene = character(), count = integer())
gene.id.prom.count <- lapply(gene.id.prom$Associated.Gene.Name, function(i) {
  for(j in 1:length(unique(gene.id.prom$Associated.Gene.Name))){
    x <- nrow(gene.id.prom[gene.id.prom$Associated.Gene.Name == unique(gene.id.prom$Associated.Gene.Name)[j],])
    return(paste(unique(gene.id.prom$Associated.Gene.Name)[j], x, sep="."))
  }                                       
},1)

### CONTINUEI AQUI
require(parallel)

uu <- epiGene[[1]]
clusters <- c(paste0("LGm", 1:6))

enrichment <- lapply(clusters, function(cluster, epiGene.ss=epiGene) {
  enrichment.group <- unlist(mclapply(epiGene.ss, function(uu, cluster.group=cluster) {
    uu$es <- uu$meth.gt0.3 & uu$exp.lt.meanUnmeth
    uu$Kclass <- uu$cluster.meth == cluster.group
    es.table <- table(uu$Kclass, uu$es)#, dnn=c("Class", "ES"))
    if((es.table[2,2]/rowSums(es.table)[2] > 0.5) & (es.table[2,2]/colSums(es.table)[2] > 0.5)) {
      return(fisher.test(table(uu$Kclass, uu$es))$p.value)
    } else {
      return(NULL)
    }
  }, mc.cores=1))
  return(enrichment.group)
})

K4enrichment <- unlist(mclapply(epiGene, function(uu) {
  uu$es <- uu$meth.gt0.3 & uu$exp.lt.meanUnmeth
  uu$K4class <- uu$cluster.meth == "LGm4"
  fisher.test(table(uu$K4class, uu$es))$p.value
}, mc.cores=1))


K3enrichment <- unlist(mclapply(epiGene, function(uu) {
  uu$es <- uu$meth.gt0.3 & uu$exp.lt.meanUnmeth
  uu$K3class <- uu$cluster.meth == "LGm3"
  fisher.test(table(uu$K3class, uu$es))$p.value
}, mc.cores=1))

K2enrichment <- unlist(mclapply(epiGene, function(uu) {
  uu$es <- uu$meth.gt0.3 & uu$exp.lt.meanUnmeth
  uu$K2class <- uu$cluster.meth == "LGm2"
  fisher.test(table(uu$K2class, uu$es))$p.value
}, mc.cores=1))

K1enrichment <- unlist(mclapply(epiGene, function(uu) {
  uu$es <- uu$meth.gt0.3 & uu$exp.lt.meanUnmeth
  uu$K1class <- uu$cluster.meth == "LGm1"
  fisher.test(table(uu$K1class, uu$es))$p.value
}, mc.cores=1))

K5enrichment <- unlist(mclapply(epiGene, function(uu) {
  uu$es <- uu$meth.gt0.3 & uu$exp.lt.meanUnmeth
  uu$K5class <- uu$cluster.meth == "LGm5"
  fisher.test(table(uu$K5class, uu$es))$p.value
}, mc.cores=1))

K6enrichment <- unlist(mclapply(epiGene, function(uu) {
  uu$es <- uu$meth.gt0.3 & uu$exp.lt.meanUnmeth
  uu$K6class <- uu$cluster.meth == "LGm6"
  fisher.test(table(uu$K6class, uu$es))$p.value
}, mc.cores=1))

enrichmentList <- list(LGm1=K1enrichment, LGm2=K2enrichment, LGm3=K3enrichment, LGm4=K4enrichment, LGm5=K5enrichment, LGm6=K6enrichment)
names(enrichmentList) <- c("LGm1", "LGm2", "LGm3", "LGm4", "LGm5", "LGm6")
enrichmentList.adj <- lapply(enrichmentList, p.adjust, method="BH")
enrichmentList.sig <- lapply(enrichmentList.adj, function(x) {
  x <- x[x < 0.05]
})


#### tentativa numero 2
names(enrichment) <- c("LGm1", "LGm2", "LGm3", "LGm4", "LGm5", "LGm6")
enrichment <- lapply(enrichment, p.adjust, method="BH")
tent.enrichmentList.sig <- lapply(enrichment, function(x) {
  x <- x[x < 0.05]
})



info.lgm3 <- setdiff(names(tent.enrichmentList.sig$LGm3),c(names(tent.enrichmentList.sig$LGm4),names(tent.enrichmentList.sig$LGm5)))
probes.lgm3 <- unlist(strsplit(info.lgm3,"[.]"))
odd_indexes<-seq(1,length(probes.lgm3),2)
probes.lgm3 <- probes.lgm3[odd_indexes]
info.lgm4 <- setdiff(names(tent.enrichmentList.sig$LGm4),c(names(tent.enrichmentList.sig$LGm3),names(tent.enrichmentList.sig$LGm5)))
probes.lgm4 <- unlist(strsplit(info.lgm4,"[.]"))
odd_indexes<-seq(1,length(probes.lgm4),2)
probes.lgm4 <- probes.lgm4[odd_indexes]
probes.lgm4 <- probes.lgm4[-79]
info.lgm5 <- setdiff(names(tent.enrichmentList.sig$LGm5),c(names(tent.enrichmentList.sig$LGm3),names(tent.enrichmentList.sig$LGm4)))
probes.lgm5 <- unlist(strsplit(info.lgm5,"[.]"))
odd_indexes<-seq(1,length(probes.lgm5),2)
probes.lgm5 <- probes.lgm5[odd_indexes]
#Hclust in each probe set
hc.lgm3 <- hclust(dist(samples[probes.lgm3,]))
probes.lgm3 <- probes.lgm3[hc.lgm3$order]
hc.lgm4 <- hclust(dist(samples[probes.lgm4,]))
probes.lgm4 <- probes.lgm4[hc.lgm4$order]
hc.lgm5<- hclust(dist(samples[probes.lgm5,]))
probes.lgm5 <- probes.lgm5[hc.lgm5$order]

rlab <- matrix("yellow", nrow = length(probes.lgm5), ncol = 1)
rlab <- c(rlab,matrix("orange", nrow = length(probes.lgm4), ncol = 1))
rlab <- c(rlab,matrix("purple", nrow = length(probes.lgm3), ncol = 1))
probes.es <- c(probes.lgm5,probes.lgm4,probes.lgm3)
heatmap.es.order <- samples[probes.es,cc.out.purity.lgg.gbm.27.450[[6]][[2]][[3]]]

####### Pilocytic Astrocytoma GSE44684
pilocytic.astro.GSE44684 <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Pilocytic_astrocytoma_GSE44684/GSE44684_beta.txt")
#remove the p value columns
index <- grep("Pval",colnames(pilocytic.astro.GSE44684))
pilocytic.astro.GSE44684 <- pilocytic.astro.GSE44684[,-index]
rownames(pilocytic.astro.GSE44684) <- pilocytic.astro.GSE44684[,1]
pilocytic.astro.GSE44684 <- pilocytic.astro.GSE44684[,c(1:38,40:68,39)] #coloca os 6 controles pro final
pilocytic.astro.GSE44684 <- pilocytic.astro.GSE44684[,c(1:62,65:67,63:64,68)] #manter a ordem que esta no GEO para coincidir com os dados de idade
PA.age <- t(read.table("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Pilocytic_astrocytoma_GSE44684/PA_age.txt"))

rownames(normalBrain_dnameth) <- normalBrain_dnameth$ID_REF

### CpG island (downloaded from Encode Aug, 28th)
cpg <- read.table("/dados/ResearchProjects/thais/TCGA/CpGisland.txt")
library(GenomicRanges)
cpg.is <- GRanges(seqnames = cpg$V1,ranges = IRanges(start = cpg$V2, end=cpg$V3), ID = cpg$V4)
aux <- GBM.LGG.27.450k[rownames(dat.lgg.gbm.27.450.noXY.dic.oecg),]
probe <- GRanges(seqnames = paste0("chr",aux$Chromosome.x),ranges = IRanges(start = aux$Genomic_Coordinate.x, end=aux$Genomic_Coordinate.x), ID = aux$Composite.Element.REF)
overlap <- as.data.frame(findOverlaps(probe,cpg.is))
rlab <- matrix("white", nrow = nrow(aux), ncol = 2)
rlab[overlap$queryHits,] <- "green"

normal.age <- read.table("/dados/ResearchProjects/thais/TCGA/LGG.GBM/normal_age.txt",header=T)


GBM.LGG.27.450k.noXY.s <- GBM.LGG.27.450k.noXY[rownames(dat.lgg.gbm.27.450.noXY.dic.oecg) ,5:909]
#GBM.LGG.27.450k.noXY.s <- GBM.LGG.27.450k[rownames(teste),5:909]
lgg.gbm.27.450.order <- GBM.LGG.27.450k.noXY.s[,cc.out.purity.lgg.gbm.27.450[[6]][[2]][[3]]]
metadata <- metadata[colnames(lgg.gbm.27.450.order),]
#clab.6 <- colors.label(as.matrix(a$K6),as.matrix(a$GBM.DNA.Methyl.Clusters),as.matrix(a$LGG.DNA.Methyl.Clusters),as.matrix(a$TumorType),as.matrix(a$age),as.matrix(a$IDH1),as.matrix(a$Cluster),as.matrix(a$RNA.cluster),as.matrix(a$COCCluster),as.numeric(a$leukocyte),as.matrix(a$estimates))
clab.6 <- colors.label(as.matrix(metadata$K6),as.matrix(metadata$DNA.Methylation.Clusters),as.matrix(metadata$MethylationCluster),as.matrix(metadata$tumor.type),as.matrix(metadata$age),as.matrix(metadata$IDH.status),as.matrix(metadata$Cluster),as.matrix(metadata$cluster.mRNA),as.matrix(metadata$COCCluster),as.numeric(metadata$leukocyte),as.matrix(metadata$rna.estimates))
clab.6 <- clab.6[c(856:905,1:333,794:855,334:489,490:708,709:793),]
lgg.gbm.27.450.order <- lgg.gbm.27.450.order[,c(856:905,1:333,794:855,334:489,490:708,709:793)]
#lgg.gbm.27.450.order <- lgg.gbm.27.450.order[heatmap.lgg.gbm$rowInd,]
normals.sel <- norms[rownames(lgg.gbm.27.450.order),]
pilocytic.astro.GSE44684.order <- pilocytic.astro.GSE44684[rownames(lgg.gbm.27.450.order),]
normalBrain_dnameth.order <- normalBrain_dnameth[rownames(lgg.gbm.27.450.order),]
colnames(normalBrain_dnameth.order)[507] <- "UMARY.1648.TCTX"
lgg.gbm.27.450.order <- cbind(normals.sel,lgg.gbm.27.450.order)
lgg.gbm.27.450.order <- cbind(lgg.gbm.27.450.order,pilocytic.astro.GSE44684.order[,-1])
lgg.gbm.27.450.order <- cbind(normalBrain_dnameth.order[,-1],lgg.gbm.27.450.order)
norm.label <- matrix("white", nrow = 77, ncol = 11)
norm.label[,1] <- jet.colors(100)[ceiling((as.numeric(normal.age[,colnames(normals.sel)])-min(as.numeric(normal.age[,colnames(normals.sel)]),na.rm=T))/(max(as.numeric(normal.age[,colnames(normals.sel)]),na.rm=T)-min(as.numeric(normal.age[,colnames(normals.sel)]),na.rm=T))*(length(jet.colors(100))-1))+1]
region.label <- matrix("white", nrow = 506, ncol = 11)
region.label[1:121,11] <- "aquamarine" #cerebelo
region.label[122:254,11] <- "pink" #frontal cortex
region.label[255:380,11] <- "deeppink" #pons 
region.label[381:506,11] <- "brown" #temporal contex
normalBrain.info.order <- normalBrain.info[rownames(normalBrain.info) %in% colnames(normalBrain_dnameth.order)[2:507],]
aux <- data.frame(setdiff(colnames(normalBrain_dnameth.order)[2:507],rownames(normalBrain.info)))
aux$age <- NA
colnames(aux) <- c("ID","Age")
normalBrain.info.order <- rbind(normalBrain.info.order,aux)
rownames(normalBrain.info.order) <- normalBrain.info.order$ID
normalBrain.info.order <- normalBrain.info.order[colnames(normalBrain_dnameth.order)[2:507],]
normalBrain.info.order$Age <- as.numeric(normalBrain.info.order$Age)
region.label[,1] <- jet.colors(100)[ceiling((normalBrain.info.order$Age-min(normalBrain.info.order$Age,na.rm=T))/(max(normalBrain.info.order$Age,na.rm=T)-min(normalBrain.info.order$Age,na.rm=T))*(length(jet.colors(100))-1))+1]
region.label[is.na(region.label)] <- "white"
norm.label[,11] <- "darkblue"
clab.6 <- rbind(norm.label,clab.6)
pa.label <- matrix("white", nrow = 67, ncol = 11)
pa.label[1:61,11] <- "khaki3" #tumor type
pa.label[62:67,11] <- "plum" #tumor type
pa.label[,1] <- jet.colors(100)[ceiling((PA.age-min(PA.age,na.rm=T))/(max(PA.age,na.rm=T)-min(PA.age,na.rm=T))*(length(jet.colors(100))-1))+1]
clab.6 <- rbind(clab.6,pa.label)
clab.6 <- rbind(region.label,clab.6)
colnames(clab.6) <- 1:11
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_905_s/heatmap_K6_normals_leuk4_pilocytic_age.png",res=700,width=20,height=20, units="in")
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_905_s/heatmap_K6_dend.pdf")
heatmap.plus.sm(as.matrix(lgg.gbm.27.450.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                #Rowv = NA,
                ColSideColors = clab.6,
                RowSideColors = rlab
                
)
dev.off()

### New Sample (scatter plot)
newsample <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/jhu-usc.edu_LGG.HumanMethylation450.2.lvl-3.TCGA-DU-6402-01A-11D-1706-05.txt")
newsample <- newsample[,c(1,3,4,5,2)]
colnames(newsample)[5] <- "TCGA-DU-6402-01A-11D-1706-05"
rownames(newsample) <- newsample$Composite.Element.REF
newsample <- newsample[rownames(GBM.LGG.27.450k),]
lgm6 <- subset(metadata, cluster.meth == "LGm6")
lm1 <- subset(lgm6, MethylationCluster == "M1"))
lgm6.s <- GBM.LGG.27.450k[,rownames(lgm6)]
lm1.s <- GBM.LGG.27.450k[,rownames(lm1)]
meanLgm6 <- apply(lgm6.s,1,mean,na.rm=T)
meanLm1 <- apply(lm1.s,1,mean,na.rm=T)
meanLgm6 <- cbind(meanLgm6,newsample[,5])
meanLm1 <- cbind(meanLm1,newsample[,5])
meanLgm6 <- as.data.frame(meanLgm6)
meanLm1 <- as.data.frame(meanLm1)
colnames(meanLm1)[2] <- "V2"


ggplot() +
  
  stat_smooth(data = meanLm1, aes(x = meanLm1, y = V2), method = "glm") +
  geom_point(data = meanLm1, aes(x = meanLm1, y = V2,colour=meanLm1)) +
  # scale_y_log10() +
  ylab("TCGA-DU-6402-01A") +
  xlab("Lm1") 
dev.off()


### Unsupervised analysis (LGG+GBM+PA)
PA.LGG.GBM <- merge(GBM.LGG.27.450k, pilocytic.astro.GSE44684[,2:62], by = 0)
rownames(PA.LGG.GBM) <- PA.LGG.GBM$Row.names
PA.LGG.GBM <- PA.LGG.GBM[,-1]
PA.LGG.GBM.noXY <- PA.LGG.GBM[PA.LGG.GBM$Chromosome.x != "X" & PA.LGG.GBM$Chromosome.x != "Y",]
PA.LGG.GBM.noXY.dic <- PA.LGG.GBM.noXY[ ,5:970]
PA.LGG.GBM.noXY.dic <- PA.LGG.GBM.noXY.dic>0.3
storage.mode(PA.LGG.GBM.noXY.dic) <- "numeric"
## Use probes methylated more than 10% of the tumor samples
PA.LGG.GBM.noXY.dic <- PA.LGG.GBM.noXY.dic[(rowSums(PA.LGG.GBM.noXY[,5:970])/ncol(PA.LGG.GBM.noXY[,5:970]))*100 >10,]

FDb.HM450.oecg.3kb.s <- subset(FDb.HM450.oecg.3kb, oecg.3kb>1.25)
PA.LGG.GBM.noXY.dic.oecg <- PA.LGG.GBM.noXY.dic[rownames(PA.LGG.GBM.noXY.dic) %in% (FDb.HM450.oecg.3kb.s$probeid),]
avg.PA.LGG.GBM.noXY <- apply(PA.LGG.GBM.noXY[,5:970], 1, mean, na.rm=TRUE)
avg.norm.PA.LGG.GBM.noXY <- merge(as.data.frame(avg.PA.LGG.GBM.noXY), avg.normals, by = 0, all.x = T)
library(LSD)
comparisonplot(avg.norm.PA.LGG.GBM.noXY$avg.normals, avg.norm.PA.LGG.GBM.noXY$avg.PA.LGG.GBM.noXY, pimp=TRUE)
PA.LGG.GBM.tumor.specific <- (subset(avg.norm.PA.LGG.GBM.noXY, avg.normals < 0.3))

dat.PA.LGG.GBM.noXY.dic.oecg <- PA.LGG.GBM.noXY.dic.oecg[(rownames(PA.LGG.GBM.noXY.dic.oecg) %in% PA.LGG.GBM.tumor.specific$Row.names), ]

library(ConsensusClusterPlus)
title="/dados/ResearchProjects/thais/TCGA/LGG.GBM/ConsensusCluster.PA.LGG.GBM"
setwd(title)
cc.out.purity.PA.LGG.GBM <- ConsensusClusterPlus(d=dat.PA.LGG.GBM.noXY.dic.oecg, 
                                                 maxK=10, 
                                                 reps=1000, 
                                                 pItem=0.8, 
                                                 pFeature=1,
                                                 title=title, 
                                                 clusterAlg="hc",
                                                 distance="binary",
                                                 innerLinkage = "ward", 
                                                 finalLinkage = "ward",
                                                 seed=1022,
                                                 plot="pdf",
                                                 #writeTable=TRUE,
                                                 verbose = TRUE)
ccICL = calcICL(cc.out.purity.PA.LGG.GBM, plot='pdf')
cc.out.purity.PA.LGG.GBM.2 <- cbind(as.matrix(cc.out.purity.PA.LGG.GBM[[2]][[3]]), 
                                    as.matrix(cc.out.purity.PA.LGG.GBM[[3]][[3]]), 
                                    as.matrix(cc.out.purity.PA.LGG.GBM[[4]][[3]]), 
                                    as.matrix(cc.out.purity.PA.LGG.GBM[[5]][[3]]), 
                                    as.matrix(cc.out.purity.PA.LGG.GBM[[6]][[3]]),
                                    as.matrix(cc.out.purity.PA.LGG.GBM[[7]][[3]]),
                                    as.matrix(cc.out.purity.PA.LGG.GBM[[8]][[3]]),
                                    as.matrix(cc.out.purity.PA.LGG.GBM[[9]][[3]]),
                                    as.matrix(cc.out.purity.PA.LGG.GBM[[10]][[3]]))

cc.out.purity.PA.LGG.GBM.2 <- as.data.frame(cc.out.purity.PA.LGG.GBM.2)
names(cc.out.purity.PA.LGG.GBM.2) <- c("K2",  "K3",  "K4" , "K5" , "K6" , "K7",  "K8" , "K9" , "K10")

cc.out.purity.PA.LGG.GBM.2$TCGAIDnew <- substr(rownames(cc.out.purity.PA.LGG.GBM.2),1,12)
cc.out.purity.PA.LGG.GBM.2$TCGAID <- rownames(cc.out.purity.PA.LGG.GBM.2)

PA.LGG.GBM.anno <- merge(mutation.clinical.cc, brennan.class.s[,c("DNA.Methylation.Clusters","TCGAIDnew")], by = "TCGAIDnew", all.x =T)
rownames(PA.LGG.GBM.anno) <- PA.LGG.GBM.anno$TCGAID
xPA.LGG.GBM.anno <- merge(PA.LGG.GBM.anno, clinical.cluster.genetics[,c("Tumor" ,"MethylationCluster","COCCluster")], by.x = "TCGAIDnew", by.y = "Tumor", all.x = T)
rownames(xPA.LGG.GBM.anno) <- xPA.LGG.GBM.anno$TCGAID
xPA.LGG.GBM.anno <- merge(xPA.LGG.GBM.anno,gbm.molecular.subtype[,c("sample.id","Cluster")], by.x = "TCGAIDnew", by.y = "sample.id", all.x = T)
rownames(xPA.LGG.GBM.anno) <- xPA.LGG.GBM.anno$TCGAID

xxMODr <- merge(xPA.LGG.GBM.anno,leukocyte, by.x=0, by.y= "x",all.x=T)
rownames(xxMODr) <- xxMODr$Row.names
xxMODr <- xxMODr[,-1]
leukocytecat <- cut(xxMODr$leukocyte, c(0.0136,0.10340341,0.15807936,0.24441133,0.77))
xxMODr$leukocyte.strat <- leukocytecat
metadata.PA <- merge(xxMODr,rna.estimates,by.x="TCGAIDnew",by.y=0,all.x=T)
colnames(metadata.PA)[66] <- "rna.estimates"
rownames(metadata.PA) <- metadata.PA$TCGAID
aux <- rownames(subset(metadata.PA,is.na(tumor.type)))
metadata.PA[aux,"tumor.type"] <- "PA"
metadata.PA$cluster.mrna.new <- NA
metadata.PA[metadata.PA$cluster.mRNA %in% "yellow","cluster.mrna.new"] <- "LGr1"
metadata.PA[metadata.PA$cluster.mRNA %in% "gray","cluster.mrna.new"] <- "LGr2"
metadata.PA[metadata.PA$cluster.mRNA %in% "cyan","cluster.mrna.new"] <- "LGr3"
metadata.PA[metadata.PA$cluster.mRNA %in% "magenta","cluster.mrna.new"] <- "LGr4"
metadata.PA[metadata.PA$cluster.mRNA %in% "blue","cluster.mrna.new"] <- "LGr5"
metadata.PA[metadata.PA$cluster.mRNA %in% "green3","cluster.mrna.new"] <- "LGr6"
metadata.PA[metadata.PA$cluster.mRNA %in% "red","cluster.mrna.new"] <- "LGr7"
metadata.PA[metadata.PA$cluster.mRNA %in% "brown","cluster.mrna.new"] <- "LGr8"


##metadata.PA <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/ConsensusCluster.PA.LGG.GBM/metadata.txt")
rownames(normalBrain_dnameth) < normalBrain_dnameth$ID_REF

PA.LGG.GBM.noXY.s <- PA.LGG.GBM.noXY[rownames(dat.PA.LGG.GBM.noXY.dic.oecg) ,5:970]
PA.LGG.GBM.noXY.order <- PA.LGG.GBM.noXY.s[,cc.out.purity.PA.LGG.GBM[[6]][[2]][[3]]]
metadata.PA <- metadata.PA[colnames(PA.LGG.GBM.noXY.order),]
clab.6 <- colors.label(as.matrix(metadata.PA$K6),as.matrix(metadata.PA$DNA.Methylation.Clusters),as.matrix(metadata.PA$MethylationCluster),as.matrix(metadata.PA$tumor.type),as.matrix(metadata.PA$age),as.matrix(metadata.PA$IDH.status),as.matrix(metadata.PA$Cluster),as.matrix(metadata.PA$cluster.mRNA),as.matrix(metadata.PA$COCCluster),as.numeric(metadata.PA$leukocyte),as.matrix(metadata.PA$rna.estimates))
#clab.6 <- clab.6[c(856:905,1:333,794:855,334:489,490:708,709:793),]
clab.6 <- clab.6[c(664:823,379:663,202:378,1:201,887:966,824:886),]
PA.LGG.GBM.noXY.order <- PA.LGG.GBM.noXY.order[,c(664:823,379:663,202:378,1:201,887:966,824:886)]
PA.LGG.GBM.noXY.order <- PA.LGG.GBM.noXY.order[heatmap.order$rowInd,]
normalBrain_dnameth.order <- normalBrain_dnameth[rownames(PA.LGG.GBM.noXY.order),]
PA.LGG.GBM.noXY.order <- cbind(normalBrain_dnameth.order[,-1],PA.LGG.GBM.noXY.order)
region.label <- matrix("white", nrow = 506, ncol = 11)
region.label[1:121,11] <- "aquamarine" #cerebelo
region.label[122:254,11] <- "pink" #frontal cortex
region.label[255:380,11] <- "deeppink" #pons 
region.label[381:506,11] <- "brown" #temporal contex
clab.6 <- rbind(region.label,clab.6)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_905_s/heatmap_K6_PA.png",res=700,width=2000,height=2000)
heatmap.plus.sm(as.matrix(PA.LGG.GBM.noXY.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab.6
)
dev.off()


###Normal
normalBrain_dnameth.noXY <- normalBrain_dnameth[rownames(GBM.LGG.27.450k.noXY),]
normalBrain_dnameth.noXY$SD <- apply(normalBrain_dnameth.noXY[,2:506],1,sd,na.rm=T)
normalBrain_dnameth.noXY.sub <- subset(normalBrain_dnameth.noXY, SD > 0.13)
GBM.LGG.27.450k.noXY.s <- GBM.LGG.27.450k[rownames(normalBrain_dnameth.noXY.sub),5:909]
lgg.gbm.27.450.order <- GBM.LGG.27.450k.noXY.s[,cc.out.purity.lgg.gbm.27.450[[6]][[2]][[3]]]
metadata.PA <- metadata.PA[colnames(lgg.gbm.27.450.order),]
#clab.6 <- colors.label(as.matrix(a$K6),as.matrix(a$GBM.DNA.Methyl.Clusters),as.matrix(a$LGG.DNA.Methyl.Clusters),as.matrix(a$TumorType),as.matrix(a$age),as.matrix(a$IDH1),as.matrix(a$Cluster),as.matrix(a$RNA.cluster),as.matrix(a$COCCluster),as.numeric(a$leukocyte),as.matrix(a$estimates))
clab.6 <- colors.label(as.matrix(metadata.PA$K6),as.matrix(metadata.PA$DNA.Methylation.Clusters),as.matrix(metadata.PA$MethylationCluster),as.matrix(metadata.PA$tumor.type),as.matrix(metadata.PA$age),as.matrix(metadata.PA$IDH.status),as.matrix(metadata.PA$Cluster),as.matrix(metadata.PA$cluster.mRNA),as.matrix(metadata.PA$COCCluster),as.numeric(metadata.PA$leukocyte),as.matrix(metadata.PA$rna.estimates))
clab.6 <- clab.6[c(856:905,1:333,794:855,334:489,490:708,709:793),]
lgg.gbm.27.450.order <- lgg.gbm.27.450.order[,c(856:905,1:333,794:855,334:489,490:708,709:793)]
#lgg.gbm.27.450.order <- lgg.gbm.27.450.order[heatmap.lgg.gbm$rowInd,]
#normals.sel <- normals.sel[rownames(lgg.gbm.27.450.order),]
pilocytic.astro.GSE44684.order <- pilocytic.astro.GSE44684[rownames(lgg.gbm.27.450.order),]
normalBrain_dnameth.order <- normalBrain_dnameth[rownames(lgg.gbm.27.450.order),2:507]
#lgg.gbm.27.450.order <- cbind(normals.sel,lgg.gbm.27.450.order)
lgg.gbm.27.450.order <- cbind(lgg.gbm.27.450.order,pilocytic.astro.GSE44684.order[,-1])
lgg.gbm.27.450.order <- cbind(normalBrain_dnameth.order,lgg.gbm.27.450.order)
#norm.label <- matrix("white", nrow = 77, ncol = 11)
region.label <- matrix("white", nrow = 506, ncol = 11)
region.label[1:121,11] <- "aquamarine" #cerebelo
region.label[122:254,11] <- "pink" #frontal cortex
region.label[255:380,11] <- "deeppink" #pons 
region.label[381:506,11] <- "brown" #temporal contex
#norm.label[,11] <- "darkblue"
#clab.6 <- rbind(norm.label,clab.6)
pa.label <- matrix("white", nrow = 67, ncol = 11)
pa.label[1:62,11] <- "khaki3" #tumor type
pa.label[63:67,11] <- "plum" #tumor type
clab.6 <- rbind(clab.6,pa.label)
clab.6 <- rbind(region.label,clab.6)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_905_s/heatmap_K6_normals_leuk4_pilocytic_region_toshi2.png",res=600,width=2000,height=2000)
heatmap.plus.sm(as.matrix(lgg.gbm.27.450.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                #Rowv = NA,
                ColSideColors = clab.6
                
)
dev.off()

#identificar as probes 450k que estao sob os genes do PA paper
PA.paper <- read.table("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Pilocytic_astrocytoma_GSE44684/PA_paper.txt",header=T) #hg19
PA.gene <- GRanges(seqnames = paste0("chr",PA.paper$Chromosome),ranges = IRanges(start = (PA.paper$Start.coordinate-2000), end=(PA.paper$End.coordinate+2000)), symbol = PA.paper$Gene.name)
probe.info <- GRanges(seqnames = paste0("chr",GBM.LGG.27.450k$Chromosome.x), ranges = IRanges(start = GBM.LGG.27.450k$Genomic_Coordinate.x, end = GBM.LGG.27.450k$Genomic_Coordinate.x), probeID = GBM.LGG.27.450k$Composite.Element.REF)
overlap <- as.data.frame(findOverlaps(PA.gene,probe.info)) #?overlap
PA.228.genes <- as.character(GBM.LGG.27.450k[overlap$subjectHits,"Composite.Element.REF"])
normal.like.id <- metadata[metadata$MethylationCluster %in% "M1","TCGAID"]
id <- metadata[metadata$cluster.meth %nin% "LGm6","TCGAID"]
normal.like <- GBM.LGG.27.450k[PA.228.genes,normal.like.id]
pa.samples <- pilocytic.astro.GSE44684[PA.228.genes,2:62]
other.samples <- GBM.LGG.27.450k[PA.228.genes,colnames(GBM.LGG.27.450k) %in% id]
other.samples <- other.samples[,-c(1:4)]
other.samples$mean <- apply(other.samples,1,mean,na.rm=T)
crbl <- normalBrain_dnameth[PA.228.genes,2:121]
crbl$mean <- apply(crbl,1,mean,na.rm=T)
pa.samples$mean <- apply(pa.samples,1,mean,na.rm=T)
normal.like$mean <- apply(normal.like,1,mean,na.rm=T)
#lgm6 <- subset(metadata, MethylationCluster != "M1" & cluster.meth == "LGm6")
lgm6 <- subset(metadata, cluster.meth == "LGm6")
lgm6 <- subset(lgm6, MethylationCluster != "M1" | is.na(MethylationCluster))
lgm6.s <- GBM.LGG.27.450k[PA.228.genes,rownames(lgm6)]
lgm6.s$mean <- apply(lgm6.s,1,mean,na.rm=T)

library(ggplot2)
x <- data.frame(mean = c(crbl$mean,normal.like$mean,pa.samples$mean,other.samples$mean,lgm6.s$mean))
x$type <- NA
x$type[1:332] <- "Cerebellum"
x$type[333:664] <- "Lm1"
x$type[665:996] <- "PA"
x$type[997:1328] <- "LGm1+LGm2+LGm3+LGm4+LGm5"
x$type[1329:1660] <- "LGm6"
#colnames(x) <- c("Cerebellum", "Lm1","PA","LGm1+LGm2+LGm3+LGm4+LGm5")
png("/dados/ResearchProjects/thais/TCGA/LGG.GBM/meanMethyGroups.png",res=600,width=20,height=15,units = "in")
p <- ggplot(x, aes(factor(type), mean))
p + geom_boxplot(aes(fill = factor(type)), notchwidth=0.25) + 
  geom_jitter(height = 0, position = position_jitter(width = .1), size=2)  + 
  # scale_fill_manual(values=c("green","red","purple", "orange", "yellow","blue"
  #), labels = c("LGm1", "LGm2", "LGm3","LGm4", "LGm5", "LGm6"
  #)) +
  #scale_fill_manual(values = rev(c("darkorchid","darkgreen","blue","red","DARKGREY"))) +
  ylab("Mean Methylation Values") +
  xlab("Groups") + 
  labs(title = "Mean Methylation Values By Groups", fill = "Groups") +
  theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white"), panel.border = element_rect(colour = "white"), legend.position="none")
dev.off()


avg.probes.cl <- cbind(crbl$mean,normal.like$mean,pa.samples$mean,other.samples$mean,lgm6.s$mean)
colnames(avg.probes.cl) <- c("Cerebellum","Lm1","PA","LGm1+LGm2+LGm3+LGm4+LGm5","LGm6-Lm1")
png("/dados/ResearchProjects/thais/TCGA/LGG.GBM/PAxEverythingElse.png",res=200,width=2000,height=2000)
heatpairs(avg.probes.cl)
dev.off()

## Probes from LGG paper [NEJM, 2014]
load("/dados/ResearchProjects/thais/LGG.puritycorrection.HeatmapOrder.rows.ID.Rda")
aux1 <- GBM.LGG.27.450k[rownames(GBM.LGG.27.450k) %in% LGG.puritycorrection.HeatmapOrder.rows.ID,5:909]
aux2 <- pilocytic.astro.GSE44684[rownames(aux1),2:62]
pa.lgg.gbm <- cbind(aux1,aux2)
pa.lgg.gbm.dic <- pa.lgg.gbm>0.3
storage.mode(pa.lgg.gbm.dic) <- "numeric"
library(ConsensusClusterPlus)
title="/dados/ResearchProjects/thais/TCGA/LGG.GBM/ConsensusCluster.LGG.GBM.PA"
setwd(title)
cc.out.purity.lgg.gbm.pa <- ConsensusClusterPlus(d=pa.lgg.gbm.dic, 
                                                 maxK=10, 
                                                 reps=1000, 
                                                 pItem=0.8, 
                                                 pFeature=1,
                                                 title=title, 
                                                 clusterAlg="hc",
                                                 distance="binary",
                                                 innerLinkage = "ward", 
                                                 finalLinkage = "ward",
                                                 seed=1022,
                                                 plot="pdf",
                                                 #writeTable=TRUE,
                                                 verbose = TRUE)

ccICL = calcICL(cc.out.purity.lgg.gbm.pa, plot='pdf')
cc.out.purity.lgg.gbm.pa.2 <- cbind(as.matrix(cc.out.purity.lgg.gbm.pa[[2]][[3]]), 
                                    as.matrix(cc.out.purity.lgg.gbm.pa[[3]][[3]]), 
                                    as.matrix(cc.out.purity.lgg.gbm.pa[[4]][[3]]), 
                                    as.matrix(cc.out.purity.lgg.gbm.pa[[5]][[3]]), 
                                    as.matrix(cc.out.purity.lgg.gbm.pa[[6]][[3]]),
                                    as.matrix(cc.out.purity.lgg.gbm.pa[[7]][[3]]),
                                    as.matrix(cc.out.purity.lgg.gbm.pa[[8]][[3]]),
                                    as.matrix(cc.out.purity.lgg.gbm.pa[[9]][[3]]),
                                    as.matrix(cc.out.purity.lgg.gbm.pa[[10]][[3]]))

cc.out.purity.lgg.gbm.pa.2 <- as.data.frame(cc.out.purity.lgg.gbm.pa.2)
names(cc.out.purity.lgg.gbm.pa.2) <- c("PA2",  "PA3",  "PA4" , "PA5" , "PA6" , "PA7",  "PA8" , "PA9" , "PA10")

aux <- cc.out.purity.lgg.gbm.pa.2
aux[,10:276] <- NA
colnames(aux)[10:276] <- colnames(metadata)
teste <- data.frame(lapply(metadata, as.character), stringsAsFactors=FALSE)
rownames(teste) <- teste$TCGAID
aux[rownames(metadata),10:276] <- teste
aux[is.na(aux$tumor.type),"tumor.type"] <- "PA"
#aux$cluster.meth <- NA
#aux[rownames(metadata),"cluster.meth"] <- metadata$cluster.meth


#mudar a funcao para tumor type e adicionar cluster.meth
pa.lgg.gbm.order <- pa.lgg.gbm[,cc.out.purity.lgg.gbm.pa[[7]][[2]][[3]]]
aux <- aux[colnames(pa.lgg.gbm.order),]
clab.6 <- colors.label(as.matrix(aux$cluster.meth),as.matrix(aux$PA7),as.matrix(aux$DNA.Methylation.Clusters),as.matrix(aux$MethylationCluster),as.matrix(aux$tumor.type),as.matrix(aux$age),as.matrix(aux$IDH.status),as.matrix(aux$Cluster),as.matrix(aux$cluster.mRNA),as.matrix(aux$COCCluster),as.numeric(aux$leukocyte),as.matrix(aux$rna.estimates))
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_905_s/heatmap_K7_PA_LGG_GBM.png",res=200,width=2000,height=2000)
heatmap.plus.sm(as.matrix(pa.lgg.gbm.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                #Rowv = NA,
                ColSideColors = clab.6
)
dev.off()

### Polycomb
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/polycomb.Infinium.Ku.ESC.rda")
dim(Ku.ESC.polycomb)


####### mRNA for GSE15745 (control)
normalBrain_mrna <- read.delim("/dados/ResearchProjects/thais/TCGA/Normals/mRNA_Brain/07-03-14/GSE15745-GPL6104_series_matrix.txt",skip=95)
normalBrain_mrna_id <- read.delim("/dados/ResearchProjects/thais/TCGA/Normals/mRNA_Brain/07-03-14/id_mrna.txt")
colnames(normalBrain_mrna)[2:585] <- colnames(normalBrain_mrna_id)
normalBrain_dnameth <- read.delim("/dados/ResearchProjects/thais/TCGA/Normals/mRNA_Brain/07-03-14/GSE15745-GPL8490_series_matrix.txt",skip=96)
normalBrain_dnameth_id <- read.delim("/dados/ResearchProjects/thais/TCGA/Normals/mRNA_Brain/07-03-14/id_dnameth.txt") 
colnames(normalBrain_dnameth)[2:507] <- colnames(normalBrain_dnameth_id)
colnames(normalBrain_dnameth)[507] <- "UMARY.1648.TCTX"
#Selecionar as amostras em que ha informacao de gene exp e dna methylation
normalBrain_mrna.sel <- normalBrain_mrna[,intersect(colnames(normalBrain_mrna)[2:585],colnames(normalBrain_dnameth)[2:507])]
normalBrain_dnameth.sel <- normalBrain_dnameth[,intersect(colnames(normalBrain_mrna)[2:585],colnames(normalBrain_dnameth)[2:507])]
#colocar a informacao da probe como primeira coluna
normalBrain_dnameth.sel <- cbind(normalBrain_dnameth[,1],normalBrain_dnameth.sel)
normalBrain_mrna.sel <- cbind(normalBrain_mrna[,1],normalBrain_mrna.sel)
colnames(normalBrain_mrna.sel)[1] <- "Probe.id"
colnames(normalBrain_dnameth.sel)[1] <- "Probe.id"

#age
normalBrain.info <- read.table("/dados/ResearchProjects/thais/TCGA/Normals/mRNA_Brain/07-03-14/normalBrainMetadata.txt")
normalBrain.info <- t(normalBrain.info)
colnames(normalBrain.info) <- normalBrain.info[1,]
normalBrain.info <- as.data.frame(normalBrain.info[-1,])
rownames(normalBrain.info) <- normalBrain.info$ID


### to map the probe id to a gene, see file /dados/scripts/sandiego/mapToGenome.R
load("/dados/ResearchProjects/thais/TCGA/Normals/mRNA_Brain/GeneExpMap.rda") #obj final

#selecionar as mesmas probes para dna methylation
rownames(normalBrain_dnameth.sel) <- normalBrain_dnameth.sel[,1]
normalBrain_dnameth.final <- normalBrain_dnameth.sel[rownames(GBM.LGG.27.450k.nearest.gene.sel),]
normalBrain_dnameth.final <- cbind(GBM.LGG.27.450k.nearest.gene.sel$EntrezGene.ID,normalBrain_dnameth.final)
colnames(normalBrain_dnameth.final)[1] <- "EntrezGene.ID"

#selecionar para esta plataforma os mesmos genes que foram usados para fazer a analise com os dados do TCGA
normalBrain_mrna.final <- merge(final[,c(1,5,6)],normalBrain_mrna.sel,by.x="qName",by.y="Probe.id")

cpg.gene.normal <- merge(na.omit(normalBrain_dnameth.final),normalBrain_mrna.final,by.x="EntrezGene.ID",by.y="EntrezGene.ID")
dimnames(cpg.gene.normal)[[1]] <- paste(paste(cpg.gene.normal[,"Probe.id"], cpg.gene.normal[,"Associated.Gene.Name"], sep="."),cpg.gene.normal[,"qName"],sep=".") #o id da probe de exp tambem teve que ser colocado porque mais de uma probe para a expressao analisa um mesmo gene



meth.g.norm <- cpg.gene.normal[,3:493]
colnames(meth.g.norm) <- substr(colnames(meth.g.norm), 1, (nchar(names(meth.g.norm))-3))
expr.g.norm <- cpg.gene.normal[,496:986]
colnames(expr.g.norm) <- substr(colnames(expr.g.norm), 1, (nchar(names(expr.g.norm))-3))
epiGene.normal <- list()
j=1

for(i in 1:dim(cpg.gene.normal)[1]) {
  meth.gu <- names(meth.g)[meth.g.norm[i, ] < 0.3]
  meth.gm <- names(meth.g.norm)[meth.g.norm[i, ] >= 0.3]
  # if(fisher.test(table(yy))$p.value < 0.05){
  if(length(meth.gu)>1 & length(meth.gm)>1){
    if(length(meth.gu)>1){
      expr.mean.gm <- mean(as.numeric(expr.g.norm[i, meth.gm]), na.rm=T)
      expr.mean.gu <- mean(as.numeric(expr.g.norm[i, meth.gu]), na.rm = T)
      yy <- cbind(t(meth.g.norm[i, ] >= 0.3), t(expr.g.norm[i, ] < expr.mean.gu))
      yy <- as.data.frame(yy)
      gene <- colnames(yy)[1]
      colnames(yy) <- c("meth.gt0.3", "exp.lt.meanUnmeth")
      yy$geneName <- gene
      yy$id <- rownames(yy)
      yy$methValues <- as.numeric(meth.g.norm[i, ])
      yy$exprValues <- as.numeric(expr.g.norm[i, ])    
      #yy <- merge(metadata, yy, by.x = "TCGAIDnew", by.y = "id")
      epiGene.normal[[j]] <- as.data.frame(yy)
      names(epiGene.normal)[[j]] <- as.character(gene)
      j <- j+1
    }else{
      cat("no methylated group", i, "\n")
    }
  }
} 
  
  aux <- unlist(lapply(c(names(tent.enrichmentList.sig$LGm3),names(tent.enrichmentList.sig$LGm4),names(tent.enrichmentList.sig$LGm5)), function (x) { grep(x,names(epiGene.normal)) }))
  index <- unlist(lapply(c(names(tent.enrichmentList.sig$LGm3),names(tent.enrichmentList.sig$LGm4),names(tent.enrichmentList.sig$LGm5)), function (x) { grep(x,rownames(cpg.gene.normal)) }))
  gene.normal.subset <- cpg.gene.normal[index,]
  epiGene.subset <- lapply(c(names(tent.enrichmentList.sig$LGm3),names(tent.enrichmentList.sig$LGm4),names(tent.enrichmentList.sig$LGm5)), function(x) {epiGene[[x]]})
  
  #normalize the expression data
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
  
  #unique genes for each group (method 1)
  ordem <- sort(enrichmentList.sig$LGm3)
  gene.names.LGm3 <- unlist(strsplit(names(ordem),"[.]"))
  even_indexes<-seq(2,length(gene.names.LGm3),3)
  gene.names.LGm3 <- data.frame(gene.names.LGm3[even_indexes])
  gene.names.LGm3$P.value <-  as.vector(ordem)
  write.table(gene.names.LGm3,file="gene.names.LGG.txt",quote=F,row.names=F,sep="\t")
  
  ordem <- c(enrichmentList.sig$LGm4,enrichmentList.sig$LGm5)
  ordem <- sort(ordem)
  gene.names.LGm4 <- unlist(strsplit(names(ordem),"[.]"))
  even_indexes<-seq(2,length(gene.names.LGm4),3)
  gene.names.LGm4 <- data.frame(gene.names.LGm4[even_indexes])
  gene.names.LGm4$P.value <-  as.vector(ordem)
  
  write.table(gene.names.LGm4,file="gene.names.GBM.txt",quote=F,row.names=F,sep="\t")
  
  
  
  #unique genes for each group (method 2)
  gene.names.LGm2 <- unlist(strsplit(names(enrichmentList.sig$LGm6)[1:1218],"[.]"))
  even_indexes<-seq(2,length(gene.names.LGm2),2)
  gene.names.LGm2 <- gene.names.LGm2[even_indexes]
  gene.names.LGm2 <- unique(gene.names.LGm2)
  gene.names.LGm2[988] <- "AC091878.1"
  gene.names.LGm2.a <- unlist(strsplit(names(enrichmentList.sig$LGm6)[1220:1261],"[.]"))
  even_indexes<-seq(2,length(gene.names.LGm2.a),2)
  gene.names.LGm2.a <- gene.names.LGm2.a[even_indexes]
  gene.names.LGm2.a <- unique(gene.names.LGm2.a)
  write.table(c(gene.names.LGm2,gene.names.LGm2.a),file="gene.names.LGm6.txt",quote=F,row.names=F,sep="\t")
  
  
  names(enrichmentList.sig$LGm3)[!(names(enrichmentList.sig$LGm3) %in% c(names(enrichmentList.sig$LGm6),names(enrichmentList.sig$LGm1),names(enrichmentList.sig$LGm2),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm5)))]                     
  
  #unique probes for each group
  
  info.lgm3 <- setdiff(names(enrichmentList.sig$LGm3),c(names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm5)))
  probes.lgm3 <- unlist(strsplit(info.lgm3,"[.]"))
  odd_indexes<-seq(2,length(probes.lgm3),3)
  probes.lgm3 <- probes.lgm3[odd_indexes]
  info.lgm4 <- setdiff(names(enrichmentList.sig$LGm4),c(names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm3)))
  probes.lgm4 <- unlist(strsplit(info.lgm4,"[.]"))
  odd_indexes<-seq(1,length(probes.lgm4),3)
  probes.lgm4 <- probes.lgm4[odd_indexes]
  info.lgm5 <- setdiff(names(enrichmentList.sig$LGm5),c(names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm3)))
  probes.lgm5 <- unlist(strsplit(info.lgm5,"[.]"))
  odd_indexes<-seq(1,length(probes.lgm5),3)
  probes.lgm5 <- probes.lgm5[odd_indexes]
  
  #Hclust in each probe set
  samples <- GBM.LGG.27.450k[,5:909]
  hc.lgm3 <- hclust(dist(samples[probes.lgm3,]))
  probes.lgm3 <- probes.lgm3[hc.lgm3$order]
  hc.lgm4 <- hclust(dist(samples[probes.lgm4,]))
  probes.lgm4 <- probes.lgm4[hc.lgm4$order]
  
  rlab <- matrix("yellow", nrow = length(probes.lgm5), ncol = 1)
  rlab <- c(rlab,matrix("orange", nrow = length(probes.lgm4), ncol = 1))
  rlab <- c(rlab,matrix("purple", nrow = length(probes.lgm3), ncol = 1))
  probes.es <- c(probes.lgm5,probes.lgm4,probes.lgm3)
  heatmap.es.order <- samples[probes.es,cc.out.purity.lgg.gbm.27.450[[6]][[2]][[3]]]
  
  
  
  metadata <- metadata[colnames(heatmap.es.order),]
  clab.6 <- colors.label(as.matrix(metadata$K6),as.matrix(metadata$DNA.Methylation.Clusters),as.matrix(metadata$MethylationCluster),as.matrix(metadata$tumor.type),as.matrix(metadata$age),as.matrix(metadata$IDH.status),as.matrix(metadata$Cluster),as.matrix(metadata$cluster.mRNA),as.matrix(metadata$COCCluster),as.numeric(metadata$leukocyte),as.matrix(metadata$rna.estimates))
  clab.6 <- clab.6[c(856:905,1:333,794:855,334:489,490:708,709:793),]
  heatmap.es.order<- heatmap.es.order[,c(856:905,1:333,794:855,334:489,490:708,709:793)]
  #normals.sel <- normals.sel[rownames(heatmap.es.order),]
  pilocytic.astro.GSE44684.order <- pilocytic.astro.GSE44684[rownames(heatmap.es.order),]
  normalBrain_dnameth.order <- normalBrain_dnameth[rownames(heatmap.es.order),]
  #heatmap.es.order <- cbind(normals.sel,heatmap.es.order)
  heatmap.es.order <- cbind(heatmap.es.order,pilocytic.astro.GSE44684.order[,-1])
  heatmap.es.order <- cbind(normalBrain_dnameth.order[,-1],heatmap.es.order)
  #norm.label <- matrix("white", nrow = 77, ncol = 11)
  region.label <- matrix("white", nrow = 506, ncol = 11)
  region.label[1:121,11] <- "aquamarine" #cerebelo
  region.label[122:254,11] <- "pink" #frontal cortex
  region.label[255:380,11] <- "deeppink" #pons 
  region.label[381:506,11] <- "brown" #temporal contex
  norm.label[,11] <- "darkblue"
  # clab.6 <- rbind(norm.label,clab.6)
  pa.label <- matrix("white", nrow = 67, ncol = 11)
  pa.label[1:62,11] <- "khaki3" #tumor type
  pa.label[63:67,11] <- "plum" #tumor type
  clab.6 <- rbind(clab.6,pa.label)
  clab.6 <- rbind(region.label,clab.6)
  png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_905_s/heatmap_ES5_normal.png",res=400,width=2000,height=2000)
  heatmap.plus.sm(as.matrix(heatmap.es.order),
                  col = jet.colors(75),
                  scale = "none",
                  labRow = NA,
                  labCol = NA,
                  Colv = NA,
                  #Rowv = NA,
                  ColSideColors = clab.6,
                  RowSideColors = cbind(rlab,rlab)
  )
  dev.off()
  
  
  
  gene.names.LGm2.a <- unlist(strsplit(names(enrichmentList.sig$LGm6),"[.]"))
  even_indexes<-seq(2,length(gene.names.LGm2.a),2)
  gene.names.LGm2.a <- gene.names.LGm2.a[even_indexes]
  gene.names.LGm2.a <- unique(gene.names.LGm2.a)
  write.table(gene.names.LGm2.a,file="gene.names.LGm5.txt",quote=F,row.names=F,sep="\t")
  
  K4enrichment.count <- unlist(mclapply(epiGene, function(uu) {
    uu$es <- uu$meth.gt0.3 & uu$exp.lt.meanUnmeth
    uu$K5class <- uu$cluster.meth == "LGm5"
    as.data.frame(table(uu$K4class, uu$es))[4, 3] > 32 & as.data.frame(table(uu$K4class, uu$es))[3, 3] < 5
  }, mc.cores=1))
  
  ss <- K4enrichment.count[K4enrichment.count==TRUE]
  ss <- str_split(as.character(names(ss)), pattern="\\.")
  ss <- do.call(rbind.data.frame, ss)
  ss.count <- unlist(lapply(split(ss[, 2],f=ss[, 2]),length))
  ss.count <- sort(ss.count)
  dimnames(enrichmentList.geneNames[[4]])[[1]] <- paste(enrichmentList.geneNames[[4]][, 1], enrichmentList.geneNames[[4]][, 2], sep=".")
  rownames(ss) <- paste(ss[,1], ss[,2], sep=".")
  ss.m <- merge(ss, enrichmentList.geneNames[[4]], by = 0, all.x = T)
  ss.m <- na.omit(ss.m)
  
  
  ss.m.count <- unlist(lapply(split(ss.m[, 5],f=ss.m[, 5]),length))
  
  
  
  epiplot <- epiGene.2[["cg03079681.CDKN2A.ILMN_1717714"]]
  #epiplot$exprValues <- epiplot$exprValues + 1
  expr.mean.gu <- mean(as.numeric(epiplot[epiplot$meth.gt0.3==FALSE, "exprValues"]), na.rm = T)
  
  ## K4 specific
  p <- ggplot(epiplot, aes(methValues, exprValues))
  p + 
    stat_smooth(method = "lm", fullrange=FALSE, na.rm=TRUE, se=FALSE, size=1) +
    geom_vline(xintercept = 0.3, linetype = "longdash", colour="black") +
    geom_hline(yintercept = expr.mean.gu, linetype = "longdash", colour="black") +
    geom_point(size=3, aes(colour = factor(cluster), shape = factor(IDH.status))) +
    scale_colour_manual(breaks=c("LGm1","LGm2","LGm3","LGm4","LGm5","LGm6","Cerebellum", "Frontal.cortex", "Temporal.cortex", "Pons"),
                        values = c("LGm1" = "green","LGm2" = "red","LGm3" = "purple", "LGm4" = "orange", "LGm5" = "yellow", "LGm6" = "blue", "Cerebellum" = "aquamarine", "Frontal.cortex" = "pink", "Temporal.cortex" = "brown", "Pons" = "deeppink")) +
    scale_x_continuous(limits=c(0,1), breaks=c(0, 0.3, 0.5, 0.7, 1)) +
    #scale_y_log10() +
    ylab("RNA seq expression (log2)") +
    xlab(paste("DNA Methylation Values\n", "cg03079681.CDKN2A.ILMN_1717714")) + 
    labs(title = paste("DNA Methylation vs Expression", "cg03079681.CDKN2A.ILMN_1717714", sep=" - "), colour = "DNA Meth. Cluster", shape = "IDH status") 
  #ggsave(file=paste0("DNAmethylationVSexpression-", "cg18230216.NYX", ".png"))
  
  
  
  
  epirho <- mclapply(epiGene, function(ii) {
    e <- (ii[, "exprValues"])
    m <- (ii[, "methValues"])
    cor.pvalue <- cor.test(e, m, method="spearman",alternative="two.sided")
    return(c(rho=as.numeric(cor.pvalue$estimate), p.value=cor.pvalue$p.value))
  }, mc.cores=1)
  
  epirho <- do.call(rbind.data.frame, epirho)
  colnames(epirho) <- c("rho", "p.value")
  
  #require(stringr)
  epirho$gene.names <- rownames(epirho)
  
  gene.names.rle <- epirho$gene.names
  
  epirho$ismin <- NA
  start <- 1
  for(i in gene.names.rle$lengths) {
    end <- start + i - 1
    minrho <- min(epirho[c(start:end), "rho"], na.rm=T)
    epirho[c(start:end), "ismin"] <- epirho[c(start:end), "rho"] == minrho
    start <- end + 1
  }
  epirho.s <- epirho[epirho$ismin, ]
  epirho.s <- subset(epirho.s, p.value <=0.05)
  epirho.s$logp <- -log10(epirho.s$p.value)
  epirho.s[epirho.s[, "logp"] == Inf, ] <- 83
  
  epiGene.s.s <- epiGene.s[rownames(epirho.s)]
  
  x <- lapply(epiGene.s.s, function(z) {
    z.es <- z$meth.gt0.3 & z$exp.lt.meanUnmeth
    z.names <- z$new
    z.es <- data.frame(z.es)
    row.names(z.es) <- z.names
    return(z.es)
  })
  
  x <- do.call(cbind.data.frame, x)
  colnames(x) <- names(epiGene.s.s)
  x <- t(x)
  
  
  epiplot <- epiGene[["cg00027083.EPB41L3.ILMN_1720637"]]
  epiplot$exprValues <- epiplot$exprValues + 1
  expr.mean.gu <- mean(as.numeric(epiplot[epiplot$meth.gt0.3==FALSE, "exprValues"]), na.rm = T)
  
  p <- ggplot(epiplot, aes(methValues, exprValues))
  p + 
    stat_smooth(method = "lm", fullrange=FALSE, na.rm=TRUE, se=FALSE, size=1) +
    geom_vline(xintercept = 0.3, linetype = "longdash", colour="black") +
    geom_hline(yintercept = expr.mean.gu, linetype = "longdash", colour="black") +
    geom_point(size=3, aes(colour = factor(cluster.meth), shape = factor(codel1p19q))) +
    scale_colour_manual(values = c("LGm1" = "green","LGm2" = "red","LGm3" = "purple", "LGm4" = "orange", "LGm5" = "yellow", "LGm6" = "blue")) +
    scale_x_continuous(limits=c(0,1), breaks=c(0, 0.3, 0.5, 0.7, 1)) +
    scale_y_log10() +
    ylab("RNA seq expression (log2)") +
    xlab(paste("DNA Methylation Values\n")) + 
    labs(title = "DNA Methylation vs Expression", colour = "DNA Meth. Cluster", shape = "1p19qStatus") 
  ggsave(file=paste0("DNAmethylationVSexpression-", gene.of.interest, ".png"))
  
  LGG.GBM.new.noXY.s <- LGG.GBM.new.noXY[rownames(full.heatmap.orderd) ,5:936]
  LGG.GBM.new.order <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
  a <- metadata.1129.samples.20141030[substr(colnames(LGG.GBM.new.order),1,12),]
  a <- a[c(1:63,64:327,405:531,532:682,683:932,328:404),]
  p <- ggplot(a, aes(y=age, x=1:932))
  png("/dados/ResearchProjects/thais/TCGA/LGG.GBM/age2.png",res=200,width=4000,height=1000)
  p + geom_point(size=3, aes(colour = factor(cluster.meth)),na.rm=T) + 
    scale_colour_manual(values = c("LGm1" = "green","LGm2" = "red","LGm3" = "purple", "LGm4" = "orange", "LGm5" = "yellow", "LGm6" = "blue")) +
    theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white"), panel.border = element_rect(colour = "white"), legend.position="none") +
    xlab("DNA Methylation Clusters") + 
    ylab("Age") +
    labs("TCGA LGG and GBM")
  #ggsave(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/age.pdf")
  dev.off()
  
  
  #scatter plot for leukocyte score
  ggplot() +
    
    stat_smooth(data = metadata, aes(x = leukocyte, y = rna.estimates), method = "glm") +
    geom_point(data = metadata, aes(x = leukocyte, y = rna.estimates,colour = cluster.meth)) +
    scale_colour_manual(values = c("LGm1" = "green","LGm2" = "red","LGm3" = "purple", "LGm4" = "orange", "LGm5" = "yellow", "LGm6" = "blue")) +
    # scale_y_log10() +
    ylab("Leukocyte Score for Expression") +
    xlab("Leukocyte Score for DNA methylation") 
  
  #  facet_grid(. ~ cluster.meth, scales="free")
  
  
  ### cell line analysis
  
  #LGG
  setwd("/dados/ResearchProjects/thais/TCGA/LGG.GBM/TCGA_Cell_Line")
  files = list.files(pattern="txt",full.names=T) 
  for(z in 1:16)
    if(is.null(lgg.cell.line)){ #read the first file
      lgg.cell.line = read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1) 
      lgg.cell.line = subset(lgg.cell.line,select=c(Composite.Element.REF,Gene_Symbol,Chromosome,Genomic_Coordinate,Beta_value)) 
      colnames(lgg.cell.line)[5] = patientID
    }
  else{
    aux = read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1)
    colnames(aux)[2] = patientID
    lgg.cell.line = merge(lgg.cell.line, aux[,1:2], by.x ="Composite.Element.REF", by.y = "Composite.Element.REF")
    
    #GBM
    batch.gbm <- read.delim("gbm_brennan2013_info.txt")
    batch.gbm <- batch.gbm[batch.gbm$PLATFORM_NAME == "HumanMethylation450" | batch.gbm$PLATFORM_NAME == "HumanMethylation27",]  
    batch.gbm <- batch.gbm[batch.gbm$DATA_LEVEL == "3",]
    
    
    
    agecat <- cut(xxMODr$age, c(10.90,38.80,52,63,89.30))
    xxMODr$age.strat <- as.character(agecat)
    
    #Information about GBM clusters from Brennan, 2013
    brennan.class <- read.delim(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Brennan2013_GBM27K-450K/DNA.methylation.k6.txt", sep="\t", header=T)
    #cc.out.purity.lgg.gbm.27.450.2$TCGAIDnew <- substr(rownames(cc.out.purity.lgg.gbm.27.450.2),1,12)
    brennan.class$TCGAIDnew <- substr(brennan.class$TCGAID, 1, 12)
    brennan.class.s <- brennan.class[brennan.class$DNA.Methylation.Clusters != "UNKNOWN", ]
    
    cc.out.purity.lgg.gbm.27.450.2$TCGAID <- rownames(cc.out.purity.lgg.gbm.27.450.2)
    lgg.gbm.anno.27.450 <- merge(cc.out.purity.lgg.gbm.27.450.2, brennan.class.s[,c("DNA.Methylation.Clusters","TCGAIDnew")], by = "TCGAIDnew", all.x =T)
    rownames(lgg.gbm.anno.27.450) <- lgg.gbm.anno.27.450$TCGAID
    
    
    ##load in LGG clinical annotations
    ### March 14, 2014 from synapse:
    #library("synapseClient")
    #synapseLogin()
    #syn2369104 <- synGet(id='syn2369104')
    #syn2369104 <- synGet(id='syn2369104', load=T) #updated march 14
    load(file="~/.synapseCache/615/490615/LGG-Annotation.3-14-2014.RData")
    recurrent.cases <- data.frame(
      c("TCGA-DH-A669",  ## not in data freeze
        "TCGA-TM-A7CF",  ## not in data freeze
        "TCGA-FG-5963",
        "TCGA-FG-5965",
        "TCGA-FG-A4MT"),
      c(rep("Recurrent", times=5))
    )
    colnames(recurrent.cases) <- c("Tumor", "Recurrent.status")
    Clinical <- merge(Clinical, recurrent.cases, by = "Tumor", all.x = T)
    clinical.cluster <- merge(Clustering, Clinical, by = "Tumor")
    clinical.cluster.genetics <- merge(clinical.cluster, Genetics, by = "Tumor")
    oncosign <- read.delim(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/oncosign_assignments.txt", sep="\t", header=T) ## new oncosign assignments provided by giovanni (April 5, 2014)
    clinical.cluster.genetics <- merge(clinical.cluster.genetics, oncosign, by.x = "Tumor", by.y = "sample", all.x =T)
    
    xlgg.gbm.anno.27.450 <- merge(lgg.gbm.anno.27.450, clinical.cluster.genetics[,c("Tumor","MethylationCluster","COCCluster","gender","age_at_initial_pathologic","death01","daystolastordeath","ProgFreeSurvEvent01","ProgFreeSurvTime_days","IDH1")], by.x = "TCGAIDnew", by.y = "Tumor", all.x = T)
    rownames(xlgg.gbm.anno.27.450) <- xlgg.gbm.anno.27.450$TCGAID
    colnames(xlgg.gbm.anno.27.450)[13] <- "GBM.DNA.Methyl.Clusters"
    colnames(xlgg.gbm.anno.27.450)[14] <- "LGG.DNA.Methyl.Clusters"
    
    library(xlsx)
    gbm.molecular.subtype <- read.xlsx("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Molecular_subtype_classification.xlsx",1)
    rownames(gbm.molecular.subtype) <- gbm.molecular.subtype$sample.id
    aa <- merge(xlgg.gbm.anno.27.450,gbm.molecular.subtype, by.x = "TCGAIDnew", by.y = 0, all.x = T)
    
    
    ##Adding the clinical data from GBM to our dataframe
    gbm.clinical.data <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Brennan2013_GBM27K-450K/Brennan_clinicalStuff/Gcimp_call.april2013.csv")
    gbm.clinical.data.nR <- subset(gbm.clinical.data, primary.recurrent!="recurrent")
    gbm.clinical.data.nR$TCGAID <- substr(gbm.clinical.data.nR$TCGAID,1,12)
    #gbm.clinical.data.sub <- gbm.clinical.data[,c("TCGAID","Age.at.Procedure","IDH1.status","Vital.Status","OS..days.","PFS..days.")]
    gbm.clinical.data.nR$Vital.Status <- str_replace(gbm.clinical.data.nR$Vital.Status,"DECEASED","1")
    gbm.clinical.data.nR$Vital.Status <- str_replace(gbm.clinical.data.nR$Vital.Status,"LIVING","0")
    gbm.clinical.data.nR$PROGRESSION.Status <- str_replace(gbm.clinical.data.nR$PROGRESSION.Status,"PROGRESSION","1")
    gbm.clinical.data.nR$PROGRESSION.Status <- str_replace(gbm.clinical.data.nR$PROGRESSION.Status,"STABLE","0")
    gbm.clinical.data.nR$IDH1.status <- str_replace(gbm.clinical.data.nR$IDH1.status,"WT","wt")
    gbm.clinical.data.nR$IDH1.status <- str_replace(gbm.clinical.data.nR$IDH1.status,"R132C|R132G|R132H","Mutant")
    colnames(gbm.clinical.data.nR)[c(27,32,34,35,36,37)] <- c("age_at_initial_pathologic","IDH1","death01","daystolastordeath","ProgFreeSurvEvent01","ProgFreeSurvTime_days")
    xx <- merge(xlgg.gbm.anno.27.450,gbm.clinical.data.nR[,c("TCGAID","Gender","Cluster","age_at_initial_pathologic","Gender","IDH1","death01","daystolastordeath","ProgFreeSurvEvent01","ProgFreeSurvTime_days")],by.x = "TCGAIDnew", by.y = "TCGAID", all.x = TRUE)
    #write.table(xx,file="lgg_gbm_27_450_905_samples.txt",quote=F,row.names=F,sep="\t")
    ###downloaded the RNAcluster from Synapse
    #https://www.synapse.org/#!Synapse:syn2407035
    
    ###opened it in excel, and manually massaged the data
    xxMod <- read.delim(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/lgg_gbm_27_450_905_Samples_FINAL.txt")
    
    rnacluster <- read.delim(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/RNAseqClusters (syn2407035)//pooledClassification-11-4-2014.csv", sep=",", header=T)
    names(rnacluster) <- c("RNA.TCGAID", "RNA.cluster", "RNA.BrennanCluster", "RNA.platform", "RNA.tumor")
    xxMODr <- merge(xxMod, rnacluster, by.x = "TCGAIDnew", by.y = "RNA.TCGAID", all.x = T)
    rownames(xxMODr) <- xxMODr$TCGAID
    
    clinicaldata.all <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/pdata.516lgg.606gbm.txt")
    
    
    #Clinical Information for new LGG samples
    lgg.newid <- xlgg.gbm.anno.27.450[colnames(lgg.patients)[5:224],"TCGAIDnew"]
    new.lgg.patients <- clinicaldata.all[lgg.newid,]
    a <- xxMODr
    new.lgg.patients$vital <- str_replace(new.lgg.patients$vital,"alive","0")
    new.lgg.patients$vital <- str_replace(new.lgg.patients$vital,"dead","1")
    new.lgg.patients$mut.IDH1 <- str_replace(new.lgg.patients$mut.IDH1,"Mut","Mutant")
    new.lgg.patients$mut.IDH1 <- str_replace(new.lgg.patients$mut.IDH1,"WT","wt")
    xxMODr[colnames(lgg.patients)[5:224],c("age","death01","daystolastordeath","IDH1")] <- new.lgg.patients[,c("age","vital","os.days","mut.IDH1")]
    
    
    ##add in age
    agecat <- cut(xxMODr$age, c(10.90,38.80,52,63,89.30))
    xxMODr$age.strat <- as.character(agecat)
    
    leukocyte <- read.csv("/dados/ResearchProjects/thais/TCGA/LGG.GBM/leukocyte_for_houtan2.csv")
    leukocyte <- leukocyte[,-1] #remove column with row number
    quantile(leukocyte$leukocyte,probs=seq(0,1,.25), na.rm=T) #divide in 5 quantiles
    
    xxMODr <- merge(xxMODr,leukocyte, by.x=0, by.y= "x",all.x=T)
    rownames(xxMODr) <- xxMODr$Row.names
    xxMODr <- xxMODr[,-1]
    leukocytecat <- cut(xxMODr$leukocyte, c(0.0136,0.10340341,0.15807936,0.24441133,0.77))
    xxMODr$leukocyte.strat <- leukocytecat
    #save(xxMODr, file="summaryTables.Rda")
    
    estimate.affy <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/estimate.gbm529.affy.u133a.txt",skip=2)
    estimate.unc <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/estimate.gbm154.lgg463.unc.rnaseqv2.txt",skip=2)
    library(Hmisc)
    uniques.affy <- estimate.affy[,colnames(estimate.affy) %nin% colnames(estimate.unc)]
    uniques.unc <- estimate.unc[,colnames(estimate.unc) %nin% colnames(estimate.affy)]
    duplicates <- estimate.unc[,colnames(estimate.unc) %in% colnames(estimate.affy)[3:531]]
    uniques.affy <- uniques.affy[3,]
    uniques.unc <- uniques.unc[3,]
    duplicates <- duplicates[3,]
    #duplicates <- uniques.unc[colnames(estimate.unc) %nin% colnames(uniques.unc)]
    rna.estimates <- cbind(uniques.affy,uniques.unc,duplicates)
    rna.estimates <- as.matrix(t(rna.estimates))
    rownames(rna.estimates) <- str_replace_all(rownames(rna.estimates),"[.]","-")
    
    a <- merge(xxMODr,rna.estimates,by.x="TCGAIDnew",by.y=0,all.x=T)
    colnames(a)[65] <- "rna.estimates"
    rownames(a) <- a$TCGAID
    
    library(heatmap.plus)
    library(matlab)
    
    
    colors.label <- function(cluster,gbm.cluster,lgg.cluster,tumor,IDH,rna.dna.c,rna.cluster,coc.c,leukocyte,estimates,sturm,new.cl,anno){
      for(i in 1:length(cluster)){
        if(is.na(cluster[i]))
          cluster[i] <- "white"
        else{ 
          if(cluster[i] == "1")
            cluster[i] <- "red"
          else if(cluster[i] == "2")
            cluster[i] <- "yellow"
          else if(cluster[i] == "3")
            cluster[i] <- "green"
          else if(cluster[i] == "4")
            cluster[i] <- "blue"
          else if(cluster[i] == "5")
            cluster[i] <- "orange"
          else if(cluster[i] == "6")
            cluster[i] <- "purple"
          else 
            cluster[i] <- "darkblue"
        }
      }
      
      for(i in 1:length(gbm.cluster)){
        if(is.na(gbm.cluster[i]))
          gbm.cluster[i] <- "white"
        else{
          if(gbm.cluster[i] == "G-CIMP")
            gbm.cluster[i] <- "yellow"
          else if(gbm.cluster[i] == "M1")
            gbm.cluster[i] <- "green"
          else if(gbm.cluster[i] == "M2")
            gbm.cluster[i] <- "black"
          else if(gbm.cluster[i] == "M3")
            gbm.cluster[i] <- "blue"
          else if(gbm.cluster[i] == "M4")
            gbm.cluster[i] <- "red"
          else if(gbm.cluster[i] == "M6")
            gbm.cluster[i] <- "darksalmon"
          else 
            gbm.cluster[i] <- "white"
        }
      }
      
      for(i in 1:length(lgg.cluster)){
        if(is.na(lgg.cluster[i]))
          lgg.cluster[i] <- "white"
        else{
          if(lgg.cluster[i] == "M5")
            lgg.cluster[i] <- "darkorchid4"
          else if(lgg.cluster[i] == "M1")
            lgg.cluster[i] <- "black"
          else if(lgg.cluster[i] == "M2")
            lgg.cluster[i] <- "red"
          else if(lgg.cluster[i] == "M3")
            lgg.cluster[i] <- "blue"
          else if(lgg.cluster[i] == "M4")
            lgg.cluster[i] <- "forestgreen"
          else 
            lgg.cluster[i] <- "white"
        }
      }
      
     for(i in 1:length(sturm)){
        if(is.na(sturm[i]))
          sturm[i] <- "white"
        else{
          if(sturm[i] == "IDH")
            sturm[i] <- "red3"
          else if(sturm[i] == "Mesenchymal")
            sturm[i] <- "darkgoldenrod3"
          else if(sturm[i] == "RTK I 'PDGFRA'")
            sturm[i] <- "tan3"
          else if(sturm[i] == "RTK II 'Classic'")
            sturm[i] <- "firebrick4"
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
            tumor[i] <- "darkgrey"
          else
            tumor[i] <- "green"
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
            IDH[i] <- "green"
          
        }
      }
      
      #age <- jet.colors(100)[ceiling((age-min(age,na.rm=T))/(max(age,na.rm=T)-min(age,na.rm=T))*(length(jet.colors(100))-1))+1]
      #age[is.na(age)] <- "white"
      estimates <- jet.colors(100)[ceiling((estimates-min(estimates,na.rm=T))/(max(estimates,na.rm=T)-min(estimates,na.rm=T))*(length(jet.colors(100))-1))+1]
      estimates[is.na(estimates)] <- "white"
      leukocyte <- jet.colors(100)[ceiling((leukocyte-min(leukocyte,na.rm=T))/(max(leukocyte,na.rm=T)-min(leukocyte,na.rm=T))*(length(jet.colors(100))-1))+1]
      
      for(i in 1:length(rna.dna.c)){
        if(is.na(rna.dna.c[i]))
          rna.dna.c[i] <- "white"
        else{
          if(rna.dna.c[i] == "Classical")
            rna.dna.c[i] <- "red"
          else if(rna.dna.c[i] == "G-CIMP")
            rna.dna.c[i] <- "yellow"
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
            rna.cluster[i] <- "gray"
          else if(rna.cluster[i] == "LGr3")
            rna.cluster[i] <- "cyan"
          else if(rna.cluster[i] == "LGr4")
            rna.cluster[i] <- "magenta"
          else if(rna.cluster[i] == "LGr5")
            rna.cluster[i] <- "blue"
          else if(rna.cluster[i] == "LGr6")
            rna.cluster[i] <- "green3"
          else if(rna.cluster[i] == "LGr7")
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
      
      for(i in 1:length(new.cl)){
        if(is.na(new.cl[i]))
          new.cl[i] <- "white"
        else{
          if(new.cl[i] == "IDHmut")
            new.cl[i] <- "brown1"
          else if(new.cl[i] == "IDHmut.codel")
            new.cl[i] <- "darkseagreen4"
          else if(new.cl[i] == "IDHwt")
            new.cl[i] <- "chartreuse3"
          else
            new.cl[i] <- "white"
          
        }
      }
      
      
      for(i in 1:length(anno)){
        if(is.na(anno[i]))
          anno[i] <- "white"
        else{
          if(anno[i] == "Case submitted is found to be a recurrence after submission")
            anno[i] <- "yellow"
          else if(anno[i] == "Case submitted is found to be a recurrence after submission/Item does not meet study protocol")
            anno[i] <- "gray"
          else if(anno[i] == "History of unacceptable prior treatment related to a prior/other malignancy")
            anno[i] <- "cyan"
          else if(anno[i] == "Item in special subset")
            anno[i] <- "magenta"
          else if(anno[i] == "Molecular analysis outside specification")
            anno[i] <- "blue"
          else if(anno[i] == "Neoadjuvant therapy")
            anno[i] <- "green3"
          else if(anno[i] == "Neoadjuvant therapy/Case submitted is found to be a recurrence after submission")
            anno[i] <- "red"
          else if(anno[i] == "Normal tissue origin incorrect")
            anno[i] <- "brown"
          else if(anno[i] == "Item may not meet study protocol")
            anno[i] <- "pink"
          else if(anno[i] == "Pathology outside specification")
            anno[i] <- "palevioletred4"
          else if(anno[i] == "Prior malignancy")
            anno[i] <- "seagreen1"
          else if(anno[i] == "Prior malignancy/History of acceptable prior treatment related to a prior/other malignancy")
            anno[i] <- "orange"
          else if(anno[i] == "Qualified in error")
            anno[i] <- "black"
          else if(anno[i] == "Qualified in error/Item in special subset")
            anno[i] <- "royalblue1"
          else if(anno[i] == "Subject withdrew consent")
            anno[i] <- "purple"
          
          else 
            anno[i] <- "white"
        }
      }
      
      
      
      cluster <- cbind(anno,IDH,sturm,rna.dna.c,gbm.cluster,coc.c,lgg.cluster,rna.cluster,cluster,new.cl,estimates,leukocyte,tumor)
      colnames(cluster) <- c("Annotation","IDH1","Sturm et al","GBM Clusters [Brennan, 2013]","GBM DNA Meth Cl [Brennan, 2013]","COC Clusters [NEJM, unpublished]", "LGG Clusters [NEJM, unpublished]", "RNA Clusters", "Consensus Cluster","New Cluster","Estimates", "Leukocyte", "Tumor Type")
      #colnames(cluster) <- c("Age", "IDH1","GBM Clusters [Brennan, 2013]","GBM DNA Meth Cl [Brennan, 2013]","COC Clusters [NEJM, unpublished]", "LGG Clusters [NEJM, unpublished]", "RNA Clusters", "Consensus Cluster","Estimates", "Leukocyte", "Tumor Type")
      return(cluster)
    }
    
    source("/dados/ResearchProjects//thais/heatmap.plus.R")
    #Heatmap for 4 Clusters
    GBM.LGG.27.450k.noXY.s <- GBM.LGG.27.450k.noXY[rownames(dat.lgg.gbm.27.450.noXY.dic.oecg) ,5:909]
    lgg.gbm.27.450.order <- GBM.LGG.27.450k.noXY.s[,cc.out.purity.lgg.gbm.27.450[[4]][[2]][[3]]]
    xxMODr <- xxMODr[colnames(lgg.gbm.27.450.order),]
    clab.4 <- colors.label(as.matrix(xxMODr$K4),as.matrix(xxMODr$GBM.DNA.Methyl.Clusters),as.matrix(xxMODr$LGG.DNA.Methyl.Clusters),as.matrix(xxMODr$TumorType),as.matrix(xxMODr$age.strat),as.matrix(xxMODr$IDH1),as.matrix(xxMODr$Cluster),as.matrix(xxMODr$RNA.cluster),as.matrix(xxMODr$COCCluster))
    lgg.gbm.27.450.order <- lgg.gbm.27.450.order[heatmap.lgg.gbm$rowInd,]
    normals.sel <- normals.sel[rownames(lgg.gbm.27.450.order),]
    lgg.gbm.27.450.order <- cbind(normals.sel,lgg.gbm.27.450.order)
    norm.label <- matrix("white", nrow = 77, ncol = 9)
    clab.4 <- rbind(norm.label,clab.4)
    png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_905_s/heatmap_K4_normals_3.png",res=200,width=2000,height=2000)
    heatmap.lgg.gbm <- heatmap.plus.sm(as.matrix(lgg.gbm.27.450.order),
                                       col = jet.colors(75),
                                       scale = "none",
                                       labRow = NA,
                                       labCol = NA,
                                       Colv = NA,
                                       Rowv = NA,
                                       ColSideColors = clab.4
    )
    dev.off()
    
    #Heatmap for 5 Clusters
    GBM.LGG.27.450k.noXY.s <- GBM.LGG.27.450k.noXY[rownames(dat.lgg.gbm.27.450.noXY.dic.oecg) ,5:909]
    lgg.gbm.27.450.order <- GBM.LGG.27.450k.noXY.s[,cc.out.purity.lgg.gbm.27.450[[5]][[2]][[3]]]
    xxMODr <- xxMODr[colnames(lgg.gbm.27.450.order),]
    clab.5 <- colors.label(as.matrix(xxMODr$K5),as.matrix(xxMODr$GBM.DNA.Methyl.Clusters),as.matrix(xxMODr$LGG.DNA.Methyl.Clusters),as.matrix(xxMODr$TumorType),as.matrix(xxMODr$age.strat),as.matrix(xxMODr$IDH1),as.matrix(xxMODr$Cluster),as.matrix(xxMODr$RNA.cluster),as.matrix(xxMODr$COCCluster))
    clab.5 <- clab.5[c(1:361,822:905,362:530,531:744,745:821),]
    lgg.gbm.27.450.order <- lgg.gbm.27.450.order[,c(1:361,822:905,362:530,531:744,745:821)]
    lgg.gbm.27.450.order <- lgg.gbm.27.450.order[heatmap.lgg.gbm$rowInd,]
    normals.sel <- normals.sel[rownames(lgg.gbm.27.450.order),]
    lgg.gbm.27.450.order <- cbind(normals.sel,lgg.gbm.27.450.order)
    norm.label <- matrix("white", nrow = 77, ncol = 9)
    clab.5 <- rbind(norm.label,clab.5)
    png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_905_s/heatmap_K5_normals3.png",res=200,width=2000,height=2000)
    heatmap.lgg.gbm <- heatmap.plus.sm(as.matrix(lgg.gbm.27.450.order),
                                       col = jet.colors(75),
                                       scale = "none",
                                       labRow = NA,
                                       labCol = NA,
                                       Colv = NA,
                                       Rowv = NA,
                                       ColSideColors = clab.5
    )
    dev.off()
    
    #Heatmap for 6 Clusters
    GBM.LGG.27.450k.noXY.s <- GBM.LGG.27.450k.noXY[rownames(dat.lgg.gbm.27.450.noXY.dic.oecg) ,5:909]
    lgg.gbm.27.450.order <- GBM.LGG.27.450k.noXY.s[,cc.out.purity.lgg.gbm.27.450[[6]][[2]][[3]]]
    #a <- a[colnames(lgg.gbm.27.450.order),]
    metadata <- metadata[colnames(lgg.gbm.27.450.order),]
    clab.6 <- colors.label(as.matrix(metadata$K6),as.matrix(metadata$DNA.Methylation.Clusters),as.matrix(metadata$MethylationCluster),as.matrix(metadata$tumor.type),as.matrix(metadata$age),as.matrix(metadata$IDH.status),as.matrix(metadata$Cluster),as.matrix(metadata$cluster.mRNA),as.matrix(metadata$COCCluster),as.numeric(metadata$leukocyte),as.matrix(metadata$rna.estimates))
    #clab.6 <- colors.label(as.matrix(a$K6),as.matrix(a$GBM.methylation.cl),as.matrix(a$LGG.methylation.c),as.matrix(a$tumor.type),as.matrix(a$age),as.matrix(a$IDH.status),as.matrix(a$GBM.cluster),as.matrix(a$cluster.mRNA),as.matrix(a$COCCluster),as.numeric(a$leukocyte),as.matrix(a$rna.estimates))
    clab.6 <- clab.6[c(856:905,1:333,794:855,334:489,490:708,709:793),]
    lgg.gbm.27.450.order <- lgg.gbm.27.450.order[,c(856:905,1:333,794:855,334:489,490:708,709:793)]
    lgg.gbm.27.450.order <- lgg.gbm.27.450.order[heatmap.lgg.gbm$rowInd,]
    normals.sel <- normals.sel[rownames(lgg.gbm.27.450.order),]
    lgg.gbm.27.450.order <- cbind(normals.sel,lgg.gbm.27.450.order)
    norm.label <- matrix("white", nrow = 77, ncol = 11)
    clab.6 <- rbind(norm.label,clab.6)
    png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_905_s/heatmap_K6_normals_leuk4_testerow.png",res=200,width=2000,height=2000)
    heatmap.plus.sm(as.matrix(lgg.gbm.27.450.order),
                    col = jet.colors(75),
                    scale = "none",
                    labRow = NA,
                    labCol = NA,
                    Colv = NA,
                    #Rowv = NA,
                    #ColSideColors = clab.6
    )
    dev.off()
    
    ## Key
    z <- seq(min(LGG.GBM.250914[,5:936],na.rm=TRUE), max(LGG.GBM.250914[,5:936],na.rm=TRUE), length = 200) #min(lgg.gbm.27.450.order) and #max(lgg.gbm.27.450.order)
    n <- max(LGG.GBM.250914[,5:936],na.rm=TRUE)
    png("/dados/ResearchProjects/thais/TCGA/LGG.GBM/key.png",res=600,width=4000,height=2000)
    image(matrix(z, ncol = 1), col = jet.colors(75),
          xaxt = "n", yaxt = "n", main = "Overlap count Key")
    box()
    par(usr = c(0, n, 0, n))
    axis(1, at = c(min(z):max(z)))
    
    
    #Heatmap for 7 Clusters
    GBM.LGG.27.450k.noXY.s <- GBM.LGG.27.450k.noXY[rownames(dat.lgg.gbm.27.450.noXY.dic.oecg) ,5:909]
    lgg.gbm.27.450.order <- GBM.LGG.27.450k.noXY.s[,cc.out.purity.lgg.gbm.27.450[[7]][[2]][[3]]]
    xxMODr <- xxMODr[colnames(lgg.gbm.27.450.order),]
    clab.7 <- colors.label(as.matrix(xxMODr$K7),as.matrix(xxMODr$GBM.DNA.Methyl.Clusters),as.matrix(xxMODr$LGG.DNA.Methyl.Clusters),as.matrix(xxMODr$TumorType),as.matrix(xxMODr$age.strat),as.matrix(xxMODr$IDH1),as.matrix(xxMODr$Cluster),as.matrix(xxMODr$RNA.cluster),as.matrix(xxMODr$COCCluster))
    clab.7 <- clab.7[c(396:445,1:228,293:395,229:292,523:680,681:905,446:522),]
    lgg.gbm.27.450.order <- lgg.gbm.27.450.order[,c(396:445,1:228,293:395,229:292,523:680,681:905,446:522)]
    lgg.gbm.27.450.order <- lgg.gbm.27.450.order[heatmap.lgg.gbm$rowInd,]
    normals.sel <- normals.sel[rownames(lgg.gbm.27.450.order),]
    lgg.gbm.27.450.order <- cbind(normals.sel,lgg.gbm.27.450.order)
    norm.label <- matrix("white", nrow = 77, ncol = 9)
    clab.7 <- rbind(norm.label,clab.7)
    png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_905_s/heatmap_K7_normals2.png",res=200,width=2000,height=2000)
    heatmap.lgg.gbm <- heatmap.plus.sm(as.matrix(lgg.gbm.27.450.order),
                                       col = jet.colors(75),
                                       scale = "none",
                                       labRow = NA,
                                       labCol = NA,
                                       Colv = NA,
                                       Rowv = NA,
                                       ColSideColors = clab.7
    )
    dev.off()
    
    #Normal samples 
    library(Biobase)
    load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/normals.gse41826.rda")
    norms <- exprs(controls.rm)
    normals.sel <- norms[rownames(lgg.gbm.27.450.order),]
    #normals.sel <- normals.sel[,-1]
    normals.sel <- normals.sel[heatmap.lgg.gbm$rowInd,]
    
    
    ## plot age at diagnosis
    library(ggplot2)
    x <- metadata
    #x$K7 <- factor(x$K7,c("3","1","5","7","6","2","4")) #change the order
    x$K6 <- factor(x$cluster.meth) #change the order
    ages <- x[,c("K6", "age")]
    #ages <- na.omit(ages)
    p <- ggplot(ages, aes(factor(K6), age))
    png("/dados/ResearchProjects/thais/TCGA/LGG.GBM/age.png",res=600,width=20,height=15,units = "in")
    p + geom_boxplot(aes(fill = factor(K6)), notchwidth=0.25) + 
      geom_jitter(height = 0, position = position_jitter(width = .1), size=2)  + 
      scale_fill_manual(values=c("green","red","purple", "orange", "yellow","blue"
      ), labels = c("LGm1", "LGm2", "LGm3","LGm4", "LGm5", "LGm6"
      )) +
      #scale_fill_manual(values = rev(c("darkorchid","darkgreen","blue","red","DARKGREY"))) +
      ylab("Age At Diagnosis") +
      xlab("DNA Methylation Clusters") + 
      labs(title = "DNA Methylation Clusters by Age At Diagnosis", fill = "DNA Meth. Cluster") +
      theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white"), panel.border = element_rect(colour = "white"), legend.position="none")
    dev.off()
    
    
    t.test(metadata[metadata$cluster.meth == "LGm5","leukocyte"],metadata[metadata$cluster.meth == "LGm6","leukocyte"])
    
    ## File to upload at Synapse.org
    a <- xxMODr
    a <- a[,c("TCGAID","K6")]
    a$K6 <- str_replace(a$K6,"^1$","LGm2")
    a$K6 <- str_replace(a$K6,"^2$","LGm5")
    a$K6 <- str_replace(a$K6,"^3$","LGm1")
    a$K6 <- str_replace(a$K6,"^4$","LGm6")
    a$K6 <- str_replace(a$K6,"^5$","LGm4")
    a$K6 <- str_replace(a$K6,"^6$","LGm3")
    write.table(a,file="TCGAID_K6.txt",quote=F,row.names=F,sep="\t")
    
    xxMODr$order <- xxMODr$K6
    xxMODr$order <- str_replace(xxMODr$order,"^1$","LGm2")
    xxMODr$order <- str_replace(xxMODr$order,"^2$","LGm5")
    xxMODr$order <- str_replace(xxMODr$order,"^3$","LGm1")
    xxMODr$order <- str_replace(xxMODr$order,"^4$","LGm6")
    xxMODr$order <- str_replace(xxMODr$order,"^5$","LGm4")
    xxMODr$order <- str_replace(xxMODr$order,"^6$","LGm3")
    
    
    #Beta Value Diff M1 x M2
    M1 <- metadata[metadata$cluster.meth == "LGm1" | metadata$cluster.meth == "LGm2" | metadata$cluster.meth == "LGm3" ,]
    Lgg.like <-  M1[M1$tumor.type == "GBM",]
    Lgg.real <- M1[M1$tumor.type == "LGG",]
    M2 <- metadata[metadata$cluster.meth == "LGm4" | metadata$cluster.meth == "LGm5" | metadata$cluster.meth == "LGm6" ,]
    Gbm.like <- M2[M2$tumor.type == "LGG",]
    Gbm.real <- M2[M2$tumor.type == "GBM",]
    #Selecting the Samples
    data <- GBM.LGG.27.450k[,5:909]
    M1 <- data[,as.character(Gbm.like$TCGAID)]
    M2 <- data[,as.character(Gbm.real$TCGAID)]
    volcano <- cbind(M1,M2)
    
    values <- t(na.omit(volcano))
    values <- data.frame(values)
    #values.names <- colnames(values)
    
    require(exactRankTests)
    require(parallel)
    system.time(w.p.values <- unlist(mclapply(values,
                                              function(probe) {
                                                zz <- wilcox.exact(probe[1:30],probe[31:445], exact=TRUE)
                                                z <- zz$p.value
                                                return(z)
                                              }, mc.cores=8)))
    
    #lggs.p.value <- as.matrix(w.p.values)
    #gbms.p.value <- as.matrix(w.p.values)
    
    w.p.values.adj <- p.adjust(w.p.values, method = "BH")
    
    volcano$meanM1 <- apply(volcano[,1:30],1,mean,na.rm=T)
    volcano$meanM2 <- apply(volcano[,31:445],1,mean,na.rm=T)
    volcano$DiffMean <- volcano$meanM1 - volcano$meanM2
    volcano <- na.omit(volcano)
    volcano$p.value.adj <- w.p.values.adj
    volcano$p.value <- w.p.values
    
    volcano$threshold <- "1"
    a <- subset(volcano, p.value.adj < 0.05)
    b <- subset(a, DiffMean < -0.25) #hyper no da direita
    c <- subset(a, DiffMean < 0 & DiffMean > -0.25)
    volcano[rownames(b),"threshold"] <- "2"
    volcano[rownames(c),"threshold"] <- "4"
    b <- subset(a, DiffMean > 0.25) #hyper no da esquerda
    c <- subset(a, DiffMean > 0 & DiffMean < 0.25)
    volcano[rownames(b),"threshold"] <- "3"
    volcano[rownames(c),"threshold"] <- "4"
    png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_905_s/volcano_LGG.png",res=200,width=2000,height=2000)
    ggplot(data=volcano, aes(x=DiffMean, y=-1*log10(p.value.adj), colour=threshold)) +
      geom_point() +
      xlim(c(-1,1)) + ylim(c(0,35)) +
      xlab("DNA Methylation Diff") + ylab("-1 * log10 of the Significance") +
      labs(title = "Volcano Plot") +
      scale_color_manual(breaks=c("1","2","3","4"), # color scale (for points)
                         values=c("black", "green", "red","grey"),
                         labels=c("Not Significant","Hypermethylated in LGG","Hypermethylated in LGG-like","Not significant"),
                         name="Legend") 
    dev.off()
    
    
    ### GBM like, GBM, LGG and LGG like comparison
    a <- GBM.LGG.27.450k[,5:909]
    M1 <- xxMODr[xxMODr$K6 == "3" | xxMODr$K6 == "1" | xxMODr$K6 == "6" ,]
    M2 <- xxMODr[xxMODr$K6 == "5" | xxMODr$K6 == "2" | xxMODr$K6 == "4" ,]
    Lgg.like <-  M1[M1$TumorType == "GBM",]
    Gbm.real <- M2[M2$TumorType == "GBM",]
    M1 <- a[,Lgg.like$TCGAID]
    M2 <- a[,Gbm.real$TCGAID]
    volcano <- cbind(M1,M2)
    volcano$meanM1 <- apply(volcano[,1:30],1,mean,na.rm=T)
    volcano$meanM2 <- apply(volcano[,31:396],1,mean,na.rm=T)
    volcano$DiffMean <- volcano$meanM1 - volcano$meanM2
    
    volcano <- na.omit(volcano)
    volcano$p.value.adj <- w.p.values.adj
    volcano$p.value <- w.p.values
    
    volcano$threshold <- "1"
    a <- subset(volcano, p.value.adj < 0.05)
    b <- subset(a, DiffMean < 0) #hyper no da direita
    beta.v.pe[rownames(b),"threshold"] <- "2"
    b <- subset(a, DiffMean > 0) #hyper no da esquerda
    beta.v.pe[rownames(b),"threshold"] <- "3"
    png(filename="/dados/ResearchProjects/thais/LGG",res=200,width=2000,height=2000)
    ggplot(data=beta.v.pe, aes(x=DiffMean, y=-1*log10(p.value.adj), colour=threshold)) +
      geom_point() +
      xlim(c(-1,1)) + ylim(c(0,4)) +
      xlab("DNA Methylation Diff") + ylab("-1 * log10 of the Significance") +
      labs(title = "Volcano Plot") +
      scale_color_manual(breaks=c("1","2","3"), # color scale (for points)
                         values=c("black", "green", "red"),
                         labels=c("Not Significant","Hypermethylated in G0","Hypermethylated in primary"),
                         name="Legend") 
    dev.off()
    
    
    
    
    a <- xxMODr
    a$Classification <- "nothing"
    b <- a[a$K6 == "3" | a$K6 == "1" | a$K6 == "6" ,]
    c <- subset(b, TumorType == "LGG")
    a[rownames(c),"Classification"] <- "LGG"
    c <- subset(b, TumorType == "GBM")
    a[rownames(c),"Classification"] <- "LGG-like"
    b <- a[a$K6 == "5" | a$K6 == "2" | a$K6 == "4" ,]
    c <- subset(b, TumorType == "LGG")
    a[rownames(c),"Classification"] <- "GBM-like"
    c <- subset(b, TumorType == "GBM")
    a[rownames(c),"Classification"] <- "GBM"
    
    library(LSD)
    
    ###heatpairs by tumor type
    Lgg.like <- a[a$Classification == "LGG-like",]
    Gbm.real <- a[a$Classification == "GBM",]
    Gbm.like <- a[a$Classification == "GBM-like",]
    Lgg.real <- a[a$Classification == "LGG",]
    b <- GBM.LGG.27.450k[,5:909]
    normal <- as.data.frame(norms[rownames(norms) %in% rownames(lgg.gbm.27.450.order),])
    M1 <- b[rownames(normal),as.character(Lgg.like$TCGAID)]
    M2 <- b[rownames(normal),as.character(Gbm.real$TCGAID)]
    M3 <- b[rownames(normal),as.character(Gbm.like$TCGAID)]
    M4 <- b[rownames(normal),as.character(Lgg.real$TCGAID)]
    M1$meanLgg.Like <- apply(M1,1,mean,na.rm=T)
    M2$meanGbm.real <- apply(M2,1,mean,na.rm=T)
    M3$meanGbm.like <- apply(M3,1,mean,na.rm=T)
    M4$meanLgg.real <- apply(M4,1,mean,na.rm=T)
    normal$mean <- apply(normal,1,mean,na.rm=T)
    avg.probes.cl <- cbind(normal$mean,M4$meanLgg.real,M1$meanLgg.Like,M3$meanGbm.like,M2$meanGbm.real)
    rownames(avg.probes.cl) <- rownames(M1)
    colnames(avg.probes.cl) <- c("normal","LGG","LGG-like","GBM-like","GBM")
    heatpairs(avg.probes.cl)
    
    a <- xxMODr
    b <- a[a$K6 == "3" | a$K6 == "1" | a$K6 == "6" ,]
    idh1.mut <- subset(b, IDH1 == "Mutant")
    idh1.wt <- subset(b, IDH1 == "wt")
    c <- GBM.LGG.27.450k[,5:909]
    normal <- as.data.frame(norms[rownames(norms) %in% rownames(lgg.gbm.27.450.order),])
    M1 <- c[rownames(normal),as.character(idh1.mut$TCGAID)]
    M2 <- c[rownames(normal),as.character(idh1.wt$TCGAID)]
    M1$meanM1 <- apply(M1,1,mean,na.rm=T)
    M2$meanM2 <- apply(M2,1,mean,na.rm=T)
    normal$mean <- apply(normal,1,mean,na.rm=T)
    avg.probes.cl <- cbind(normal$mean,M1$meanM1,M2$meanM2)
    rownames(avg.probes.cl) <- rownames(M1)
    colnames(avg.probes.cl) <- c("normal","IDH1 Mut","IDH1 Wt")
    heatpairs(avg.probes.cl)
    
    ###heatpairs by clusters
    K1 <- a[a$LGm.clusters == "LGm1",]
    K2 <- a[a$LGm.clusters == "LGm2",]
    K3 <- a[a$LGm.clusters == "LGm3",]
    K4 <- a[a$LGm.clusters == "LGm4",]
    K5 <- a[a$LGm.clusters == "LGm5",]
    K6 <- a[a$LGm.clusters == "LGm6",]
    b <- GBM.LGG.27.450k[,5:909]
    normal <- as.data.frame(norms[rownames(norms) %in% rownames(GBM.LGG.27.450k),])
    M1 <- b[rownames(normal),as.character(K1$TCGAID)]
    M2 <- b[rownames(normal),as.character(K2$TCGAID)]
    M3 <- b[rownames(normal),as.character(K3$TCGAID)]
    M4 <- b[rownames(normal),as.character(K4$TCGAID)]
    M5 <- b[rownames(normal),as.character(K5$TCGAID)]
    M6 <- b[rownames(normal),as.character(K6$TCGAID)]
    M1$meanM1 <- apply(M1,1,mean,na.rm=T)
    M2$meanM2 <- apply(M2,1,mean,na.rm=T)
    M3$meanM3 <- apply(M3,1,mean,na.rm=T)
    M4$meanM4 <- apply(M4,1,mean,na.rm=T)
    M5$meanM5 <- apply(M5,1,mean,na.rm=T)
    M6$meanM6 <- apply(M6,1,mean,na.rm=T)
    normal$mean <- apply(normal,1,mean,na.rm=T)
    avg.probes.cl <- cbind(normal$mean,M1$meanM1,M2$meanM2,M3$meanM3,M4$meanM4,M5$meanM5,M6$meanM6)
    rownames(avg.probes.cl) <- rownames(M1)
    colnames(avg.probes.cl) <- c("normal","LGm1","LGm2","LGm3","LGm4","LGm5","LGm6")
    heatpairs(avg.probes.cl)
    
    ## Boxplots
    
    library(ggplot2)
    
    ### leukocyte plot
    x <- a
    #x$xxMODr <- factor(x$K7,c("4","1","2","7","6","3","5")) #change the order
    leukocytes <- x[,c("LGm.clusters", "estimates")]
    leukocytes <- na.omit(leukocytes) #para o caso de não haver informaçao para todas as amostras
    p <- ggplot(leukocytes, aes(factor(LGm.clusters), estimates)) #Aqui seleciona os dados. LGm.clusters no eixo x e leukocytes no eixo y
    p + geom_boxplot(aes(fill = factor(LGm.clusters)), notchwidth=0.25) +  #Aqui diz para fazer um boxplot
      geom_jitter(height = 0, position = position_jitter(width = .1), size=2)  +  #mostra os pontinhos
      scale_fill_manual(values=c("green","red","purple", "orange", "yellow","blue"
      ), labels = c("LGm1", "LGm2", "LGm3","LGm4", "LGm5", "LGm6"
      )) +
      #scale_fill_manual(values = rev(c("darkorchid","darkgreen","blue","red","DARKGREY"))) +
      ylab("Leukocyte infiltration Expression Estimate") +
      xlab("DNA Methylation Clusters") + 
      labs(title = "DNA Methylation Clusters by Leukocyte infiltration Expression Estimate", 
           fill = "DNA Meth. Cluster") 