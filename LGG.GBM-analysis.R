##Part LGG, GBM
setwd("/dados/normalized.geneExp/rename")
patients <- list.files(pattern="TCGA*")

gene.id <- read.delim(patients[1], header=TRUE, sep="\t")[,1]

require(plyr)
gene.express <- sapply(patients, function(x) {
  gene.expres <- read.delim(x, header=TRUE, sep="\t")[,2]
})

require(stringr)
gene.id <- str_split(as.character(gene.id), pattern="\\|")
gene.id <- do.call(rbind.data.frame, gene.id)
colnames(gene.id)  <- c("SYM", "GENEID")
gene.matrix <- cbind(gene.id, gene.express)
colnames(gene.matrix)[-c(1:2)] <- substr(colnames(gene.matrix)[-c(1:2)], 1, 12)
crap.genes.map <- read.delim("crapgenesmap.txt", header=T, sep="\t") ## inputted GeneID into pubmed and extracted this table.
x <- gene.id[1:29, ]
xx <- merge(x, crap.genes.map, by.x="GENEID", by.y="GeneID")[, c("Symbol", "GENEID")]
gene.matrix[1:29, c("SYM", "GENEID")]
gene.matrix$SYM <- as.character(gene.matrix$SYM)
gene.matrix$SYM[1:29] <- as.character(xx[, 1]); rm(x); rm(xx)



load("/dados/methylation.LGG.289.rda")
##updated annotation
require(synapseClient)
synapseLogin()
# Load the annotationUpdate
syn2369104 <- synGet(id='syn2369104', load=T) #updated march 14
#syn2366689 <- synGet(id='syn2366689', load=T) #updated march 10
#syn2359167 <- synGet(id='syn2359167', load=T)
#load(file="~/.synapseCache/172/479172/LGG-Annotation.2-28-2014.RData")
#load(file="~/.synapseCache/313/488313//LGG-Annotation.3-10-2014.RData") #updated march 10
load(file="~/.synapseCache/615/490615/LGG-Annotation.3-14-2014.RData") #updated march 14
##creates 4 files: Annotations, Genetics, Clustering, Clinical
## pull in the 7 miRNA cluster
miR7 <- read.delim("~/.synapseCache/713/484713//LGG.miRNAseq.NMF-7groups.cluster-calls.BCGSC.20140115.txt")
miR7$new <- substr(miR7$ID, 1, 12)
names(miR7)[4] <- "miRNACluster7"
Clustering <- merge(Clustering, miR7[,c(7,4)], by.x="Tumor", by.y = "new")

##recurrent status
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

cpg.matrix <- dat
rm(dat)
colnames(cpg.matrix)[-c(1:4)] <- substr(colnames(cpg.matrix)[-c(1:4)], 1, 12)

rsid <- grepl("rs", cpg.matrix$Composite.Element.REF)
cpg.matrix.s <- cpg.matrix[!rsid,]
cpg.matrix.s$Gene_Symbol[cpg.matrix.s$Gene_Symbol==""] <- NA
require(GenomicRanges)
cpg.granges <- GRanges(seqnames=paste0("chr", as.character(cpg.matrix.s$Chromosome)),
                       ranges=IRanges(start=cpg.matrix.s$Genomic_Coordinate,
                                      end=cpg.matrix.s$Genomic_Coordinate,
                                      names=cpg.matrix.s$Composite.Element.REF))

#Aquired from ENSEMBL biomart website - 2014-03-03
## issue: gene symbols are not updated. and in order to get chrom:start:end, we used bioMart website.
gene.annotation <- read.delim("~/Galaxy1-[Homo_sapiens_genes_(GRCh37.p13)].tabular", header=T, sep="\t")
dups <- duplicated(gene.annotation$EntrezGene.ID)
uniques <- gene.annotation[!dups, ]
gene.matrix.loc <- merge(gene.matrix[,1:2], uniques, by.x = "GENEID", by.y = "EntrezGene.ID", all.x=TRUE, all.y=FALSE, incomparables = NA, sort=FALSE)
gene.matrix.loc <- gene.matrix.loc[, -13]
gene.matrix.loc <- na.omit(gene.matrix.loc)
#dropping NA's


gene.loc.granges <- GRanges(seqnames=paste0("chr", gene.matrix.loc$Chromosome.Name),
                            ranges=IRanges(start=gene.matrix.loc$Gene.Start..bp.,
                                           end=gene.matrix.loc$Gene.End..bp.,
                                           names=gene.matrix.loc$GENEID),
                            strand=gene.matrix.loc$Strand,
                            EnsemblName=gene.matrix.loc$Associated.Gene.Name,
                            TCGA=gene.matrix.loc$SYM)

require(ChIPpeakAnno)
cpg.gene.loc.granges <- annotatePeakInBatch(RangedData(cpg.granges), AnnotationData=RangedData(gene.loc.granges), output="nearestStart")

#Please Don't convert like this - sorry - i only realized afterwards...
cpg.gene.loc.df <- as.data.frame(cpg.gene.loc.granges)
dimnames(gene.matrix)[[1]] <- gene.matrix$GENEID

cpg.matrix.sm <- merge(cpg.gene.loc.df, cpg.matrix.s, by.x = "peak", by.y = "Composite.Element.REF")

gene.names <- names(gene.matrix)[3:267]
cpg.names <- names(cpg.matrix.sm)[18:306]
id.common <- intersect(cpg.names, gene.names)
cpg.names.loc <- c(names(cpg.matrix.sm)[1:17], id.common)
gene.names.loc <- c(names(gene.matrix)[1:2], id.common)
cpg.matrix.sm <- cpg.matrix.sm[,cpg.names.loc]
gene.matrix.sm <- gene.matrix[,gene.names.loc]
dimnames(cpg.matrix.sm)[[1]] <- paste(cpg.matrix.sm$peak, cpg.matrix.sm$feature, sep=".")

##remove all NAs in DNA methylation data (SNPs).
tmp <- cpg.matrix.sm[,id.common]
tmp <- na.omit(tmp)
clean.rows <- row.names(tmp)
cpg.matrix.sm <- cpg.matrix.sm[clean.rows, ]

## merge gene and dna methylation
cpg.gene <- merge(cpg.matrix.sm, gene.matrix.sm, by.x = "feature", by.y = "GENEID", all.x = T)
cpg.gene <- cpg.gene[, c(1:17,279,18:278,280:540)] # reorder columns to put annotation in the first 18 columns.
dimnames(cpg.gene)[[1]] <- paste(cpg.gene[,"peak"], cpg.gene[,"SYM"], sep=".")

dnamethylation.samples <- c(19:279)
geneexpression.samples <- c(280:540)
cpg.gene$corValue <- NA
cpg.gene$corpValue <- NA
tmp <- matrix(NA, nrow = dim(cpg.gene)[1], ncol=2)
dimnames(tmp)[[1]] <- dimnames(cpg.gene)[[1]]
dimnames(tmp)[[2]] <- c("corValue", "corpValue")


## need to use mclapply...
for ( i in 1:dim(cpg.gene)[1] ){
  e <- t(cpg.gene[i, geneexpression.samples])
  m <- t(cpg.gene[i, dnamethylation.samples])
  cor.value <- cor(e, m, method="spearman",use="complete.obs")
  cor.pvalue <- cor.test(e, m, method="spearman",alternative="two.sided")
  tmp[i, "corValue"] <- cor.value[1, 1]
  tmp[i, "corpValue"] <- cor.pvalue$p.value
  cat("completed", i, "\n")
  rm(e, m, cor.value, cor.pvalue)
}

cpg.gene$corpValue <- tmp[,"corpValue"]
cpg.gene$corValue <- tmp[,"corValue"]
cpg.gene$corpValue.adj <- p.adjust(cpg.gene$corpValue, method="BH")


meth.g <- cpg.gene[, dnamethylation.samples]
names(meth.g) <- substr(names(meth.g), 1, 12)
expr.g <- cpg.gene[, geneexpression.samples]
names(expr.g) <- substr(names(expr.g), 1, 12)


#tmp <- matrix(NA, nrow = dim(cpg.gene)[1], ncol=2)
#dimnames(tmp)[[1]] <- dimnames(cpg.gene)[[1]]
#dimnames(tmp)[[2]] <- c("corValue", "corpValue")

epiGene <- list()
j=1

###mclapply would make this much faster...currently can take up to 3 hours to complete 350K rows.
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
        yy <- merge(samps.purity, yy, by.x = "new", by.y = "id", all.x=T)
        epiGene[[j]] <- yy
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

### start from here
setwd("/dados/normalized.geneExp/rename")
load("LGG-meth.expr.analysis_05032014.Rda")


require(stringr)
promoter.probes <- rownames(subset(cpg.gene, shortestDistance <= 5000))
gene.id.prom <- str_split(as.character(intersect(names(epiGene), promoter.probes)), pattern="\\.")
gene.id.prom <- do.call(rbind.data.frame, gene.id.prom)
colnames(gene.id.prom)  <- c("cpgID1", "GENEID1", "cpgID2", "GENEID2")

all.probes <- rownames(subset(cpg.gene, shortestDistance <= 5000))
all.probes.id <- str_split(as.character(all.probes), pattern="\\.")
all.probes.id <- do.call(rbind.data.frame, all.probes.id)
colnames(all.probes.id)  <- c("cpgID1", "GENEID1", "cpgID2", "GENEID2")

gene.id.prom.count <- unlist(lapply(split(gene.id.prom$GENEID1,f=gene.id.prom$GENEID1),length))
all.probes.id.count <- unlist(lapply(split(all.probes.id$GENEID1,f=all.probes.id$GENEID1),length))
gene.id.df <- as.data.frame(gene.id.prom.count)
all.id.df <- as.data.frame(all.probes.id.count)


all.gene.id <- merge(gene.id.df, all.id.df, by = 0, all.x =T)
all.gene.id$prop <- all.gene.id$gene.id.prom.count / all.gene.id$all.probes.id.count



tmp <- subset(all.gene.id, prop > 0.5)
tmp$type <- TRUE
tmp <- tmp[,c(1,5)]
tmps <- merge(gene.id.prom, tmp, by.x="GENEID1", by.y="Row.names", all.x =T)
tmps <- subset(tmps, type==TRUE)
tmps$newID <- paste(tmps$cpgID1, tmps$GENEID1,  sep=".")

epiGene.s <- epiGene[tmps$newID]
save(epiGene.s, file="epiGene.s.Rda")

load(file="epiGene.s.Rda")

uu <- epiGene.s[[1]]

require(parallel)

clusters <- c(paste0("K", 1:5))

enrichment <- lapply(clusters, function(cluster, epiGene.ss=epiGene.s) {
  enrichment.group <- unlist(mclapply(epiGene.ss, function(uu, cluster.group=cluster) {
    uu$es <- uu$meth.gt0.3 & uu$exp.lt.meanUnmeth
    uu$K4class <- uu$K5.purity == cluster.group
    es.table <- table(uu$K4class, uu$es)#, dnn=c("Class", "ES"))
    if((es.table[2,2]/rowSums(es.table)[2] > 0.5) & (es.table[2,2]/colSums(es.table)[2] > 0.5)) {
      return(fisher.test(table(uu$K4class, uu$es))$p.value)
    } else {
      return(NULL)
    }
  }, mc.cores=1))
  return(enrichment.group)
})

K4enrichment <- unlist(mclapply(epiGene.s, function(uu) {
  uu$es <- uu$meth.gt0.3 & uu$exp.lt.meanUnmeth
  uu$K4class <- uu$K5.purity == "K4"
  es.table <- table(uu$K4class, uu$es, dnn=c("Class", "ES"))
  if(es.table[2,2]/rowSums(es.table)[2] > 0.5) {
    if(es.table[2,2]/colSums(es.table)[2] > 0.5) {
      return(fisher.test(table(uu$K4class, uu$es))$p.value)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}, mc.cores=4))


K3enrichment <- unlist(mclapply(epiGene.s, function(uu) {
  uu$es <- uu$meth.gt0.3 & uu$exp.lt.meanUnmeth
  uu$K3class <- uu$K5.purity == "K3"
  fisher.test(table(uu$K3class, uu$es))$p.value
}, mc.cores=4))

K2enrichment <- unlist(mclapply(epiGene.s, function(uu) {
  uu$es <- uu$meth.gt0.3 & uu$exp.lt.meanUnmeth
  uu$K2class <- uu$K5.purity == "K2"
  fisher.test(table(uu$K2class, uu$es))$p.value
}, mc.cores=4))

K1enrichment <- unlist(mclapply(epiGene.s, function(uu) {
  uu$es <- uu$meth.gt0.3 & uu$exp.lt.meanUnmeth
  uu$K1class <- uu$K5.purity == "K1"
  fisher.test(table(uu$K1class, uu$es))$p.value
}, mc.cores=4))

K5enrichment <- unlist(mclapply(epiGene.s, function(uu) {
  uu$es <- uu$meth.gt0.3 & uu$exp.lt.meanUnmeth
  uu$K5class <- uu$K5.purity == "K5"
  fisher.test(table(uu$K5class, uu$es))$p.value
}, mc.cores=4))

enrichmentList <- list(K1=K1enrichment, K2=K2enrichment, K3=K3enrichment, K4=K4enrichment, K5=K5enrichment)
names(enrichment) <- c("K1", "K2", "K3", "K4", "K5")
enrichment <- lapply(enrichment, p.adjust, method="BH")
enrichmentList.sig <- lapply(enrichment, function(x) {
  x <- x[x < 0.05]
})

# enrichmentList.count <- lapply(K4enrichment.count, function(x) {
#   x <- x == TRUE
# })



require(stringr)
enrichmentList.geneNames <- lapply(enrichmentList.sig, function(x) {
  x <- str_split(as.character(names(x)), pattern="\\.")
  x <- do.call(rbind.data.frame, x)
  names(x) <- c("cpgID", "GeneName")
  return(x)
})

enrichmentList.geneNames2 <- do.call(rbind, enrichmentList.geneNames)
names <- str_split(row.names(enrichmentList.geneNames2), pattern="\\.")
names <- sapply(names, "[", 1)
enrichmentList.geneNames2$cluster<- names
save(enrichmentList.geneNames2, file="enrichment.cluster.cpg.gene.Rda")

require(venneuler)
require(VennDiagram)
K4enrichment.sig <- K4enrichment[K4enrichment < 0.05]

### finding the lowest correlation value to define ES.
require(parallel)
epirho <- mclapply(epiGene.s, function(ii) {
  e <- (ii[, "exprValues"])
  m <- (ii[, "methValues"])
  cor.pvalue <- cor.test(e, m, method="spearman",alternative="two.sided")
  return(c(rho=as.numeric(cor.pvalue$estimate), p.value=cor.pvalue$p.value))
}, mc.cores=4)

epirho <- do.call(rbind.data.frame, epirho)
colnames(epirho) <- c("rho", "p.value")

require(stringr)
gene.names <- str_split(rownames(epirho), pattern="\\.")
gene.names <- sapply(gene.names, '[', 2)
epirho$gene.names <- gene.names

gene.names.rle <- rle(gene.names)

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


color.map.meth.cl <- function(mol.biol) {
  if (mol.biol=="M5") "darkorchid" #darkorchid
  else if (mol.biol=="M3") "blue" #blue
  else if (mol.biol=="M4") "darkgreen" #darkgreen
  else if (mol.biol=="M1") "black" #black
  else if (mol.biol=="M2") "red" #red
  else if (mol.biol==6) "yellow"
  else "#FFFFFF"
} #cluster.consensus.clustering
color.map.batch <- function(mol.biol) {
  if (mol.biol=="Batch_078") "#9488f9" #lightpurpleTOdarkpurple
  else if (mol.biol=="Batch_112") "#857ae0" #lightpurpleTOdarkpurple
  else if (mol.biol=="Batch_146") "#766cc7" #lightpurpleTOdarkpurple
  else if (mol.biol=="Batch_163") "#675fae" #lightpurpleTOdarkpurple
  else if (mol.biol=="Batch_189") "#585195" #lightpurpleTOdarkpurple
  else if (mol.biol=="Batch_219") "#4a447c" #lightpurpleTOdarkpurple
  else if (mol.biol=="Batch_245") "#3b3663" #lightpurpleTOdarkpurple
  else if (mol.biol=="Batch_282") "#2c284a" #lightpurpleTOdarkpurple
  else if (mol.biol=="Batch_292") "#1d1b31" #lightpurpleTOdarkpurple
  else if (mol.biol=="Batch_295") "#0e0d18" #lightpurpleTOdarkpurple
  else if (mol.biol=="Batch_306") "#000000" #lightpurpleTOdarkpurple
  else "#800000" #moroon
} #Batch
color.map.histo<-function(mol.biol) {
  if (mol.biol=="Astrocytoma") "yellow3"
  else if (mol.biol=="Oligoastrocytoma") "tan4"
  else if (mol.biol=="Oligodendroglioma") "olivedrab3"
  else "#FFFFFF"
}#Gene Expression Cluster
color.map.grade<-function(mol.biol) {
  if (mol.biol=="G2") "White"
  else if (mol.biol=="G3") "Black"
  else "Gray"
}#grade
color.map.mutation<-function(mol.biol) {
  if (mol.biol=="Mutant") "black"
  else if (mol.biol=="wt") "Gray"
  else "White"
}#
color.map.IDH1p19q<-function(mol.biol) {
  if (mol.biol=="TRUE") "black"
  else if (mol.biol=="FALSE") "Gray"
  else "White"
}#
color.map.RNAseq<-function(mol.biol) {
  if (mol.biol==1) "blue" #ok
  else if (mol.biol==2) "red" #ok
  else if (mol.biol==3) "darkorchid" #ok
  else if (mol.biol==4) "darkgreen" # decided these were not real clusters
  else if (mol.biol==5) "#FF8C00" #ok
  else if (mol.biol==6) "yellow" # decided these were not real clusters
  else if (mol.biol==7) "#DC143C"
  else "#FFFFFF" 
}#

color.map.Oncosign<-function(mol.biol) {
  if (mol.biol=="OSC1") "#0D9344" #ok
  else if (mol.biol=="OSC2") "#EE1D23" #ok
  else if (mol.biol=="OSC3") "#24AAE1" #ok
  else if (mol.biol=="Unclassified") "lightgrey"
  else "#FFFFFF" 
}#
color.map.gender<-function(mol.biol) {
  if (mol.biol=="MALE") "lightblue" #ok
  else if (mol.biol=="FEMALE") "pink" #ok
  else "#FFFFFF" 
}#
color.map.family<-function(mol.biol) {
  if (mol.biol=="NO") "grey" #ok
  else if (mol.biol=="YES") "black" #ok
  else "#FFFFFF" 
}#
color.map.location<-function(mol.biol) {
  if (mol.biol=="Posterior Fossa, Cerebellum") "#857ae0" #ok
  else if (mol.biol=="Supratentorial, Frontal Lobe") "#EFCDB8" #ok
  else if (mol.biol=="Supratentorial, Not Otherwise Sp") "#FDDB6D" #ok
  else if (mol.biol=="Supratentorial, Occipital Lobe") "#585195" #ok
  else if (mol.biol=="Supratentorial, Parietal Lobe") "#4a447c" #ok
  else if (mol.biol=="Supratentorial, Temporal Lobe") "#9F8170" #ok
  else "#FFFFFF" 
}#
color.map.symptom<-function(mol.biol) {
  if (mol.biol=="Headaches") "#A8E4A0" #ok
  else if (mol.biol=="Mental Status Changes") "#0000CD" #ok
  else if (mol.biol=="Motor/Movement Changes") "#00BFFF" #ok
  else if (mol.biol=="Seizures") "#4682B4" #ok
  else if (mol.biol=="Sensory Changes") "#00CED1" #ok
  else if (mol.biol=="Visual Changes") "#AFEEEE" #ok
  else "#FFFFFF" 
}#
color.map.ages<-function(mol.biol) {
  if (mol.biol=="(13,31.2]") "red" #ok
  else if (mol.biol=="(31.2,40]") "orange" #ok
  else if (mol.biol=="(40,54]") "green" #ok
  else if (mol.biol=="(54,75]") "blue" #ok
  else "#FFFFFF" 
}#
color.map.type<-function(mol.biol) {
  if (mol.biol=="GBM") "red" #ok
    else if (mol.biol=="(31.2,40]") "orange" #ok
    else if (mol.biol=="LGG") "blue" #ok
  else "#FFFFFF" 
}#
source("Func_List.r")

load(file="LGG.PurityCorrection-POSTCLUSTER.RData")
load(file="LGG.puritycorrection.HeatmapOrder.Rda")
## consensus cluster
temp <- x[, LGG.puritycorrection.HeatmapOrder]
temp <- as.data.frame(temp)
##to remove 28 samples not profiled due to no expression data
temp.t <- t(temp)
temp.t <- na.omit(temp.t)
temp <- t(temp.t)
temp <- as.data.frame(temp)

samps.purity <- epiGene.s[[1]]
##updated Annotations
#synapse
require(synapseClient)
synapseLogin()
# Load the annotationUpdate
syn2369104 <- synGet(id='syn2369104', load=T) #updated march 14
#syn2366689 <- synGet(id='syn2366689', load=T) #updated march 10
#syn2359167 <- synGet(id='syn2359167', load=T)
#load(file="~/.synapseCache/172/479172/LGG-Annotation.2-28-2014.RData")
#load(file="~/.synapseCache/313/488313//LGG-Annotation.3-10-2014.RData") #updated march 10
load(file="~/.synapseCache/615/490615/LGG-Annotation.3-14-2014.RData") #updated march 14
## pull in the 7 miRNA cluster
syn2362918 <- synGet(id='syn2362918', load=T)
miR7 <- read.delim("~/.synapseCache/713/484713//LGG.miRNAseq.NMF-7groups.cluster-calls.BCGSC.20140115.txt")
miR7$new <- substr(miR7$ID, 1, 12)
names(miR7)[4] <- "miRNACluster7"
Clustering <- merge(Clustering, miR7[,c(7,4)], by.x="Tumor", by.y = "new")

##recurrent status
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

i=2
names(clinical.cluster.genetics)[i];
table(as.factor(clinical.cluster.genetics[,i]), clinical.cluster.genetics$MethylationCluster); 
i <- i+1
require(reshape)
library(scales) #important for percent
##methylation by somatic mutation
ccg.m <- melt(clinical.cluster.genetics[,c(4, 100:165)], id.vars="MethylationCluster")
ccg.m <- na.omit(ccg.m)
c <- ggplot(ccg.m, aes(factor(MethylationCluster), fill=value, colour=MethylationCluster))
c +
  scale_colour_manual(values = c("M1" = "Black","M2" = "Red","M3" = "Blue", "M4" = "Darkgreen", "M5" = "Darkorchid")) +
  scale_fill_manual(values = c("Mutant" = "black", "wt" = "gray")) +
  #geom_bar() +
  geom_bar(position="fill") +
  labs(title = "DNA Methylation Cluster ~ Somatic Mutation", colour = "DNA Meth. Cluster", fill = "Mutation") +
  xlab("DNA Methylation Cluster") +
  ylab("Percent Total") +
  scale_y_continuous(labels = percent) +
  facet_wrap(~ variable, nrow=7) 
ggsave(width=15, height=15, units="in",  dpi=200, file="DNAmethylation.SomaticMutation-fill.png")
#ggsave(file="DNAmethylation.SomaticMutation-fill.png")

##ID+1p19q
ccg.m <- melt(clinical.cluster.genetics[,c(4, 2)], id.vars="MethylationCluster")
ccg.m <- na.omit(ccg.m)
c <- ggplot(ccg.m, aes(factor(MethylationCluster), fill=value
                       #, colour=MethylationCluster
                       ))

c +
  #scale_colour_manual(values = c("K1" = "blue","K2" = "pink","K3" = "red", "K4" = "black", "K5" = "green")) +
  scale_fill_manual(values = c("IDHmut-codel" = "purple4",
                               "IDHmut-non-codel" = "darkcyan",
                               "IDHwt" = "red4")) +
  #geom_bar() +
  geom_bar(position="fill") +
  labs(title = "DNA Methylation Cluster ~ IDH+1p19q Mutation", colour = "DNA Meth. Cluster", fill = "Mutation") +
  xlab("DNA Methylation Cluster") +
  ylab("Percent Total") +
  scale_y_continuous(labels = percent)
ggsave(width=15,height=15, units="in",  dpi=200, file="DNAmethylation.IDH1p19qMutation-fill.png")

##methylation by clusters
clst <- (clinical.cluster.genetics[,c("MethylationCluster", 
                                      "RNASeqCluster", 
                                      "CNCluster", 
                                      "miRNACluster", 
                                      "RPPACluster", 
                                      "COCCluster", 
                                      "Oncosign")])
# clst[clst==1] <- "c1"
# clst[clst==2] <- "c2"
# clst[clst==3] <- "c3"
# clst[clst==4] <- "c4"
# clst[clst==5] <- "c5"
# clst[clst==6] <- "c6"
# clst[clst==7] <- "c7"
# clst$OncosignCluster <- as.character(clst$OncosignCluster)
# clst[clst=="OSC1"] <- "c1"
# clst[clst=="OSC2"] <- "c2"
# clst[clst=="OSC3"] <- "c3"
ccg.m <- melt(clst, id.vars="MethylationCluster")
ccg.m <- na.omit(ccg.m)
c <- ggplot(ccg.m, aes(factor(MethylationCluster), fill=value))
c +
  #scale_colour_manual(values = c("K1" = "blue","K2" = "pink","K3" = "red", "K4" = "black", "K5" = "green")) +
  scale_fill_manual(values = c("R1" = "red","R2" = "#0D9344","R3" = "#24AAE1", "R4" = "#223F99", 
                               "C1" = "red","C2" = "#0D9344","C3" = "#24AAE1", "C4" = "#223F99",
                               "mi1" = "red","mi2" = "#0D9344","mi3" = "#24AAE1", "mi4" = "#223F99",
                               "P1" = "red","P2" = "#0D9344","P3" = "#24AAE1", "P4" = "#223F99",
                               "coc1" = "red","coc2" = "#0D9344","coc3" = "#24AAE1", 
                               "OSC1" = "red","OSC2" = "#0D9344","OSC3" = "#24AAE1", "Unclassified" = "lightgrey" 
                               )) +
  #geom_bar() +
  geom_bar(position="fill") +
  labs(title = "DNA Methylation Cluster ~ Other Clusters", colour = "DNA Meth. Cluster", fill = "Other Clusters") +
  xlab("DNA Methylation Cluster") +
  ylab("Percent Total") +
  scale_y_continuous(labels = percent) +
  guides(fill = guide_legend(ncol = 2)) +
  facet_grid(.~variable) 
ggsave(width=15,height=15, units="in", dpi=200, file="DNAmethylation.Clusters-fill.png")

##methylation by clinical features
yr <- c("MethylationCluster", 
        "gender","neoplasm_histologic_grade","histological_type",
        "family_history_of_cancer", "tumor_location", "first_presenting_symptom")
ccg.m <- melt(clinical.cluster.genetics[, yr], id.vars="MethylationCluster")
ccg.m <- na.omit(ccg.m)
c <- ggplot(ccg.m, aes(factor(MethylationCluster), fill=value))
c +
  #scale_colour_manual(values = c("K1" = "blue","K2" = "pink","K3" = "red", "K4" = "black", "K5" = "green")) +
   scale_fill_manual(values = c("FEMALE" = "pink",
                                "MALE" = "lightblue",
                                "G2" = "white", 
                                "G3" = "black", 
                                "Astrocytoma" = "yellow3",
                                "Oligoastrocytoma" = "tan4", 
                                "Oligodendroglioma" = "olivedrab3",
                                "NO" = "grey", 
                                "YES" = "black", 
                                "Posterior Fossa, Cerebellum"    =     "#857ae0",
                                "Supratentorial, Frontal Lobe"  =    "#EFCDB8",
                                "Supratentorial, Not Otherwise Sp"   = "#FDDB6D",
                                "Supratentorial, Occipital Lobe"    =  "#585195",
                                "Supratentorial, Parietal Lobe"    =  "#4a447c",
                                "Supratentorial, Temporal Lobe"   =   "#9F8170",
                                "Headaches"         = "#A8E4A0",
                                "Mental Status Changes" = "#0000CD",
                                "Motor/Movement Changes" ="#00BFFF",
                                "Seizures"          = "#4682B4",
                                "Sensory Changes" = "#00CED1",
                                "Visual Changes" = "#AFEEEE"
                                ),
                     labels = c("Gender: FEMALE", "Gendere: MALE",
                                "Grade: 2", "Grade: 3", 
                                "Histology: Astrocytoma", "Histology: Oligoastrocytoma", "Histology: Oligodendroglioma",
                                "History: No", "History: Yes",
                                "Location: Cerebellum", "Location: Fronal Lobe", "Location: UNK", "Location: Occipital Lobe", 
                                "Location: Parietal Lobe", "Location: Temporal Lobe",
                                "1st Symptom: Headaches", "1st Symptom: Mental Change", "1st Symptom: Moter/Move. Change", "1st Symptom: Seizures", "1st Symptom: Sensory Change", "1st Symptom: Visual Change"

)
                     ) +
  #geom_bar() +
  geom_bar(position="fill") +
  labs(title = "DNA Methylation Cluster ~ Clinical Features", colour = "DNA Meth. Cluster", fill = "Clinical") +
  xlab("DNA Methylation Cluster") +
  ylab("Percent Total") +
  scale_y_continuous(labels = percent) +
  guides(fill = guide_legend(ncol = 2)) +
  facet_wrap(~ variable) 
ggsave(width=15,height=10, units="in", dpi=200, file="DNAmethylation.Clinical-fill.png")
##interesting annotations associated with K4

##add in age
agecat <- cut(clinical.cluster.genetics$age_at_initial_pathologic, c(13,31.25,40,54,75))
clinical.cluster.genetics$age.strat <- as.character(agecat)


conclus.5<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"MethylationCluster")
conclus.5<-as.character(conclus.5)
conclus.5[is.na(conclus.5)]<-c("0")
cc.conclus.5<-unlist(lapply(conclus.5, color.map.meth.cl))

#age
ages<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"age.strat")
ages<-as.character(ages)
ages[is.na(ages)]<-c("0")
cc.ages<-unlist(lapply(ages, color.map.ages))

#rnacluster
rnacluster<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"RNASeqCluster")
rnacluster<-as.character(rnacluster)
rnacluster[is.na(rnacluster)]<-c("0")
cc.rnacluster<-unlist(lapply(rnacluster, color.map.RNAseq))

#CNCcluster
CNCcluster<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics, "CNCluster")
CNCcluster<-as.character(CNCcluster)
CNCcluster[is.na(CNCcluster)]<-c("0")
cc.CNCcluster<-unlist(lapply(CNCcluster, color.map.RNAseq))

#miRNAcluster
miRNAcluster<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics, "miRNACluster")
miRNAcluster<-as.character(miRNAcluster)
miRNAcluster[is.na(miRNAcluster)]<-c("0")
cc.miRNAcluster<-unlist(lapply(miRNAcluster, color.map.RNAseq))

#miRNAcluster7
miRNAcluster7<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics, "miRNACluster7")
miRNAcluster7<-as.character(miRNAcluster7)
miRNAcluster7[is.na(miRNAcluster7)]<-c("0")
cc.miRNAcluster7<-unlist(lapply(miRNAcluster7, color.map.RNAseq))

#RPPAcluster
RPPAcluster<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics, "RPPACluster")
RPPAcluster<-as.character(RPPAcluster)
RPPAcluster[is.na(RPPAcluster)]<-c("0")
cc.RPPAcluster<-unlist(lapply(RPPAcluster, color.map.RNAseq))

#Oncosigncluster
OSCcluster<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics, "Oncosign")
OSCcluster<-as.character(OSCcluster)
OSCcluster[is.na(OSCcluster)]<-c("0")
cc.OSCcluster<-unlist(lapply(OSCcluster, color.map.Oncosign))

#histo
histo<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"histological_type")
histo<-as.character(histo)
histo[is.na(histo)]<-c("0")
cc.histo<-unlist(lapply(histo, color.map.histo))

#grade
grade<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"neoplasm_histologic_grade")
grade<-as.character(grade)
grade[is.na(grade)]<-c("0")
cc.grade<-unlist(lapply(grade, color.map.grade))

#gender
gender<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"gender")
gender<-as.character(gender)
gender[is.na(gender)]<-c("0")
cc.gender<-unlist(lapply(gender, color.map.gender))

#history
family<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"family_history_of_cancer")
family<-as.character(family)
family[is.na(family)]<-c("0")
cc.family<-unlist(lapply(family, color.map.family))

#location
location<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"tumor_location")
location<-as.character(location)
location[is.na(location)]<-c("0")
cc.location<-unlist(lapply(location, color.map.location))

#symptom
symptom<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"first_presenting_symptom")
symptom<-as.character(symptom)
symptom[is.na(symptom)]<-c("0")
cc.symptom<-unlist(lapply(symptom, color.map.symptom))

#batch
batch<-apply(as.matrix(names(temp)),1,func.list$vlookup,subset(samps.purity, HumanMethylation450==1),"batch")
batch<-as.character(batch)
batch[is.na(batch)]<-c("0")
cc.batch<-unlist(lapply(batch, color.map.batch))

#idh1
idh1<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"IDH1")
idh1<-as.character(idh1)
idh1[is.na(idh1)]<-c("0")
cc.idh1<-unlist(lapply(idh1, color.map.mutation))

#idh2
idh2<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"IDH2")
idh2<-as.character(idh2)
idh2[is.na(idh2)]<-c("0")
cc.idh2<-unlist(lapply(idh2, color.map.mutation))

#1p19q
Chr1p19q<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"X1p_19q_co_del_status")
Chr1p19q<-as.character(Chr1p19q)
Chr1p19q[is.na(Chr1p19q)]<-c("0")
cc.Chr1p19q<-unlist(lapply(Chr1p19q, color.map.IDH1p19q))

#atrx
atrx<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"ATRX")
atrx<-as.character(atrx)
atrx[is.na(atrx)]<-c("0")
cc.atrx<-unlist(lapply(atrx, color.map.mutation))

#cic
cic<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"CIC")
cic<-as.character(cic)
cic[is.na(cic)]<-c("0")
cc.cic<-unlist(lapply(cic, color.map.mutation))

#EGFR
egfr<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"EGFR")
egfr<-as.character(egfr)
egfr[is.na(egfr)]<-c("0")
cc.egfr<-unlist(lapply(egfr, color.map.mutation))

#NF1
nf1<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"NF1")
nf1<-as.character(nf1)
nf1[is.na(nf1)]<-c("0")
cc.nf1<-unlist(lapply(nf1, color.map.mutation))

#NOTCH1
notch1<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"NOTCH1")
notch1<-as.character(notch1)
notch1[is.na(notch1)]<-c("0")
cc.notch1<-unlist(lapply(notch1, color.map.mutation))

#FUBP1
fubp1<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"FUBP1")
fubp1<-as.character(fubp1)
fubp1[is.na(fubp1)]<-c("0")
cc.fubp1<-unlist(lapply(fubp1, color.map.mutation))

#PTEN
pten<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"PTEN")
pten<-as.character(pten)
pten[is.na(pten)]<-c("0")
cc.pten<-unlist(lapply(pten, color.map.mutation))

#p53
p53<-apply(as.matrix(names(temp)),1,func.list$vlookup,clinical.cluster.genetics,"TP53")
p53<-as.character(p53)
p53[is.na(p53)]<-c("0")
cc.p53<-unlist(lapply(p53, color.map.mutation))

#cc.col<-matrix(as.character(c(cc.batch,cc.p53,cc.pten,cc.nf1,cc.idh1,cc.egfr,cc.expclus,cc.gcimp,cc.conclus)), nrow=288,ncol=9)
cc.col<-matrix(as.character(c(cc.batch,
                              cc.ages,
                              cc.location,
                              cc.family,
                              cc.gender,
                              cc.symptom,
                              cc.histo, 
                              cc.grade, 
                              cc.p53, 
                              cc.pten, 
                              cc.notch1, 
                              cc.nf1, 
                              cc.fubp1, 
                              cc.egfr, 
                              cc.cic, 
                              cc.atrx, 
                              cc.idh2, 
                              cc.idh1, 
                              cc.Chr1p19q, 
                              cc.OSCcluster,
                              cc.RPPAcluster,
                              cc.miRNAcluster7,
                              cc.miRNAcluster,
                              cc.CNCcluster,
                              cc.rnacluster,
                              cc.conclus.5)), 
               nrow=289,
               ncol=26)
#colnames(cc.col) = c("BATCH","TP53","PTEN","NF1","IDH1","EGFR","GeneExp CC", "Noushmehr et al - GCIMP", "DNA MethCC")
colnames(cc.col) = c("BATCH",
                     "ages",
                     "Location",
                     "Family History",
                     "Gender",
                     "Symptom",
                     "Histo", 
                     "Grade",
                     "p53",
                     "PTEN",
                     "NOTCH1",
                     "NF1",
                     "FUBP1", 
                     "EGFR", 
                     "CIC", 
                     "ATRX",
                     "IDH2",
                     "IDH1",
                     "1p19q", 
                     "OSC CC", 
                     "RPPA CC", 
                     "miRNA CC (7)" ,
                     "miRNA CC (4)", 
                     "CN CC", 
                     "RNA CC (4/6)", 
                     "DNA CC (5)")

#library(RColorBrewer)
#hmColors <- colorRampPalette(c("darkblue", "yellow"))(256)
library(heatmap.plus)
library(matlab)
source(file="heatmap.plus.R")
temp.1 <- temp * 1
temp.1 <- temp.1 + 1
temp.1[is.na(temp.1)] <- 0

## redo the heatmap. Load in the RData file to include datT
dimnames(datT)[[2]] <- substr(dimnames(datT)[[2]], 1, 12)
load(file="LGG.puritycorrection.HeatmapOrder.rows.Rda")
datT <- datT[LGG.puritycorrection.HeatmapOrder.rows, LGG.puritycorrection.HeatmapOrder]

#png(filename = "esheatmap-ordered-Ks_ver2d.png", bg="white", res=200, width=2000, height=2000) ## ES calls (temp.1)
png(filename = "K5-sd0.3,LGG,batches78,112,146,163,189,219,245,282,292,295,306-ccmatrix_using723_probes_PURITY_ver2d.png", bg="white", res=200, width=2000, height=2000)
pdf(file = "K5-sd0.3,LGG,batches78,112,146,163,189,219,245,282,292,295,306-ccmatrix_using723_probes_PURITY_ver2i.pdf", width=10, height=10)
pdf(file="oncosignTrack.pdf", width=10, height=10)
hv2.lgg<-heatmap.plus.sm(
  #as.matrix(temp.1),
  as.matrix(datT[]),
  na.rm=TRUE,
  scale="none",
  #RowSideColor=probe.cc,
  #ColSideColors=cc.col[,1:5],
  #ColSideColors=cc.col[,6:10],
  #ColSideColors=cc.col[,11:15],
  ColSideColors=cc.col[,c(20,26)],
  #ColSideColors=cc.col[,1:6], ##added in additional clinical
  #col=c("white", "#A0A0A0", "black"),
  col=jet.colors(75),
  #key=TRUE,
  symkey=FALSE,
  density.info="none",
  trace="none",
  Rowv=TRUE,
  Colv=NA,
  cexRow=1,
  cexCol=0.6,
  keysize=1,
  dendrogram=c("none"),
  #main = paste("LGG DNA meth Clustering, K5 ",dim(temp)[1],"P.; ",dim(as.matrix(temp[,(as.character(cl.5.sort[,1]))]))[2],"S."),
  #labCol=NA
  labRow=NA
)
dev.off()



##get probe colors
## rerun above using heatmap.plus and then you have hv2.lgg
require(stringr)
temp.names <- str_split(rownames(temp.1), pattern="\\.")
temp.names <- sapply(temp.names, '[', 2)
temp.2 <- temp.1
dimnames(temp.2)[[1]] <- temp.names


library(FDb.InfiniumMethylation.hg19)
hm450.hg19 <- getPlatform(platform='HM450', genome='hg19')
data(hg19.islands)
CGI.probes <- subsetByOverlaps(hm450.hg19, hg19.islands)
head(CGI.probes)
tail(CGI.probes)
hg19.shores <- c(flank(hg19.islands, 2000, start=TRUE),
                 flank(hg19.islands, 2000, start=FALSE))
shore.probes <- subsetByOverlaps(hm450.hg19, hg19.shores)
head(shore.probes)
tail(shore.probes)

shore.probes <- as.data.frame(shore.probes)
shore.probes$ID <- rownames(shore.probes)
shore.probes$Type <- "Shores"
shore.probes <- shore.probes[,c("ID","Type")]
CGI.probes <- as.data.frame(CGI.probes)
CGI.probes$ID <- rownames(CGI.probes)
CGI.probes$Type <- "Islands"
CGI.probes <- CGI.probes[,c("ID","Type")]

shore.probes[intersect(shore.probes[,1], CGI.probes[,1]),"Type"] <- "Islands"
shore.probes <- subset(shore.probes, Type!="Islands")
cgi.id <- rbind(CGI.probes, shore.probes)
cgi.id <- as.data.frame(cgi.id)

probe.labelsbyK <- read.delim(file="venn/probe_labels.uniqueKs.txt", sep="\t")
dimnames(probe.labelsbyK)[[1]] <- probe.labelsbyK$Gene

color.map.probe.cl <- function(mol.biol) {
  if (mol.biol=="Islands") "darkorchid" #darkorchid
  else if (mol.biol=="Shores") "blue" #blue
  else if (mol.biol=="K4") "darkgreen" #darkgreen
  else if (mol.biol=="K1") "black" #black
  else if (mol.biol=="K2") "red" #red
  else if (mol.biol==6) "yellow"
  else "#FFFFFF"
} #cluster.consensus.clustering

kprobes<-apply(as.matrix(dimnames(temp.2[hv2.lgg[[1]], ])[[1]]),1,func.list$vlookup,probe.labelsbyK,"Clust")
kprobes<-as.character(kprobes)
kprobes[is.na(kprobes)]<-c("0")
cc.kprobes<-unlist(lapply(kprobes, color.map.probe.cl))

probe.col<-matrix(as.character(c(cc.kprobes, cc.kprobes)), nrow=1257,ncol=2)
colnames(probe.col) = c("Kcluster", "Enrichment")
png(filename = "esheatmap-ordered-Ks-withprobes.png", bg="white", res=200, width=2000, height=2000)
pdf(file = "K5-sd0.3,LGG,batches78,112,146,163,189,219,245,282,292,295,306-ccmatrix_using723_probes_PURITY_ver2g.pdf", width=10, height=10)

hv3.lgg<-heatmap.plus.sm(
  #as.matrix(temp[hv1[[1]],(as.character(cl.sort[,1]))]),
  #as.matrix(temp[,(as.character(cl.sort[,1]))]),
  as.matrix(temp.2[hv2.lgg[[1]], ]), 
  #temp.cc,
  na.rm=TRUE,
  scale="none",
  RowSideColor=probe.col,
  ColSideColors=cc.col,
  col=c("white", "#A0A0A0", "black"),
  #key=TRUE,
  symkey=FALSE,
  density.info="none",
  trace="none",
  Rowv=NA,
  Colv=NA,
  cexRow=1,
  cexCol=0.6,
  keysize=1,
  dendrogram=c("none"),
  #main = paste("Clustering, K6 ",dim(temp)[1],"Probes; ",dim(temp)[2],"Samples"),
  #labCol=NA
  labRow=NA
)
dev.off()






### plot specific genes
require(ggplot2)




## fishing for enrichments with more than 32 elemens in K4
K4enrichment.count <- unlist(mclapply(epiGene.s, function(uu) {
  uu$es <- uu$meth.gt0.3 & uu$exp.lt.meanUnmeth
  uu$K4class <- uu$K5.purity == "K4"
  as.data.frame(table(uu$K4class, uu$es))[4, 3] > 32 & as.data.frame(table(uu$K4class, uu$es))[3, 3] < 5
}, mc.cores=4))

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




##to find a gene/probe set
gene.of.interest <- "ZNF492"
Kclust.of.interest <- "K4"

#enrichmentList.geneNames[[Kclust.of.interest]][enrichmentList.geneNames[[Kclust.of.interest]][, 2] == gene.of.interest,]
#epirho.s[epirho.s$gene.names==gene.of.interest,]
gg <- ss.m[ss.m[,5]==gene.of.interest, "Row.names"]
gg <- gg[1]
gg <- rownames(epirho.s[epirho.s$gene.names==gene.of.interest,])

epiplot <- epiGene.s[[gg]]
epiplot$exprValues <- epiplot$exprValues + 1
expr.mean.gu <- mean(as.numeric(epiplot[epiplot$meth.gt0.3==FALSE, "exprValues"]), na.rm = T)

## K4 specific
p <- ggplot(epiplot, aes(methValues, exprValues))
p + 
  stat_smooth(method = "lm", fullrange=FALSE, na.rm=TRUE, se=FALSE, size=1) +
  geom_vline(xintercept = 0.3, linetype = "longdash", colour="black") +
  geom_hline(yintercept = expr.mean.gu, linetype = "longdash", colour="black") +
  geom_point(size=3, aes(colour = factor(K5.purity), shape = factor(X1p_19q_co_del_status))) +
  scale_colour_manual(values = c("M1" = "black","M2" = "red","M3" = "blue", "M4" = "darkgreen", "M5" = "darkorchid")) +
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.3, 0.5, 0.7, 1)) +
  scale_y_log10() +
  ylab("RNA seq expression (log2)") +
  xlab(paste("DNA Methylation Values\n", gg)) + 
  labs(title = paste("DNA Methylation vs Expression", gg, sep=" - "), colour = "DNA Meth. Cluster", shape = "1p19qStatus") 
ggsave(file=paste0("DNAmethylationVSexpression-", gene.of.interest, ".png"))



library(survival);
### progress free survival time

#gb.surv.cl <- survfit(Surv(daystolastordeath, death01) ~ factor(MethylationCluster), data = clinical.cluster.genetics)
#gg<-survdiff(Surv(ProgFreeSurvTime_days, ProgFreeSurvEvent01) ~ factor(K5.purity), data = samps.purity)
#x11()
png(filename = "K5puritycorrection-survival_plots-ProgFreeSurvTime_days.png", bg="white", res=300, width=2000, height=2000)
plot(gb.surv.cl, 
     lwd=5,
     col=c("black","red","blue","arkgreen","darkorchid"),
     main="Kaplan-Meier Survival Curves\nTCGA LGG SAMPLES (DNA Methylation Clusters)", 
     xlab="TIME SINCE DIAGNOSIS (YEARS)", 
     ylab="PERCENT PROBABILITY OF SURVIVAL", 
     yscale=100,
     xscale=365, 
     bg="black"
)
box(col="black", lwd=3);
abline(v=5*52, lwd=5, lty=2); text(4.9*52, 0.75, "5 YEARS", srt=90)
abline(h=0, lwd=5, lty=1)
legend("topright", legend=c(paste("K1 (n=",gg$n[1],")",sep=""),paste("K2 (n=",gg$n[2],")",sep=""),paste("K3 (n=",gg$n[3],")",sep=""),paste("K4 (n=",gg$n[4],")",sep=""),paste("K5 (n=",gg$n[5],")",sep="")),col=rev(c("GREEN","BLACK","RED","PINK","BLUE")),lwd=3,title="DNA METHYLATION Clusters",box.lwd=3, bg="white")
text(2,0.05,paste("Log-Rank p-value:",round(1 - pchisq(gg$chisq, length(gg$n) - 1),10)), col="black")
dev.off()

## days to last or death
## Survival Curves
## days to last or death
#adjust
agecat <- cut(clinical.cluster.genetics$age_at_initial_pathologic, c(13,31.25,40,54,75))
fdata <- clinical.cluster.genetics
fdata$sex <- fdata$gender
fdata$age <- fdata$age_at_initial_pathologic
fdata$age2 <- agecat
fdata$group <- as.factor(fdata$MethylationCluster)
data(uspop2)
refpop <- uspop2[as.character(14:75), , "2000"]
temp <- as.numeric(cut(14:75, c(13,31.25,40,54,75) + 0.5))
pi.us<- tapply(refpop, list(temp[row(refpop)], col(refpop)), sum)/sum(refpop)
tab2 <- with(fdata, table(age2, sex, group))/ nrow(fdata)
us.wt <- rep(pi.us, 5)/ tab2
index <- with(fdata, cbind(as.numeric(age2), as.numeric(sex), as.numeric(group)))
fdata$uswt <- us.wt[index]
sfit1a.d <- survfit(Surv(daystolastordeath, death01) ~ group, fdata)
sfit3a.d <-survfit(Surv(daystolastordeath, death01) ~ group, data=fdata, weight=uswt)

#sfit1a.p <- survfit(Surv(ProgFreeSurvTime_days, ProgFreeSurvEvent01) ~ group, fdata)
sfit3a.p <-survfit(Surv(ProgFreeSurvTime_days, ProgFreeSurvEvent01) ~ group, data=fdata, weight=uswt)

#plot(sfit3a, mark.time=F, col=c("red", "orange", "yellow", "green", "blue"), lty=1, lwd=2, xscale=365.25, xlab="Years from Sample", ylab="Survival")
#
#x11()
survLGG <- sfit3a.d

png(filename = "K5puritycorrection-survival_plots-OverallSurvival-adjustedbyAge.png", bg="white", res=300, width=2000, height=2000)
#par(mfrow=c(1,2))
plot(survLGG, 
     lwd=5,
     col=c("black","red","blue","darkgreen","darkorchid"),
     main="Overall Survival\nKaplan-Meier Survival Curves\nTCGA LGG SAMPLES (DNA Methylation Clusters)", 
     #main="Progression Free\nKaplan-Meier Survival Curves\nTCGA LGG SAMPLES (DNA Methylation Clusters)", 
     xlab="TIME SINCE DIAGNOSIS (YEARS)", 
     ylab="PERCENT PROBABILITY OF SURVIVAL", 
     yscale=100,
     xscale=365.25, 
     bg="black"
)
#lines(sfit1a.d, mark.time=T, col=c("black","red","blue","darkgreen","darkorchid"), lty=2, lwd=1, xscale=365.25)
box(col="black", lwd=3);
#abline(v=5*52, lwd=5, lty=2); text(4.9*52, 0.75, "5 YEARS", srt=90)
abline(h=0, lwd=5, lty=1)
legend("topright", 
       legend=c(paste("M1 (n=",survLGG$n[1],")",sep=""),
                            paste("M2 (n=",survLGG$n[2],")",sep=""),
                            paste("M3 (n=",survLGG$n[3],")",sep=""),
                            paste("M4 (n=",survLGG$n[4],")",sep=""),
                            paste("M5 (n=",survLGG$n[5],")",sep="")),
       col=c("black","red","blue","darkgreen","darkorchid"),
       lwd=3,
       title="DNA METHYLATION Clusters",
       box.lwd=3,
       bg="white")
#text(2,0.05,paste("Log-Rank p-value:",round(1 - pchisq(gg$chisq, length(gg$n) - 1),10)), col="black")
dev.off()

### 


## plot age at diagnosis
library(ggplot2)
ages <- clinical.cluster.genetics[,c("MethylationCluster", "age_at_initial_pathologic")]
ages <- na.omit(ages)
p <- ggplot(ages, aes(factor(MethylationCluster), age_at_initial_pathologic))
p + geom_boxplot(aes(fill = factor(MethylationCluster)), notchwidth=0.25) + 
  geom_jitter(height = 0, position = position_jitter(width = .1), size=2)  + 
  scale_fill_manual(values = rev(c("darkorchid","darkgreen","blue","red","DARKGREY"))) +
  ylab("Age At Diagnosis") +
  xlab("DNA Methylation Clusters") + 
  labs(title = "DNA Methylation Clusters by Age at Diagnosis", fill = "DNA Meth. Cluster")
ggsave(file="age_diagnosis_K5purity_boxplot.png")

## plot purity
library(ggplot2)
purity <- clinical.cluster.genetics[,c("MethylationCluster", "absolute_extract_purity")]
purity <- na.omit(purity)
p <- ggplot(purity, aes(factor(MethylationCluster), absolute_extract_purity))
p + geom_boxplot(aes(fill = factor(MethylationCluster)), notchwidth=0.25) + 
  geom_jitter(height = 0, position = position_jitter(width = .1), size=2)  + 
  scale_fill_manual(values = rev(c("darkorchid","darkgreen","blue","red","DARKGREY"))) +
  ylab("Purity Estimates") +
  xlab("DNA Methylation Clusters") + 
  labs(title = "DNA Methylation Clusters by Purity Estimates", fill = "DNA Meth. Cluster")
ggsave(file="purityestimates_K5purity_boxplot.png")


## plot ploidy
library(ggplot2)
ploidy <- clinical.cluster.genetics[,c("MethylationCluster", "absolute_extract_ploidy")]
ploidy <- na.omit(ploidy)
p <- ggplot(ploidy, aes(factor(MethylationCluster), absolute_extract_ploidy))
p + geom_boxplot(aes(fill = factor(MethylationCluster)), notchwidth=0.25) + 
  geom_jitter(height = 0, position = position_jitter(width = .1), size=2)  + 
  scale_fill_manual(values = rev(c("darkorchid","darkgreen","blue","red","DARKGREY"))) +
  ylab("Ploidy Values") +
  xlab("DNA Methylation Clusters") + 
  labs(title = "DNA Methylation Clusters by Ploidy", fill = "DNA Meth. Cluster")
ggsave(file="ploidy_K5purity_boxplot.png")
###normals from GEO GSE15745 (Methylation)

## get 11,977 CpG probe Id for figure.
load(file="LGG.puritycorrection.HeatmapOrder.rows.ID.Rda")
dimnames(gse15745)[[1]] <- gse15745[,1]
gse15745 <- gse15745[,-1]
## data
gse15745 <- read.delim(file="GSE15745.brainMethylation_methylation.txt", sep="\t")
## sample manifest
gse15745.m <- read.delim(file="GSE15745-GPL8490_series_matrix_samples.txt", sep="\t")
gse15745.m <- t(gse15745.m)
dimnames(gse15745.m)[[2]] <- gse15745.m[1,]
gse15745.m <- gse15745.m[-1,]
gse15745.m <- as.data.frame(gse15745.m)
gse15745.m$SampleID <- dimnames(gse15745.m)[[1]]
gse15745.m <- gse15745.m[,c(6:7)]
dimnames(gse15745.m)[[2]] <- c("source_name_ch1", "SampleID")
head(gse15745.m)


##average by brain region
gse15745.crblm <- apply(gse15745[,as.character(subset(gse15745.m, source_name_ch1=="Human Brain Tissue: CRBLM")[,2])], 1, mean, na.rm=T)
gse15745.pons <- apply(gse15745[,as.character(subset(gse15745.m, source_name_ch1=="Human Brain Tissue: PONS")[,2])], 1, mean, na.rm=T)
gse15745.fctx <- apply(gse15745[,as.character(subset(gse15745.m, source_name_ch1=="Human Brain Tissue: FCTX")[,2])], 1, mean, na.rm=T)
gse15745.tctx <- apply(gse15745[,as.character(subset(gse15745.m, source_name_ch1=="Human Brain Tissue: TCTX")[,2])], 1, mean, na.rm=T)
brain.regions <- cbind(gse15745.crblm, gse15745.fctx,gse15745.pons,gse15745.tctx)
dimnames(brain.regions)[[1]] <- dimnames(gse15745)[[1]]

png(filename="brainregions.png", bg="white", res=300, width=2000, height=2000)
heatmapairs(brain.regions)
dev.off()


#######################
### GBM vs LGG Analysis
#######################
load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/TCGA.AWG.October.2013.Rda")
load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/normals.gse41826.rda")
load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/20130413_HM450.oecg.3kb.rda")
avg.normals <- as.data.frame(avg.normals)
dimnames(lgg.gbm)[[1]] <- lgg.gbm[,1]
lgg.gbm.sub <- lgg.gbm[rownames(dat.s.rm),]
lgg.gbm.sub.x <- subset(lgg.gbm.sub, Chromosome!="X")
lgg.gbm.sub.xy <- subset(lgg.gbm.sub.x, Chromosome!="Y")
temps <- na.omit(lgg.gbm.sub.xy[,5:403])
lgg.gbm.sub.xy <- lgg.gbm.sub.xy[rownames(temps),]
save(avg.normals, lgg.gbm, clinical.cluster.genetics, cc.DNAmethylation, FDb.HM450.oecg.3kb, lgg.gbm.sub.xy, file="/dados/LGG.GBM/lgg.gbm.normals.plusAnno.rda")

##start here if interested in 450K only - GBM.LGG merge###
load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG.GBM.losangeles/lgg.gbm.normals.plusAnno.rda")
load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Brennan2013_GBM27K-450K//gbm.27.450.Brennan2013.rda")
GBM.LGG.27.450k <- merge(dat, GBM27.450.brennan2013, by.x = "Composite.Element.REF", by.y = "TargetID", all.x = F)
avg.lgg.gbm <- apply(lgg.gbm.sub.xy[,5:403], 1, mean, na.rm=TRUE)
avg.lgg.gbm.norm <- cbind(avg.normals[rownames(lgg.gbm.sub.xy),], avg.lgg.gbm)
colnames(avg.lgg.gbm.norm) <- c("avg.normals", "avg.lgg.gbm")
avg.lgg.gbm.norm <- as.data.frame(avg.lgg.gbm.norm)
avg.lgg.gbm.norm.s <- (subset(avg.lgg.gbm.norm, avg.normals <0.2))
## Dichotomize the data ##TODO: redo the consensus clustering without
## dichotomizing the data
avg.lgg.gbm.s <- lgg.gbm.sub.xy[rownames(avg.lgg.gbm.norm.s),5:403]
datT.lgg.gbm <- avg.lgg.gbm.s>0.3
storage.mode(datT.lgg.gbm) <- "numeric"
## Use probes methylated more than 10% of the tumor samples
datT.lgg.gbm <- datT.lgg.gbm[(rowSums(datT.lgg.gbm)/ncol(datT.lgg.gbm))*100 >10,]
load(file="/dados/20130413_HM450.oecg.3kb.rda")
head(FDb.HM450.oecg.3kb)
FDb.HM450.oecg.3kb.s <- subset(FDb.HM450.oecg.3kb, oecg.3kb>1.5)
datT.lgg.gbm.oecg <- datT.lgg.gbm[rownames(datT.lgg.gbm) %in% (FDb.HM450.oecg.3kb.s$probeid),]
## consensus cluster
library(ConsensusClusterPlus)
title="/dados/LGG.GBM/lgg.gbm.c.clu.out"
cc.out.purity.lgg.gbm.1 <-  ConsensusClusterPlus(d=datT.lgg.gbm.oecg, 
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
ccICL = calcICL(cc.out.purity.lgg.gbm.1,plot='pdf')
cc.out.purity.lgg.gbm.2 <- cbind(as.matrix(cc.out.purity.lgg.gbm.1[[2]][[3]]), 
                       as.matrix(cc.out.purity.lgg.gbm.1[[3]][[3]]), 
                       as.matrix(cc.out.purity.lgg.gbm.1[[4]][[3]]), 
                       as.matrix(cc.out.purity.lgg.gbm.1[[5]][[3]]), 
                       as.matrix(cc.out.purity.lgg.gbm.1[[6]][[3]]),
                       as.matrix(cc.out.purity.lgg.gbm.1[[7]][[3]]),
                       as.matrix(cc.out.purity.lgg.gbm.1[[8]][[3]]),
                       as.matrix(cc.out.purity.lgg.gbm.1[[9]][[3]]),
                       as.matrix(cc.out.purity.lgg.gbm.1[[10]][[3]]))

cc.out.purity.lgg.gbm.2 <- as.data.frame(cc.out.purity.lgg.gbm.2)
names(cc.out.purity.lgg.gbm.2) <- c("K2",  "K3",  "K4" , "K5" , "K6" , "K7",  "K8" , "K9" , "K10")

cc.out.purity.lgg.gbm.2$TCGAID <- dimnames(cc.out.purity.lgg.gbm.2)[[1]]
cc.out.purity.lgg.gbm.2$TCGAID <- substr(cc.out.purity.lgg.gbm.2$TCGAID, 1, 12)
cc.DNAmethylation$Tumor <- substr(cc.DNAmethylation$TCGAID, 1, 12)
cc.DNAmethylation$type <- "GBM"
clinical.cluster.genetics$type <- "LGG"
lgg.gbm.anno.purity <- merge(cc.out.purity.lgg.gbm.2, clinical.cluster.genetics, by.x = "TCGAID", by.y = "Tumor", all.x = T)
lgg.gbm.anno.purity <- merge(lgg.gbm.anno.purity, cc.DNAmethylation, by.x = "TCGAID", by.y = "Tumor", all.x = T)
write.table(lgg.gbm.anno.purity, file="/dados/LGG.GBM/lgg.gbm.anno.purity.txt", sep="\t", quote=F, row.names=F) ### used excel
##using excel, merged columns (age_at_initial_pathologic, IDH1, type.x)
lgg.gbm.anno.purity <- read.delim(file="/dados/LGG.GBM/lgg.gbm.anno.purity-modifiedinExcel.txt", sep="\t", header=T)
lgg.gbm.anno.purity.k7 <- read.delim(file="/dados/LGG.GBM/lgg.gbm.anno.purity-modifiedinExcel.K7.txt", sep="\t", header=T)
##added by excel (vlookup) this file: https://www.synapse.org/#!Synapse:syn2369027
##updated gender, daystolastordeath, death01
lgg.gbm.anno.purity.k7.exp <- read.delim(file="/dados/LGG.GBM/lgg.gbm.anno.purity-modifiedinExcel.K7-addedExpressionCalls.txt", sep="\t", header=T)
temp <- avg.lgg.gbm.s[rownames(datT.lgg.gbm.oecg),cc.out.purity.lgg.gbm.1[[7]][[2]][[3]]]
names(temp) <- substr(names(temp), 1, 12)
norms <- exprs(controls.rm)
norms <- norms[rownames(datT.lgg.gbm.oecg),]
temp <- cbind(norms,temp)




#LGG.GBM purity clusters
  color.map.lgg.gbm<-function(mol.biol) {
    if (mol.biol==1) "blue" #ok
    else if (mol.biol==2) "red" #ok
    else if (mol.biol==3) "darkorchid" #ok
    else if (mol.biol==4) "darkgreen" # decided these were not real clusters
    else if (mol.biol==5) "#FF8C00" #ok
    else if (mol.biol==6) "yellow" # decided these were not real clusters
    else if (mol.biol==7) "#DC143C"
    else "#FFFFFF" 
  }#
  
  
conclus.6<-apply(as.matrix(names(temp)),1,func.list$vlookup,lgg.gbm.anno.purity ,"K6.x")
conclus.6<-as.character(conclus.6)
conclus.6[is.na(conclus.6)]<-c("0")
cc.conclus.6<-unlist(lapply(conclus.6, color.map.lgg.gbm))

conclus.7<-apply(as.matrix(names(temp)),1,func.list$vlookup,lgg.gbm.anno.purity ,"K7.x")
conclus.7<-as.character(conclus.7)
conclus.7[is.na(conclus.7)]<-c("0")
cc.conclus.7<-unlist(lapply(conclus.7, color.map.lgg.gbm))

##LGG clusters
color.map.meth.cl <- function(mol.biol) {
  if (mol.biol=="M5") "darkorchid" #darkorchid
  else if (mol.biol=="M3") "blue" #blue
  else if (mol.biol=="M4") "darkgreen" #darkgreen
  else if (mol.biol=="M1") "black" #black
  else if (mol.biol=="M2") "red" #red
  else if (mol.biol==6) "yellow"
  else "#FFFFFF"
} #cluster.consensus.clustering
conclus.5<-apply(as.matrix(names(temp)),1,func.list$vlookup,lgg.gbm.anno.purity ,"MethylationCluster")
conclus.5<-as.character(conclus.5)
conclus.5[is.na(conclus.5)]<-c("0")
cc.conclus.5<-unlist(lapply(conclus.5, color.map.meth.cl))

#rnacluster
color.map.lgg.gbmRNAseq<-function(mol.biol) {
  if (mol.biol=="R1") "red" #ok
  else if (mol.biol=="R2") "#0D9344" #ok
  else if (mol.biol=="R3") "#24AAE1" #ok
  else if (mol.biol=="R4") "#223F99" # decided these were not real clusters
  else "#FFFFFF" 
}#
rnacluster<-apply(as.matrix(names(temp)),1,func.list$vlookup,lgg.gbm.anno.purity ,"RNASeqCluster")
rnacluster<-as.character(rnacluster)
rnacluster[is.na(rnacluster)]<-c("0")
cc.rnacluster<-unlist(lapply(rnacluster, color.map.lgg.gbmRNAseq))


#rnaLGG-GBMcluster (pased on https://www.synapse.org/#!Synapse:syn2369027)
rnaclusterLGGGBM<-apply(as.matrix(names(temp)),1,func.list$vlookup,lgg.gbm.anno.purity.k7.exp ,"consensusClustering")
rnaclusterLGGGBM<-as.character(rnaclusterLGGGBM)
rnaclusterLGGGBM[is.na(rnaclusterLGGGBM)]<-c("#FFFFFF")
cc.rnaclusterLGGGBM<-rnaclusterLGGGBM

#CNCcluster


color.map.lgg.gbmCNseq<-function(mol.biol) {
  if (mol.biol=="C1") "red" #ok
  else if (mol.biol=="C2") "#0D9344" #ok
  else if (mol.biol=="C3") "#24AAE1" #ok
  else if (mol.biol=="C4") "#223F99" # decided these were not real clusters
  else "#FFFFFF" 
}#
CNCcluster<-apply(as.matrix(names(temp)),1,func.list$vlookup,lgg.gbm.anno.purity , "CNCluster")
CNCcluster<-as.character(CNCcluster)
CNCcluster[is.na(CNCcluster)]<-c("0")
cc.CNCcluster<-unlist(lapply(CNCcluster, color.map.lgg.gbmCNseq))

#miRNAcluster
color.map.lgg.gbmmiRseq<-function(mol.biol) {
  if (mol.biol=="mi1") "red" #ok
  else if (mol.biol=="mi2") "#0D9344" #ok
  else if (mol.biol=="mi3") "#24AAE1" #ok
  else if (mol.biol=="mi4") "#223F99" # decided these were not real clusters
  else "#FFFFFF" 
}#
miRNAcluster<-apply(as.matrix(names(temp)),1,func.list$vlookup,lgg.gbm.anno.purity , "miRNACluster")
miRNAcluster<-as.character(miRNAcluster)
miRNAcluster[is.na(miRNAcluster)]<-c("0")
cc.miRNAcluster<-unlist(lapply(miRNAcluster, color.map.lgg.gbmmiRseq))

#RPPAcluster
color.map.lgg.gbmRPPA<-function(mol.biol) {
  if (mol.biol=="P1") "red" #ok
  else if (mol.biol=="P2") "#0D9344" #ok
  else if (mol.biol=="P3") "#24AAE1" #ok
  else if (mol.biol=="P4") "#223F99" # decided these were not real clusters
  else "#FFFFFF" 
}#
RPPAcluster<-apply(as.matrix(names(temp)),1,func.list$vlookup,lgg.gbm.anno.purity , "RPPACluster")
RPPAcluster<-as.character(RPPAcluster)
RPPAcluster[is.na(RPPAcluster)]<-c("0")
cc.RPPAcluster<-unlist(lapply(RPPAcluster, color.map.lgg.gbmRPPA))

#Oncosigncluster
color.map.lgg.gbmOSC<-function(mol.biol) {
  if (mol.biol=="OSC1") "red" #ok
  else if (mol.biol=="OSC2") "#0D9344" #ok
  else if (mol.biol=="OSC3") "#24AAE1" #ok
  else if (mol.biol=="Unclassified") "lightgrey" # decided these were not real clusters
  else "#FFFFFF" 
}#

OSCcluster<-apply(as.matrix(names(temp)),1,func.list$vlookup,lgg.gbm.anno.purity , "Oncosign")
OSCcluster<-as.character(OSCcluster)
OSCcluster[is.na(OSCcluster)]<-c("0")
cc.OSCcluster<-unlist(lapply(OSCcluster, color.map.lgg.gbmOSC))

#GBMcluster
color.map.lgg.gbmCluster<-function(mol.biol) {
  if (mol.biol=="G-CIMP") "red" #ok
  else if (mol.biol=="Classical") "#0D9344" #ok
  else if (mol.biol=="Mesenchymal") "#24AAE1" #ok
  else if (mol.biol=="Neural") "lightgrey" # decided these were not real clusters
  else if (mol.biol=="Proneural") "blue" # decided these were not real clusters
  else "#FFFFFF" 
}#
gbmcluster<-apply(as.matrix(names(temp)),1,func.list$vlookup,lgg.gbm.anno.purity , "Cluster")
gbmcluster<-as.character(gbmcluster)
gbmcluster[is.na(gbmcluster)]<-c("0")
cc.gbmcluster<-unlist(lapply(gbmcluster, color.map.lgg.gbmCluster))

#idh1
color.map.mutation<-function(mol.biol) {
  if (mol.biol=="Mutant") "black"
  else if (mol.biol=="wt") "Gray"
  else "White"
}#
idh1<-apply(as.matrix(names(temp)),1,func.list$vlookup,lgg.gbm.anno.purity ,"IDH1")
idh1<-as.character(idh1)
idh1[is.na(idh1)]<-c("0")
cc.idh1<-unlist(lapply(idh1, color.map.mutation))

#type
color.map.type<-function(mol.biol) {
  if (mol.biol=="GBM") "red" #ok
  else if (mol.biol=="(31.2,40]") "orange" #ok
  else if (mol.biol=="LGG") "blue" #ok
  else "#FFFFFF" 
}#
type<-apply(as.matrix(names(temp)),1,func.list$vlookup,lgg.gbm.anno.purity ,"type.x")
type<-as.character(type)
type[is.na(type)]<-c("0")
cc.type<-unlist(lapply(type, color.map.type))

#1p19q
color.map.IDH1p19q<-function(mol.biol) {
  if (mol.biol=="TRUE") "black"
  else if (mol.biol=="FALSE") "Gray"
  else "White"
}#
codel<-apply(as.matrix(names(temp)),1,func.list$vlookup,lgg.gbm.anno.purity ,"X1p_19q_co_del_status")
codel<-as.character(codel)
codel[is.na(codel)]<-c("0")
cc.codel<-unlist(lapply(codel, color.map.IDH1p19q))

#cc.col<-matrix(as.character(c(cc.batch,cc.p53,cc.pten,cc.nf1,cc.idh1,cc.egfr,cc.expclus,cc.gcimp,cc.conclus)), nrow=288,ncol=9)
cc.col<-matrix(as.character(c(
                              cc.type,
                              cc.codel,
                              cc.idh1, 
                              cc.OSCcluster,
                              cc.RPPAcluster,
                             # cc.miRNAcluster7,
                              cc.miRNAcluster,
                              cc.CNCcluster,
                              cc.rnacluster,
                              cc.gbmcluster,
                              cc.conclus.5,
                              cc.rnaclusterLGGGBM,
                              cc.conclus.7)), 
               nrow=476,
               ncol=12)
#colnames(cc.col) = c("BATCH","TP53","PTEN","NF1","IDH1","EGFR","GeneExp CC", "Noushmehr et al - GCIMP", "DNA MethCC")
colnames(cc.col) = c(
                    "tumor type",
                    "1p19q codel",
                     "IDH1",
                     "OSC CC", 
                     "RPPA CC", 
                     "miRNA CC (4)", 
                     "CN CC", 
                     "RNA CC (4/6)", 
                     "GBM cluster",
                     "LGG cluster",
                    "LGG.GBM RNA cluster",
                     "LGG.GBM DNAme cluster")

#library(RColorBrewer)
#hmColors <- colorRampPalette(c("darkblue", "yellow"))(256)
library(heatmap.plus)
library(matlab)
source(file="heatmap.plus.R")


pdf(file="/dados/LGG.GBM/heatmapLGGGBM.purity.ver2-k7.pdf", width=10, height=10)
hv2.lgg.gbm<-heatmap.plus.sm(
  #as.matrix(temp.1),
  as.matrix(temp),
  na.rm=TRUE,
  scale="none",
  #RowSideColor=probe.cc,
  ColSideColors=cc.col,
  #col=c("white", "#A0A0A0", "black"),
  col=jet.colors(75),
  #key=TRUE,
  symkey=FALSE,
  density.info="none",
  trace="none",
  Rowv=TRUE,
  Colv=NA,
  cexRow=1,
  cexCol=0.6,
  keysize=1,
  dendrogram=c("none"),
  #main = paste("LGG DNA meth Clustering, K5 ",dim(temp)[1],"P.; ",dim(as.matrix(temp[,(as.character(cl.5.sort[,1]))]))[2],"S."),
  #labCol=NA
  labRow=NA
)
dev.off()

if (mol.biol==1) "blue" #ok
else if (mol.biol==2) "red" #ok
else if (mol.biol==3) "darkorchid" #ok
else if (mol.biol==4) "darkgreen" # decided these were not real clusters
else if (mol.biol==5) "#FF8C00" #ok
else if (mol.biol==6) "yellow" # decided these were not real clusters
else if (mol.biol==7) "#DC143C"
## plot age at diagnosis
library(ggplot2)
ages <- lgg.gbm.anno.purity.k7.exp[,c("K7.x", "age_at_initial_pathologic","type.x")]
ages <- na.omit(ages)
p <- ggplot(ages, aes(factor(K7.x), age_at_initial_pathologic))
p + geom_boxplot(aes(fill = factor(K7.x)), notchwidth=0.25) + 
  geom_jitter(aes(colour = type.x),height = 0, position = position_jitter(width = .1), size=2)  + 
  scale_fill_manual(values = rev(c("darkgreen",
                                   "#DC143C",
                                   "darkorchid",
                                   "#FF8C00",
                                   "yellow",
                                   "red",
                                   "blue"))) +
  scale_colour_manual(values = c("red", "blue")) +
  ylab("Age At Diagnosis") +
  xlab("DNA Methylation Clusters") + 
  labs(title = "DNA Methylation Clusters by Age at Diagnosis", fill = "DNA Meth. Cluster") +
  facet_grid(type.x ~ .)
ggsave(file="/dados/LGG.GBM/age_diagnosis_K7purity_boxplot_byTYPE.png", width=10, height=10)


library(survival);
## days to last or death
## Survival Curves
## days to last or death
#adjust
agecat <- cut(lgg.gbm.anno.purity.k7.exp$age_at_initial_pathologic, c(13,31.25,40,54,75))
fdata <- lgg.gbm.anno.purity.k7.exp
fdata$sex <- fdata$gender
fdata$age <- fdata$age_at_initial_pathologic
fdata$age2 <- agecat
fdata$group <- as.factor(fdata$K7.x)
data(uspop2)
refpop <- uspop2[as.character(14:75), , "2000"]
temp <- as.numeric(cut(14:75, c(13,31.25,40,54,75) + 0.5))
pi.us<- tapply(refpop, list(temp[row(refpop)], col(refpop)), sum)/sum(refpop)
tab2 <- with(fdata, table(age2, sex, group))/ nrow(fdata)
us.wt <- rep(pi.us, 7)/ tab2
index <- with(fdata, cbind(as.numeric(age2), as.numeric(sex), as.numeric(group)))
fdata$uswt <- us.wt[index]
sfit1a.d <- survfit(Surv(daystolastordeath, death01) ~ group, fdata)
sfit3a.d <-survfit(Surv(daystolastordeath, death01) ~ group, data=fdata, weight=uswt)

#sfit1a.p <- survfit(Surv(ProgFreeSurvTime_days, ProgFreeSurvEvent01) ~ group, fdata)
#sfit3a.p <-survfit(Surv(ProgFreeSurvTime_days, ProgFreeSurvEvent01) ~ group, data=fdata, weight=uswt)

#plot(sfit3a, mark.time=F, col=c("red", "orange", "yellow", "green", "blue"), lty=1, lwd=2, xscale=365.25, xlab="Years from Sample", ylab="Survival")
#
#x11()
survLGGGBM <- sfit3a.d

png(filename = "/dados/LGG.GBM/K7puritycorrection-survival_plots-OverallSurvival-adjustedbyAge.png", bg="white", res=300, width=2000, height=2000)
#par(mfrow=c(1,2))
plot(survLGGGBM, 
     lwd=5,
     col=rev(c("darkgreen",
     "#DC143C",
     "darkorchid",
     "#FF8C00",
     "yellow",
     "red",
     "blue")),
     #col=c("black","red","blue","darkgreen","darkorchid"),
     main="Overall Survival\nKaplan-Meier Survival Curves\nTCGA LGG.GBM SAMPLES (DNA Methylation Clusters)", 
     #main="Progression Free\nKaplan-Meier Survival Curves\nTCGA LGG SAMPLES (DNA Methylation Clusters)", 
     xlab="TIME SINCE DIAGNOSIS (YEARS)", 
     ylab="PERCENT PROBABILITY OF SURVIVAL", 
     yscale=100,
     xscale=365.25, 
     bg="black"
)
#lines(sfit1a.d, mark.time=T, col=c("black","red","blue","darkgreen","darkorchid"), lty=2, lwd=1, xscale=365.25)
box(col="black", lwd=3);
#abline(v=5*52, lwd=5, lty=2); text(4.9*52, 0.75, "5 YEARS", srt=90)
abline(h=0, lwd=5, lty=1)
legend("topright", 
       legend=c(paste("M1 (n=",survLGGGBM$n[1],")",sep=""),
                paste("M2 (n=",survLGGGBM$n[2],")",sep=""),
                paste("M3 (n=",survLGGGBM$n[3],")",sep=""),
                paste("M4 (n=",survLGGGBM$n[4],")",sep=""),
                paste("M5 (n=",survLGGGBM$n[5],")",sep=""),
                paste("M6 (n=",survLGGGBM$n[6],")",sep=""),
                paste("M7 (n=",survLGGGBM$n[7],")",sep="")
                ),
       col=rev(c("darkgreen",
             "#DC143C",
             "darkorchid",
             "#FF8C00",
             "yellow",
             "red",
             "blue")),
       lwd=3,
       title="DNA METHYLATION Clusters",
       box.lwd=3,
       bg="white")
#text(2,0.05,paste("Log-Rank p-value:",round(1 - pchisq(gg$chisq, length(gg$n) - 1),10)), col="black")
dev.off()

### 
##start here to work with GBM + LGG merge using 27k.450k merged GBM platforms ###
##start here###
load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG.GBM.losangeles/lgg.gbm.normals.plusAnno.rda")
load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Brennan2013_GBM27K-450K//gbm.27.450.Brennan2013.rda") ##SNPs have already been removed
load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG.GBM.losangeles/lgg_289_samples.Rda")
GBM.LGG.27.450k <- merge(dat, GBM27.450.brennan2013, by.x = "Composite.Element.REF", by.y = "TargetID", all.x = F)
GBM.LGG.27.450k <- GBM.LGG.27.450k[,-c(294:296)]
dimnames(GBM.LGG.27.450k)[[1]] <- GBM.LGG.27.450k$Composite.Element.REF
#GBM.LGG.27.450k.noX <- subset(GBM.LGG.27.450k, Chromosome.x != "X")
#GBM.LGG.27.450k.noXY <- subset(GBM.LGG.27.450k.noX, Chromosome.x != "Y")
GBM.LGG.27.450k.noXY <- GBM.LGG.27.450k[GBM.LGG.27.450k$Chromosome.x != "X" & GBM.LGG.27.450k$Chromosome.x != "Y",]

## dichotomizing the data
GBM.LGG.27.450k.noXY.dic <- GBM.LGG.27.450k.noXY[ ,5:689]
GBM.LGG.27.450k.noXY.dic <- GBM.LGG.27.450k.noXY.dic>0.3
storage.mode(GBM.LGG.27.450k.noXY.dic) <- "numeric"
## Use probes methylated more than 10% of the tumor samples
GBM.LGG.27.450k.noXY.dic <- GBM.LGG.27.450k.noXY.dic[(rowSums(GBM.LGG.27.450k.noXY[,5:689])/ncol(GBM.LGG.27.450k.noXY[,5:689]))*100 >10,]

head(FDb.HM450.oecg.3kb)
FDb.HM450.oecg.3kb.s <- subset(FDb.HM450.oecg.3kb, oecg.3kb>1.25)
GBM.LGG.27.450k.noXY.dic.oecg <- GBM.LGG.27.450k.noXY.dic[rownames(GBM.LGG.27.450k.noXY.dic) %in% (FDb.HM450.oecg.3kb.s$probeid),]

#TODO:
## need to select tissue specific probes (>0.3)
avg.lgg.gbm.27.450.noXY <- apply(GBM.LGG.27.450k.noXY[,5:689], 1, mean, na.rm=TRUE)
avg.norm.lgg.gbm.27.450.noXY <- merge(as.data.frame(avg.lgg.gbm.27.450.noXY), avg.normals, by = 0, all.x = T)
library(LSD)
comparisonplot(avg.norm.lgg.gbm.27.450.noXY$avg.normals, avg.norm.lgg.gbm.27.450.noXY$avg.lgg.gbm.27.450.noXY, pimp=TRUE)
lgg.gbm.27.450.tumor.specific <- (subset(avg.norm.lgg.gbm.27.450.noXY, avg.normals < 0.3))

##selecting tumor specific probes for lgg.gbm
dat.lgg.gbm.27.450.noXY.dic.oecg <- GBM.LGG.27.450k.noXY.dic.oecg[(rownames(GBM.LGG.27.450k.noXY.dic.oecg) %in% lgg.gbm.27.450.tumor.specific$Row.names), ]

## consensus cluster
library(ConsensusClusterPlus)
title="/dados/ResearchProjects/thais/TCGA/LGG.GBM/ConsensusCluster.LGG.GBM.27.450"
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

cc.out.purity.lgg.gbm.27.450.2$TumorType <- "UNK"
cc.out.purity.lgg.gbm.27.450.2[1:289,"TumorType"] <- "LGG"
cc.out.purity.lgg.gbm.27.450.2[290:685,"TumorType"] <- "GBM"

brennan.class <- read.delim(file="../Brennan2013_GBM27K-450K/DNA.methylation.k6.txt", sep="\t", header=T)
cc.out.purity.lgg.gbm.27.450.2$TCGAIDnew <- substr(rownames(cc.out.purity.lgg.gbm.27.450.2),1,12)
brennan.class$TCGAIDnew <- substr(brennan.class$TCGAID, 1, 12)
brennan.class.s <- brennan.class[brennan.class$DNA.Methylation.Clusters != "UNKNOWN", ]

cc.out.purity.lgg.gbm.27.450.2$TCGAID <- rownames(cc.out.purity.lgg.gbm.27.450.2)
lgg.gbm.anno.27.450 <- merge(cc.out.purity.lgg.gbm.27.450.2, brennan.class.s[,c("DNA.Methylation.Clusters","TCGAIDnew")], by = "TCGAIDnew", all.x =T)
rownames(lgg.gbm.anno.27.450) <- lgg.gbm.anno.27.450$TCGAID

##load in LGG clinical annotations
### March 14, 2014 from synapse:
syn2369104 <- synGet(id='syn2369104')
synapseLogin()
syn2369104 <- synGet(id='syn2369104', load=T) #updated march 14
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
oncosign <- read.delim(file="/dados/oncosign_assignments.txt", sep="\t", header=T) ## new oncosign assignments provided by giovanni (April 5, 2014)
oncosign <- read.delim(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/oncosign_assignments.txt", sep="\t", header=T) ## new oncosign assignments provided by giovanni (April 5, 2014)
clinical.cluster.genetics <- merge(clinical.cluster.genetics, oncosign, by.x = "Tumor", by.y = "sample", all.x =T)

xlgg.gbm.anno.27.450 <- merge(lgg.gbm.anno.27.450, clinical.cluster.genetics, by.x = "TCGAIDnew", by.y = "Tumor", all = T)
#We should have done -> all.x = T
# there are 4 extra lgg patients 
# setdiff(xlgg.gbm.anno.27.450$TCGAIDnew,lgg.gbm.anno.27.450$TCGAIDnew)

##add in age
agecat <- cut(xlgg.gbm.anno.27.450$age_at_initial_pathologic, c(13,31.25,40,54,75))
xlgg.gbm.anno.27.450$age.strat <- as.character(agecat)

###TODO:
Need to decide clusters. K5 or K6.  
To determine clusters we need additional clinical data from GBM.  Pull in the clinical data, massage it into the current object (xlgg.gbm.anno.27.450)
library(downloader)
download("http://tcga-data.nci.nih.gov/docs/publications/gbm_2013/supplement/Molecular_subtype_classification.xlsx", "/dados/ResearchProjects/thais/TCGA/LGG.GBM/Molecular_subtype_classification.xlsx")
#I removed manually the first 2 rows and deleted the row 316 because it was duplicated
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
xx <- merge(xlgg.gbm.anno.27.450,gbm.clinical.data.nR,by.x = "TCGAIDnew", by.y = "TCGAID", all.x = TRUE)
write.table(xx,file="lgg_gbm_27_450.txt",quote=F,row.names=F,sep="\t")
###opened it in excel, and manually massaged the data
xxMod <- read.delim(file="lgg_gbm_27_450_xlsMod.csv", sep="\t", header=T)
###downloaded the RNAcluster from Synapse
#https://www.synapse.org/#!Synapse:syn2407035
rnacluster <- read.delim(file="RNAseqClusters (syn2407035)//pooledClassification-11-4-2014.csv", sep=",", header=T)
names(rnacluster) <- c("RNA.TCGAID", "RNA.cluster", "RNA.BrennanCluster", "RNA.platform", "RNA.tumor")
xxMODr <- merge(xxMod, rnacluster, by.x = "TCGAIDnew", by.y = "RNA.TCGAID", all = T)


##add in age
agecat.x <- cut(xxMod$age_at_initial_pathologic.x, c(13,31.25,40,54,75))
xxMod$age.strat.x <- as.character(agecat.x)

## plot age at diagnosis
library
x <- xxMod
x$K4.x <- factor(x$K4.x,c("1","4","3","2")) #change the order
ages <- x[,c("K4.x", "age_at_initial_pathologic.x")]
ages <- na.omit(ages)
p <- ggplot(ages, aes(factor(K4.x), age_at_initial_pathologic.x))
p + geom_boxplot(aes(fill = factor(K4.x)), notchwidth=0.25) + 
  geom_jitter(height = 0, position = position_jitter(width = .1), size=2)  + 
  scale_fill_manual(values=c("red","blue","green", "yellow"#, "green"
                             ), labels = c("M1", "M2", "M3","M4"#, "M5"#, "M6"
                                           )) +
  #scale_fill_manual(values = rev(c("darkorchid","darkgreen","blue","red","DARKGREY"))) +
  ylab("Age At Diagnosis") +
  xlab("DNA Methylation Clusters") + 
  labs(title = "DNA Methylation Clusters by Age at Diagnosis", fill = "DNA Meth. Cluster")
ggsave(file="age_diagnosis_K6_LGG_GBM_boxplot.png")




library(survival)
## days to last or death
## Survival Curves
## days to last or death
#adjust
#agecat <- cut(xxMod$age_at_initial_pathologic.x, c(13,31.25,40,54,75))
agecat <- xxMod$age.strat.x
fdata <- xxMod
fdata$sex <- fdata$gender
fdata$age <- fdata$age_at_initial_pathologic.x
fdata$age2 <- agecat
fdata$group <- as.factor(fdata$K4.x)
data(uspop2)
refpop <- uspop2[as.character(14:75), , "2000"]
temp <- as.numeric(cut(14:75, c(13,31.2,40,54,75) + 0.5))
pi.us<- tapply(refpop, list(temp[row(refpop)], col(refpop)), sum)/sum(refpop)
tab2 <- with(fdata, table(age2, sex, group))/ nrow(fdata)
us.wt <- rep(pi.us, 4)/ tab2
index <- with(fdata, cbind(as.numeric(age2), as.numeric(sex), as.numeric(group)))
fdata$uswt <- us.wt[index]
sfit1a.d <- survfit(Surv(daystolastordeath.x, death01.x) ~ group, fdata)
#sfit3a.d <- survfit(Surv(daystolastordeath.x, death01.x) ~ group, data=fdata, weight=as.character(fdata$uswt))

#sfit1a.p <- survfit(Surv(ProgFreeSurvTime_days, ProgFreeSurvEvent01) ~ group, fdata)
#sfit3a.p <-survfit(Surv(ProgFreeSurvTime_days, ProgFreeSurvEvent01) ~ group, data=fdata, weight=uswt)

#plot(sfit3a, mark.time=F, col=c("red", "orange", "yellow", "green", "blue"), lty=1, lwd=2, xscale=365.25, xlab="Years from Sample", ylab="Survival")
#
#x11()
survLGGGBM.x <- sfit1a.d

png(filename = "/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_GBM_K6.png", bg="white", res=300, width=2000, height=2000)
#par(mfrow=c(1,2))
plot(survLGGGBM.x, 
     lwd=4,
     col=c("red",
               "yellow",
               "green",
               #"#FF8C00",
              # "yellow",
               #"red",
               "blue"
                #"orange",
                #"purple"
           ),
     #col=c("black","red","blue","darkgreen","darkorchid"),
     main="Overall Survival\nKaplan-Meier Survival Curves\nTCGA LGG.GBM SAMPLES (DNA Methylation Clusters)", 
     #main="Progression Free\nKaplan-Meier Survival Curves\nTCGA LGG SAMPLES (DNA Methylation Clusters)", 
     xlab="TIME SINCE DIAGNOSIS (YEARS)", 
     ylab="PERCENT PROBABILITY OF SURVIVAL", 
     yscale=100,
     xscale=365.25, 
     bg="black"
)
#lines(sfit1a.d, mark.time=T, col=c("black","red","blue","darkgreen","darkorchid"), lty=2, lwd=1, xscale=365.25)
box(col="black", lwd=3);
#abline(v=5*52, lwd=5, lty=2); text(4.9*52, 0.75, "5 YEARS", srt=90)
abline(h=0, lwd=4, lty=1)
legend("topright", 
       legend=c(paste("M1 (n=",survLGGGBM.x$n[1],")",sep=""),
                paste("M2 (n=",survLGGGBM.x$n[4],")",sep=""),
                paste("M3 (n=",survLGGGBM.x$n[3],")",sep=""),
                paste("M4 (n=",survLGGGBM.x$n[2],")",sep="")
               # paste("M5 (n=",survLGGGBM.x$n[3],")",sep="")
              #  paste("M6 (n=",survLGGGBM.x$n[3],")",sep="")
              #  paste("M7 (n=",survLGGGBM$n[7],")",sep="")
       ),
       col=c("red",
                  "blue",
                  #"orange",
                  #"#FF8C00",
                  # "yellow",
                  #"red",
                  #"purple",
                  "green",
                  "yellow"),
       lwd=3,
       title="DNA METHYLATION Clusters",
       box.lwd=3,
       bg="white")
#text(2,0.05,paste("Log-Rank p-value:",round(1 - pchisq(gg$chisq, length(gg$n) - 1),10)), col="black")
dev.off()




##HEATMAPS
library(heatmap.plus)
library(matlab)

#Function to label the samples (heatmap)
colors <- function(cluster,gbm.cluster,lgg.cluster,tumor,age,IDH,rna.dna.c,rna.cluster,coc.c){
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
                else 
                  cluster[i] <- "purple"
    }
  }
  
  for(i in 1:length(gbm.cluster)){
    if(is.na(gbm.cluster[i]))
     aux[i] <- "white"
    else{
      if(gbm.cluster[i] == "G-CIMP")
        aux[i] <- "yellow"
      else if(gbm.cluster[i] == "M1")
        aux[i] <- "green"
        else if(gbm.cluster[i] == "M2")
          aux[i] <- "black"
          else if(gbm.cluster[i] == "M3")
            aux[i] <- "blue"
            else if(gbm.cluster[i] == "M4")
              aux[i] <- "red"
              else if(gbm.cluster[i] == "M6")
                aux[i] <- "darksalmon"
                else 
                  aux[i] <- "white"
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
  for(i in 1:length(tumor)){
    if(is.na(tumor[i]))
      tumor[i] <- "white"
    else{ 
    if(tumor[i] == "LGG")
      tumor[i] <- "black"
    else
      tumor[i] <- "grey"
    }
    
  }
  
  for(i in 1:length(IDH)){
    if(is.na(IDH[i]))
      IDH[i] <- "white"
    else{
      if(IDH[i] == "Mutant")
        IDH[i] <- "black"
      else if(IDH[i] == "wt")
          IDH[i] <- "grey"
        else
          IDH[i] <- "white"
    
    }
  }

  for(i in 1:length(age)){
    if(is.na(age[i]))
      age[i] <- "white"
    else{
      if(age[i] == "(13,31.2]")
        age[i] <- "red"
      else if(age[i] == "(31.2,40]")
        age[i] <- "yellow"
      else if(age[i] == "(40,54]")
        age[i] <- "green"
      else if(age[i] == "(54,75]")
        age[i] <- "blue"
      else 
        age[i] <- "white"
    }
  }
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
  rna.cluster[is.na(rna.cluster)] <- "white"
  
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
  
  cluster <- cbind(age,IDH,rna.dna.c,aux,coc.c,lgg.cluster,rna.cluster,cluster,tumor)
  colnames(cluster) <- c("Age", "IDH1","GBM Clusters [Brennan, 2013]","GBM DNA Meth Cl [Brennan, 2013]","COC Clusters [NEJM, unpublished]", "LGG Clusters [NEJM, unpublished]", "RNA Clusters", "Consensus Cluster", "Tumor Type")
  return(cluster)
}

#Add IDH1, AGE-CUT, LGG CLUSTERS information
#bb <- lgg.gbm.anno.27.450 #dim 685 13
lgg.gbm.anno.27.450 <- merge(lgg.gbm.anno.27.450,xxMod[,c("TCGAIDnew","MethylationCluster","IDH1.x","age.strat.x")],by.x="TCGAIDnew",by.y="TCGAIDnew",all.x=T)
rownames(lgg.gbm.anno.27.450) <- lgg.gbm.anno.27.450$TCGAID
lgg.gbm.anno.27.450.rna <- merge(lgg.gbm.anno.27.450,xxMODr[,c("TCGAIDnew","RNA.cluster","Cluster","COCCluster")],by.x="TCGAIDnew",by.y="TCGAIDnew",all.x=T)
rownames(lgg.gbm.anno.27.450.rna) <- lgg.gbm.anno.27.450.rna$TCGAID
       
#Heatmap for 4 Clusters
GBM.LGG.27.450k.noXY.s <- GBM.LGG.27.450k.noXY[rownames(dat.lgg.gbm.27.450.noXY.dic.oecg) ,5:689]
lgg.gbm.27.450.order <- GBM.LGG.27.450k.noXY.s[,cc.out.purity.lgg.gbm.27.450[[4]][[2]][[3]]]
lgg.gbm.anno.27.450.rna <- lgg.gbm.anno.27.450.rna[colnames(lgg.gbm.27.450.order),]
clab.4 <- colors(as.matrix(lgg.gbm.anno.27.450.rna$K4),as.matrix(lgg.gbm.anno.27.450.rna$DNA.Methylation.Clusters),as.matrix(lgg.gbm.anno.27.450.rna$MethylationCluster),as.matrix(lgg.gbm.anno.27.450.rna$TumorType),as.matrix(lgg.gbm.anno.27.450.rna$age.strat.x),as.matrix(lgg.gbm.anno.27.450.rna$IDH1.x),as.matrix(lgg.gbm.anno.27.450.rna$Cluster),as.matrix(lgg.gbm.anno.27.450.rna$RNA.cluster),as.matrix(lgg.gbm.anno.27.450.rna$COCCluster))
lgg.gbm.27.450.order <- lgg.gbm.27.450.order[heatmap.lgg.gbm$rowInd,]
normals.sel <- normals.sel[rownames(lgg.gbm.27.450.order),]
lgg.gbm.27.450.order <- cbind(normals.sel,lgg.gbm.27.450.order)
norm.label <- matrix("white", nrow = 80, ncol = 9)
clab.4 <- rbind(norm.label,clab.4)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps/heatmap_K4_normals.png",res=200,width=2000,height=2000)
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

lgg.gbm.27.450.order <- GBM.LGG.27.450k.noXY.s[,cc.out.purity.lgg.gbm.27.450[[5]][[2]][[3]]]
lgg.gbm.anno.27.450.rna <- lgg.gbm.anno.27.450.rna[colnames(lgg.gbm.27.450.order),]
clab.5 <- colors(as.matrix(lgg.gbm.anno.27.450.rna$K5),as.matrix(lgg.gbm.anno.27.450.rna$DNA.Methylation.Clusters),as.matrix(lgg.gbm.anno.27.450.rna$MethylationCluster),as.matrix(lgg.gbm.anno.27.450.rna$TumorType),as.matrix(lgg.gbm.anno.27.450.rna$age.strat.x),as.matrix(lgg.gbm.anno.27.450.rna$IDH1.x),as.matrix(lgg.gbm.anno.27.450.rna$Cluster),as.matrix(lgg.gbm.anno.27.450.rna$RNA.cluster),as.matrix(lgg.gbm.anno.27.450.rna$COCCluster))
clab.5 <- clab.5[c(1:215,554:599,216:364,365:553,600:685),]
lgg.gbm.27.450.order <- lgg.gbm.27.450.order[,c(1:215,554:599,216:364,365:553,600:685)]
lgg.gbm.27.450.order <- lgg.gbm.27.450.order[heatmap.lgg.gbm$rowInd,]
normals.sel <- normals.sel[rownames(lgg.gbm.27.450.order),]
lgg.gbm.27.450.order <- cbind(normals.sel,lgg.gbm.27.450.order)
norm.label <- matrix("white", nrow = 80, ncol = 9)
clab.5 <- rbind(norm.label,clab.5)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps/heatmap_K5_normals.png",res=200,width=2000,height=2000)
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
lgg.gbm.27.450.order <- GBM.LGG.27.450k.noXY.s[,cc.out.purity.lgg.gbm.27.450[[6]][[2]][[3]]]
lgg.gbm.anno.27.450.rna <- lgg.gbm.anno.27.450.rna[colnames(lgg.gbm.27.450.order),]
clab.6 <- colors(as.matrix(lgg.gbm.anno.27.450.rna$K6),as.matrix(lgg.gbm.anno.27.450.rna$DNA.Methylation.Clusters),as.matrix(lgg.gbm.anno.27.450.rna$MethylationCluster),as.matrix(lgg.gbm.anno.27.450.rna$TumorType),as.matrix(lgg.gbm.anno.27.450.rna$age.strat.x),as.matrix(lgg.gbm.anno.27.450.rna$IDH1.x),as.matrix(lgg.gbm.anno.27.450.rna$Cluster),as.matrix(lgg.gbm.anno.27.450.rna$RNA.cluster),as.matrix(lgg.gbm.anno.27.450.rna$COCCluster))
clab.6 <- clab.6[c(1:215,580:625,216:295,296:364,365:579,626:685),]
lgg.gbm.27.450.order <- lgg.gbm.27.450.order[,c(1:215,580:625,216:295,296:364,365:579,626:685)]
lgg.gbm.27.450.order <- lgg.gbm.27.450.order[heatmap.lgg.gbm$rowInd,]
normals.sel <- normals.sel[rownames(lgg.gbm.27.450.order),]
lgg.gbm.27.450.order <- cbind(normals.sel,lgg.gbm.27.450.order)
norm.label <- matrix("white", nrow = 80, ncol = 9)
clab.6 <- rbind(norm.label,clab.6)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps/heatmap_K6_normals.png",res=200,width=2000,height=2000)
heatmap.lgg.gbm <- heatmap.plus.sm(as.matrix(lgg.gbm.27.450.order),
                                col = jet.colors(75),
                                scale = "none",
                                labRow = NA,
                                labCol = NA,
                                Colv = NA,
                                Rowv = NA,
                                ColSideColors = clab.6
)
dev.off()
       
       
#Beta Value Diff M1 x M2
M1 <- xxMod[xxMod$K5.x == "1" & xxMod$IDH1.x == "Mutant",]
M1 <- subset(M1, IDH1.x != "NA")
M2 <-  xxMod[xxMod$K5.x == "4" & xxMod$IDH1.x == "Mutant",]
M2 <- subset(M2, IDH1.x != "NA")
#Selecting the Samples
M1 <- GBM.LGG.27.450k[, M1$TCGAID]
M2 <- GBM.LGG.27.450k[, M2$TCGAID]
volcano <- cbind(M1,M2)
volcano$meanM1 <- apply(volcano[,1:204],1,mean,na.rm=T)
volcano$meanM2 <- apply(volcano[,205:244],1,mean,na.rm=T)
volcano$DiffMean <- volcano$meanM1 - volcano$meanM2
       
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

###normals
load(file="normals.gse41826.rda")
norms <- exprs(controls.rm)
#dat.normals
normals <- merge(norms,dat.normals[,c("BRCA","KIRC","LUSC")], by=0)
rownames(normals) <- normals$Row.names
normals <- normals[,-1]
normals.sel <- normals[rownames(lgg.gbm.27.450.order),]
normals.sel <- normals.sel[,-1]
normals.sel <- normals.sel[heatmap.lgg.gbm$rowInd,]
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps/normals.png",res=200,width=2000,height=2000)
heatmap.plus.sm(as.matrix(normals.sel),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA
                #ColSideColors = clab.6
)
dev.off()

