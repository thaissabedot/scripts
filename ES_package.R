####### START mapping CpG probes

library("openxlsx")
library(rtracklayer)
### Hg19
#Data from: "Validation of a DNA methylation microarray for 850,000 CpG sites of the human genome enriched in enhancer sequences"
a <- read.xlsx("/home/thais/Downloads/850kand450k.xlsx",sheet=1)
b <- read.xlsx("/home/thais/Downloads/850k_only.xlsx",sheet=1)
load("/home/thais/Cloud_SanDiego/TCGA/LGG.GBM/NewSamples_16-09/gbm450.Rda")
CpG.probes <- data.frame(Composite.Element.REF=c(a$IlmnID,b$IlmnID),Chromosome=c(a$CHR,b$CHR),Genomic_Coordinate=c(a$MAPINFO,b$MAPINFO))
gbm.450k <- gbm.450k[(!(as.character(gbm.450k$Composite.Element.REF) %in% as.character(CpG.probes$Composite.Element.REF)) & substr(gbm.450k$Composite.Element.REF,1,2) %in% "cg"),] #get 450k probes that are not in 850k
CpG.probes <- rbind(CpG.probes,gbm.450k[,c(1,3,4)])
CpG.probes$Chromosome <- paste0("chr",CpG.probes$Chromosome)
hg19 <- CpG.probes
colnames(hg19) <- c("probeID","chrom","start")
rownames(hg19) <- as.character(hg19$probeID)

CpG.probes$end <- CpG.probes$Genomic_Coordinate

colnames(CpG.probes) <- c("name","chrom","chromStart","chromEnd")
CpG.probes <- CpG.probes[,c(2,3,4,1)]
CpG.probes$chromStart <- as.integer(CpG.probes$chromStart)
CpG.probes$chromEnd <- as.integer(CpG.probes$chromEnd) #896.166 (853.307 HM850/HM450 + 42.859 HM450)

### Hg18
chain <- import.chain("/home/thais/Package/hg19ToHg18.over.chain")
hg18 <- unlist(liftOver(makeGRangesFromDataFrame(CpG.probes,TRUE),chain))
hg18 <- as.data.frame(hg18)
hg18 <- hg18[,c(6,1,2)]
colnames(hg18) <- c("probeID","chrom","start")
rownames(hg18) <- as.character(hg18$probeID)

### Hg38
chain <- import.chain("/home/thais/Package/hg19ToHg38.over.chain")
hg38 <- unlist(liftOver(makeGRangesFromDataFrame(CpG.probes,TRUE),chain))
hg38 <- as.data.frame(hg38)
hg38 <- hg38[,c(6,1,2)]
colnames(hg38) <- c("probeID","chrom","start")
rownames(hg38) <- as.character(hg38$probeID)

CpG.probes <- list(hg18=hg18,hg19=hg19,hg38=hg38)
rm(a,b,chain,gbm.450k,hg18,hg19,hg38); gc();

####### END mapping CpG probes

####### START mapping genes
library(org.Hs.eg.db)

### Hg19
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
genes <- as.data.frame(genes)
genes$symbol <- as.matrix(unlist(mapIds(org.Hs.eg.db, genes$gene_id, "SYMBOL", "ENTREZID", multiVals = "first")))
colnames(genes) <- c("chrom","start","end","width","strand","gene.id","gene.symbol")
aux <- as.data.frame(distanceToNearest(makeGRangesFromDataFrame(CpG.probes$hg19,end="start"), makeGRangesFromDataFrame(genes)))
CpG.probes$hg19 <- cbind(CpG.probes$hg19,genes[aux$subjectHits,c("gene.id","gene.symbol")],distance=aux$distance)

### Hg18
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
genes <- genes(TxDb.Hsapiens.UCSC.hg18.knownGene)
genes <- as.data.frame(genes)
genes$symbol <- as.matrix(unlist(mapIds(org.Hs.eg.db, genes$gene_id, "SYMBOL", "ENTREZID", multiVals = "first")))
colnames(genes) <- c("chrom","start","end","width","strand","gene.id","gene.symbol")
aux <- as.data.frame(distanceToNearest(makeGRangesFromDataFrame(CpG.probes$hg18,end="start"), makeGRangesFromDataFrame(genes)))
CpG.probes$hg18 <- cbind(CpG.probes$hg18,genes[aux$subjectHits,c("gene.id","gene.symbol")],distance=aux$distance)

### Hg38
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
genes <- as.data.frame(genes)
genes$symbol <- as.matrix(unlist(mapIds(org.Hs.eg.db, genes$gene_id, "SYMBOL", "ENTREZID", multiVals = "first")))
colnames(genes) <- c("chrom","start","end","width","strand","gene.id","gene.symbol")
aux <- as.data.frame(distanceToNearest(makeGRangesFromDataFrame(CpG.probes$hg38,end="start"), makeGRangesFromDataFrame(genes)))
CpG.probes$hg38 <- cbind(CpG.probes$hg38,genes[aux$subjectHits,c("gene.id","gene.symbol")],distance=aux$distance)

rm(genes,aux); gc();

####### END mapping genes


####### START getting data to test

library(stringr)
load("/home/thais/Cloud_SanDiego/TCGA/LGG.GBM/RNAseq_05-11-14/mRNA.Rda")
mrna <- merge(lgg.mrna,gbm.mrna,by="gene_id")
colnames(mrna) <- substr(colnames(mrna),1,12)
mrna <- mrna[,-which(names(mrna) %in% c("TCGA-HT-A61A"))] #remove this sample
entrezID <- str_extract(mrna$gene_id,"[^|]*$")
mrna$entrezID <- as.matrix(entrezID)
rownames(mrna) <- as.character(mrna$entrezID)
mrna <- mrna[,-c(1,669)]
load("/home/thais/Package/control_test.rda")
rm(lgg.mrna,gbm.mrna,entrezID); gc();

####### END getting data to test

####### START getting percent GC by gene
host <- 'ensembl.org'
ensembl <- useMart(host=host, biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
dframe <- getBM(attributes=c("entrezgene", "percentage_gc_content"), filters=c("entrezgene"), values=rownames(mrna), mart=ensembl)
hg38 <- dframe[!duplicated(dframe$entrezgene),]
rownames(hg38) <- as.character(hg38$entrezgene)

host <- 'grch37.ensembl.org'
ensembl <- useMart(host=host, biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
dframe <- getBM(attributes=c("entrezgene", "percentage_gc_content"), filters=c("entrezgene"), values=rownames(mrna), mart=ensembl)
hg19 <- dframe[!duplicated(dframe$entrezgene),]
rownames(hg19) <- as.character(hg19$entrezgene)

host <- 'may2009.archive.ensembl.org'
ensembl <- useMart(host=host, biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
dframe <- getBM(attributes=c("entrezgene", "percentage_gc_content"), filters=c("entrezgene"), values=rownames(mrna), mart=ensembl)
hg18 <- dframe[!duplicated(dframe$entrezgene),]
rownames(hg18) <- as.character(hg18$entrezgene)

percent.GC <- list(hg18=hg18,hg19=hg19,hg38=hg38)
save(percent.GC,file="/home/thais/Package/EScall/R/percent_GC.rda")
rm(hg18,hg19,hg38,dframe,ensembl,host); gc();
####### END getting percent GC by gene

### Libraries
library(EDASeq)
library(steFunctions)
library(biomaRt)
library(stringr)


############ START TESTING DATA
load("/home/thais/Cloud_SanDiego/TCGA/LGG.GBM/LGG_GBM_obj_20150127.Rda")
load("/home/thais/LGG-GBM/LGG-GBM.merged.data.FB20150930v3.Rdata") #pdata from Sep 2015 (LGG-GBM project). Created by Floris Barthel
rownames(pd) <- as.character(pd$case.id)
genomeVersion <- "hg19"
meth <- LGG.GBM.250914[,-c(1:4)]
cut.off <- 0.5
control.cutoff <- 0.5
rm(cc.out.purity.lgg.gbm.new,full.heatmap.orderd,metadata.1129.samples.20150112,norms,enrichmentList.sig,dat.lgg.gbm.new.noXY.dic.oecg,cpg.gene,cpg.gene.normal, epiGene.RNAseq,LGG.GBM.250914); gc();
############ END TESTING DATA

meth <- meth[,intersect(colnames(meth),colnames(mrna))]
mrna <- mrna[,intersect(colnames(meth),colnames(mrna))]

############ START TESTING GROUPS
groups <- as.character(pd[colnames(meth),"clustM.olddiscovery"]); rm(pd); gc();
############ END TESTING GROUPS

normalization <- function(mrna,genomeVersion) {
  load("/home/thais/Package/EScall/R/percent_GC.rda")
  percent.GC <- percent.GC[[genomeVersion]]
  percent.GC <- percent.GC[rownames(percent.GC) %in% rownames(mrna),]
  percent.GC <-
  mrna <- mrna[rownames(percent.GC),]
  mrna <- floor(as.matrix(mrna))
  tmp <- newSeqExpressionSet(mrna, featureData = percent.GC)
  tmp <- withinLaneNormalization(tmp, "percentage_gc_content", which = "upper", offset = TRUE)
  tmp <- betweenLaneNormalization(tmp, which = "upper", offset = TRUE)
  normCounts <- log(mrna + .1) + offst(tmp)
  normCounts <- floor(exp(normCounts) - .1)
  tmp <- t(quantileNormalization(t(normCounts)))
  mrna <- floor(tmp)
  mrna <- as.data.frame(mrna)
  return(mrna)
}

#' @param meth A DNA matrix...
#' @param mrna
#' @param groups
#' @param genomeVersion
#' @param control.meth
#' @param control.mrna
#' @param normalize
#' @param cores
#' @param cutoff
#' @param control.cutoff
#' @import biomaRt EDASeq steFunctions
ES_call <- function(meth,
                    mrna,
                    groups=NULL,
                    genomeVersion="hg38",
                    control.meth=NULL,
                    control.mrna=NULL,
                    normalize=FALSE,
                    cores=1,
                    cutoff=0.5,
                    control.cutoff=NULL) {

  if((!genomeVersion %in% c("hg18","hg19","hg38")) | length(genomeVersion) > 1)
    stop("ERROR: Genome version not supported")

  if(cutoff > 1 | cutoff < 0)
    stop("ERROR: cutoffs must be between 0 and 1")

  if(!is.null(control.cutoff)){
    if(control.cutoff > 1 | control.cutoff < 0)
      stop("ERROR: cutoffs must be between 0 and 1")
  }

  if(cores > detectCores() | cores == 0)
    cores <- detectCores()

  if(Sys.info()["sysname"] == "Windows") cores <- 1;

  ### START select probes and genes based on genome version
  meth.info <- CpG.probes[[genomeVersion]]
  meth.info <- meth.info[(meth.info$probeID %in% rownames(meth)) & (meth.info$gene.id %in% rownames(mrna)),]
  ### END select probes and genes based on genome version

  ### START prepare the data and sort the genes according to CpG mapping
  meth <- meth[as.character(meth.info$probeID),intersect(colnames(meth),colnames(mrna))]
  mrna <- mrna[as.character(meth.info$gene.id),intersect(colnames(meth),colnames(mrna))]
  ### END prepare the data and sort the genes according to CpG mapping

  ### START normalization
  if(normalize==TRUE){
    if(!is.null(control.mrna)) { #Normalize data with control samples
      mrna <- merge(mrna,control.mrna,by=0,all.x=TRUE,sort=F)
      rownames(mrna) <- as.character(mrna$Row.names)
      mrna <- mrna[,-1]
    }
    mrna <- normalization(mrna,genomeVersion)
    if(!is.null(control.mrna)) { #Separate control samples from the dataset
      control.mrna <- mrna[,colnames(control.mrna)]
      mrna <- mrna[,!(colnames(mrna) %in% colnames(control.mrna))]
    }

  }
  ### END normalization


  if(length(groups) != ncol(meth))
    stop("Group assignments don't match with number of samples")

  if(ncol(meth)<2 | nrow(meth)<2)
    stop("Not enough observations")

  if(!is.null(control.mrna) & !is.null(control.meth)) {
    control.meth <- control.meth[rownames(meth),intersect(colnames(control.meth),colnames(control.mrna))]
    control.mrna <- control.mrna[rownames(mrna),intersect(colnames(control.meth),colnames(control.mrna))]
    if(ncol(control.meth)<2 | nrow(control.meth)<2 | nrow(control.mrna)<2) {
      writeLines("Not enough observations for control group \n Setting to NULL")
      control.meth <- NULL
      control.mrna <- NULL
    }
    else if(is.null(control.cutoff)){
      print("Setting cut off for control to 50%")
      control.cutoff <- 0.5
    }

  }


  if(!is.null(control.mrna) & !is.null(control.meth)) {
    epiGene <- mclapply(1:nrow(meth),
                        function(i){
                          meth.gu <- names(meth)[meth[i, ] < 0.3]
                          meth.gm <- names(meth)[meth[i, ] >= 0.3]
                          if(length(meth.gu)>1 & length(meth.gm)>1){
                            expr.mean.gm <- mean(as.numeric(mrna[i, meth.gm]), na.rm=T)
                            if(expr.mean.gm < 1.28*sd(as.numeric(mrna[i, meth.gu]), na.rm=T)){
                              expr.mean.gu <- mean(as.numeric(mrna[i, meth.gu]), na.rm = T)
                              yy <- cbind(t(meth[i, ] >= 0.3), t(mrna[i, ] < expr.mean.gu))
                              colnames(yy) <- c("meth.gt0.3", "exp.lt.meanUnmeth")
                              aux <- cbind(t(control.meth[rownames(meth[i,]), ] >= 0.3), t(control.mrna[rownames(mrna[i,]), ] < expr.mean.gu))
                              colnames(aux) <- c("meth.gt0.3", "exp.lt.meanUnmeth")
                              yy <- rbind(yy, aux)
                              yy <- as.data.frame(yy)
                              probeID <- rownames(meth[i,])
                              if(dim(table(yy$meth.gt0.3[1:ncol(meth)])) == 2  & dim(table(yy$exp.lt.meanUnmeth[1:ncol(meth)])) == 2){
                                if(as.matrix(table(yy$meth.gt0.3[1:ncol(meth)], yy$exp.lt.meanUnmeth[1:ncol(meth)]))[2,2]/as.numeric(table(yy$meth.gt0.3[1:nrow(meth)])[2])>0.8){
                                  yy$probeID <- probeID
                                  yy$id <- rownames(yy)
                                  yy$methValues <- NA
                                  yy$exprValues <- NA
                                  yy$methValues[1:ncol(meth)] <- as.numeric(meth[i, ])
                                  yy$exprValues[1:ncol(meth)] <- as.numeric(mrna[i, ])
                                  yy$methValues[(ncol(meth)+1):nrow(yy)] <- as.vector(as.matrix(t(control.meth[rownames(meth[i,]), ])))
                                  yy$exprValues[(ncol(meth)+1):nrow(yy)] <- as.vector(as.matrix(t(control.mrna[rownames(mrna[i,]), ])))
                                  yy$groups <- NA
                                  if(!is.null(groups)) yy$groups[1:ncol(meth)] <- as.matrix(groups);
                                  return(yy)
                                }
                              }
                            }
                          }
                        }
                        , mc.cores=cores)
    names(epiGene) <- rownames(meth)
    epiGene <- epiGene[!sapply(epiGene, is.null)]

    uu.nl <- epiGene[[1]]
    clusters <- levels(factor(groups))


    enrichment.normal <- lapply(clusters, function(cluster, epiGene.ss=epiGene) {
      enrichment.group <- unlist(mclapply(epiGene.ss, function(uu.nl, cluster.group=cluster) {
        control <- (!uu.nl$meth.gt0.3[(ncol(meth)+1):nrow(uu.nl)] & !uu.nl$exp.lt.meanUnmeth[(ncol(meth)+1):nrow(uu.nl)])
        uu.nl$es <- NA
        uu.nl$es[1:ncol(meth)] <- uu.nl$meth.gt0.3[1:ncol(meth)] & uu.nl$exp.lt.meanUnmeth[1:ncol(meth)]
        uu.nl$Kclass <- NA
        uu.nl$Kclass <- uu.nl$groups == cluster.group
        es.table <- table(uu.nl$Kclass, uu.nl$es)
        if(nrow(es.table) == 2 & ncol(es.table) == 2){
          if((es.table[2,2]/rowSums(es.table)[2] > cutoff ) & (es.table[2,2]/colSums(es.table)[2] > cutoff ) & (table(control)[2]/(table(control)[2] + table(control)[1]) > control.cutoff  | is.na(table(control)[2]/(table(control)[2] + table(control)[1])))) {
            return(fisher.test(table(uu.nl$Kclass, uu.nl$es))$p.value)
          }
        } else {
          return(NULL)
        }
      }, mc.cores=cores))
      return(enrichment.group)
    })

    names(enrichment.normal) <- clusters
    enrichment.normal.adj <- lapply(enrichment.normal, p.adjust, method="BH")
    enrichmentList.sig <- lapply(enrichment.normal.adj, function(x) {
      x <- x[x < 0.05]
    })

  } else {

    epiGene <- mclapply(1:nrow(meth),
                        function(i){
                          meth.gu <- names(meth)[meth[i, ] < 0.3]
                          meth.gm <- names(meth)[meth[i, ] >= 0.3]
                          if(length(meth.gu)>1 & length(meth.gm)>1){
                            expr.mean.gm <- mean(as.numeric(mrna[i, meth.gm]), na.rm=T)
                            if(expr.mean.gm < 1.28*sd(as.numeric(mrna[i, meth.gu]), na.rm=T)){
                              expr.mean.gu <- mean(as.numeric(mrna[i, meth.gu]), na.rm = T)
                              yy <- cbind(t(meth[i, ] >= 0.3), t(mrna[i, ] < expr.mean.gu))
                              colnames(yy) <- c("meth.gt0.3", "exp.lt.meanUnmeth")
                              yy <- as.data.frame(yy)
                              probeID <- rownames(meth[i,])
                              if(dim(table(yy$meth.gt0.3[1:ncol(meth)])) == 2  & dim(table(yy$exp.lt.meanUnmeth[1:ncol(meth)])) == 2){
                                if(as.matrix(table(yy$meth.gt0.3[1:ncol(meth)], yy$exp.lt.meanUnmeth[1:ncol(meth)]))[2,2]/as.numeric(table(yy$meth.gt0.3[1:nrow(meth)])[2])>0.8){
                                  yy$probeID <- probeID
                                  yy$id <- rownames(yy)
                                  yy$methValues <- as.numeric(meth[i, ])
                                  yy$exprValues <- as.numeric(mrna[i, ])
                                  if(!is.null(groups)) yy$groups <- as.matrix(groups);
                                  return(yy)
                                }
                              }
                            }
                          }
                        }
                        ,mc.cores=cores)
    names(epiGene) <- rownames(meth)
    epiGene <- epiGene[!sapply(epiGene, is.null)]

    uu.nl <- epiGene[[1]]
    clusters <- levels(factor(groups))

    enrichment.normal <- lapply(clusters, function(cluster, epiGene.ss=epiGene) {
      enrichment.group <- unlist(mclapply(epiGene.ss, function(uu.nl, cluster.group=cluster) {
        uu.nl$es <- NA
        uu.nl$es <- uu.nl$meth.gt0.3 & uu.nl$exp.lt.meanUnmeth
        uu.nl$Kclass <- NA
        uu.nl$Kclass <- uu.nl$groups == cluster.group
        es.table <- table(uu.nl$Kclass, uu.nl$es)
        if(nrow(es.table) == 2 & ncol(es.table) == 2){
          if((es.table[2,2]/rowSums(es.table)[2] > cutoff) & (es.table[2,2]/colSums(es.table)[2] > cutoff)) {
            return(fisher.test(table(uu.nl$Kclass, uu.nl$es))$p.value)
          }
        } else {
          return(NULL)
        }
      }, mc.cores=cores)) 
      return(enrichment.group)
    })

    names(enrichment.normal) <- clusters
    enrichment.normal.adj <- lapply(enrichment.normal, p.adjust, method="BH")
    enrichmentList.sig <- lapply(enrichment.normal.adj, function(x) {
      x <- x[x < 0.05]
    })
  }
  return(enrichmentList.sig)
}

