library(downloader)
download("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE15745&format=file&file=GSE15745%5FmRNA%5FPONS%5Fno%5Fnormalization%2Etxt%2Egz", "GSE15745_mRNA_PONS_no_normalization.txt.gz")
download("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE15745&format=file&file=GSE15745%5FmRNA%5FCRBLM%5Fno%5Fnormalization%2Etxt%2Egz", "GSE15745_mRNA_CRBLM_no_normalization.txt.gz")
download("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE15745&format=file&file=GSE15745%5FmRNA%5FFCTX%5Fno%5Fnormalization%2Etxt%2Egz","GSE15745_mRNA_FCTX_no_normalization.txt.gz")
download("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE15745&format=file&file=GSE15745%5FmRNA%5FTCTX%5Fno%5Fnormalization%2Etxt%2Egz", "GSE15745_mRNA_TCTX_no_normalization.txt.gz")

#metadata
metadata.mrna <- read.delim("GPL6104-11576.txt",skip=26) #primeiras linhas sao comentarios
library(seqRFLP)
dataframe2fas(metadata.mrna[,c("ID","SEQUENCE")],file="GPL6104.fa") #converts a data frame to a fasta file

################## MAPPING THROUGH BLAT

reverseProbes <- metadata.mrna[metadata.mrna$Probe_Chr_Orientation == "-","ID"]
forwardProbes <- metadata.mrna[metadata.mrna$Probe_Chr_Orientation == "+","ID"]
noInfoProbes <- metadata.mrna[metadata.mrna$Probe_Chr_Orientation == "","ID"]
#hg18
resulthg18 <- read.table("result_Hg18.psl",skip=2,header=T)
hg18resultForReverseProbes <- resulthg18[resulthg18$qName %in% reverseProbes,]
hg18resultForReverseProbes <- resultForReverseProbes[resultForReverseProbes$strand == "-",]
hg18resultForForwardProbes <- resulthg18[resulthg18$qName %in% forwardProbes,]
hg18resultForForwardProbes <- resultForForwardProbes[resultForForwardProbes$strand == "+",]
hg18resultForNoInfoProbes <- resulthg18[resulthg18$qName %in% noInfoProbes,]

#hg19
resulthg19 <- read.table("result_Hg19.psl",skip=2,header=T)
resultForReverseProbes <- resulthg19[resulthg19$qName %in% reverseProbes,]
hg19resultForReverseProbes <- resultForReverseProbes[resultForReverseProbes$strand == "-",]
resultForForwardProbes <- resulthg19[resulthg19$qName %in% forwardProbes,]
hg19resultForForwardProbes <- resultForForwardProbes[resultForForwardProbes$strand == "+",]
hg19resultForNoInfoProbes <- resulthg19[resulthg19$qName %in% noInfoProbes,]

#Remove duplicates
hg19resultForNoInfoProbes <- hg19resultForNoInfoProbes[!(duplicated(hg19resultForNoInfoProbes$qName) | duplicated(hg19resultForNoInfoProbes$qName, fromLast = TRUE)),]
hg19resultForReverseProbes <- hg19resultForReverseProbes[!(duplicated(hg19resultForReverseProbes$qName) | duplicated(hg19resultForReverseProbes$qName, fromLast = TRUE)),]
hg19resultForForwardProbes <- hg19resultForForwardProbes[!(duplicated(hg19resultForForwardProbes$qName) | duplicated(hg19resultForForwardProbes$qName, fromLast = TRUE)),]

#Select the highest scores among the duplicates
a <- hg19resultForNoInfoProbes[(duplicated(hg19resultForNoInfoProbes$qName) | duplicated(hg19resultForNoInfoProbes$qName, fromLast = TRUE)),]
b <- hg19resultForReverseProbes[(duplicated(hg19resultForReverseProbes$qName) | duplicated(hg19resultForReverseProbes$qName, fromLast = TRUE)),]
c <- hg19resultForForwardProbes[(duplicated(hg19resultForForwardProbes$qName) | duplicated(hg19resultForForwardProbes$qName, fromLast = TRUE)),]
HighestScore <- function(obj, duplicates){
  for(i in 1:length(levels(factor(duplicates$qName)))){
    aux <- duplicates[duplicates$qName %in% unique(duplicates$qName)[i],]
    if(nrow(aux[aux$matches == max(aux$matches),])==1)
      obj <- rbind(obj,aux[aux$matches == max(aux$matches),])
  }
    return(obj)
}

hg19resultForNoInfoProbes <- HighestScore(hg19resultForNoInfoProbes,a)
hg19resultForReverseProbes <- HighestScore(hg19resultForReverseProbes,b)
hg19resultForForwardProbes <- HighestScore(hg19resultForForwardProbes,c)

################## MAPPING THROUGH SIMON'S FILE (GENE ANNOTATION)

gene.annotation <- read.delim("Galaxy1-[Homo_sapiens_genes_(GRCh37.p13)].tabular", header=T, sep="\t")
dups <- duplicated(gene.annotation$EntrezGene.ID)
unique.genes <- gene.annotation[!dups, ]
a <- str_detect(unique.genes$Chromosome.Name, "^H")
b <- unique.genes[!a,]
a <- str_detect(b$Chromosome.Name, "^L")
b <- b[!a,]
a <- str_detect(b$Chromosome.Name, "^G")
b <- b[!a,]
library(GenomicRanges)
strand <- regmatches(b$Strand, regexpr("^[^[:digit:]]*", b$Strand))
for (i in 1:length(strand)){
  if(strand[i] == "")
    strand[i] <- "+"
}
b$Strand <- as.factor(strand)
gene.coord <- GRanges(seqnames = paste0("chr",b$Chromosome.Name),
                      ranges = IRanges(start = b$Gene.Start..bp., end=b$Gene.End..bp.), 
                      strand = b$Strand, 
                      symbol = b$Associated.Gene.Name)
result.blat <- GRanges(seqnames = resulthg19$tName,ranges = IRanges(start = resulthg19$tStart, end=resulthg19$tEnd), strand = resulthg19$strand, ProbeID = resulthg19$qName)
overlap.gene <- findOverlaps(result.blat,gene.coord)
overlap.gene <- as.data.frame(overlap.gene)
overlap.gene <- overlap.gene[unique(overlap.gene$queryHits),]
c <- resulthg19[overlap.gene$queryHits,]
gene.info <- b[overlap.gene$subjectHits,]

final <- cbind(c[,c("qName","tStart", "tEnd", "tName")],gene.info$Associated.Gene.Name)
save(final, file="GeneExpMap.rda")



load("epigenerowId.rda")
epiGene.probe.ids <- str_split(rownames(x), "[.]")
for(i in 1:length(epiGene.probe.ids)){
  epiGene.probe.ids[[i]] <- epiGene.probe.ids[[i]][[1]] 
}
epiGene.probe.ids <- as.matrix(epiGene.probe.ids)
library(FDb.InfiniumMethylation.hg19)
hm27.hg19 <- getPlatform(platform='HM27', genome='hg19')
hm27.hg19 <- as.data.frame(hm27.hg19)
overlap <- epiGene.probe.ids %in% rownames(hm27.hg19)


########### MAPPING 27K TO THE NEAREST GENES (SIMON'S FILE)

library(FDb.InfiniumMethylation.hg19)
hm27.hg19 <- getPlatform(platform='HM27', genome='hg19')
gene.annotation <- read.delim("Galaxy1-[Homo_sapiens_genes_(GRCh37.p13)].tabular", header=T, sep="\t")
dups <- duplicated(gene.annotation$EntrezGene.ID)
unique.genes <- gene.annotation[!dups, ]
a <- str_detect(unique.genes$Chromosome.Name, "^H")
b <- unique.genes[!a,]
a <- str_detect(b$Chromosome.Name, "^L")
b <- b[!a,]
a <- str_detect(b$Chromosome.Name, "^G")
b <- b[!a,]
library(GenomicRanges)
strand <- regmatches(b$Strand, regexpr("^[^[:digit:]]*", b$Strand))
for (i in 1:length(strand)){
  if(strand[i] == "")
    strand[i] <- "+"
}
b$Strand <- as.factor(strand)
gene.coord <- GRanges(seqnames = paste0("chr",b$Chromosome.Name),ranges = IRanges(start = b$Gene.Start..bp., end=b$Gene.End..bp.), strand = b$Strand, symbol = b$Associated.Gene.Name)
require(ChIPpeakAnno)
cpg.gene.loc.granges <- annotatePeakInBatch(RangedData(hm27.hg19), AnnotationData=RangedData(gene.coord, output="nearestStart")
cpg.gene.loc.granges <- as.data.frame(cpg.gene.loc.granges)
gene.coord <- as.data.frame(gene.coord)
for(i in 1:nrow(cpg.gene.loc.granges)){
    cpg.gene.loc.granges$Symbol[i] <- as.character(gene.coord[gene.coord$start == cpg.gene.loc.granges$start_position[i] & gene.coord$end == cpg.gene.loc.granges$end_position[i] & gene.coord$seqnames %in% cpg.gene.loc.granges$space[i], "symbol"])
}
nearest.gene <- cpg.gene.loc.granges[,c("peak","space","start","end","strand","distancetoFeature","Symbol")]
colnames(nearest.gene) <- c("ProbeID","Chromossome","Start","End", "Strand", "DistanceToGene","GeneSymbol")
save(nearest.gene,file="NearestGene27k.rda")

################# Gene expression data
Crblm <- read.delim("GSE15745_mRNA_CRBLM_no_normalization.txt")
Fctx <- read.delim("GSE15745_mRNA_FCTX_no_normalization.txt")
Pons <- read.delim("GSE15745_mRNA_PONS_no_normalization.txt")
Tctx <- read.delim("GSE15745_mRNA_TCTX_no_normalization.txt")
#Leave only the columns with signal and p value information
Columns <- function(obj){
  col.names <- str_detect(colnames(obj), "(Signal)$|(Pval)$")
  obj <- cbind(obj[,c("PROBE_ID", "SYMBOL")],obj[,col.names])
  return(obj)
}
Crblm <- Columns(Crblm)
Fctx <- Columns(Fctx)
Pons <- Columns(Pons)
Tctx <- Columns(Tctx)
#Mask the probes that have p.value > 0.05
P.value.mask <- function(obj){
   p.cols <- str_detect(colnames(obj), "(Pval)$")
   s.cols <- str_detect(colnames(obj), "(Signal)$")
   a <- obj[,p.cols]
   b <- obj[,s.cols]
   a <- a <= 0.05
   a[a==FALSE] <- NA 
   a[a==TRUE] <- FALSE
   obj[,s.cols] <- a + b
   obj$mean <- apply(obj[,s.cols], 1, mean,na.rm=TRUE)
   return(obj)
   
}
Crblm <- P.value.mask(Crblm)
Fctx <- P.value.mask(Fctx)
Pons <- P.value.mask(Pons)
Tctx <- P.value.mask(Tctx)
brain.regions.exp <- cbind(Crblm$mean, Fctx$mean, Pons$mean, Tctx$mean)
colnames(brain.regions.exp) <- c("Crblm","Fctx","Pons","Tctx")
brain.regions.exp <- as.data.frame(brain.regions.exp)
brain.regions.exp$ProbeID <- Crblm$PROBE_ID
brain.regions.exp$Symbol.Original <- Crblm$SYMBOL
brain.regions.exp<- merge(brain.regions.exp, final[,c(1,5)], by.x='ProbeID', by.y = 'qName', all.x =T)                                      
library(LSD)
colnames(brain.regions.exp)[7] <- "Symbol.New"
rownames(brain.regions.exp) <- brain.regions.exp$ProbeID
heatpairs(as.matrix(log10(brain.regions.exp[intersect(final$qName,brain.regions.exp$ProbeID),1:4])))
                                            

                                            
                                            
                                            
                                            
                                            
### download data using getGEO and using the default Symbols derived from GEO.
gse15745.e <- getGEO(GEO="GSE15745", destdir="./")
exp.brain <- gse15745.e$`GSE15745-GPL6104_series_matrix.txt.gz`
exp.brain.norm <- lumiN(exp.brain)
exp.brain.pdata <- pData(exp.brain.norm)
Crblm.e <- apply(exprs(exp.brain.norm)[,rownames(subset(exp.brain.pdata, description=="Human Brain Tissue: CRBLM"))], 1, mean, na.rm=T)
Fctx.e <- apply(exprs(exp.brain.norm)[,rownames(subset(exp.brain.pdata, description=="Human Brain Tissue: FCTX"))], 1, mean, na.rm=T)
Pons.e <- apply(exprs(exp.brain.norm)[,rownames(subset(exp.brain.pdata, description=="Human Brain Tissue: PONS"))], 1, mean, na.rm=T)
Tctx.e <- apply(exprs(exp.brain.norm)[,rownames(subset(exp.brain.pdata, description=="Human Brain Tissue: TCTX"))], 1, mean, na.rm=T)
brain.regions.exp.N <- cbind(Crblm.e, Fctx.e, Pons.e, Tctx.e)
brain.regions.exp.N <- as.data.frame(brain.regions.exp.N)
brain.regions.exp.N$Symbols <- as.character(fData(exp.brain.norm)$Symbol)
save(brain.regions.exp.N, file="brain.regions.exp.N.rda")                                            
                                            
                                            
                                            
                                            
                                            