##Genomic annotation of 450k probes
##Get 450k data
setwd("/dados/camila/TCGA_GBM_Recurrent_Project/")
dim(LGGrecurrent)  ##485577
table(LGGrecurrent$Chromosome)  ##chr X,Y,1:22
HM450k <- LGGrecurrent[,c(1:4)]  ##485577
HM450k$Genomic_Coordinate_end <- HM450k$Genomic_Coordinate
HM450k$probetype <- as.character(rownames(HM450k))  ##485577
library(stringr)
HM450k$probetype <- substr(HM450k$probetype,1,2)
table(HM450k$probetype)  ##cg: 482421  ch: 3091  rs: 65
HM450k.cg <- subset(HM450k,probetype %in% paste("cg"))  ##482421
HM450k.cg$type <- "NA"  ##482421
##Download CpG islands data
##UCSC Genome Bioinformatics (genome.ucsc.edu)
##Table Browser
##Assembly: Feb.2009 (GRCh37/hg19)
##Group: Regulation
##Track: CpG Islands
##island <- read.delim("/dados/camila/Recurrent.gliomas/Genomic_Files/island")  ##data.frame  ##28691
##table(island$chrom)
##island <- subset(island, chrom %in% paste("chr",c(1:22,"X","Y"),sep=""))  ##take into consideration probes located at chromosomes X,Y,1:22  ##27718
##table(island$chrom)
##island$chrom <- factor(island$chrom)
##table(island$chrom)
##Create GenomicRanges objects
library(GenomicRanges)
HM450k.cg_GR <- GRanges(seqnames = 
                paste0("chr",HM450k.cg$Chromosome),ranges = IRanges(start =
                HM450k.cg$Genomic_Coordinate,end = HM450k.cg$Genomic_Coordinate_end),
                cgID = HM450k.cg$Composite.Element.REF,geneID = HM450k.cg$Gene_Symbol,site = HM450k.cg$type)
length(HM450k.cg_GR)  ##482421
##island_GR <- GRanges(seqnames = island$chrom,ranges = IRanges(start =
             ##island$chromStart,end = island$chromEnd),
             ##id = island$name)
##length(island_GR)  ##27718
##Shores
##shore_GR <- flank(island_GR, width = 2000, both = TRUE)
##length(shore_GR)  ##27718
##Find overlaps between GenomicRanges objects
HM450k.cg.island <- findOverlaps(HM450k.cg_GR,island_GR)
length(HM450k.cg.island)  ##150253
HM450k.cg.shore <- findOverlaps(HM450k.cg_GR,shore_GR)
length(HM450k.cg.shore)  ##278003
table(HM450k.cg$type)  ##NA: 482421
##Identify if there is common probes between HM450k.cg.island and HM450k.cg.shore. If positive, consider these probes as located at CpG islands
HM450k.cg.island
HM450k.cg.island_df <- as.data.frame(HM450k.cg.island)  ##150253
HM450k.cg.shore_df <- as.data.frame(HM450k.cg.shore)  ##278003
removeprobe450k <- (intersect(rownames(HM450k.cg[HM450k.cg.island_df$queryHits,]), rownames(HM450k.cg[HM450k.cg.shore_df$queryHits,])))
retainprobe450k <- setdiff(rownames(HM450k.cg[unique(HM450k.cg.shore_df$queryHits),]), removeprobe450k)
##Catalogue the genomic location of probes
HM450k.cg[retainprobe450k,"type"] <- "Shores (flanking CpG Island - 2000bp)"
HM450k.cg[HM450k.cg.island_df$queryHits,"type"] <- "CpG Island (Irizarry's HMM)"
HM450k.cg[HM450k.cg$type == "NA","type"] <- "Open Sea"
table(HM450k.cg$type)  ##CpG Island (Irizarry's HMM): 150253  ##Open Sea: 229269  ##Shores (flanking CpG Island - 2000bp): 102899
##Identify which probes are located at coding regions
library(VariantAnnotation)
library(AnnotationHub)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
loc_hg19_450k <- locateVariants(HM450k.cg_GR, txdb_hg19, CodingVariants())
table(loc_hg19_450k$LOCATION)  ##coding: 110669
loc_hg19_450k <- as.data.frame(loc_hg19_450k)  ##110669
probes450k <- merge(HM450k.cg,loc_hg19_450k,by.x=c("Chromosome","Genomic_Coordinate","Genomic_Coordinate_end"),by.y=c("seqnames","start","end"),sort=FALSE)
HM450k.cg$location <- NA
HM450k.cg[unique(probes450k$Composite.Element.REF),"location"] <- "Coding"
dim(HM450k.cg)  ##482421
table(HM450k.cg$location)  ##Coding: está dando 0 (deveria dar 110669)
length(unique(biofeatures$name))  
110669/482421  ##23% of 450k cg probes are located in coding regions