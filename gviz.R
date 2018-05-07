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


mean.Y.p <- apply(overlap.cpg[,substr(colnames(LGGrecurrent),1,16) %in% id.p.Y],1,mean, na.rm=TRUE)
mean.N.p <- apply(overlap.cpg[,substr(colnames(LGGrecurrent),1,16) %in% id.p.N],1,mean, na.rm=TRUE)
mean.Y.r <- apply(overlap.cpg[,substr(colnames(LGGrecurrent),1,16) %in% id.r.Y],1,mean, na.rm=TRUE)
mean.N.r <- apply(overlap.cpg[,substr(colnames(LGGrecurrent),1,16) %in% id.r.N],1,mean, na.rm=TRUE)


probes.df <- overlap.cpg[,1:4]
probes.GR <- GRanges(seqnames = paste0("chr",probes.df$Chromosome),ranges = IRanges(start = probes.df$Genomic_Coordinate, end=probes.df$Genomic_Coordinate),mean.Y.p=mean.Y.p,mean.N.p=mean.N.p,mean.Y.r=mean.Y.r,mean.N.r=mean.N.r)


genome(probes.GR) <- "hg19"
atrack <- AnnotationTrack(probes.GR, name = "CpG",stacking="dense")
#plotTracks(atrack)
gtrack <- GenomeAxisTrack()
#plotTracks(list(gtrack, atrack))
gen <- genome(probes.GR)
chr <- as.character(unique(seqnames(probes.GR)))
itrack <- IdeogramTrack(genome = gen, chromosome = chr)
#plotTracks(list(itrack, gtrack, atrack))
probes.GR.mean <- GRanges(seqnames = paste0("chr",probes.df$Chromosome),ranges = IRanges(start = probes.df$Genomic_Coordinate, end=probes.df$Genomic_Coordinate),mean.Y.p=mean.Y.p,mean.N.p=mean.N.p,mean.Y.r=mean.Y.r,mean.N.r=mean.N.r)

genome(probes.GR.mean) <- "hg19"


dTrack <- DataTrack(probes.GR.mean, name = "DNA methylation")
#plotTracks(list(itrack, gtrack, atrack, grtrack,dTrack), groups = c("Rad.p","NonRad.p","Rad.r","NonRad,r"), type = c("a", "p", "confint"))


knownGenes <- UcscTrack(genome = "hg19", chromosome = "chr2",track = "knownGene", from = min(probes.df$Genomic_Coordinate), to = max(probes.df$Genomic_Coordinate), trackType = "GeneRegionTrack",rstarts = "exonStarts", rends = "exonEnds", gene = "name", symbol = "name", transcript = "name", strand = "strand", name = "UCSC Genes")

ucsf <- UcscTrack(genome = "hg19", chromosome = "chr2",track = "ucsfBrainMethyl", from = min(probes.df$Genomic_Coordinate), to = max(probes.df$Genomic_Coordinate),trackType = "AnnotationTrack", start = "chromStart",end = "chromEnd", id = "name", shape = "box",fill = "#006400", name = "UCSF Brain Methyl")

miRNA <- UcscTrack(genome = "hg19", chromosome = "chr2",track = "wgRna", from = min(probes.df$Genomic_Coordinate), to = max(probes.df$Genomic_Coordinate),trackType = "GeneRegionTrack", start = "chromStart",end = "chromEnd", id = "name", shape = "box",fill = "#006400", name = "miRNA")

pdf("teste.pdf")
plotTracks(list(itrack, gtrack, atrack,dTrack,knownGenes,ucsf), groups = c("Rad.p","NonRad.p","Rad.r","NonRad.r"), type = c("horiz"),showSampleNames = TRUE,ylim=c(0,1))
dev.off()

