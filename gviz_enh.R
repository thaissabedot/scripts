https://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.pdf
library(Gviz)

#chromosome
a <- AnnotationTrack(start=174000000, end=175000000,chromosome="chr4", genome="hg19", name="AnnotationTrack")
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg19", chromosome = "chr4")

#genes
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
txTranscripts_v1 <- GeneRegionTrack(txdb_hg19, genome="hg19", chromosome="chr4", showId=TRUE, geneSymbol=TRUE, name="UCSC")
symbols <- unlist(mapIds(org.Hs.eg.db, gene(txTranscripts_v1), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(txTranscripts_v1) <- symbols[gene(txTranscripts_v1)]

#se você quiser colocar se o gene está up ou downregulated com cores diferentes, rodar essa parte. Senão pode pular
overlaps <- findOverlaps(makeGRangesFromDataFrame(volcano.map), aux.GR)
b <- volcano.map[queryHits(overlaps),]
b <- b[with(b,order(seqnames,start)),]
txTranscripts_v1 <- GeneRegionTrack(makeGRangesFromDataFrame(b[,c("seqnames","start","end")]),symbol=as.character(b$symbol), genome="hg19", chromosome="chr4", showId=TRUE, geneSymbol=TRUE, name="Genes",fill=as.character(b$sig))
#fim cor do gene up and down

aux.GR <- GRanges(seqnames = "chr4",ranges = IRanges(start = 174000000, end=175000000))

#aqui são seus enhancers 
overlaps <- findOverlaps(makeGRangesFromDataFrame(LaminB1), aux.GR)
b <- LaminB1[queryHits(overlaps),]
foo2 <- makeGRangesFromDataFrame(b,TRUE) 
genome(foo2) <- "hg19"
foo2 <- AnnotationTrack(foo2,name="LaminB1",col="black",fill="black")

#se quiser colocar dados de expressão gênica, metilação do DNA, etc.
a <- tudao[tudao$chrom %in% "chr4",]
a.GR <- makeGRangesFromDataFrame(a,TRUE)
overlaps <- findOverlaps(a.GR, aux.GR)
a <- a[queryHits(overlaps),]
meth <- makeGRangesFromDataFrame(a[,c(1:6,10,11)],TRUE) 
genome(meth) <- "hg19"
meth <- DataTrack(meth,name="DNA methylation",type="smooth",groups=c("GCIMP.low",rep("GCIMP.high",2),rep("Non-tumor",2)),col=c("firebrick","darkgreen","grey"))
#fim dados de metilação do DNA

plotTracks(list(itrack,gtrack,txTranscripts_v1,foo2,meth),legend=TRUE,from=174000000,to=175000000)

