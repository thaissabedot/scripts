library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(DiffBind)
load("/home/thais/sirt1_act.RData")
sirt_act <- dba.report(test_rmdup_sirt_act,th=1,method=DBA_DESEQ)

txdb_mm10 <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb_mm10)
symbols <- unlist(mapIds(org.Mm.eg.db, genes$gene_id, "SYMBOL", "ENTREZID", multiVals = "first"))
genes <- as.data.frame(genes)
genes$symbol <- as.matrix(symbols)
genes <- makeGRangesFromDataFrame(genes,TRUE)

a <- distanceToNearest(sirt_act,genes)
sirt_act.DF <- as.data.frame(sirt_act)
rownames(sirt_act.DF) <- NULL

sirt_act.DF$Gene_symbol <- NA
sirt_act.DF$Gene_EntrezID <- NA
sirt_act.DF$Distance <- NA

aux <- genes[subjectHits(a),]
sirt_act.DF[queryHits(a),"Gene_symbol"] <- mcols(aux)[,2]
sirt_act.DF[queryHits(a),"Gene_EntrezID"] <- mcols(aux)[,1]
sirt_act.DF[queryHits(a),"Distance"] <- mcols(a)[,1]

write.table(sirt_act.DF,file="/home/thais/NearestGene.txt",quote=F,row.names=F,sep="\t")