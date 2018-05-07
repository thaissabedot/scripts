#### find the top -10/+10 genes that are around an enhancer
#get gene coordinates
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- genes(txdb_hg19)
symbols <- unlist(mapIds(org.Hs.eg.db, genes$gene_id, "SYMBOL", "ENTREZID", multiVals = "first"))
genes <- as.data.frame(genes)
genes$symbol <- as.matrix(symbols)

dados.RNAseq <- merge(genes,dados.RNAseq,by="gene_id")
dados.RNAseq.GR <- makeGRangesFromDataFrame(dados.RNAseq[,1:7],TRUE) #nao precisa por os dados de expressão aqui, só as coordenadas dos genes 

enhancer.GR <- makeGRangesFromDataFrame(enhancer,TRUE)

genes.precede <- matrix(data=NA,nrow=nrow(enhancer),ncol=10,dimnames=list(1:nrow(enhancer),paste("Precede.round",1:10,sep=".")))

for(i in 1:10){
  dados.RNAseq.GR <- makeGRangesFromDataFrame(dados.RNAseq.GR,TRUE)
  aux <- precede(enhancer.GR,dados.RNAseq.GR)
  genes.precede[which(!is.na(aux)),i] <- as.character(as.data.frame(mcols(dados.RNAseq.GR[na.omit(aux),])[1])[,1]) #a sua primeira coluna de metadado no GRanges "dados.RNAseq.GR" deve ter o entrezID
  dados.RNAseq.GR <- dados.RNAseq.GR[-(na.omit(aux)),]
  
}

dados.RNAseq.GR <- makeGRangesFromDataFrame(dados.RNAseq[,1:7],TRUE) 
enhancer.GR <- makeGRangesFromDataFrame(enhancer,TRUE)

genes.follow <- matrix(data=NA,nrow=nrow(enhancer),ncol=10,dimnames=list(1:nrow(enhancer),paste("Follow.round",1:10,sep=".")))

for(i in 1:10){
  dados.RNAseq.GR <- makeGRangesFromDataFrame(dados.RNAseq.GR,TRUE)
  aux <- follow(enhancer.GR,dados.RNAseq.GR)
  genes.precede[which(!is.na(aux)),i] <- as.character(as.data.frame(mcols(dados.RNAseq.GR[na.omit(aux),])[1])[,1])
  dados.RNAseq.GR <- dados.RNAseq.GR[-aux,]
  
}

genes.all <- cbind(genes.precede,genes.follow)
genes.all <- as.data.frame(genes.all)
dados.RNAseq.GR <- makeGRangesFromDataFrame(dados.RNAseq[,1:7],TRUE) 
aux <- findOverlaps(enhancer.GR,dados.RNAseq.GR)
genes.all$overlap <- NA
genes.all[queryHits(aux),"overlap"] <- dados.RNAseq[subjectHits(aux),"gene_id"]

#### ordenar o p-value dos 20 genes (10 antes e 10 depois) e pegar o menor 
