#After running "/home/thais/Tathi_PanCan/get-data.R"
setwd("~/Tathi_PanCan")
pcbc_rnaseq <- read.table("pcbc1.tsv",header=TRUE)

load("~/Tathi_PanCan/PCBC.groups.Rda")
load("~/Tathi_PanCan/PCBC.p62.Rda")

id <- as.matrix(colnames(pcbc_rnaseq)[-1])

library(org.Hs.eg.db)

a <- read.table("/home/thais/Tathi_PanCan/neighbor_genes.txt",header=T)
a$probeID <- rownames(a)
nearest.genes <- melt(a,id.vars = "probeID")
ensembl <- as.matrix(unlist(mapIds(org.Hs.eg.db, as.character(nearest.genes$value), "ENSEMBL", "ENTREZID", multiVals = "first")))
identical(rownames(ensembl),as.character(nearest.genes$value)) #tem que dar TRUE, senao me avisa e nao continua aqui ainda
nearest.genes$ensembl_id <- ensembl #59 entrez ID nao foram mapeados pra ensembl ID
rownames(pcbc_rnaseq) <- as.character(pcbc_rnaseq$GeneID)
pcbc_rnaseq <- pcbc_rnaseq[,-1]

nearest.genes <- nearest.genes[as.character(nearest.genes$ensembl_id) %in% rownames(pcbc_rnaseq),] #373 genes don't have gene expression data
gene_expression <- pcbc_rnaseq[as.character(nearest.genes$ensembl_id),]
dna_meth <- PCBC.p62[as.character(nearest.genes$probeID),]

#### IMPORTANTE: aqui precisa selecionar as amostras que tem ambos os dados e deixa-las na mesma ordem
cpg.gene <- cbind(dna_meth[,1:20],gene_expression[,1:20]) 
n.samples <- 0 #colocar aqui quantas amostras pareadas você tem

nearest.genes$correlation <- apply(cpg.gene,1, function(i) {
  zz <- cor(i[1:n.samples],i[(n.samples+1):ncol(cpg.gene)],method="spearman",use="complete")
  return(zz)
})

nearest.genes <- nearest.genes[nearest.genes$correlation < -0.2,] #pensei em ja filtrar os genes que tem uma correlacao positiva ou perto de zero .. mas nao eh obrigatorio
nearest.genes$probeID <- as.factor(nearest.genes$probeID)
nearest.genes <- do.call(rbind, unname(by(nearest.genes, nearest.genes$probeID, function(x) x[x$value == min(x$value),])))
