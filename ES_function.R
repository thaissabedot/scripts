
#O objeto Tumor.DNAmethylation já deve conter uma coluna com o gene mais próximo mapeado
cpg.gene <- merge(Tumor.DNAmethylation,Tumor.GeneExpression,by="GeneID") 
#O nome da linha vai ser o nome da probe + . + gene. Ex: cg123456.ABC
dimnames(cpg.gene)[[1]] <- paste(cpg.gene[,"ProbeID"], cpg.gene[,"GeneID"], sep=".")

meth.g <- cpg.gene[,1:100] #Colunas que contenham SOMENTE os dados de metilação do DNA 
names(meth.g) <- substr(names(meth.g), 1, 12) #Barcode precisa ser cortado dos caracteres de 1 à 12 para ficarem iguais entre experimentos diferentes.
expr.g <- cpg.gene[,101:200] #Colunas que contenham SOMENTE os dados de expressão gênica
names(expr.g) <- substr(names(expr.g), 1, 12) #Barcode precisa ser cortado dos caracteres de 1 à 12 para ficarem iguais entre experimentos diferentes.

#As amostras precisam estar na mesma ordem para metilação do DNA e expressão gênica, ou seja, colnames(meth.g)[5] == colnames(expr.g)[5]


cpg.gene.normal <- merge(NonTumor.DNAmethylation,NonTumor.GeneExpression,by="GeneID")
dimnames(cpg.gene.normal)[[1]] <- paste(cpg.gene[,"ProbeID"], cpg.gene[,"GeneID"], sep=".")


#As linhas de cpg.gene e cpg.gene.normal precisam estar exatamente na mesma ordem. 
#identical(rownames(cpg.gene), rownames(cpg.gene.normal)) == TRUE

meth.g.norm <- cpg.gene.normal[,1:100]
names(meth.g.norm) <- substr(names(meth.g.norm), 1, 12)
expr.g.norm <- cpg.gene.normal[,101:200]
names(expr.g.norm) <- substr(names(expr.g.norm), 1, 12)
#As amostras precisam estar na mesma ordem para metilação do DNA e expressão gênica, ou seja, colnames(meth.g.norm)[5] = colnames(expr.g.norm)[5]

epiGene <- list()
j=1
metadata <- metadata[colnames(meth.g),] 
#Objeto de metadata (com informação de subgrupos, idade, etc.) precisa estar na mesma ordem que as amostras em meth.g ou expr

for(i in 1:dim(cpg.gene.normal)[1]) {
  
  meth.gu <- names(meth.g)[meth.g[i, ] < 0.3]
  meth.gm <- names(meth.g)[meth.g[i, ] >= 0.3]
  if(length(meth.gu)>1 & length(meth.gm)>1){
    expr.mean.gm <- mean(as.numeric(expr.g[i, meth.gm]), na.rm=T)
    if(expr.mean.gm < 1.28*sd(as.numeric(expr.g[i, meth.gu]), na.rm=T)){
      expr.mean.gu <- mean(as.numeric(expr.g[i, meth.gu]), na.rm = T)
      yy <- cbind(t(meth.g[i, ] >= 0.3), t(expr.g[i, ] < expr.mean.gu))
      colnames(yy) <- c("meth.gt0.3", "exp.lt.meanUnmeth")
      aux <- cbind(t(meth.g.norm[rownames(meth.g[i,]), ] >= 0.3), t(expr.g.norm[rownames(meth.g[i,]), ] < expr.mean.gu))
      colnames(aux) <- c("meth.gt0.3", "exp.lt.meanUnmeth")
      yy <- rbind(yy, aux)
      yy <- as.data.frame(yy)
      gene <- rownames(meth.g.norm[i,])
      if(as.matrix(table(yy$meth.gt0.3[1:ncol(meth.g)], yy$exp.lt.meanUnmeth[1:ncol(meth.g)]))[2,2]/as.numeric(table(yy$meth.gt0.3[1:nrow(meth.g)])[2])>0.8){
        yy$geneName <- gene
        yy$id <- rownames(yy)
        yy$methValues <- NA
        yy$exprValues <- NA
        yy$methValues[1:ncol(meth.g)] <- as.numeric(meth.g[i, ])
        yy$exprValues[1:ncol(meth.g)] <- as.numeric(expr.g[i, ])  
        yy$methValues[(ncol(meth.g)+1):nrow(yy)] <- as.vector(as.matrix(t(meth.g.norm[rownames(meth.g[i,]), ])))
        yy$exprValues[(ncol(meth.g)+1):nrow(yy)] <- as.vector(as.matrix(t(expr.g.norm[rownames(meth.g[i,]), ])))
        yy <- merge(yy, metadata, by=0, all.x=TRUE)
        epiGene.RNAseq[[j]] <- as.data.frame(yy)
        names(epiGene.RNAseq)[[j]] <- as.character(gene)
        j <- j+1
      }   
    }
  }
}  


uu.nl <- epiGene.RNAseq[[1]]
clusters <- as.character(metadata$subgroups)


enrichment.normal <- lapply(clusters, function(cluster, epiGene.ss=epiGene.RNAseq) {
  enrichment.group <- unlist(mclapply(epiGene.ss, function(uu.nl, cluster.group=cluster) {
    control <- (!uu.nl$meth.gt0.3[(ncol(meth.g)+1):length(uu.nl$meth.gt0.3)] & !uu.nl$exp.lt.meanUnmeth[(ncol(meth.g)+1):length(uu.nl$exp.lt.meanUnmeth)])
    uu.nl$es <- NA
    uu.nl$es[1:ncol(meth.g)] <- uu.nl$meth.gt0.3[1:ncol(meth.g)] & uu.nl$exp.lt.meanUnmeth[1:ncol(meth.g)]
    uu.nl$Kclass <- NA
    uu.nl$Kclass <- uu.nl$cluster == cluster.group
    es.table <- table(uu.nl$Kclass, uu.nl$es)
    if((es.table[2,2]/rowSums(es.table)[2] > 0.5) & (es.table[2,2]/colSums(es.table)[2] > 0.5)  & sum(control) > 0 &(table(control)[2]/(table(control)[2] + table(control)[1]) > 0.5 | is.na(table(control)[2]/(table(control)[2] + table(control)[1])))) {
      return(fisher.test(table(uu.nl$Kclass, uu.nl$es))$p.value)
    } else {
      return(NULL)
    }
  }, mc.cores=1))
  return(enrichment.group)
})

names(enrichment.normal) <- as.character(metadata$subgroups)
enrichment.normal.adj <- lapply(enrichment.normal, p.adjust, method="BH")
enrichmentList.sig <- lapply(enrichment.normal.adj, function(x) {
  x <- x[x < 0.05]
})

