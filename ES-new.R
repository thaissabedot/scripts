load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG-GBM.merged.data.FB20150918.Rdata")
pdata <- subset(pd,!is.na(cluster.meth2))


library(GenomicRanges)
load("LGG_GBM_obj_20150127.Rda")

gene.GR <- GRanges(seqnames = paste0("chr",gene.location$Chromosome.Name),ranges = IRanges(start = gene.location$Gene.Start..bp., end=gene.location$Gene.End..bp.), strand = gene.location$Strand, symbol = gene.location$Associated.Gene.Name, EntrezID = gene.location$EntrezGene.ID)
probe.info <- GRanges(seqnames = paste0("chr",LGG.GBM.new$Chromosome), ranges = IRanges(start = LGG.GBM.new$Genomic_Coordinate, end = LGG.GBM.new$Genomic_Coordinate), probeID = LGG.GBM.new$Composite.Element.REF)
distance <- as.data.frame(distanceToNearest(probe.info,gene.GR)) #closest gene to each 450k probe
rownames(gene.location) <- NULL
gene.order.by.distance <- gene.location[distance$subjectHits,]
gene.order.by.distance$distance <- as.matrix(distance$distance)
LGG.GBM <- LGG.GBM.new[,colnames(dat.lgg.gbm.new.noXY.dic.oecg)]
rownames(LGG.GBM) <- as.character(LGG.GBM.new$Composite.Element.REF)
GBM.LGG.27.450k.nearest.gene <- cbind(LGG.GBM.new[,1:4],LGG.GBM)
rownames(GBM.LGG.27.450k.nearest.gene) <- as.character(GBM.LGG.27.450k.nearest.gene$Composite.Element.REF)
GBM.LGG.27.450k.nearest.gene <- cbind(GBM.LGG.27.450k.nearest.gene, gene.order.by.distance[,c("Associated.Gene.Name","EntrezGene.ID","distance")])
info <- GBM.LGG.27.450k.nearest.gene[,c(1:4,937:939)]
colnames(GBM.LGG.27.450k.nearest.gene) <- substr(colnames(GBM.LGG.27.450k.nearest.gene),1,12)
GBM.LGG.27.450k.nearest.gene <- GBM.LGG.27.450k.nearest.gene[,colnames(GBM.LGG.27.450k.nearest.gene) %in% as.character(pdata$case.id)]
GBM.LGG.27.450k.nearest.gene <- GBM.LGG.27.450k.nearest.gene[,colnames(GBM.LGG.27.450k.nearest.gene) %in% substr(colnames(rna.seq.norm)[1:667],1,12)]

GBM.LGG.27.450k.nearest.gene <- cbind(info,GBM.LGG.27.450k.nearest.gene)


cpg.gene <- merge(GBM.LGG.27.450k.nearest.gene,rna.seq.norm[,colnames(GBM.LGG.27.450k.nearest.gene)[8:636]],by.x="Associated.Gene.Name",by.y=0) 
#23774 probes - 1258 samples - 7 metadata columns
dimnames(cpg.gene)[[1]] <- paste(cpg.gene[,"Composite.Element.REF"], cpg.gene[,"Associated.Gene.Name"], sep=".")


meth.g <- cpg.gene[,8:636] #only tumor
names(meth.g) <- substr(names(meth.g), 1, 12)
expr.g <- cpg.gene[,637:1265] #only tumor
names(expr.g) <- substr(names(expr.g), 1, 12)
#identical(colnames(meth.g),colnames(expr.g))


normalTCGA.mrna <- rna.seq.norm[,colnames(rna.seq.norm) %in% colnames(normal.mrna)]
normalTCGA.mrna <- normalTCGA.mrna[,-111]
rownames(normal.dnameth) <- as.character(normal.dnameth$Composite.El)
#identical(rownames(normal.dnameth),rownames(GBM.LGG.27.450k.nearest.gene)
normal.dnameth <- cbind(gene.order.by.distance[,c("Associated.Gene.Name","EntrezGene.ID","distance")],normal.dnameth)
cpg.gene.normal <- merge(normal.dnameth,normalTCGA.mrna[,colnames(normal.dnameth)[8:117]],by.x="Associated.Gene.Name",by.y=0) 
#23774 probes - 1258 samples - 7 metadata columns
dimnames(cpg.gene.normal)[[1]] <- paste(cpg.gene.normal[,"Composite.El"], cpg.gene.normal[,"Associated.Gene.Name"], sep=".")


meth.g.norm <- cpg.gene.normal[,8:117]
names(meth.g.norm) <- substr(names(meth.g.norm), 1, 12)
expr.g.norm <- cpg.gene.normal[,118:227]
names(expr.g.norm) <- substr(names(expr.g.norm), 1, 12)

epiGene.RNAseq <- list()
j=1
a <- subset(pd,case.id %in% colnames(meth.g))
rownames(a) <- as.character(a$case.id)

for(i in 1:dim(cpg.gene.normal)[1]) {
  #for(i in 1:20) {
  
  meth.gu <- names(meth.g)[meth.g[i, ] < 0.3]
  if(sum(is.na(meth.gu)) == length(meth.g))
    meth.gu <- NA
  meth.gm <- names(meth.g)[meth.g[i, ] >= 0.3]
  if(sum(is.na(meth.gm)) == length(meth.g))
    meth.gm <- NA
  meth.gu <- na.omit(meth.gu)
  meth.gm <- na.omit(meth.gm)
  # if(fisher.test(table(yy))$p.value < 0.05){
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
      if(as.matrix(table(yy$meth.gt0.3[1:629], yy$exp.lt.meanUnmeth[1:629]))[2,2]/as.numeric(table(yy$meth.gt0.3[1:629])[2])>0.8){
        yy$geneName <- gene
        yy$id <- rownames(yy)
        yy$methValues <- NA
        yy$exprValues <- NA
        yy$methValues[1:629] <- as.numeric(meth.g[i, ])
        yy$exprValues[1:629] <- as.numeric(expr.g[i, ])  
        #ic <-248+(491*(k-1))
        #fim <- 247+(491*k)
        yy$methValues[630:nrow(yy)] <- as.vector(as.matrix(t(meth.g.norm[rownames(meth.g[i,]), ])))
        yy$exprValues[630:nrow(yy)] <- as.vector(as.matrix(t(expr.g.norm[rownames(meth.g[i,]), ])))
        #yy$methValues[248:nrow(yy)] <- as.numeric(meth.g.norm[index[k], ])
        #yy$exprValues[248:nrow(yy)] <- as.numeric(expr.g.norm[index[k], ])
        info <- a[rownames(yy)[1:629],c("cluster.meth","IDH.status","Codel.1p.19q","cluster.meth.gbm","cluster.meth.lgg","tumor.type","cluster.meth2")]
        yy$IDH.status <- "WT"
        yy$codel1p19q <- "non-codel"
        yy$cluster <- NA
        yy$cluster.gbm <- NA
        yy$cluster.lgg <- NA
        yy$tumor.type <- NA
        yy$cluster.meth2 <- NA
        yy$IDH.status[1:629] <- as.character(info$IDH.status)
        yy$cluster.gbm[1:629] <- as.character(info$cluster.meth.gbm)
        yy$cluster.lgg[1:629] <- as.character(info$cluster.meth.lgg)
        yy$tumor.type[1:629] <- as.character(info$tumor.type)
        yy$codel1p19q[1:629] <- as.character(info$Codel.1p.19q)
        yy$cluster[1:629] <- as.character(info$cluster.meth)
        yy$cluster.meth2[1:629] <- as.character(info$cluster.meth2)
        #yy <- merge(metadata, yy, by.x = "TCGAIDnew", by.y = "id")
        epiGene.RNAseq[[j]] <- as.data.frame(yy)
        names(epiGene.RNAseq)[[j]] <- as.character(gene)
        j <- j+1
        cat("found a ES", j, "\n")
      }
      else{
        cat("no step3", "\n")
      }
      
    }else{
      cat("no ES group", "\n")
    }
  }else{
    cat("no methylated group", i, "\n")
  }
}  


uu.nl <- epiGene.RNAseq[[1]]
#clusters <- c(paste0("LGm", 1:6))
clusters <- c("IDHmut-K1","IDHmut-K2","IDHmut-K3","IDHwt-K1","IDHwt-K2","IDHwt-K3")
#clusters <- "LGm1"


enrichment.normal <- lapply(clusters, function(cluster, epiGene.ss=epiGene.RNAseq) {
  enrichment.group <- unlist(mclapply(epiGene.ss, function(uu.nl, cluster.group=cluster) {
    control <- (!uu.nl$meth.gt0.3[630:length(uu.nl$meth.gt0.3)] & !uu.nl$exp.lt.meanUnmeth[630:length(uu.nl$exp.lt.meanUnmeth)])
    uu.nl$es <- NA
    uu.nl$es[1:629] <- uu.nl$meth.gt0.3[1:629] & uu.nl$exp.lt.meanUnmeth[1:629]
    uu.nl$Kclass <- NA
    uu.nl$Kclass <- uu.nl$cluster.meth2 == cluster.group
    es.table <- table(uu.nl$Kclass, uu.nl$es)#, dnn=c("Class", "ES"))
    #print()
    if((es.table[2,2]/rowSums(es.table)[2] > 0.5) & (es.table[2,2]/colSums(es.table)[2] > 0.5)  & (table(control)[2]/(table(control)[2] + table(control)[1]) > 0.5 | is.na(table(control)[2]/(table(control)[2] + table(control)[1])))) {
      return(fisher.test(table(uu.nl$Kclass, uu.nl$es))$p.value)
    } else {
      return(NULL)
    }
    #}
    #}  
  }, mc.cores=1))
  return(enrichment.group)
})

names(enrichment.normal) <- c("IDHmut-K1","IDHmut-K2","IDHmut-K3","IDHwt-K1","IDHwt-K2","IDHwt-K3")
enrichment.normal.adj <- lapply(enrichment.normal, p.adjust, method="BH")
enrichmentList.sig <- lapply(enrichment.normal.adj, function(x) {
  x <- x[x < 0.05]
})

aux <- c(names(enrichmentList.sig$`IDHmut-K2`),names(enrichmentList.sig$`IDHmut-K3`),names(enrichmentList.sig$`IDHwt-K2`))
probe.meth <- str_extract(aux,"^*[^.]*")
c <- cpg.gene
d <- cpg.gene.normal
rownames(c) <- as.character(c$Composite.Element.REF)
rownames(d) <- as.character(d$Composite.El)
LGG.GBM.new.order <- c[as.character(probe.meth),8:636]
aux <- as.character(c[as.character(probe.meth),"Associated.Gene.Name"])
#rownames(LGG.GBM.new.order) <- as.character(aux)
normals.sel <- d[as.character(probe.meth),8:117]
#rownames(normals.sel) <- as.character(aux)
a <- subset(pd, !is.na(cluster.meth2))

ES.probes <- LGG.GBM[as.character(probe.meth),]
colnames(ES.probes) <- substr(colnames(ES.probes),1,12)
ES.probes <- ES.probes[,as.character(a$case.id)]
IDHmutK1 <- as.character(subset(a,cluster.meth2 %in% "IDHmut-K1")$case.id)
IDHmutK1 <- ES.probes[,IDHmutK1]
IDHmut.K1.s <- hclust(dist(t(IDHmutK1)))
IDHmutK2 <- as.character(subset(a,cluster.meth2 %in% "IDHmut-K2")$case.id)
IDHmutK2 <- ES.probes[,IDHmutK2]
IDHmut.K2.s <- hclust(dist(t(IDHmutK2)))
IDHmutK3 <- as.character(subset(a,cluster.meth2 %in% "IDHmut-K3")$case.id)
IDHmutK3 <- ES.probes[,IDHmutK3]
IDHmut.K3.s <- hclust(dist(t(IDHmutK3)))
IDHwtK1 <- as.character(subset(a,cluster.meth2 %in% "IDHwt-K1")$case.id)
IDHwtK1 <- ES.probes[,IDHwtK1]
IDHwt.K1.s <- hclust(dist(t(IDHwtK1)))
IDHwtK2 <- as.character(subset(a,cluster.meth2 %in% "IDHwt-K2")$case.id)
IDHwtK2 <- ES.probes[,IDHwtK2]
IDHwt.K2.s <- hclust(dist(t(IDHwtK2)))
IDHwtK3 <- as.character(subset(a,cluster.meth2 %in% "IDHwt-K3")$case.id)
IDHwtK3 <- ES.probes[,IDHwtK3]
IDHwt.K3.s <- hclust(dist(t(IDHwtK3)))


#order <- rbind(ESg1[ESg1.p$order,],ESg2[ESg2.p$order,],ESg3[ESg3.p$order,],ESg4[ESg4.p$order,],ESg5[ESg5.p$order,])
#a <- metadata.1129.samples.20150312[colnames(order),]
order.s <- cbind(IDHmutK1[,IDHmut.K1.s$order],IDHmutK3[,IDHmut.K3.s$order],IDHmutK2[,IDHmut.K2.s$order],IDHwtK2[,IDHwt.K2.s$order],IDHwtK1[,IDHwt.K1.s$order],IDHwtK3[,IDHwt.K3.s$order])

probe.meth <- str_extract(names(sort(enrichmentList.sig$`IDHmut-K2`)),"^*[^.]*")
IDHmutK2 <- ES.probes[as.character(probe.meth),]
IDHmutK2.p <- hclust(dist(IDHmutK2))
probe.meth <- str_extract(names(sort(enrichmentList.sig$`IDHmut-K3`)),"^*[^.]*")
IDHmutK3 <- ES.probes[as.character(probe.meth),]
IDHmutK3.p <- hclust(dist(IDHmutK3))
probe.meth <- str_extract(names(sort(enrichmentList.sig$`IDHwt-K2`)),"^*[^.]*")
IDHwtK2 <- ES.probes[as.character(probe.meth),]
IDHwtK2.p <- hclust(dist(IDHwtK2))

order.p <- rbind(IDHwtK2[IDHwtK2.p$order,],IDHmutK2[IDHmutK2.p$order,],IDHmutK3[IDHmutK3.p$order,])

LGG.GBM.normal <- cbind(normal.dnameth[rownames(order.p),8:117],order.s[rownames(order.p),])
a$cluster.meth <- as.factor(a$cluster.meth)
a$cluster.meth2 <- as.factor(a$cluster.meth2)
levels(a$cluster.meth) <- c("green","red","purple","orange","yellow","blue")
levels(a$cluster.meth2) <- c("green","purple","red","yellow","orange","blue")
rownames(a) <- as.character(a$case.id)
a <- a[colnames(order.s),]
norm.label <- data.frame(cluster.meth=rep("white",110), cluster.meth2=rep("white",110))

clab <- rbind(norm.label,a[,c("cluster.meth","cluster.meth2")])
rlab <- c(rep("pink",nrow(IDHwtK2)),rep("cyan",nrow(IDHmutK2)),rep("violet",nrow(IDHmutK3)))

pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_RNAseqBoth2.pdf")
heatmap.plus.sm(as.matrix(LGG.GBM.normal),
                col=jet.colors(75),
                scale = "none",
                #trace = "none",
                #labRow = as.character(aux),
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = as.matrix(clab),
                RowSideColors = cbind(rlab,rlab)
                
)
dev.off()

c <- cpg.gene
rownames(c) <- as.character(c$Composite.Element.REF)
c <- c[rownames(LGG.GBM.normal),637:1265]
genes <- as.character(cpg.gene[rownames(LGG.GBM.normal),"Associated.Gene.Name"])
colnames(c) <- substr(colnames(c),1,12)
a <- subset(pd, !is.na(cluster.meth2))
IDHmutK1 <- as.character(subset(a,cluster.meth2 %in% "IDHmut-K1")$case.id)
IDHmutK1 <- c[,colnames(c) %in% IDHmutK1]
IDHmutK1.mean <- apply(IDHmutK1,1,mean,na.rm=TRUE)
IDHmutK2 <- as.character(subset(a,cluster.meth2 %in% "IDHmut-K2")$case.id)
IDHmutK2 <- c[,colnames(c) %in% IDHmutK2]
IDHmutK2.mean <- apply(IDHmutK2,1,mean,na.rm=TRUE)
IDHmutK3 <- as.character(subset(a,cluster.meth2 %in% "IDHmut-K3")$case.id)
IDHmutK3 <- c[,colnames(c) %in% IDHmutK3]
IDHmutK3.mean <- apply(IDHmutK3,1,mean,na.rm=TRUE)
IDHwtK1 <- as.character(subset(a,cluster.meth2 %in% "IDHwt-K1")$case.id)
IDHwtK1 <- c[,colnames(c) %in% IDHwtK1]
IDHwtK1.mean <- apply(IDHwtK1,1,mean,na.rm=TRUE)
IDHwtK2 <- as.character(subset(a,cluster.meth2 %in% "IDHwt-K2")$case.id)
IDHwtK2 <- c[,colnames(c) %in% IDHwtK2]
IDHwtK2.mean <- apply(IDHwtK2,1,mean,na.rm=TRUE)
IDHwtK3 <- as.character(subset(a,cluster.meth2 %in% "IDHwt-K3")$case.id)
IDHwtK3 <- c[,colnames(c) %in% IDHwtK3]
IDHwtK3.mean <- apply(IDHwtK3,1,mean,na.rm=TRUE)



mean <- cbind(IDHmutK1.mean,IDHmutK2.mean,IDHmutK3.mean,IDHwtK1.mean,IDHwtK2.mean,IDHwtK3.mean)
mean <- log2(mean)

pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_RNAseqBoth2.pdf")
heatmap.plus.sm(as.matrix(mean),
                col=rev(redgreen(1000)),
                scale = "row",
                #trace = "none",
                labRow = as.character(genes),
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = cbind(c("green","red","purple","orange","yellow","blue"),c("green","red","purple","orange","yellow","blue"))
                #RowSideColors = rlab
                
)
dev.off()
