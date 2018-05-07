load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG-GBM.merged.data.FB20150924.Rdata")
levels(pd$clustM.supervised) <- c("IDHmut-codel","IDHmut-non-codel","IDHmut-hypo","IDHmut-non-codel","IDHwt-Mesenchymal","IDHwt-Classic","IDHwt-hypo")
pd$clustM.supervised <- as.character(pd$clustM.supervised)
pd[pd$clustM.supervised %in% "IDHwt-hypo" & pd$tumor.type == "GBM", "clustM.supervised"] <- "IDHwt-hypo-GBM"
pd[pd$clustM.supervised %in% "IDHwt-hypo" & pd$tumor.type == "LGG", "clustM.supervised"] <- "IDHwt-hypo-LGG"

meth.g <- cpg.gene[,8:643] #only tumor
names(meth.g) <- substr(names(meth.g), 1, 12)
expr.g <- cpg.gene[,644:1279] #only tumor
names(expr.g) <- substr(names(expr.g), 1, 12)

aux <- subset(pd, clustM.supervised %in% c("IDHwt-Mesenchymal","IDHwt-Classic","IDHwt-hypo-GBM","IDHwt-hypo-LGG"))
#aux <- subset(pd, clustM.supervised %in% c("IDHmut-codel","IDHmut-hypo","IDHmut-non-codel"))

meth.g <- meth.g[,colnames(meth.g) %in% as.character(aux$case.id)] #204 samples IDHwt - 424 samples IDHmut
expr.g <- expr.g[,colnames(meth.g)]


meth.g.norm <- cpg.gene.normal[,8:117]
names(meth.g.norm) <- substr(names(meth.g.norm), 1, 12)
expr.g.norm <- cpg.gene.normal[,118:227]
names(expr.g.norm) <- substr(names(expr.g.norm), 1, 12)
epiGene.RNAseq <- list()
j=1
#metadata <- subset(metadata.1129.samples.20141030,Study == "TCGA")
a <- subset(pd, clustM.supervised %in% c("IDHwt-Mesenchymal","IDHwt-Classic","IDHwt-hypo-GBM","IDHwt-hypo-LGG"))
rownames(a) <- as.character(a$case.id)

for(i in 1:dim(cpg.gene.normal)[1]) {
 #for(i in 138:250) {
  
  meth.gu <- names(meth.g)[meth.g[i, ] < 0.3]
  meth.gm <- names(meth.g)[meth.g[i, ] >= 0.3]
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
      if(as.matrix(table(yy$meth.gt0.3[1:204], yy$exp.lt.meanUnmeth[1:204]))[2,2]/as.numeric(table(yy$meth.gt0.3[1:204])[2])>0.8){
        yy$geneName <- gene
        yy$id <- rownames(yy)
        yy$methValues <- NA
        yy$exprValues <- NA
        yy$methValues[1:204] <- as.numeric(meth.g[i, ])
        yy$exprValues[1:204] <- as.numeric(expr.g[i, ])  
        #ic <-248+(491*(k-1))
        #fim <- 247+(491*k)
        yy$methValues[205:nrow(yy)] <- as.vector(as.matrix(t(meth.g.norm[rownames(meth.g[i,]), ])))
        yy$exprValues[205:nrow(yy)] <- as.vector(as.matrix(t(expr.g.norm[rownames(meth.g[i,]), ])))
        #yy$methValues[248:nrow(yy)] <- as.numeric(meth.g.norm[index[k], ])
        #yy$exprValues[248:nrow(yy)] <- as.numeric(expr.g.norm[index[k], ])
        info <- a[rownames(yy)[1:204],c("clustM.supervised","IDH.status","Codel.1p.19q","clustM.GBM","clustM.LGG","tumor.type")]
        yy$IDH.status <- "WT"
        yy$codel1p19q <- "non-codel"
        yy$clustM.supervised <- NA
        yy$clustM.GBM <- NA
        yy$clustM.LGG <- NA
        yy$tumor.type <- NA
        yy$IDH.status[1:204] <- as.character(info$IDH.status)
        yy$clustM.GBM[1:204] <- as.character(info$clustM.GBM)
        yy$clustM.LGG[1:204] <- as.character(info$clustM.LGG)
        yy$tumor.type[1:204] <- as.character(info$tumor.type)
        yy$codel1p19q[1:204] <- as.character(info$Codel.1p.19q)
        yy$clustM.supervised[1:204] <- as.character(info$clustM.supervised)
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
clusters <- c("IDHwt-Mesenchymal","IDHwt-Classic","IDHwt-hypo-GBM","IDHwt-hypo-LGG")
#clusters <- "LGm1"


enrichment.normal <- lapply(clusters, function(cluster, epiGene.ss=epiGene.RNAseq) {
  enrichment.group <- unlist(mclapply(epiGene.ss, function(uu.nl, cluster.group=cluster) {
    control <- (!uu.nl$meth.gt0.3[205:length(uu.nl$meth.gt0.3)] & !uu.nl$exp.lt.meanUnmeth[205:length(uu.nl$exp.lt.meanUnmeth)])
    uu.nl$es <- NA
    uu.nl$es[1:204] <- uu.nl$meth.gt0.3[1:204] & uu.nl$exp.lt.meanUnmeth[1:204]
    uu.nl$Kclass <- NA
    uu.nl$Kclass <- uu.nl$clustM.supervised == cluster.group
    es.table <- table(uu.nl$Kclass, uu.nl$es)#, dnn=c("Class", "ES"))
    #print()
    if((es.table[2,2]/rowSums(es.table)[2] > 0.5) & (es.table[2,2]/colSums(es.table)[2] > 0.5)  & sum(control) > 0 &(table(control)[2]/(table(control)[2] + table(control)[1]) > 0.5 | is.na(table(control)[2]/(table(control)[2] + table(control)[1])))) {
      return(fisher.test(table(uu.nl$Kclass, uu.nl$es))$p.value)
    } else {
      return(NULL)
    }
    #}
    #}  
  }, mc.cores=7))
  return(enrichment.group)
})

names(enrichment.normal) <- c("IDHwt-Mesenchymal","IDHwt-Classic","IDHwt-hypo-GBM","IDHwt-hypo-LGG")
enrichment.normal.adj <- lapply(enrichment.normal, p.adjust, method="BH")
enrichmentList.IDHwt2 <- lapply(enrichment.normal.adj, function(x) {
  x <- x[x < 0.05]
})


#a <- c(names(sort(enrichmentList.IDHwt2$`IDHwt-Mesenchymal`)[1:10]),names(sort(enrichmentList.IDHwt2$`IDHwt-Classic`)[1:10]),names(sort(enrichmentList.IDHmut2$`IDHmut-codel`)[1:10]),names(sort(enrichmentList.IDHmut2$`IDHmut-non-codel`)[1:10]))
a <- c(names(sort(enrichmentList.IDHwt2$`IDHwt-Classic`)[1:10]),names(sort(enrichmentList.IDHmut2$`IDHmut-codel`)[1:10]),names(sort(enrichmentList.IDHmut2$`IDHmut-non-codel`)[1:10]))
probe.meth <- str_extract(a,"^*[^.]*")

rownames(LGG.GBM.new) <- as.character(LGG.GBM.new$Composite.Element.REF)
LGG.GBM.new <- LGG.GBM.new[,colnames(dat.lgg.gbm.new.noXY.dic.oecg)]
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/order_Heatmap_IDHmut_K3.Rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/order_Heatmap_IDHwt_K3.Rda")

rownames(pd) <- as.character(pd$case.id)
LGG.GBM.new.noXY.s <- LGG.GBM.new[ as.character(probe.meth),]
colnames(LGG.GBM.new.noXY.s) <- substr(colnames(LGG.GBM.new.noXY.s),1,12)
LGG.GBM.new.order2 <- LGG.GBM.new.noXY.s[,c(order_Heatmap_IDHmut_K3,order_Heatmap_IDHwt_K3)]


a <- pd[colnames(LGG.GBM.new.order2),c("clustM.supervised","histology","case.id")]

hypo <- subset(a, clustM.supervised %in% "IDHmut-hypo")
hypo <- LGG.GBM.new.order2[,as.character(hypo$case.id)]
hypo.s <- hclust(dist(t(hypo)))

codel <- subset(a, clustM.supervised %in% "IDHmut-codel")

noncod <- subset(a, clustM.supervised %in% "IDHmut-non-codel")
noncod <- LGG.GBM.new.order2[,as.character(noncod$case.id)]
noncod.s <- hclust(dist(t(noncod)))

classic <- subset(a, clustM.supervised %in% "IDHwt-Classic")

mesen <- subset(a, clustM.supervised %in% "IDHwt-Mesenchymal")

lgg <- subset(a, clustM.supervised %in% "IDHwt-hypo-LGG")
lgg <- LGG.GBM.new.order2[,as.character(lgg$case.id)]
lgg.s <- hclust(dist(t(lgg)))

gbm <- subset(a, clustM.supervised %in% "IDHwt-hypo-GBM")
gbm <- LGG.GBM.new.order2[,as.character(gbm$case.id)]
gbm.s <- hclust(dist(t(gbm)))

col.order <- c(colnames(hypo[,hypo.s$order]),colnames(noncod[,noncod.s$order]), as.character(codel$case.id), as.character(classic$case.id), as.character(mesen$case.id), colnames(gbm[,gbm.s$order]),colnames(lgg[,lgg.s$order]))


ESg1 <- LGG.GBM.new.order2[str_extract(names(sort(enrichmentList.IDHmut2$`IDHmut-non-codel`)[1:10]),"^*[^.]*"),]
ESg1.p <- hclust(dist(ESg1))
probe.meth <- str_extract(names(sort(enrichmentList.IDHmut2$`IDHmut-codel`)[1:10]),"^*[^.]*")
ESg2 <- LGG.GBM.new.order2[as.character(probe.meth),]
ESg2.p <- hclust(dist(ESg2))
probe.meth <- str_extract(names(sort(enrichmentList.IDHwt2$`IDHwt-Classic`)[1:10]),"^*[^.]*")
ESg3 <- LGG.GBM.new.order2[as.character(probe.meth),]
ESg3.p <- hclust(dist(ESg3))
#probe.meth <- str_extract(names(sort(enrichmentList.IDHwt2$`IDHwt-Mesenchymal`)[1:10]),"^*[^.]*")
#ESg4 <- LGG.GBM.new.order2[as.character(probe.meth),]
#ESg4.p <- hclust(dist(ESg4))



order <- rbind(ESg3[ESg3.p$order,],ESg2[ESg2.p$order,],ESg1[ESg1.p$order,])

order <- order[,col.order]

#a <- colnames(LGG.GBM.new.order) %in% colnames(GBM.LGG.27.450k.nearest.gene)
#teste <- LGG.GBM.new.order[,a]
#LGG.GBM.new.order <- subset(LGG.GBM.new.order,select=colnames(GBM.LGG.27.450k.nearest.gene)[8:643])
a <- pd[colnames(order),c("clustM.supervised","histology")]
#LGG.GBM.new.order2 <- LGG.GBM.new.order2[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
#a <- a[c(1:63,64:327,405:531,532:682,683:932,328:404),]
a$clustM.supervised <- as.factor(a$clustM.supervised)
levels(a$clustM.supervised) <- c("purple","darkgreen","red","orange","blue","lightblue","yellow")
levels(a$histology) <- c("red","purple","cyan","green")
rownames(normal.dnameth) <- as.character(normal.dnameth$Composite.Element.REF)
normals.sel <- meth.g.norm
rownames(normals.sel) <- as.character(cpg.gene.normal$Composite.El)
normals.sel <- normals.sel[rownames(order),]
norm.label <- matrix("white", nrow = 110, ncol = 2)
clab <- rbind(norm.label,setNames(as.matrix(a),names(norm.label)))
order <- cbind(normals.sel,order)

#rlab <- c(rep("lightpink",10),rep("hotpink1",10),rep("hotpink3",10),rep("hotpink4",10))

rlab <- c(rep("hotpink1",10),rep("hotpink3",10),rep("hotpink4",10))
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_DNAMeth_IDHKn4.png",res=1200,width=8000,height=8000)
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_DNAMeth_tert.pdf")
heatmap.plus.sm(as.matrix(order),
                col = jet.colors(75),
                scale = "none",
                #labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab,
		RowSideColors = cbind(rlab,rlab),
		cexRow = 0.8
                #margins = c(1, 6)
                #RowSideColors = rlab
                
)
dev.off()




###### IDHwt Wilcoxon Tests
### hypo: LGG x GBM (450k)
load("col.order")
LGG.GBM.new <- LGG.GBM.250914[,5:936]
colnames(LGG.GBM.new) <- substr(colnames(LGG.GBM.new),1,12)
LGG.GBM.new.order <- LGG.GBM.new[,col.order]
a <- colnames(LGG.GBM.new.order)[colnames(LGG.GBM.new.order) %in% substr(colnames(all.450k),1,12)]
colnames(all.450k) <- substr(colnames(all.450k),1,12)
LGG.GBM.new.order <- all.450k[,a]

clab <- c(rep("orange",73),rep("yellow",101),rep("blue",12),rep("lightblue",26))
#clab <- c(rep("blue",12),rep("lightblue",26))
clab <- cbind(clab,clab)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/heatmap_2.png",res=300,width=2000,height=2000)

heatmap.plus(as.matrix(LGG.GBM.new.order[rownames(subset(IDHwt.LGGxGBM, p.value.adj < 0.05 & abs(DiffMean) > 0.25)),426:637]),
                col = jet.colors(75),
			#trace = "none",
                scale = "none",
                #labRow = NA,
                labCol = NA,
                Colv = NA,
               # Rowv = NA,
                ColSideColors = as.matrix(clab),
		#RowSideColors = cbind(as.character(rlab$cor),as.character(rlab$cor)),
                margins = c(1, 6)
                #RowSideColors = rlab
                
)
dev.off()


### Mesenchymal x Classic
LGG.GBM.new <- LGG.GBM[,5:936]
colnames(LGG.GBM.new) <- substr(colnames(LGG.GBM.new),1,12)
LGG.GBM.new <- LGG.GBM.new[,col.order]
LGG.GBM.new <- LGG.GBM.new[c(str_extract(names(sort(enrichmentList.IDHwt2$`IDHwt-Classic`))[1:10],"^*[^.]*"),rownames(subset(IDHwt.ClaxMes, p.value.adj < 0.05 & (DiffMean < -0.1 | DiffMean > 0.3))),rownames(subset(IDHwt.GBMxCl.Ms.lgg, p.value.adj < 0.05 & abs(DiffMean) > 0.3))),449:878]

clab <- c(rep("orange",148),rep("yellow",215),rep("blue",41),rep("lightblue",26))
#clab <- c(rep("blue",12),rep("lightblue",26))
clab <- cbind(clab,clab)

rlab <- c(rep("violet",30),rep("pink",16))

#png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/heatmap_IDHwt-GBMxClass_e_Mesen.png",res=300,width=2000,height=2000)

heatmap.plus(as.matrix(LGG.GBM.new),
                col = jet.colors(75),
                scale = "none",
                #labRow = NA,
                labCol = NA,
                Colv = NA,
               # Rowv = NA,
                ColSideColors = as.matrix(clab),
		RowSideColors = cbind(rlab,rlab),
                margins = c(1, 6)
                #RowSideColors = rlab
                
)
dev.off()

