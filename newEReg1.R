load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/EReg1_new_12genes.Rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/tudao.27.450.rda")
load("./LGG.GBM/LGG_GBM_obj_20150127.Rda")

rownames(pd) <- as.character(pd$case.id)

aux <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/EReg5.txt")
a <- c(as.character(aux$id),names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm3),names(sort(enrichmentList.sig$LGm2)[1:15]),as.character(EReg1_new$Row.names))
#a <- as.character(EReg1_new$Row.names)
probe.meth <- str_extract(a,"^*[^.]*")

rownames(LGG.GBM.new) <- as.character(LGG.GBM.new$Composite.Element.REF)
LGG.GBM.new <- LGG.GBM.new[,colnames(dat.lgg.gbm.new.noXY.dic.oecg)]
LGG.GBM.new.noXY.s <- LGG.GBM.new[ as.character(probe.meth),]
LGG.GBM.new.order2 <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
colnames(LGG.GBM.new.order2) <- substr(colnames(LGG.GBM.new.order2),1,12)
#a <- colnames(LGG.GBM.new.order) %in% colnames(GBM.LGG.27.450k.nearest.gene)
#teste <- LGG.GBM.new.order[,a]
#LGG.GBM.new.order <- subset(LGG.GBM.new.order,select=colnames(GBM.LGG.27.450k.nearest.gene)[8:643])
a <- pd[substr(colnames(LGG.GBM.new.order2),1,12),c("cluster.meth","histology")]
LGG.GBM.new.order2 <- LGG.GBM.new.order2[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
a <- a[c(1:63,64:327,405:531,532:682,683:932,328:404),]

levels(a$cluster.meth) <- c("green","red","purple","orange","yellow","blue")
levels(a$histology) <- c("red","purple","cyan","green")
rownames(normal.dnameth) <- as.character(normal.dnameth$Composite.Element.REF)
normals.sel <- normal.dnameth[rownames(LGG.GBM.new.order2),5:114]
norm.label <- matrix("white", nrow = 110, ncol = 2)
clab <- rbind(norm.label,setNames(as.matrix(a),names(norm.label)))

LGG.GBM.new.order2 <- cbind(normals.sel,LGG.GBM.new.order2)


ESg1 <- LGG.GBM.new.order2[as.character(EReg1_new$Row.names),]
ESg1.p <- hclust(dist(ESg1))
probe.meth <- str_extract(names(sort(enrichmentList.sig$LGm2)[1:15]),"^*[^.]*")
ESg2 <- LGG.GBM.new.order2[as.character(probe.meth),]
ESg2.p <- hclust(dist(ESg2))
probe.meth <- str_extract(names(enrichmentList.sig$LGm3),"^*[^.]*")
ESg3 <- LGG.GBM.new.order2[as.character(probe.meth),]
ESg3.p <- hclust(dist(ESg3))
probe.meth <- str_extract(c(names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4)),"^*[^.]*")
ESg4 <- LGG.GBM.new.order2[as.character(probe.meth),]
ESg4.p <- hclust(dist(ESg4))
probe.meth <- str_extract(as.character(aux$id),"^*[^.]*")
ESg5 <- LGG.GBM.new.order2[as.character(probe.meth),]
ESg5.p <- hclust(dist(ESg5))

order <- rbind(ESg5[ESg5.p$order,],ESg4[ESg4.p$order,],ESg3[ESg3.p$order,],ESg2[ESg2.p$order,],ESg1[ESg1.p$order,])

png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_DNAMeth_new.png",res=1200,width=8000,height=8000)
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_DNAMeth_tert.pdf")
heatmap.plus.sm(as.matrix(order),
                col = jet.colors(75),
                scale = "none",
                #labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab,
		cexRow = 0.8
                #margins = c(1, 6)
                #RowSideColors = rlab
                
)
dev.off()

load("./LGG.GBM/RNAseq_05-11-14/rnaSeqNorm.Rda")


aux <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/EReg5.txt")
a <- c(as.character(aux$id),names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm3),names(sort(enrichmentList.sig$LGm2)[1:15]))
#a <- as.character(EReg1_new$Row.names)
probe.meth <- str_extract(a,"^*[^.]*")
c <- cpg.gene
rownames(c) <- as.character(c$Composite.El)
LGG.GBM.new.order <- c[as.character(probe.meth),644:1279] #57 genes
colnames(LGG.GBM.new.order) <- substr(colnames(LGG.GBM.new.order),1,12)
aux <- rna.seq.norm[as.character(EReg1_new$gene),] 
aux <- aux[,colnames(LGG.GBM.new.order)]
LGG.GBM.new.order <- rbind(LGG.GBM.new.order,aux)

a <- pd[pd$case.id %in% colnames(LGG.GBM.new.order),]
lgm1 <- subset(a, cluster.meth == "LGm1")
lgm1.s <- LGG.GBM.new.order[,as.character(lgm1$case.id)]
mean <- apply(lgm1.s, 1, mean, na.rm=TRUE)

lgm2 <- subset(a, cluster.meth == "LGm2")
lgm2.s <- LGG.GBM.new.order[,as.character(lgm2$case.id)]
mean <- cbind(mean,apply(lgm2.s, 1, mean, na.rm=TRUE))

lgm3 <- subset(a, cluster.meth == "LGm3")
lgm3.s <- LGG.GBM.new.order[,as.character(lgm3$case.id)]
mean <- cbind(mean,apply(lgm3.s, 1, mean, na.rm=TRUE))

lgm4 <- subset(a, cluster.meth == "LGm4")
lgm4.s <- LGG.GBM.new.order[,as.character(lgm4$case.id)]
mean <- cbind(mean,apply(lgm4.s, 1, mean, na.rm=TRUE))

lgm5 <- subset(a, cluster.meth == "LGm5")
lgm5.s <- LGG.GBM.new.order[,as.character(lgm5$case.id)]
mean <- cbind(mean,apply(lgm5.s, 1, mean, na.rm=TRUE))

lgm6 <- subset(a, cluster.meth == "LGm6")
lgm6.s <- LGG.GBM.new.order[,as.character(lgm6$case.id)]
mean <- cbind(mean,apply(lgm6.s, 1, mean, na.rm=TRUE))
mean <- as.data.frame(mean)
colnames(mean) <- c("LGm1-hypo","LGm1","LGm2","LGm3","LGm4","LGm5","LGm6")
