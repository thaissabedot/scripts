load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/EReg1.new3_12genes.Rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Control_EScalls/normal_DNAMeth.rda")
load("LGG-GBM.merged.data_20150714.Rdata")
b <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/probesLGm2.txt",header=F)
aux <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/EReg5.txt")
a <- c(as.character(aux$id),names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm3),names(sort(enrichmentList.sig$LGm2)[1:15]),as.character(EReg1$cgID))
probe.meth <- str_extract(a,"^*[^.]*")

LGG.GBM.new.noXY.s <- LGG.GBM.250914[ as.character(probe.meth),5:936]
LGG.GBM.new.order2 <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
colnames(LGG.GBM.new.order2) <- substr(colnames(LGG.GBM.new.order2),1,12)
LGG.GBM.new.order2 <- LGG.GBM.new.order2[,c(1:63,64:327,405:531,532:682,683:932,328:404)]

ESg1 <- LGG.GBM.new.order2[as.character(EReg1$cgID),]
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


rownames(pd) <- as.character(pd$case.id)
#order <- rbind(ESg5[ESg5.p$order,],ESg4[ESg4.p$order,],ESg3[ESg3.p$order,],ESg2[ESg2.p$order,],ESg1[ESg1.p$order,])
order <- ESg1[ESg1.p$order,]
d <- normal.dnameth
rownames(d) <- as.character(d$Composite.El)
d <- d[rownames(order),5:114]


clab <- pd[colnames(order),]
levels(clab$TERT.expression) <- c("purple","black")
clab$TERT.expression <- as.character(clab$TERT.expression)
clab[is.na(clab$TERT.expression),"TERT.expression"] <- "white"
rownames(clab) <- as.character(clab$case.id)
clab <- clab[colnames(order),]
levels(clab$cluster.meth) <- c("green","red","purple","orange","yellow","blue")
clab$cluster.meth <- as.character(clab$cluster.meth)
clab[is.na(clab$cluster.meth),"cluster.meth"] <- "white"
clab$histology <- as.factor(clab$histology)
levels(clab$histology) <- c("red","purple","cyan","green")
clab$histology <- as.character(clab$histology)
clab[is.na(clab$histology),"histology"] <- "white"
clab <- clab[,c("histology","TERT.expression","cluster.meth")]
norm.label <- matrix("white", nrow = 110, ncol = 3,dimnames = list(c(1:110), c("histology","TERT.expression","cluster.meth")))

clab <- rbind(norm.label,clab)
order <- cbind(d,order)

png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_newEReg1_v2.png",res=1200,width=8000,height=8000)
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_newTracks.pdf")
heatmap.plus.sm(as.matrix(order),
                col = jet.colors(75),
                scale = "none",
                #labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = as.matrix(clab),
                margins = c(1, 6)
                #RowSideColors = rlab
                
)
dev.off()

########## VALIDATION

validation.set.p <- validation.set[rownames(order),] #136 DKFZ, 61 PA, 46 Mur, 81 Turcan

png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_Validation_newEReg1.png",res=1200,width=8000,height=8000)
heatmap.plus.sm(as.matrix(validation.set.p[,rownames(validation.info)]),
                col = jet.colors(75),
                scale = "none",
                #labRow = NA,
                #cexRow = 0.4,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = cbind(as.character(validation.info[,"cluster.meth.all"]),as.character(validation.info[,"cluster.meth.all"]))
                
)
dev.off()


c <- cpg.gene
rownames(c) <- as.character(c$Composite.El)
LGG.GBM.new.order <- c[rownames(order),644:1279] #tem que ter 12 probes so
aux <- as.character(c[rownames(order),"Associated.Gene.Name"])

colnames(LGG.GBM.new.order) <- substr(colnames(LGG.GBM.new.order),1,12)
a <- metadata.1129.samples.20150312[metadata.1129.samples.20150312$id %in% substr(colnames(LGG.GBM.new.order),1,12),]
lgm1 <- subset(a, cluster.meth == "LGm1")
lgm1.s <- LGG.GBM.new.order[,as.character(lgm1$id)]
mean <- apply(lgm1.s, 1, mean, na.rm=TRUE)
lgm2 <- subset(a, cluster.meth == "LGm2")
lgm2.s <- LGG.GBM.new.order[,as.character(lgm2$id)]
mean <- cbind(mean,apply(lgm2.s, 1, mean, na.rm=TRUE))
lgm3 <- subset(a, cluster.meth == "LGm3")
lgm3.s <- LGG.GBM.new.order[,as.character(lgm3$id)]
mean <- cbind(mean,apply(lgm3.s, 1, mean, na.rm=TRUE))
lgm4 <- subset(a, cluster.meth == "LGm4")
lgm4.s <- LGG.GBM.new.order[,as.character(lgm4$id)]
mean <- cbind(mean,apply(lgm4.s, 1, mean, na.rm=TRUE))
lgm5 <- subset(a, cluster.meth == "LGm5")
lgm5.s <- LGG.GBM.new.order[,as.character(lgm5$id)]
mean <- cbind(mean,apply(lgm5.s, 1, mean, na.rm=TRUE))
lgm6 <- subset(a, cluster.meth == "LGm6")
lgm6.s <- LGG.GBM.new.order[,as.character(lgm6$id)]
mean <- cbind(mean,apply(lgm6.s, 1, mean, na.rm=TRUE))
mean <- as.data.frame(mean)
colnames(mean) <- c("LGm1","LGm2","LGm3","LGm4","LGm5","LGm6")

clab.2 <- c("green","red","purple","orange","yellow","blue")
#png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_MeanRNA3.png",res=1200,width=8000,height=8000)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_MeanRNA3.pdf")
heatmap.plus.sm(as.matrix(mean),
                col=rev(redgreen(1000)),
                #scale = "none",
                #trace = "none",
                labRow = as.character(aux),
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                #ColSideColors = clab.6
                ColSideColors = cbind(clab.2,clab.2)
                #RowSideColors = rlab
                
)
dev.off()

