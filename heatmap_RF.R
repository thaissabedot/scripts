colnames(LGG.GBM.149) <- substr(colnames(LGG.GBM.149),1,12)
rownames(pd) <- as.character(pd$case.id)
aux <- subset(pd,clustM.supervised2 %in% c("G-CIMP-low","G-CIMP-high","IDHmut-codel"))
a <- as.character(subset(aux, clustM.supervised2 %in% c("G-CIMP-low"))$case.id)
a1 <- LGG.GBM.149[,a]
a1.s <- hclust(dist(t(a1)))
a <- as.character(subset(aux, clustM.supervised2 %in% c("G-CIMP-high"))$case.id)
a2 <- LGG.GBM.149[,a]
a2.s <- hclust(dist(t(a2)))
a <- as.character(subset(aux, clustM.supervised2 %in% c("IDHmut-codel"))$case.id)
a3 <- LGG.GBM.149[,a]
a3.s <- hclust(dist(t(a3)))

order <- cbind(a1[,a1.s$order],a2[,a2.s$order],a3[,a3.s$order])


clab <- aux[colnames(order),c("clustM.supervised2","histology")]
clab$histology <- as.factor(clab$histology)
clab$clustM.supervised2 <- factor(clab$clustM.supervised2)
levels(clab$histology) <- c("red","purple","cyan","green")
levels(clab$clustM.supervised2) <- c("firebrick","darkgreen","purple")
clab <- as.matrix(clab)

a1 <- order[rownames(order) %in% GCIMPlow.probes,]
a1.s <- hclust(dist(a1))
a2 <- order[rownames(order) %in% Noncodels.probes,]
a2.s <- hclust(dist(a2))

order <- rbind(a2[a2.s$order,],a1[a1.s$order,])

rlab <- c(rep("brown",nrow(a2)),rep("pink",nrow(a1)))

#png("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmap_TCGA_GCIMPlow_LGm.png",res=500,width=2000,height=2000)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmap_TCGA_GCIMPlow_LGmLABELS.pdf")
heatmap.plus.sm(as.matrix(order[1:2,]),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab,
		#RowSideColors = cbind(rlab,rlab),	
                margins = c(1,6)
                #RowSideColors = rlab
                
)
dev.off()

a <- as.character(subset(ids, clustM.supervised %in% c("G-CIMP-low"))$id)
a1 <- all.149[,a]
a1.s <- hclust(dist(t(a1)))
a <- as.character(subset(ids, clustM.supervised %in% c("G-CIMP-high"))$id)
a2 <- all.149[,a]
a2.s <- hclust(dist(t(a2)))
a <- as.character(subset(ids, clustM.supervised %in% c("IDHmut-codel"))$id)
a3 <- all.149[,a]
a3.s <- hclust(dist(t(a3)))

order2 <- cbind(a1[,a1.s$order],a2[,a2.s$order],a3[,a3.s$order])


clab <- ids[colnames(order2),c("clustM.supervised","original.cluster","Study","Histology")]
clab$original.cluster <- as.factor(clab$original.cluster)
clab$clustM.supervised <- factor(clab$clustM.supervised)
clab$Study <- as.factor(clab$Study)
clab$Histology <- factor(clab$Histology)
levels(clab$clustM.supervised) <- c("purple","darkgreen","firebrick")
levels(clab$original.cluster) <- c("cornflowerblue","chartreuse4","yellow","lightgreen","black","blue4","orchid","green3","gold3","firebrick4","darkorange2","lightblue4","red") 
levels(clab$Study) <- c("lightpink4","magenta","salmon","tomato4")
levels(clab$Histology) <- c("red","purple","cyan","green","darkblue")
clab <- as.matrix(clab)
order2 <- order2[rownames(order),]

#png("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmap_Validation_GCIMPlow_LGmLABELS.png",res=500,width=2000,height=2000)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmap_Validation_GCIMPlow_LGmLABELS.pdf")
heatmap.plus.sm(as.matrix(order2[1:2,]),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab,
		#RowSideColors = cbind(rlab,rlab),	
                margins = c(1,6)
                #RowSideColors = rlab
                
)
dev.off()

