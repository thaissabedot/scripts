######## Proteomics
load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/LGG-GBM.merged.data.FB20150930v3.Rdata") #pdata from Sep 2015 (LGG-GBM project). Created by Floris Barthel
rownames(pd) <- as.character(pd$case.id)

#TCGA PanCan
setwd("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/Proteomics")
Prot.PanCan <- read.csv("PanCan/TCGA-PANCAN19-L4.csv")
#LGG 257
#GBM 205
Prot.PanCan <- Prot.PanCan[Prot.PanCan$Cancer_Type %in% c("LGG","GBM"),]
rownames(Prot.PanCan) <- substr(as.character(Prot.PanCan$Sample_ID),1,12)
Prot.PanCan <- Prot.PanCan[,-c(1:4)]
Prot.PanCan <- merge(Prot.PanCan,pd[,c("case.id","cartoon")],by.x=0,by.y="case.id")
Prot.PanCan.GCIMP <- rbind(Prot.PanCan[Prot.PanCan$clustM.supervised2 %in% c("G-CIMP-low"),],Prot.PanCan[Prot.PanCan$clustM.supervised2 %in% c("G-CIMP-high"),])
Prot.PanCan.GCIMP$clustM.supervised2 <- factor(Prot.PanCan.GCIMP$clustM.supervised2)
rownames(Prot.PanCan.GCIMP) <- as.character(Prot.PanCan.GCIMP$Row.names)
Prot.PanCan.GCIMP <- Prot.PanCan.GCIMP[,-1]
#GCIMP-low: 8
#GCIMP-high: 128
#Proteinas: 258

p.value <- NA
for(i in 1:(ncol(Prot.PanCan.GCIMP)-1)){
  aux <- Prot.PanCan.GCIMP[,c(colnames(Prot.PanCan.GCIMP)[i],"clustM.supervised2")]
  colnames(aux)[1] <- "protein"
  if(nrow(na.omit(aux)) > 0)
    p.value <- c(p.value,pvalue(wilcox_test(protein ~ clustM.supervised2,aux,exact=TRUE)))
  else
    p.value <- c(p.value,NA)
  
}
p.value <- p.value[-1] #57
p.value.adj <- p.adjust(p.value,method = "fdr") #2, linhas: 42 (CYCLINB1) e 68 (IGFBP2)

ggplot(Prot.PanCan.GCIMP, aes(clustM.supervised2, CYCLINB1)) +
  geom_boxplot(aes(fill = clustM.supervised2)) +
  geom_jitter() +  scale_fill_manual(values = c("firebrick","darkgreen"))

### LGG
Prot.LGG <- read.csv("LGG/TCGA-LGG-L4.csv")
#LGG 427
rownames(Prot.LGG) <- substr(as.character(Prot.LGG$Sample_ID),1,12)
Prot.LGG <- Prot.LGG[,-c(1:4)]
Prot.LGG <- merge(Prot.LGG,pd[,c("case.id","cartoon")],by.x=0,by.y="case.id")
Prot.LGG.CIMP <- rbind(Prot.LGG[Prot.LGG$cartoon %in% c("GCIMP-low"),],Prot.LGG[Prot.LGG$cartoon %in% c("GCIMP-high"),])
Prot.LGG.CIMP$cartoon <- factor(Prot.LGG.CIMP$cartoon)
rownames(Prot.LGG.CIMP) <- as.character(Prot.LGG.CIMP$Row.names)
Prot.LGG.CIMP <- Prot.LGG.CIMP[,-1]
#GCIMP-low: 10
#GCIMP-high: 191
#Proteinas: 219

p.value <- NA
for(i in 1:(ncol(Prot.LGG.CIMP)-1)){
  aux <- Prot.LGG.CIMP[,c(colnames(Prot.LGG.CIMP)[i],"clustM.supervised2")]
  colnames(aux)[1] <- "protein"
  if(nrow(na.omit(aux)) > 0)
    p.value <- c(p.value,pvalue(wilcox_test(protein ~ clustM.supervised2,aux,exact=TRUE)))
  else
    p.value <- c(p.value,NA)
  
}
p.value <- p.value[-1] #110
sum(p.value < 0.05,na.rm=TRUE)
p.value.adj <- p.adjust(p.value,method = "fdr") #80, linhas: 42 (CYCLINB1) e 68 (IGFBP2)
sum(p.value.adj < 0.05,na.rm=TRUE)

Prot.LGG.CIMP.diff <- Prot.LGG.CIMP[,p.value.adj < 0.05]
a <- Prot.LGG.CIMP.diff[Prot.LGG.CIMP.diff$clustM.supervised2 %in% c("G-CIMP-low"),]
aux <- hclust(dist(a[,-c(81)]))
ordem <- a[aux$order,]
clab <- rep("darkgreen",nrow(a))
a <- Prot.LGG.CIMP.diff[Prot.LGG.CIMP.diff$clustM.supervised2 %in% c("G-CIMP-high"),]
aux <- hclust(dist(a[,-c(81)]))
ordem <- rbind(ordem,a[aux$order,])
clab <- c(clab,rep("firebrick",nrow(a)))

aux <- volcano[volcano$threshold %in% c(2), 1:201]
a <- ordem[rownames(ordem) %in% rownames(aux),]
aux <- hclust(dist(a))
rows <- rownames(a[aux$order,])
rlab <- rep("darkblue",nrow(a))
aux <- volcano[volcano$threshold %in% c(3), 1:201]
a <- ordem[rownames(ordem) %in% rownames(aux),]
aux <- hclust(dist(a))
rows <- c(rows,rownames(a[aux$order,]))
rlab <- c(rlab,rep("darkblue",nrow(a)))

#ha1 = HeatmapAnnotation(df = data.frame(Subtype=Prot.LGG.CIMP.diff$clustM.supervised2), 
#                        col = list(Subtype = c("G-CIMP-low" =  "darkgreen", "G-CIMP-high" = "firebrick")))

#Heatmap(t(ordem[,-c(81)]), 
#        name = "Protein values", 
#        col=greenred(75),
#        cluster_rows = TRUE, 
#        cluster_columns = FALSE,
#        top_annotation = ha1,
#        show_column_names = FALSE,
#        row_names_gp = gpar(fontsize = 5))

heatmap.2(t(ordem[,-c(81)]),
             Colv=FALSE,
             col=greenred(75),
             ColSideColors = clab,
             dendrogram = "row",
             #labRow = NA,
             labCol = NA,
             scale = "row",
             trace = "none")

volcano <- as.data.frame(t(Prot.LGG.CIMP[,-c(220)]))
volcano$meanM1 <- apply(volcano[,1:10],1,mean,na.rm=T)
volcano$meanM2 <- apply(volcano[,11:201],1,mean,na.rm=T)
volcano$DiffMean <- volcano$meanM1 - volcano$meanM2
volcano$p.value.adj <- p.value.adj
volcano$p.value <- p.value

volcano$threshold <- "1"
volcano[volcano$p.value.adj < 0.05 & volcano$DiffMean < -0.25, "threshold"] <- "2" #hypo GCIMP-low
volcano[volcano$p.value.adj < 0.05 & volcano$DiffMean > 0.5, "threshold"] <- "3" #hyper GCIMP-low
volcano$protein_name <- rownames(volcano)

ggplot(data=volcano,aes(x=DiffMean, y=-1*log10(p.value.adj),colour=threshold)) +
  geom_point() +
  xlab("Protein Expression (fold-change)") + ylab("-1 * log10 of the Significance") +
  labs(title = "Volcano Plot") +
  scale_color_manual(breaks=c("1","2","3"), # color scale (for points)
                     values=c("grey", "green", "red"),
                     labels=c("Not Significant","Downregulated in GCIMP-low","Upregulated in GCIMP-low"),name="Legend")  +
  geom_label_repel(data=subset(volcano,threshold %in% c("2","3")),aes(label=protein_name),size = 5,box.padding = unit(1, "lines"),point.padding = unit(1, "lines"))

volcano2 <- volcano[volcano$threshold %in% c(2,3), 1:201]
aux <- as.character(subset(pd, clustM.supervised2 %in% c("G-CIMP-low"))$case.id)
a <- volcano2[,colnames(volcano2) %in% aux]
aux <- hclust(dist(t(a)))
ordem <- a[,aux$order]
clab <- rep("darkgreen",ncol(a))
aux <- as.character(subset(pd, clustM.supervised2 %in% c("G-CIMP-high"))$case.id)
a <- volcano2[,colnames(volcano2) %in% aux]
aux <- hclust(dist(t(a)))
ordem <- cbind(ordem,a[,aux$order])
clab <- c(clab,rep("firebrick",ncol(a)))

aux <- volcano[volcano$threshold %in% c(2), 1:201]
a <- ordem[rownames(ordem) %in% rownames(aux),]
aux <- hclust(dist(a))
rows <- rownames(a[aux$order,])
rlab <- rep("darkblue",nrow(a))
aux <- volcano[volcano$threshold %in% c(3), 1:201]
a <- ordem[rownames(ordem) %in% rownames(aux),]
aux <- hclust(dist(a))
rows <- c(rows,rownames(a[aux$order,]))
rlab <- c(rlab,rep("darkblue",nrow(a)))

heatmap.2(as.matrix(ordem[rows,]),
          Colv=FALSE,
          Rowv=FALSE,
          col=greenred(75),
          ColSideColors = clab,
          RowSideColors = rlab,
          dendrogram = "none",
          #labRow = NA,
          labCol = NA,
          scale = "row",
          trace = "none")

### GBM
Prot.GBM <- read.csv("GBM/TCGA-GBM-L4.csv")
#GBM 205
rownames(Prot.GBM) <- substr(as.character(Prot.GBM$Sample_ID),1,12)
Prot.GBM <- Prot.GBM[,-c(1:4)]
Prot.GBM <- merge(Prot.GBM,pd[,c("case.id","clustM.supervised2")],by.x=0,by.y="case.id")
Prot.GBM.CIMP <- rbind(Prot.GBM[Prot.GBM$clustM.supervised2 %in% c("G-CIMP-low"),],Prot.GBM[Prot.GBM$clustM.supervised2 %in% c("G-CIMP-high"),])
Prot.GBM.CIMP$clustM.supervised2 <- factor(Prot.GBM.CIMP$clustM.supervised2)
rownames(Prot.GBM.CIMP) <- as.character(Prot.GBM.CIMP$Row.names)
Prot.GBM.CIMP <- Prot.GBM.CIMP[,-1]
#GCIMP-low: 4
#GCIMP-high: 4
#Proteinas: 223

p.value <- NA
for(i in 1:(ncol(Prot.GBM.CIMP)-1)){
  aux <- Prot.GBM.CIMP[,c(colnames(Prot.GBM.CIMP)[i],"clustM.supervised2")]
  colnames(aux)[1] <- "protein"
  if(nrow(na.omit(aux)) > 0)
    p.value <- c(p.value,pvalue(wilcox_test(protein ~ clustM.supervised2,aux,exact=TRUE)))
  else
    p.value <- c(p.value,NA)
  
}
p.value <- p.value[-1] #110
sum(p.value < 0.05,na.rm=TRUE)
p.value.adj <- p.adjust(p.value,method = "fdr") #80, linhas: 42 (CYCLINB1) e 68 (IGFBP2)
sum(p.value.adj < 0.05,na.rm=TRUE)

Prot.GBM.CIMP.diff <- Prot.GBM.CIMP[,p.value < 0.05]
Prot.GBM.CIMP.diff$clustM.supervised2 <- Prot.GBM.CIMP$clustM.supervised2
a <- Prot.GBM.CIMP.diff[Prot.GBM.CIMP.diff$clustM.supervised2 %in% c("G-CIMP-low"),]
aux <- hclust(dist(a[,-15]))
ordem <- a[aux$order,]
clab <- rep("darkgreen",nrow(a))
a <- Prot.GBM.CIMP.diff[Prot.GBM.CIMP.diff$clustM.supervised2 %in% c("G-CIMP-high"),]
aux <- hclust(dist(a[,-15]))
ordem <- rbind(ordem,a[aux$order,])
clab <- c(clab,rep("firebrick",nrow(a)))


volcano <- as.data.frame(t(Prot.GBM.CIMP[,-c(224)]))
volcano$meanM1 <- apply(volcano[,1:4],1,mean,na.rm=T)
volcano$meanM2 <- apply(volcano[,5:8],1,mean,na.rm=T)
volcano$DiffMean <- volcano$meanM1 - volcano$meanM2
volcano$p.value.adj <- p.value.adj
volcano$p.value <- p.value

volcano$threshold <- "1"
volcano[volcano$p.value.adj < 0.05 & volcano$DiffMean < -0.25, "threshold"] <- "2" #hypo GCIMP-low
volcano[volcano$p.value.adj < 0.05 & volcano$DiffMean > 0.5, "threshold"] <- "3" #hyper GCIMP-low
volcano$protein_name <- rownames(volcano)

ggplot(data=volcano,aes(x=DiffMean, y=-1*log10(p.value.adj),colour=threshold)) +
  geom_point() +
  xlab("Protein Expression (fold-change)") + ylab("-1 * log10 of the Significance") +
  labs(title = "Volcano Plot") +
  scale_color_manual(breaks=c("1","2","3"), # color scale (for points)
                     values=c("grey", "green", "red"),
                     labels=c("Not Significant","Downregulated in GCIMP-low","Upregulated in GCIMP-low"),name="Legend")  +
  geom_label_repel(data=subset(volcano,threshold %in% c("2","3")),aes(label=protein_name),size = 5,box.padding = unit(1, "lines"),point.padding = unit(1, "lines"))

volcano2 <- volcano[volcano$threshold %in% c(2,3), 1:201]
aux <- as.character(subset(pd, clustM.supervised2 %in% c("G-CIMP-low"))$case.id)
a <- volcano2[,colnames(volcano2) %in% aux]
aux <- hclust(dist(t(a)))
ordem <- a[,aux$order]
clab <- rep("darkgreen",ncol(a))
aux <- as.character(subset(pd, clustM.supervised2 %in% c("G-CIMP-high"))$case.id)
a <- volcano2[,colnames(volcano2) %in% aux]
aux <- hclust(dist(t(a)))
ordem <- cbind(ordem,a[,aux$order])
clab <- c(clab,rep("firebrick",ncol(a)))

aux <- volcano[volcano$threshold %in% c(2), 1:201]
a <- ordem[rownames(ordem) %in% rownames(aux),]
aux <- hclust(dist(a))
rows <- rownames(a[aux$order,])
rlab <- rep("darkblue",nrow(a))
aux <- volcano[volcano$threshold %in% c(3), 1:201]
a <- ordem[rownames(ordem) %in% rownames(aux),]
aux <- hclust(dist(a))
rows <- c(rows,rownames(a[aux$order,]))
rlab <- c(rlab,rep("darkblue",nrow(a)))

heatmap.2(as.matrix(ordem[rows,]),
          Colv=FALSE,
          Rowv=FALSE,
          col=greenred(75),
          ColSideColors = clab,
          RowSideColors = rlab,
          dendrogram = "none",
          #labRow = NA,
          labCol = NA,
          scale = "row",
          trace = "none")