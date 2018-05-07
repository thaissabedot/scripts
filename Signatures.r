load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG-GBM.merged.data.FB20150930v3.Rdata")
rownames(pd) <- as.character(pd$case.id)
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/EReg1_PATENT.Rda")
newEReg1 <- volcano
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/GCIMPhigh_codel.rda")
EReg5.p <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/EReg5.txt")

w.p.values.adj <- p.adjust(w.p.values, method = "BH")

volcano$meanM1 <- apply(volcano[,1:100],1,mean,na.rm=T)
volcano$meanM2 <- apply(volcano[,101:200],1,mean,na.rm=T)
volcano$DiffMean <- volcano$meanM1 - volcano$meanM2
volcano <- na.omit(volcano)
volcano$p.value.adj <- w.p.values.adj
volcano$p.value <- w.p.values

volcano$threshold <- "1"
a <- subset(volcano, p.value.adj < 0.05)
b <- subset(a, DiffMean < -0.5) #hyper no da direita
c <- subset(a, DiffMean < 0 & DiffMean > -0.5)
volcano[rownames(b),"threshold"] <- "2"
volcano[rownames(c),"threshold"] <- "4"
b <- subset(a, DiffMean > 0.5) #hyper no da esquerda
c <- subset(a, DiffMean > 0 & DiffMean < 0.5)
volcano[rownames(b),"threshold"] <- "3"
volcano[rownames(c),"threshold"] <- "4"

probes <- rownames(subset(volcano, threshold %in% c(2,3)))
probes <- c("cg16313343","cg23555120","cg10521852","cg15439862")
LGG.GBM.new <- LGG.GBM[as.character(EReg2),5:936]
colnames(LGG.GBM.new) <- substr(colnames(LGG.GBM.new),1,12)
#LGG.GBM.new <- LGG.GBM.new > 0.1


a <- pd[colnames(LGG.GBM.new),c("clustM.supervised2","histology","case.id")]

hypo <- subset(a, clustM.supervised2 %in% "G-CIMP-low")
hypo <- LGG.GBM.new[,as.character(hypo$case.id)]
hypo.s <- hclust(dist(t(hypo)))


noncod <- subset(a, clustM.supervised2 %in% "G-CIMP-high")
noncod <- LGG.GBM.new[,as.character(noncod$case.id)]
noncod.s <- hclust(dist(t(noncod)))

codel <- subset(a, clustM.supervised2 %in% "IDHmut-codel")
codel <- LGG.GBM.new[,as.character(codel$case.id)]
codel.s <- hclust(dist(t(codel)))

classic <- subset(a, clustM.supervised2 %in% "Classic-like")
classic <- LGG.GBM.new[,as.character(classic$case.id)]
classic.s <- hclust(dist(t(classic)))

mesen <- subset(a, clustM.supervised2 %in% "Mesenchymal-like")
mesen <- LGG.GBM.new[,as.character(mesen$case.id)]
mesen.s <- hclust(dist(t(mesen)))

lgg <- subset(a, clustM.supervised2 %in% "PA-like-LGG")
lgg <- LGG.GBM.new[,as.character(lgg$case.id)]
lgg.s <- hclust(dist(t(lgg)))

gbm <- subset(a, clustM.supervised2 %in% "PA-like-GBM")
gbm <- LGG.GBM.new[,as.character(gbm$case.id)]
gbm.s <- hclust(dist(t(gbm)))

col.order1 <- c(colnames(hypo[,hypo.s$order]),colnames(noncod[,noncod.s$order]), colnames(codel[,codel.s$order]))
col.order2 <- c(colnames(classic[,classic.s$order]),  colnames(mesen[,mesen.s$order]), colnames(gbm[,gbm.s$order]),colnames(lgg[,lgg.s$order]))

EReg1 <- LGG.GBM.new[rownames(newEReg1[newEReg1$threshold %in% c(2,3),]),]
EReg1.s <- hclust(dist(EReg1))

EReg2 <- rev(rownames(volcano[volcano$threshold %in% c(2,3),]))

EReg5 <- LGG.GBM.new[as.character(EReg5.p$id),]
EReg5.s <- hclust(dist(EReg5))

row.order1 <- c(EReg2)
row.order2 <- as.character(probes)


LGG.GBM.new2 <- LGG.GBM.new[rev(row.order1),as.character(col.order1)]
#a <- as.data.frame(lapply(LGG.GBM.new2,function(z){z <- as.character(z); z[is.na(z)] <- "NA"; return(z)}))
#indx <- sapply(a, is.factor)
#a[indx] <- lapply(a[indx], function(x) as.numeric(as.character(x)))
#rownames(a) <- rownames(LGG.GBM.new)
#a <- a > 0.1


a <- LGG.GBM.new2 > 0.5

clab <- pd[col.order1,c("case.id","clustM.supervised2","histology")]
#LGG.GBM.new.order2 <- LGG.GBM.new.order2[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
#a <- a[c(1:63,64:327,405:531,532:682,683:932,328:404),]
clab$clustM.supervised2 <- factor(clab$clustM.supervised2)
clab$histology <- factor(clab$histology)
#levels(clab$clustM.supervised2) <- c("orange","red","green3","purple","yellow","darkblue","lightblue")
#levels(clab$histology) <- c("red","purple","cyan","green")


teste <- as.data.frame(t(a))
teste$id <- rownames(teste)
rownames(teste) <- gsub("\\.", "\\-" ,rownames(teste))
teste1 <- merge(teste,clab,by.x=0,by.y="case.id")
rownames(teste1) <- as.character(teste1$Row.names)
teste1 <- teste1[,-1]
teste1 <- melt(teste1,id="id")
teste1$id <- factor(teste1$id)
teste1$id <- gsub("\\.", "\\-" ,teste1$id)
teste1$id <- factor(teste1$id, levels = col.order1)

ggplot(teste1, aes(y=factor(variable),x=id)) + 
  geom_tile(aes(fill=factor(value))) +
  #scale_fill_manual(values=c("red","orange","white","brown4","green3","purple","purple","yellow","blue","green","darkblue","lightblue","black"),guide = guide_legend(title = "Legend")) +
  #scale_fill_manual(values=c("red","orange","white","purple","yellow","blue","green","lightblue","darkblue","black"),guide = guide_legend(title = "Legend")) +
  scale_fill_manual(values=c("red","white","brown4","green3","purple","purple","blue","green","black"),guide = guide_legend(title = "Legend")) +  
  theme_bw() +
  xlab("TCGA LGG and GBM") + 
  ylab("Probe ID") +
  guides(col = guide_legend(nrow = 8))


load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/p131_roworder_Heatmap_G-CIMP-low.Rda")
probes <- c("cg19154438","cg20579480","cg11719784","cg22792910","cg18123677","cg13904968","cg22367264","cg09748960")
LGG.GBM.new <- LGG.GBM[probes,5:936]
colnames(LGG.GBM.new) <- substr(colnames(LGG.GBM.new),1,12)

a <- pd[colnames(LGG.GBM.new),c("clustM.supervised2","histology","case.id")]
hypo <- subset(a, clustM.supervised2 %in% "G-CIMP-low")
hypo <- LGG.GBM.new[,as.character(hypo$case.id)]
hypo.s <- hclust(dist(t(hypo)))


noncod <- subset(a, clustM.supervised2 %in% "G-CIMP-high")
noncod <- LGG.GBM.new[,as.character(noncod$case.id)]
noncod.s <- hclust(dist(t(noncod)))

order <- LGG.GBM.new[,c(colnames(hypo[,hypo.s$order]),colnames(noncod[,noncod.s$order]))]
heatmap.plus.sm(as.matrix(validation),
                col = jet.colors(75),
                scale = "none",
                #labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
		cexRow = 0.8
                #margins = c(1, 6)
                #RowSideColors = rlab
                
)
order <- LGG.GBM.new2
rownames(order) <- paste(rownames(order),c("GPR156","FGFR2"),sep=":")
rownames(order) <- paste(rownames(order),c("Island","Island"),sep="-")
a <- order > 0.45

clab <- pd[colnames(a),c("case.id","clustM.supervised2","histology")]
#LGG.GBM.new.order2 <- LGG.GBM.new.order2[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
#a <- a[c(1:63,64:327,405:531,532:682,683:932,328:404),]
clab$clustM.supervised2 <- factor(clab$clustM.supervised2)
clab$histology <- factor(clab$histology)
#levels(clab$clustM.supervised2) <- c("orange","red","green3","purple","yellow","darkblue","lightblue")
#levels(clab$histology) <- c("red","purple","cyan","green")


teste <- as.data.frame(t(a))
teste$id <- rownames(teste)
rownames(teste) <- gsub("\\.", "\\-" ,rownames(teste))
teste1 <- merge(teste,clab,by.x=0,by.y="case.id")
rownames(teste1) <- as.character(teste1$Row.names)
teste1 <- teste1[,-1]
teste1 <- melt(teste1,id="id")
teste1$id <- factor(teste1$id)
teste1$id <- gsub("\\.", "\\-" ,teste1$id)
teste1$id <- factor(teste1$id, levels = colnames(a))

ggplot(teste1, aes(y=factor(variable),x=id)) + 
  geom_tile(aes(fill=factor(value))) +
  #scale_fill_manual(values=c("red","orange","white","brown4","green3","purple","purple","yellow","blue","green","darkblue","lightblue","black"),guide = guide_legend(title = "Legend")) +
  #scale_fill_manual(values=c("red","orange","white","purple","yellow","blue","green","lightblue","darkblue","black"),guide = guide_legend(title = "Legend")) +
  scale_fill_manual(values=c("red","white","brown4","green3","purple","purple","blue","green","black"),guide = guide_legend(title = "Legend")) + 
  #scale_fill_manual(values=c("red","white","brown4","green3","purple","blue","green","black"),guide = guide_legend(title = "Legend")) + 
  theme_bw() +
  xlab("TCGA LGG and GBM") + 
  ylab("Probe ID") +
  guides(col = guide_legend(nrow = 8))


geneAnnot <- txs(TSS=list(upstream=0, downstream=0))
probe <- GRanges(seqnames = paste0("chr",as.character(LGG.GBM[EReg2,"Chromosome"])),
range=IRanges(start = as.numeric(LGG.GBM[EReg2,"Genomic_Coordinate"]), end= as.numeric(LGG.GBM[EReg2,"Genomic_Coordinate"])),
name= as.character(LGG.GBM[EReg2,"Composite.Element.REF"]))
NearbyGenes <- GetNearGenes(geneNum=1,geneAnnot=geneAnnot,TRange=probe)


######### Validation
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/all_validationDNAmeth.rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/clinical_validation3.rda")

probes <- rownames(subset(volcano, threshold %in% c(2,3)))
probes <- c("cg16313343","cg23555120","cg10521852","cg15439862")
#probes <- c("cg19154438","cg20579480","cg11719784","cg22792910","cg18123677","cg13904968","cg22367264","cg09748960")
validation <-all.validation[as.character(probes),]
#LGG.GBM.new <- LGG.GBM.new > 0.1


a <- ids[colnames(validation),c("clustM.supervised","original.cluster","Study","id")]

hypo <- subset(a, clustM.supervised %in% "G-CIMP-low")
hypo <- validation[,as.character(hypo$id)]
hypo.s <- hclust(dist(t(hypo)))


noncod <- subset(a, clustM.supervised %in% "G-CIMP-high")
noncod <- validation[,as.character(noncod$id)]
noncod.s <- hclust(dist(t(noncod)))

codel <- subset(a, clustM.supervised %in% "IDHmut-codel")
codel <- validation[,as.character(codel$id)]
codel.s <- hclust(dist(t(codel)))

classic <- subset(a, clustM.supervised %in% "Classic-like")
classic <- validation[,as.character(classic$id)]
classic.s <- hclust(dist(t(classic)))

mesen <- subset(a, clustM.supervised %in% "Mesenchymal-like I")
mesen <- validation[,as.character(mesen$id)]
mesen.s <- hclust(dist(t(mesen)))

lgg <- subset(a, clustM.supervised %in% "PA-like")
lgg <- validation[,as.character(lgg$id)]
lgg.s <- hclust(dist(t(lgg)))

gbm <- subset(a, clustM.supervised %in% "Mesenchymal-like II")
gbm <- validation[,as.character(gbm$id)]
gbm.s <- hclust(dist(t(gbm)))

col.order1 <- c(colnames(hypo[,hypo.s$order]),colnames(noncod[,noncod.s$order]))
#col.order1 <- c(colnames(hypo[,hypo.s$order]),colnames(noncod[,noncod.s$order]), colnames(codel[,codel.s$order]))
col.order2 <- c(colnames(classic[,classic.s$order]),  colnames(mesen[,mesen.s$order]), colnames(gbm[,gbm.s$order]),colnames(lgg[,lgg.s$order]))

EReg1 <- validation[rownames(newEReg1[newEReg1$threshold %in% c(2,3),]),]
EReg1.s <- hclust(dist(EReg1))

EReg2 <- rev(rownames(volcano[volcano$threshold %in% c(2,3),]))

EReg5 <- validation[as.character(EReg5.p$id),]
EReg5.s <- hclust(dist(EReg5))

row.order1 <- c(probes)
row.order2 <- as.character(probes)


validation2 <- validation[rev(row.order2),as.character(col.order2)]
#a <- as.data.frame(lapply(LGG.GBM.new2,function(z){z <- as.character(z); z[is.na(z)] <- "NA"; return(z)}))
#indx <- sapply(a, is.factor)
#a[indx] <- lapply(a[indx], function(x) as.numeric(as.character(x)))
#rownames(a) <- rownames(LGG.GBM.new)
#a <- a > 0.1


a <- validation2 > 0.45

clab <- ids[col.order2,c("clustM.supervised","original.cluster","Study","id")]
#LGG.GBM.new.order2 <- LGG.GBM.new.order2[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
#a <- a[c(1:63,64:327,405:531,532:682,683:932,328:404),]
clab$clustM.supervised <- factor(clab$clustM.supervised)
clab$original.cluster <- factor(clab$original.cluster)
clab$Study <- factor(clab$Study)


teste <- as.data.frame(t(a))
teste$id <- rownames(teste)
teste1 <- merge(teste,clab,by.x=0,by.y="id")
rownames(teste1) <- as.character(teste1$Row.names)
teste1 <- teste1[,-1]
teste1 <- melt(teste1,id="id")
teste1$id <- factor(teste1$id)
teste1$id <- factor(teste1$id, levels = col.order2)

#levels(validation.info$Sturm.PA.merged.clusters) <- c("cornflowerblue","chartreuse4","yellow","black","lightgreen","blue4","orchid","green3","gold3","darkorange2","firebrick4","red","lightblue4") #PA1 PA2 CD-CIMP+ CIMP+ CIMP- G34 IDH K27 Mesenc RTKI RTKII CIMP+ CIMP-



ggplot(teste1, aes(y=factor(variable),x=id)) + 
  geom_tile(aes(fill=factor(value))) +
  scale_fill_manual(values=c("cornflowerblue","chartreuse4","lightgreen","orange","white","blue4","green3","lightpink4","gold3","yellow","darkblue","magenta","lightblue","firebrick4","darkorange2","salmon","black","turquoise","lightblue4"),guide = guide_legend(title = "Legend")) +
  #scale_fill_manual(values=c("black","white","brown4","green3","blue4","magenta","darkorange2","salmon","black","turquoise","red"),guide = guide_legend(title = "Legend")) +
  #scale_fill_manual(values=c("yellow","black","white","brown4","green3","blue4","purple","gold3","magenta","darkorange2","salmon","black","turquoise","red"),guide = guide_legend(title = "Legend")) +  
  theme_bw() +
  xlab("Validation Set") + 
  ylab("Probe ID") +
  guides(col = guide_legend(nrow = 8))

