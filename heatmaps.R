################# IDHwt
aux <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/EReg5.txt")
a <- c(as.character(aux$id),names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4))
probe.meth <- str_extract(a,"^*[^.]*")


LGG.GBM.new.noXY.s <- LGG.GBM.250914[ as.character(probe.meth),5:936]
LGG.GBM.new.order2 <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
colnames(LGG.GBM.new.order2) <- substr(colnames(LGG.GBM.new.order2),1,12)
LGG.GBM.new.order2 <- LGG.GBM.new.order2[,c(1:63,64:327,405:531,532:682,683:932,328:404)]

probe.meth <- str_extract(as.character(c(names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4))),"^*[^.]*")
ESg4 <- LGG.GBM.new.order2[as.character(probe.meth),]
ESg4.p <- hclust(dist(ESg4))
probe.meth <- str_extract(as.character(aux$id),"^*[^.]*")
ESg5 <- LGG.GBM.new.order2[as.character(probe.meth),]
ESg5.p <- hclust(dist(ESg5))

order <- rbind(ESg5[ESg5.p$order,],ESg4[ESg4.p$order,])

d <- normal.dnameth
rownames(d) <- as.character(d$Composite.El)
d <- d[rownames(order),5:114]

clab <- pd[colnames(order),c("TERT.expression","histology","clustM.supervised2","clustM.olddiscovery")]
clab$histology <- as.factor(clab$histology)
levels(clab$histology) <- c("red","purple","cyan","green")
clab$TERT.expression <- as.factor(clab$TERT.expression)
levels(clab$TERT.expression) <- c("black","purple")
clab$clustM.olddiscovery <- as.factor(clab$clustM.olddiscovery)
levels(clab$clustM.olddiscovery) <- c("green","red","purple","orange","yellow","blue")
clab$clustM.supervised2 <- as.factor(clab$clustM.supervised2)
levels(clab$clustM.supervised2) <- c("orange","#cf0000","#006400","purple","yellow","cornflowerblue","lightblue")
norm.label <- matrix("white", nrow = 110, ncol = 4,dimnames=list(NULL,colnames(clab)))
order <- cbind(d,order)
clab <- rbind(norm.label,clab)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_DNAMeth_Cartoon_LABELSwt.pdf")
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_DNAMeth_new.pdf")
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


c <- cpg.gene
rownames(c) <- as.character(c$Composite.El)
LGG.GBM.new.order <- c[rownames(order),644:1279]
aux <- c[rownames(order),c("Associated.Gene.Name")]
#b <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/probesLGm2.txt",header=F)
#a <- c(rownames(aux),names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm3),names(sort(enrichmentList.sig$LGm2)[1:15]),as.character(b$V1))
#probe.meth <- str_extract(a,"*[^.]*")
colnames(LGG.GBM.new.order) <- substr(colnames(LGG.GBM.new.order),1,12)
a <- pd[pd$case.id %in% colnames(LGG.GBM.new.order),]
lgm1 <- subset(a, clustM.olddiscovery %in% "LGm1")
lgm1.s <- LGG.GBM.new.order[,as.character(lgm1$case.id)]
mean <- apply(lgm1.s, 1, mean, na.rm=TRUE)
lgm2 <- subset(a, clustM.olddiscovery %in% "LGm2")
lgm2.s <- LGG.GBM.new.order[,as.character(lgm2$case.id)]
mean <- cbind(mean,apply(lgm2.s, 1, mean, na.rm=TRUE))
lgm3 <- subset(a, clustM.olddiscovery %in% "LGm3")
lgm3.s <- LGG.GBM.new.order[,as.character(lgm3$case.id)]
mean <- cbind(mean,apply(lgm3.s, 1, mean, na.rm=TRUE))
lgm4 <- subset(a, clustM.olddiscovery %in% "LGm4")
lgm4.s <- LGG.GBM.new.order[,as.character(lgm4$case.id)]
mean <- cbind(mean,apply(lgm4.s, 1, mean, na.rm=TRUE))
lgm5 <- subset(a, clustM.olddiscovery %in% "LGm5")
lgm5.s <- LGG.GBM.new.order[,as.character(lgm5$case.id)]
mean <- cbind(mean,apply(lgm5.s, 1, mean, na.rm=TRUE))
lgm6 <- subset(a, clustM.olddiscovery %in% "LGm6")
lgm6.s <- LGG.GBM.new.order[,as.character(lgm6$case.id)]
mean <- cbind(mean,apply(lgm6.s, 1, mean, na.rm=TRUE))
mean <- as.data.frame(mean)
colnames(mean) <- c("LGm1","LGm2","LGm3","LGm4","LGm5","LGm6")
clab <-c("green","red","purple","orange","yellow","blue")
mean <- mean + 1
mean <- log2(mean)
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
                ColSideColors = cbind(clab,clab)
                #RowSideColors = rlab
                
)
dev.off()




###Random Forest (mur,Sturm,PA,turcan)
validation.set.p <- all[rownames(order),]

aux <- subset(ids, LGm == "LGm1")
lgm1.p <- validation.set.p[,as.character(aux$id)]
hc.lgm1.p <- hclust(dist(t(lgm1.p)))
aux <- subset(ids, LGm == "LGm2") 
lgm2.p <- validation.set.p[,as.character(aux$id)]
hc.lgm2.p <- hclust(dist(t(lgm2.p)))
aux <- subset(ids, LGm == "LGm3") 
lgm3.p <- validation.set.p[,as.character(aux$id)]
hc.lgm3.p <- hclust(dist(t(lgm3.p)))
aux <- subset(ids, LGm == "LGm4")
lgm4.p <- validation.set.p[,as.character(aux$id)]
hc.lgm4.p <- hclust(dist(t(lgm4.p)))
aux <- subset(ids, LGm == "LGm5")
lgm5.p <- validation.set.p[,as.character(aux$id)]
hc.lgm5.p <- hclust(dist(t(lgm5.p)))
aux <- subset(ids, LGm == "LGm6")
lgm6.p <- validation.set.p[,as.character(aux$id)]
hc.lgm6.p <- hclust(dist(t(lgm6.p)))

order <- cbind(lgm1.p[,hc.lgm1.p$order],lgm2.p[,hc.lgm2.p$order],lgm3.p[,hc.lgm3.p$order],lgm4.p[,hc.lgm4.p$order],lgm5.p[,hc.lgm5.p$order],lgm6.p[,hc.lgm6.p$order])
validation.info <- ids[colnames(order),c("original.cluster","Histology","clustM.supervised","LGm")]
validation.info$Histology <- as.factor(validation.info$Histology)
validation.info$clustM.supervised <- as.factor(validation.info$clustM.supervised)
levels(validation.info$Histology) <- c("red","purple","cyan","green","gray","darkblue","orange")
validation.info$LGm <- as.factor(validation.info$LGm)
levels(validation.info$LGm) <- c("green","red","purple","orange","yellow","blue")
validation.info$original.cluster <- as.factor(validation.info$original.cluster)
levels(validation.info$original.cluster) <-  c("cornflowerblue","chartreuse4","yellow","lightgreen","black","blue4","orchid","green3","gold3","firebrick4","darkorange2","lightblue4","red") 
levels(validation.info$clustM.supervised) <- c("purple","#006400","#cf0000","orange","yellow","cornflowerblue","lightblue")

#png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_SturmPAMurTurcan5.png",res=1600,width=8000,height=8000)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_Validation_LABELSwt.pdf")
heatmap.plus.sm(as.matrix(order[1:2,]),
                col = jet.colors(75),
                scale = "none",
                #labRow = NA,
                #cexRow = 0.4,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = as.matrix(validation.info)
                #RowSideColors = rlab
                
)
dev.off()






########### IDHmut
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/EReg1_lowXhigh_021015_15genes.Rda")
a <- c(names(enrichmentList.sig$LGm3),names(sort(enrichmentList.sig$LGm2)[1:15]),as.character(EReg1_lowXhigh_021015$Row.names))
probe.meth <- str_extract(a,"^*[^.]*")


LGG.GBM.new.noXY.s <- LGG.GBM[ as.character(probe.meth),5:936]
LGG.GBM.new.order2 <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
colnames(LGG.GBM.new.order2) <- substr(colnames(LGG.GBM.new.order2),1,12)
LGG.GBM.new.order2 <- LGG.GBM.new.order2[,c(1:63,64:327,405:531,532:682,683:932,328:404)]

probe.meth <- str_extract(as.character(EReg1_lowXhigh_021015$Row.names),"^*[^.]*")
ESg1 <- LGG.GBM.new.order2[as.character(probe.meth),]
ESg1.p <- hclust(dist(ESg1))
probe.meth <- str_extract(names(sort(enrichmentList.sig$LGm2)[1:15]),"^*[^.]*")
ESg2 <- LGG.GBM.new.order2[as.character(probe.meth),]
ESg2.p <- hclust(dist(ESg2))
probe.meth <- str_extract(names(enrichmentList.sig$LGm3),"^*[^.]*")
ESg3 <- LGG.GBM.new.order2[as.character(probe.meth),]
ESg3.p <- hclust(dist(ESg3))

order <- rbind(ESg3[ESg3.p$order,],ESg2[ESg2.p$order,],ESg1[ESg1.p$order,])

d <- normal.dnameth
rownames(d) <- as.character(d$Composite.El)
d <- d[rownames(order),5:114]

a <- pd[pd$case.id %in% colnames(LGG.GBM.new.order2),]
lgm1.gcimpLow <- subset(a, clustM.supervised2 %in% "G-CIMP-low")
lgm1.gcimpLow  <- LGG.GBM.new.order2[,as.character(lgm1.gcimpLow$case.id)]
lgm1.gcimpLow.p <- hclust(dist(t(lgm1.gcimpLow)))
lgm1 <- subset(a, clustM.olddiscovery %in% "LGm1" & (clustM.supervised2 != "G-CIMP-low" | is.na(clustM.supervised2)))
lgm1 <- LGG.GBM.new.order2[,as.character(lgm1$case.id)]
lgm1.p <- hclust(dist(t(lgm1)))

order <- cbind(lgm1.gcimpLow[,lgm1.gcimpLow.p$order],lgm1[,lgm1.p$order],order[,64:932])


clab <- pd[colnames(order),c("TERT.expression","histology","clustM.supervised2","clustM.olddiscovery")]
clab$histology <- as.factor(clab$histology)
levels(clab$histology) <- c("red","purple","cyan","green")
clab$TERT.expression <- as.factor(clab$TERT.expression)
levels(clab$TERT.expression) <- c("black","purple")
clab$clustM.olddiscovery <- as.factor(clab$clustM.olddiscovery)
levels(clab$clustM.olddiscovery) <- c("green","red","purple","orange","yellow","blue")
clab$clustM.supervised2 <- as.factor(clab$clustM.supervised2)
levels(clab$clustM.supervised2) <- c("orange","red","green","purple","yellow","cornflowerblue","lightblue")
norm.label <- matrix("white", nrow = 110, ncol = 4,dimnames=list(NULL,colnames(clab)))
order <- cbind(d,order)
clab <- rbind(norm.label,clab)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_DNAMeth_Cartoon_LABELS3.pdf")
#png("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_DNAMeth_new2.png",res=1200,width=8000,height=8000)
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


c <- cpg.gene
rownames(c) <- as.character(c$Composite.El)
LGG.GBM.new.order <- c[rownames(order),644:1279]
aux <- as.character(c[rownames(order),c("Associated.Gene.Name")])
aux[32] <- "PAPPA2"
aux[33] <- "NNMT"
aux[36] <- "PLAT"
aux[43] <- "PTGFRN"

d <- mrna[aux,"entrezID"]
c <- rna.seq.norm[aux[c(32,33,36,43)],substr(colnames(LGG.GBM.new.order),1,12)] #substr(colnames(LGG.GBM.new.order),1,12)
LGG.GBM.new.order[43,] <- c[4,]
LGG.GBM.new.order[36,] <- c[3,]
LGG.GBM.new.order[33,] <- c[2,]
LGG.GBM.new.order[32,] <- c[1,]
#b <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/probesLGm2.txt",header=F)
#a <- c(rownames(aux),names(enrichmentList.sig$LGm5),names(enrichmentList.sig$LGm4),names(enrichmentList.sig$LGm3),names(sort(enrichmentList.sig$LGm2)[1:15]),as.character(b$V1))
#probe.meth <- str_extract(a,"*[^.]*")
colnames(LGG.GBM.new.order) <- substr(colnames(LGG.GBM.new.order),1,12)
a <- pd[pd$case.id %in% colnames(LGG.GBM.new.order),]
lgm1.gcimpLow <- subset(a, clustM.supervised2 %in% "G-CIMP-low")
lgm1.s.gcimpLow <- LGG.GBM.new.order[,as.character(lgm1.gcimpLow$case.id)]
mean <- apply(lgm1.s.gcimpLow, 1, mean, na.rm=TRUE)
lgm1 <- subset(a, clustM.olddiscovery %in% "LGm1" & clustM.supervised2 != "G-CIMP-low")
lgm1.s <- LGG.GBM.new.order[,as.character(lgm1$case.id)]
mean <- cbind(mean,apply(lgm1.s, 1, mean, na.rm=TRUE))
lgm2 <- subset(a, clustM.olddiscovery %in% "LGm2")
lgm2.s <- LGG.GBM.new.order[,as.character(lgm2$case.id)]
mean <- cbind(mean,apply(lgm2.s, 1, mean, na.rm=TRUE))
lgm3 <- subset(a, clustM.olddiscovery %in% "LGm3")
lgm3.s <- LGG.GBM.new.order[,as.character(lgm3$case.id)]
mean <- cbind(mean,apply(lgm3.s, 1, mean, na.rm=TRUE))
lgm4 <- subset(a, clustM.olddiscovery %in% "LGm4")
lgm4.s <- LGG.GBM.new.order[,as.character(lgm4$case.id)]
mean <- cbind(mean,apply(lgm4.s, 1, mean, na.rm=TRUE))
lgm5 <- subset(a, clustM.olddiscovery %in% "LGm5")
lgm5.s <- LGG.GBM.new.order[,as.character(lgm5$case.id)]
mean <- cbind(mean,apply(lgm5.s, 1, mean, na.rm=TRUE))
lgm6 <- subset(a, clustM.olddiscovery %in% "LGm6")
lgm6.s <- LGG.GBM.new.order[,as.character(lgm6$case.id)]
mean <- cbind(mean,apply(lgm6.s, 1, mean, na.rm=TRUE))
mean <- as.data.frame(mean)
colnames(mean) <- c("LGm1-GCIMPlow","LGm1","LGm2","LGm3","LGm4","LGm5","LGm6")
clab <-c("darkgreen","green","red","purple","orange","yellow","blue")
mean <- mean + 1
mean <- log2(mean)
#png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_MeanRNA3.png",res=1200,width=8000,height=8000)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_MeanRNA_2.pdf")
heatmap.plus.sm(as.matrix(mean),
                col=rev(redgreen(1000)),
                #scale = "none",
                #trace = "none",
                labRow = as.character(aux),
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                #ColSideColors = clab.6
                ColSideColors = cbind(clab,clab)
                #RowSideColors = rlab
                
)
dev.off()




###Random Forest (mur,Sturm,PA,turcan)
validation.set.p <- all[rownames(order),]

aux <- subset(ids, LGm == "LGm1")
lgm1.p <- validation.set.p[,as.character(aux$id)]
hc.lgm1.p <- hclust(dist(t(lgm1.p)))
aux <- subset(ids, LGm == "LGm2") 
lgm2.p <- validation.set.p[,as.character(aux$id)]
hc.lgm2.p <- hclust(dist(t(lgm2.p)))
aux <- subset(ids, LGm == "LGm3") 
lgm3.p <- validation.set.p[,as.character(aux$id)]
hc.lgm3.p <- hclust(dist(t(lgm3.p)))
aux <- subset(ids, LGm == "LGm4")
lgm4.p <- validation.set.p[,as.character(aux$id)]
hc.lgm4.p <- hclust(dist(t(lgm4.p)))
aux <- subset(ids, LGm == "LGm5")
lgm5.p <- validation.set.p[,as.character(aux$id)]
hc.lgm5.p <- hclust(dist(t(lgm5.p)))
aux <- subset(ids, LGm == "LGm6")
lgm6.p <- validation.set.p[,as.character(aux$id)]
hc.lgm6.p <- hclust(dist(t(lgm6.p)))

order <- cbind(lgm1.p[,hc.lgm1.p$order],lgm2.p[,hc.lgm2.p$order],lgm3.p[,hc.lgm3.p$order],lgm4.p[,hc.lgm4.p$order],lgm5.p[,hc.lgm5.p$order],lgm6.p[,hc.lgm6.p$order])
validation.info <- ids[colnames(order),c("Histology","original.cluster","clustM.supervised","LGm")]
validation.info$Histology <- as.factor(validation.info$Histology)
levels(validation.info$Histology) <- c("red","purple","cyan","green","gray","darkblue","orange")
validation.info$LGm <- as.factor(validation.info$LGm)
levels(validation.info$LGm) <- c("green","red","purple","orange","yellow","blue")
validation.info$original.cluster <- as.factor(validation.info$original.cluster)
levels(validation.info$original.cluster) <-  c("cornflowerblue","chartreuse4","yellow","lightgreen","black","blue4","orchid","green3","gold3","firebrick4","darkorange2","lightblue4","red") 
levels(validation.info$clustM.supervised) <- c("#9f20f0","#006400","#cf0000","#ffa500","#FFFF00","#006680","#00CCFF")


#png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_SturmPAMurTurcan_2.png",res=1600,width=8000,height=8000)
pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_Validation_LABELS3.pdf")
heatmap.plus.sm(as.matrix(order),
                col = jet.colors(75),
                scale = "none",
                #labRow = NA,
                #cexRow = 0.4,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = as.matrix(validation.info)
                #RowSideColors = rlab
                
)
dev.off()



a <- subset(ids, LGm %in% c("LGm4","LGm5","LGm6"))
a$LGm <- factor(a$LGm)
a$status <- as.numeric(as.character(a$status))
a$OS.months <- as.numeric(as.character(a$OS.months))
rownames(a) <- as.character(a$id)
aux <- subset(a, OS.months > 60)
a[rownames(aux),"OS.months"] <- 60 #colocar a coluna com a data do ultimo contato com o paciente
a[rownames(aux), "status"] <- 0 #colocar a coluna com informacao sobre o estado do paciente 
f.m = formula(Surv(OS.months,status) ~ LGm)
fit.m = survfit(f.m, data=a)

fit <- coxph(Surv(pd$os.months, pd$vital)~pd$cartoon.cluster)
summary(fit)

pdf(file = "/dados/ResearchProjects/thais/TCGA/LGG.GBM/survival_validation.pdf")
plot(fit.m, 
     lwd=4,
     #lty = c(2,1,2,1,1),
     col=c("orange",
       "yellow", #colocar a cor de acordo com a ordem de 1 a 6
           "blue"
           
           
     ),
     main="Kaplan-Meier Overall Survival Curves\nValidation Set\nDKFZ GBM, Mur and Turcan Samples ", 
     xlab="TIME SINCE DIAGNOSIS (MONTHS)", 
     ylab="PERCENT PROBABILITY OF SURVIVAL", 
     yscale=100,
     # xscale=365.25, 
     bg="black"
)
box(col="black", lwd=3);
#abline(h=0, lwd=4, lty=1)
legend("bottomleft", 
       legend=c(
         paste("LGm4 (n=",fit.m$n[1],")",sep=""),
         paste("LGm5 (n=",fit.m$n[2],")",sep=""),
         paste("LGm6 (n=",fit.m$n[3],")",sep="")
         
       ),
       #lty = c(1,3,1,3,1),
       col=c("orange",
       "yellow", #colocar a cor de acordo com a ordem de 1 a 6
           "blue"
             
       ),
       lwd=3,
       title="Legend",
       box.lwd=3,
       bg="white")
dev.off()


