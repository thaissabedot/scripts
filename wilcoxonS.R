####### Volcano LGm4+LGm5 x LGm6
#Beta Value Diff M1 x M2
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG-GBM.merged.data.FB20150924.Rdata")
levels(pd$clustM.supervised) <- c("IDHmut-codel","IDHmut-non-codel","IDHmut-hypo","IDHmut-non-codel","IDHwt-Mesenchymal","IDHwt-Classic","IDHwt-hypo")
pd$clustM.supervised <- as.character(pd$clustM.supervised)
pd[pd$clustM.supervised %in% "IDHwt-hypo" & pd$tumor.type == "GBM", "clustM.supervised"] <- "IDHwt-hypo-GBM"
pd[pd$clustM.supervised %in% "IDHwt-hypo" & pd$tumor.type == "LGG", "clustM.supervised"] <- "IDHwt-hypo-LGG"
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/tudao.27.450.rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_GBM_obj_20150127.Rda")


rownames(pd) <- as.character(pd$case.id)
aux <- pd[pd$clustM.supervised %in% "IDHwt-hypo-GBM" & !is.na(pd$clustM.supervised) & pd$has.RNAseq == "Yes",]
M1 <- LGG.GBM[,substr(colnames(LGG.GBM),1,12) %in% as.character(aux$case.id)]
aux <- pd[pd$clustM.supervised %in% "IDHwt-hypo-LGG" & !is.na(pd$clustM.supervised) & pd$has.RNAseq == "Yes",]
M2 <- LGG.GBM[,substr(colnames(LGG.GBM),1,12) %in% as.character(aux$case.id)]

volcano <- cbind(M1,M2)


IDHwt.GBMxLGG <- t(na.omit(volcano))
IDHwt.GBMxLGG <- data.frame(IDHwt.GBMxLGG)
#values.names <- colnames(values)

save(IDHwt.GBMxLGG,file='IDHwt.GBMxLGG.rda')




require(exactRankTests)
require(parallel)
dim(IDHwt.GBMxLGG)
system.time(w.p.values <- unlist(mclapply(IDHwt.GBMxLGG,
                                          function(probe) {
                                            zz <- wilcox.exact(probe[1:12],probe[13:37], exact=TRUE)
                                            z <- zz$p.value
                                            return(z)
                                          }, mc.cores=8)))

#lggs.p.value <- as.matrix(w.p.values)



#LGm4 x LGm5 - 1:136 137:362
#LGm4 x LGm6-LGG - 1:136 137:162
#LGm4 x LGm6-GBM - 1:136 137:177
#LGm5 x LGm6-GBM - 1:226 227:267
#LGm5 x LGm6-LGG - 1:226 227:252
#LGm6-GBM x LGm6-LGG - 1:41 42:67


IDHwt.GBMxLGG <- t(IDHwt.GBMxLGG)
IDHwt.GBMxLGG <- as.data.frame(IDHwt.GBMxLGG)
IDHwt.GBMxLGG$p.value <- w.p.values
IDHwt.GBMxLGG$p.value.adj <- p.adjust(IDHwt.GBMxLGG$p.value, method="BH")
IDHwt.GBMxLGG$meanM1 <- apply(IDHwt.GBMxLGG[,1:12],1,mean,na.rm=T)
IDHwt.GBMxLGG$meanM2 <- apply(IDHwt.GBMxLGG[,13:37],1,mean,na.rm=T)
IDHwt.GBMxLGG$DiffMean <- IDHwt.GBMxLGG$meanM1 - IDHwt.GBMxLGG$meanM2
dim(subset(IDHwt.GBMxLGG, p.value.adj < 0.05 & abs(DiffMean) > 0.3))  
save(IDHwt.GBMxLGG,file="IDHwt.GBMxLGG_pvalue.rda")
       
my.wilcox.test.p.value <- function(...) {
  obj<-try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error"))
    return(NA)
  else
    return(obj$p.value)
}

wilcoxt <- apply(volcano, 1, function(i) my.wilcox.test.p.value(i[1:204],i[205:244]) )
volcano$p.value = wilcoxt
volcano$p.value.adj = p.adjust(volcano$p.value, method = "fdr")

p4x5 <- rownames(subset(LGm4xLGm5, p.value.adj < 0.05 & abs(DiffMean) > 0.3))
p4x6l <- rownames(subset(LGm4xLGm6.LGG, p.value.adj < 0.05 & abs(DiffMean) > 0.3))
p4x6g <- rownames(subset(LGm4xLGm6.GBM, p.value.adj < 0.05 & abs(DiffMean) > 0.3))
p5x6 <- rownames(subset(LGm5xLGm6.GBM, p.value.adj < 0.05 & abs(DiffMean) > 0.3))
p6lx6g <- rownames(subset(LGm6.LGGxLGm6.GBM, p.value.adj < 0.05 & abs(DiffMean) > 0.3))



p <- c(p4x5,p4x6g,p4x6l,p5x6,p6lx6g,p5x6l)
rlab <- data.frame(cor=c(rep("pink",length(p4x5)),rep("pink2",length(p4x6g)),rep("pink4",length(p4x6l)),rep("purple",length(p5x6)),rep("palevioletred3",length(p6lx6g)),rep("palevioletred1",length(p5x6l))),id=p)
p <- unique(p)
p <- p6lx6g
rlab <- rlab[!duplicated(rlab$id),]
LGG.GBM.new.noXY.s <- LGG.GBM.new
rownames(LGG.GBM.new.noXY.s) <- as.character(LGG.GBM.new.noXY.s$Composite.Element.REF)
LGG.GBM.new.noXY.s <- LGG.GBM.new.noXY.s[,colnames(dat.lgg.gbm.new.noXY.dic.oecg)]
LGG.GBM.new.noXY.s <- LGG.GBM.new.noXY.s[ p,]
LGG.GBM.new.order2 <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
LGG.GBM.new.order2 <- LGG.GBM.new.order2[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
colnames(LGG.GBM.new.order2) <- substr(colnames(LGG.GBM.new.order2),1,12)
a <- pd
rownames(a) <- as.character(a$case.id)
a <- a[colnames(LGG.GBM.new.order2),]
levels(a$cluster.meth) <- c("green","red","purple","orange","yellow","blue")
a$tumor.type <- as.factor(a$tumor.type)
levels(a$tumor.type) <- c("gray","black") #LGG black
a$has.HM450 <- as.factor(a$has.HM450)
levels(a$has.HM450) <- c("gray","black") #black yes
#a <- colnames(LGG.GBM.new.order) %in% colnames(GBM.LGG.27.450k.nearest.gene)
#teste <- LGG.GBM.new.order[,a]
a <- a[856:932,c("cluster.meth","tumor.type","age","has.HM450")]
a$age <- jet.colors(100)[ceiling((a$age-min(a$age,na.rm=T))/(max(a$age,na.rm=T)-min(a$age,na.rm=T))*(length(jet.colors(100))-1))+1]
LGG.GBM.new.order2 <- as.matrix(LGG.GBM.new.order2[,856:932])
aux <- subset(a, has.HM450 == "gray")
aux2 <- subset(a, has.HM450 == "black")
a <- a[c(rownames(aux),rownames(aux2)),]

png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/heatmap_IDHwtProbes2.png",res=300,width=2000,height=2000)
#pdf("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_ES_DNAMeth_tert.pdf")
heatmap.plus.sm(LGG.GBM.new.order2[,c(rownames(aux),rownames(aux2))],
                col = jet.colors(75),
                scale = "none",
                #labRow = NA,
                labCol = NA,
                Colv = NA,
               # Rowv = NA,
                ColSideColors = as.matrix(a),
		#RowSideColors = cbind(as.character(rlab$cor),as.character(rlab$cor)),
                margins = c(1, 6)
                #RowSideColors = rlab
                
)
dev.off()






