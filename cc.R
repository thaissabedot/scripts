### LGG.GBM.new.noXY is a matrix without chromosomes X and Y and without SNPs (NAs) 


## dichotomizing the data
LGG.GBM.new.noXY.dic <- LGG.GBM.new.noXY[ ,5:936] #Remove the annotation columns. Keep only tumor samples
LGG.GBM.new.noXY.dic <- LGG.GBM.new.noXY.dic>0.3
storage.mode(LGG.GBM.new.noXY.dic) <- "numeric"
## Use probes methylated more than 10% of the tumor samples
LGG.GBM.new.noXY.dic <- LGG.GBM.new.noXY.dic[(rowSums(LGG.GBM.new.noXY.dic)/ncol(LGG.GBM.new.noXY.dic))*100 >10,]
#LGG.GBM.new.noXY.dic <- LGG.GBM.new.noXY.dic[(rowSums(LGG.GBM.new.noXY[,5:936])/ncol(LGG.GBM.new.noXY[,5:936]))*100 >10,] #Old way

#head(FDb.HM450.oecg.3kb)
FDb.HM450.oecg.3kb.s <- subset(FDb.HM450.oecg.3kb, oecg.3kb>1.25)
LGG.GBM.new.noXY.dic.oecg <- LGG.GBM.new.noXY.dic[rownames(LGG.GBM.new.noXY.dic) %in% (FDb.HM450.oecg.3kb.s$probeid),]

## need to select tissue specific probes (>0.3)
avg.lgg.gbm.new.noXY <- apply(LGG.GBM.new.noXY[,5:936], 1, mean, na.rm=TRUE)
avg.norm.lgg.gbm.new.noXY <- merge(as.data.frame(avg.lgg.gbm.new.noXY), avg.normals, by = 0, all.x = T)
#library(LSD)
#comparisonplot(avg.norm.lgg.gbm.27.450.noXY$avg.normals, avg.norm.lgg.gbm.27.450.noXY$avg.lgg.gbm.27.450.noXY, pimp=TRUE)
lgg.gbm.new.tumor.specific <- (subset(avg.norm.lgg.gbm.new.noXY, avg.normals < 0.3))

##selecting tumor specific probes for lgg.gbm
dat.lgg.gbm.new.noXY.dic.oecg <- LGG.GBM.new.noXY.dic.oecg[(rownames(LGG.GBM.new.noXY.dic.oecg) %in% lgg.gbm.new.tumor.specific$Row.names), ]

#write.table(colnames(dat.lgg.gbm.new.noXY.dic.oecg), file="TCGA_IDs.txt",quote=F,row.names=F,sep="\t")

library(ConsensusClusterPlus)
title="/folder"
setwd(title)
cc.out.purity.lgg.gbm.new <- ConsensusClusterPlus(d=dat.lgg.gbm.new.noXY.dic.oecg, 
                                                     maxK=10, 
                                                     reps=1000, 
                                                     pItem=0.8, 
                                                     pFeature=1,
                                                     title=title, 
                                                     clusterAlg="hc",
                                                     distance="binary",
                                                     innerLinkage = "ward", 
                                                     finalLinkage = "ward",
                                                     seed=1022,
                                                     plot="pdf",
                                                     #writeTable=TRUE,
                                                     verbose = TRUE)
ccICL = calcICL(cc.out.purity.lgg.gbm.new, plot='pdf')
cc.out.purity.lgg.gbm.new.2 <- cbind(as.matrix(cc.out.purity.lgg.gbm.new[[2]][[3]]), 
                                        as.matrix(cc.out.purity.lgg.gbm.new[[3]][[3]]), 
                                        as.matrix(cc.out.purity.lgg.gbm.new[[4]][[3]]), 
                                        as.matrix(cc.out.purity.lgg.gbm.new[[5]][[3]]), 
                                        as.matrix(cc.out.purity.lgg.gbm.new[[6]][[3]]),
                                        as.matrix(cc.out.purity.lgg.gbm.new[[7]][[3]]),
                                        as.matrix(cc.out.purity.lgg.gbm.new[[8]][[3]]),
                                        as.matrix(cc.out.purity.lgg.gbm.new[[9]][[3]]),
                                        as.matrix(cc.out.purity.lgg.gbm.new[[10]][[3]]))

cc.out.purity.lgg.gbm.new.2 <- as.data.frame(cc.out.purity.lgg.gbm.new.2)
names(cc.out.purity.lgg.gbm.new.2) <- c("K2.new",  "K3.new",  "K4.new" , "K5.new" , "K6.new" , "K7.new",  "K8.new" , "K9.new" , "K10.new")
cc.out.purity.lgg.gbm.new.2$TCGAIDnew <- substr(rownames(cc.out.purity.lgg.gbm.new.2),1,12)
cc.out.purity.lgg.gbm.new.2$TCGAID <- rownames(cc.out.purity.lgg.gbm.new.2)
