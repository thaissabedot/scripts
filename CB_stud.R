info <- LGG.GBM.new[,c(1:4)]
LGG.GBM.new <- LGG.GBM.new[, colnames(LGG.GBM.new) %in% colnames(dat.lgg.gbm.new.noXY.dic.oecg)]

a <- pd[,c("cluster.meth","cluster.meth.gbm","case.id","tumor.type")]
a <- na.omit(a)
a <- subset(a, tumor.type %in% "GBM" & cluster.meth.gbm != "G-CIMP")
a <- sample(a$case.id,67)
c <- LGG.GBM.new[, substr(colnames(LGG.GBM.new),1,12) %in% as.character(a)]


a <- pd[,c("cluster.meth","cluster.meth.gbm","case.id","tumor.type")]
a <- na.omit(a)
a <- subset(a, cluster.meth.gbm %in% "G-CIMP")
a <- sample(a$case.id,33)
b <- LGG.GBM.new[, substr(colnames(LGG.GBM.new),1,12) %in% as.character(a)]


b <- cbind(b,c)
b <- cbind(info,b)

a <- pd
rownames(a) <- as.character(a$case.id)
a <- a[substr(colnames(GBMs),1,12)[5:104],c("case.id","age")]
a$case.id <- colnames(GBMs)[5:104]

