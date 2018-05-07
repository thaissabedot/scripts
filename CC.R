tudao$sd <- apply(tudao,1,sd,na.rm=TRUE)
hist(tudao$sd)
dados <- subset(tudao, sd > 0.XX)
#verificar se tem NA


library(ConsensusClusterPlus)
title="/dados/ConsensusCluster"
setwd(title)
cc.out <- ConsensusClusterPlus(d=as.matrix(dados[,c(2:15)]), 
                                                     maxK=10, 
                                                     reps=1000, 
                                                     title=title, 
                                                     clusterAlg="km",
                                                     distance="euclidean",
                                                     seed=1022,
                                                     plot="pdf",
                                                     verbose = TRUE)
ccICL = calcICL(cc.out, plot='pdf')
cc.out.2 <- cbind(as.matrix(cc.out[[2]][[3]]), 
                                        as.matrix(cc.out[[3]][[3]]), 
                                        as.matrix(cc.out[[4]][[3]]), 
                                        as.matrix(cc.out[[5]][[3]]), 
                                        as.matrix(cc.out[[6]][[3]]),
                                        as.matrix(cc.out[[7]][[3]]),
                                        as.matrix(cc.out[[8]][[3]]),
                                        as.matrix(cc.out[[9]][[3]]),
                                        as.matrix(cc.out[[10]][[3]]))

cc.out.2 <- as.data.frame(cc.out.2)
names(cc.out.2) <- c("K2",  "K3",  "K4" , "K5" , "K6" , "K7",  "K8" , "K9" , "K10")

