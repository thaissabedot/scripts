####### Volcano LGm1-non.codel x LGm2-non.codel
#Beta Value Diff M1 x M2
load("pd.new.discovery_20150712.Rda")
load("all450.rda")

M1 <- pd.new.discovery[pd.new.discovery$new.discovery %in% "LGm1.Hypogroup",] ### 12 samples
M2 <- pd.new.discovery[pd.new.discovery$new.discovery %in% "LGm1.Hypergroup",] ### 35 samples

#Selecting the Samples
data <- all450
M1 <- data[,substr(colnames(data),1,12) %in% as.character(M1$case.id)]
M2 <- data[,substr(colnames(data),1,12) %in% as.character(M2$case.id)]
volcano <- cbind(M1,M2)

values <- t(na.omit(volcano))
values <- data.frame(values)


require(exactRankTests)
require(parallel)
system.time(w.p.values <- unlist(mclapply(values,
                                          function(probe) {
                                            zz <- wilcox.exact(probe[1:12],probe[13:47], exact=TRUE)
                                            z <- zz$p.value
                                            return(z)
                                          }, mc.cores=8)))


w.p.values.adj <- p.adjust(w.p.values)


volcano.450 <- volcano
volcano.450$meanHypo <- apply(volcano.450[,1:12],1,mean,na.rm=T)
volcano.450$meanHyper <- apply(volcano.450[,13:47],1,mean,na.rm=T)
volcano.450$DiffMean <- volcano.450$meanHyper - volcano.450$meanHypo #Hyper - Hypo
volcano.450 <- na.omit(volcano.450)
volcano.450$p.value.adj <- w.p.values.adj
volcano.450$p.value <- w.p.values

volcano.450 <- volcano.450[,c("meanHypo","meanHyper","DiffMean","p.value.adj","p.value")]
volcano.450$threshold <- "1"
a <- subset(volcano.450, p.value.adj < 0.05) ##10^-4 and meandiff at 0.4 to get 28 probes (16 in seas)
b <- subset(a, DiffMean < -0.25) #hyper no LGm1-hypo
volcano.450[rownames(b),"threshold"] <- "2"
c <- subset(a, DiffMean > 0.25) #hyper no LGm1-hyper
volcano.450[rownames(c),"threshold"] <- "3"

png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/volcano_LGm1hypo.png",res=200,width=2000,height=2000)
ggplot(data=volcano.450, aes(x=DiffMean, y= -1*log10(p.value.adj), colour=threshold)) +
  geom_point() +
  #scale_color_manual(values = c("black", "red", "green")) + 
  xlim(c(-1, 1)) + 
  #ylim(c(0,3)) +
  xlab("DNA Methylation Diff") + ylab("-1 * log10 of the Significance") + 
  labs(title = "DNA Methylation\nLGm1-hyper vs LGm1-hypo") +
  scale_color_manual(breaks=c("1","2","3"), # color scale (for points) 
                     values=c("black", "blue", "red"),
                     labels=c("Not Significant","Hypomethylated in LGm1-hyper","Hypermethylated in LGm1-hyper"),
                     name="Legend") + geom_hline(aes(yintercept= -1*log10(10^-4))) + geom_vline(aes(xintercept= -0.25)) + geom_vline(aes(xintercept= 0.25))
dev.off()


##homer analysis
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/gbm450.Rda")
load("Hm450k.cgprobes_June19_2015.Rda")
rownames(gbm.450k) <- gbm.450k$Composite.Element.REF
volcano.450m <- cbind(gbm.450k[rownames(volcano.450),c(2:4)], volcano.450)
volcano.450m <- cbind(volcano.450m, Hm450k.cgprobes[rownames(volcano.450m),5:6])
volcano.450m$type <- substr(volcano.450m$cgType, 1,2)

homer.450m2 <- (subset(volcano.450m, threshold !=1)[,c(2:3,12)])
homer.450m2$End <- homer.450m2$Genomic_Coordinate + 50
homer.450m2$Start <- homer.450m2$Genomic_Coordinate - 50
homer.450m2$cgID <- rownames(homer.450m2)
homer.450m2$cgID <- paste(homer.450m2$type, homer.450m2$cgID, sep=".")
homer.450m2 <- homer.450m2[,c(1,5,4,3,6)]
homer.450m2$Chromosome <- paste("chr",homer.450m2$Chromosome, sep="")
homer.450m2$strand <- "0"

###hypo e hyper

##each file was downloaded and used in my work office computer (Genetics.Office). Stored in a folder in my desktop called 'lgg.gbm450.opeanseas' (Houtan Noushmehr, Mar 27, 2015)
setw("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Homer")
write.table(subset(homer.450m2)[,c(1:4,6)], file = "homer.450m2.All.txt", sep="\t", quote=F, row.names =F)
write.table(subset(homer.450m2, type == "Op")[,c(1:4,6)], file = "homer.450m2.Op.100bp.16probes.txt", sep="\t", quote=F, row.names =F)
write.table(subset(homer.450m2, type == "Cp")[,c(1:4,6)], file = "homer.450m2.Cp.txt", sep="\t", quote=F, row.names =F)
write.table(subset(homer.450m2, type == "Sh")[,c(1:4,6)], file = "homer.450m2.Sh.txt", sep="\t", quote=F, row.names =F)


