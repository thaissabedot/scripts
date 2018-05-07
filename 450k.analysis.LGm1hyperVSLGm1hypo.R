##analysis on 450K only
##select 450K platforms within the new discovery (lgm1-hypo and lgm1-hyper)
##received platform info from Thais.
#'27_450k.txt'

tmp.plat <- pdata.new #object that contains new.discovery
platform <- read.delim(file = "/dados/ResearchProjects/tathi/27_450k.txt")
platform$ids <- as.character(substr(platform$longid, 1, 12))
tmp.plat <- merge(tmp.plat, platform, by.x = "id", by.y = "ids", all.x = T)
table(tmp.plat$new.discovery, tmp.plat$platform)
table(tmp.plat$new.discovery, tmp.plat$platform, tmp.plat$tumor.type)
tmp.plat.450 <- subset(tmp.plat, platform == "450K")
ids.newdiscovery.450 <- c(as.character(subset(tmp.plat, platform == "450K" & new.discovery == "LGm1.Hypergroup")$longid), as.character(subset(tmp.plat, platform == "450K" & new.discovery == "LGm1.Hypogroup")$longid), as.character(subset(tmp.plat, platform == "450K" & new.discovery == "LGm2.Hypergroup")$longid))
save(ids.newdiscovery.450, file = "ids.newdiscovery.450.rda") ## sent object to Thais to get the 450K data.
load(file = "/dados/ResearchProjects/tathi/LGG.GBM_Tathi/LGG.GBM.newdiscovey.450.rda")
load(file = "/dados/ResearchProjects/tathi/LGG.GBM_Tathi/ids.newdiscovery.450.rda")
library(exactRankTests)
values <- LGG.GBM.newdiscovey.450[,c(ids.newdiscovery.450)]
values <- t(na.omit(values))
values <- data.frame(values)
length(as.character(subset(tmp.plat, platform == "450K" & new.discovery == "LGm1.Hypergroup")$longid))
length(as.character(subset(tmp.plat, platform == "450K" & new.discovery == "LGm1.Hypogroup")$longid))
w.p.values <- unlist(mclapply(values,
                              function(probe) {
                                #LGm1 hyper vs LGm1 Hypo
                                zz <- wilcox.exact(as.matrix(probe[1:34]),as.matrix(probe[35:47]), exact=TRUE, paired=FALSE) 
                                z <- zz$p.value
                                return(z)
                              }, mc.cores=8))
save(w.p.values, file = "/dados/ResearchProjects/tathi/LGG.GBM_Tathi/lgm1hyperVSlgm1hypo.pvalues.rda")
load(file = "/dados/ResearchProjects/tathi/LGG.GBM_Tathi/lgm1hyperVSlgm1hypo.pvalues.rda")

### following confirms using a larger platform (450k), the changes observed between GBM/LGG within each new discovery is still the same.
lgm1.hyper.lgg <- apply(LGG.GBM.newdiscovey.450[,as.character(subset(tmp.plat, platform == "450K" & new.discovery == "LGm1.Hypergroup" & tumor.type == "LGG")$longid)], 2, mean, na.rm =T)
lgm1.hypo.lgg <- apply(LGG.GBM.newdiscovey.450[,as.character(subset(tmp.plat, platform == "450K" & new.discovery == "LGm1.Hypogroup" & tumor.type == "LGG")$longid)], 2, mean, na.rm =T)
lgm1.hypo.gbm <- mean(LGG.GBM.newdiscovey.450[,as.character(subset(tmp.plat, platform == "450K" & new.discovery == "LGm1.Hypogroup" & tumor.type == "GBM")$longid)], na.rm =T); names(lgm1.hypo.gbm) <- as.character(subset(tmp.plat, platform == "450K" & new.discovery == "LGm1.Hypogroup" & tumor.type == "GBM")$longid)
lgm1.hyper.gbm <- mean(LGG.GBM.newdiscovey.450[,as.character(subset(tmp.plat, platform == "450K" & new.discovery == "LGm1.Hypergroup" & tumor.type == "GBM")$longid)], na.rm =T); names(lgm1.hyper.gbm) <- as.character(subset(tmp.plat, platform == "450K" & new.discovery == "LGm1.Hypergroup" & tumor.type == "GBM")$longid)
boxplot(lgm1.hyper.lgg, lgm1.hyper.gbm, lgm1.hypo.lgg, lgm1.hypo.gbm)
t.test(c(lgm1.hypo.lgg, lgm1.hypo.gbm), c(lgm1.hyper.lgg, lgm1.hyper.gbm))
###end confirmation

##
lgggbm.newdis.450.order <- LGG.GBM.newdiscovey.450[,c(ids.newdiscovery.450)]
meanlgm1.hypo <- apply(lgggbm.newdis.450.order[,35:47], 1, mean, na.rm =T)
meanlgm1.hyper <- apply(lgggbm.newdis.450.order[,1:34], 1, mean, na.rm =T)
meandiff.newdis <- meanlgm1.hyper - meanlgm1.hypo
padj.newdis <- p.adjust(w.p.values)

volcano <- cbind(meandiff.newdis[colnames(values)], padj.newdis, w.p.values)
volcano <- as.data.frame(volcano); names(volcano) <- c("meandiff.newdis", "padj.newdis", "p.newdis")
volcano$threshold <- "1"
a <- subset(volcano, padj.newdis < 0.05) ##10^-4 and meandiff at 0.4 to get 28 probes (16 in seas)
b <- subset(a, meandiff.newdis < -0.25) #hyper em "b" 
volcano[rownames(b),"threshold"] <- "2"
c <- subset(a, meandiff.newdis > 0.25) #hyper em "a"
volcano[rownames(c),"threshold"] <- "3"

ggplot(data=volcano, aes(x=meandiff.newdis, y= -1*log10(padj.newdis), colour=threshold)) +
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
table(volcano$threshold)
volcano.450 <- volcano
##got volcano.27 from file 'Mutation_NonCodel.R'
volcano.27$diffmean.LGm1_hyper_LGm1_hypo <- volcano.27$mean.LGm1_hyper - volcano.27$mean.LGm1_hypo #flipped difference.
plot(volcano.27[subset(volcano.27, threshold != 1)$ID, "diffmean.LGm1_hyper_LGm1_hypo"], volcano.450[subset(volcano.27, threshold != 1)$ID,"meandiff.newdis"])

#cpgmeta450.Island <- subset(cpgmeta450, type == "CpG Island (Irizarry's HMM)") #210054
#cpgmeta450.Shore <- subset(cpgmeta450, type == "Shores (flanking CpG Island - 2000bp)") #118896
#cpgmeta450.Opensea <- subset(cpgmeta450, type == "Open Sea") #155067
#volcano.450.m <- merge(volcano.450, cpgmeta450, by.x = 0, by.y = "cgID", all.x = T)
cpgmeta450.n <- cpgmeta450[rownames(volcano.450),]
volcano.450m <- cbind(volcano.450, cpgmeta450.n)
##enrichement/fold across known features.
c(round(table(subset(volcano.450m, threshold!=1)$type)/dim(subset(volcano.450m, threshold!=1))[[1]], 2) *100) / c(round(table(volcano.450m$type)/dim(volcano.450m)[[1]], 2) *100)
heatmap.2(as.matrix(lgggbm.newdis.450.order[rownames(subset(volcano.450m, threshold != 1 & type == "Open Sea")),]), trace = "none", scale = "none", col = jet.colors(75), Colv = NA)
levels(volcano.450m$threshold) <- c("Background", "Hypermethylated in LGm1-Hyper", "Hypomethylated in LGm1-Hypo")

ggplot(subset(volcano.450m, threshold != "Hypermethylated in LGm1-Hyper"), aes(factor(threshold), fill = type, color = type)) + 
  geom_bar(position = "fill") + 
  labs(title = "LGm1-Hypo Differentially Methylated\nCpG probes across Known Genomic Features", y = "Percent Total", x = "") + 
  scale_y_continuous(labels=percent) + scale_color_manual(values = c("black","black","black")) +
  #theme_bw() +
  theme(plot.title = element_text(size = rel(2), colour = "black", face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(colour = "black",size=15,face = "bold"),
        axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="bold"),
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="bold")) +
  annotate("text", x = 1, y = .9, label = "25%\n(n=91,572)") +
  annotate("text", x = 1, y = .6, label = "29%\n(n=107,544)") +
  annotate("text", x = 1, y = .25, label = "47%\n(n=178,424)") +
  annotate("text", x = 2, y = .9, label = "18%\n(n=280)") +
  annotate("text", x = 2, y = .5, label = "67%\n(n=1,031)") +
  annotate("text", x = 2, y = .07, label = "14%\n(n=215)") 
ggsave(file = "LGm1-hypo.cpgprobes.genomic.location.pdf", scale = 2)

##homer analysis
volcano.450m2 <- cbind(LGG.GBM.newdiscovey.450[rownames(volcano.450m),c(2:4)], volcano.450m)
volcano.450m2$type2 <- substr(volcano.450m2$type, 1,2)
homer.450m2 <- (subset(volcano.450m2, threshold !=1)[,c(2:3,8,10)])


homer.450m2$End <- homer.450m2$Genomic_Coordinate + 50
homer.450m2$Start <- homer.450m2$Genomic_Coordinate - 50
homer.450m2$cgID <- paste(homer.450m2$type2, homer.450m2$cgID, sep=".")
homer.450m2 <- homer.450m2[,c(1,6,5,3,4)]
homer.450m2$Chromosome <- paste("chr",homer.450m2$Chromosome, sep="")
homer.450m2$strand <- "0"
##each file was downloaded and used in my work office computer (Genetics.Office). Stored in a folder in my desktop called 'lgg.gbm450.opeanseas' (Houtan Noushmehr, Mar 27, 2015)
write.table(subset(homer.450m2)[,c(1:4,6)], file = "homer.450m2.All.txt", sep="\t", quote=F, row.names =F)
write.table(subset(homer.450m2, type2 == "Op")[,c(1:4,6)], file = "homer.450m2.Op.100bp.16probes.txt", sep="\t", quote=F, row.names =F)
write.table(subset(homer.450m2, type2 == "Cp")[,c(1:4,6)], file = "homer.450m2.Cp.txt", sep="\t", quote=F, row.names =F)
write.table(subset(homer.450m2, type2 == "Sh")[,c(1:4,6)], file = "homer.450m2.Sh.txt", sep="\t", quote=F, row.names =F)


exp.obs.chrom <- cbind(as.matrix(round(table(volcano.450m2$Chromosome)/dim(volcano.450m2)[[1]], 3)*100), as.matrix(round(table(subset(volcano.450m2, threshold!=1)$Chromosome)/dim(subset(volcano.450m2, threshold!=1))[[1]], 3)*100))
exp.obs.chrom <- as.data.frame(exp.obs.chrom)
names(exp.obs.chrom) <- c("exp", "obs")
exp.obs.chrom$fold <- round(exp.obs.chrom[,2]/exp.obs.chrom[,1], 3)


#Pasted from untitled file (31/03/15 by Tathi)
##analysis on 450K only
##select 450K platforms within the new discovery (lgm1-hypo and lgm1-hyper)
##received platform info from Thais.
#'27_450k.txt'

tmp.plat <- pdata.new #object that contains new.discovery
platform <- read.delim(file = "/dados/ResearchProjects/tathi/27_450k.txt")
platform$ids <- as.character(substr(platform$longid, 1, 12))
tmp.plat <- merge(tmp.plat, platform, by.x = "id", by.y = "ids", all.x = T)
table(tmp.plat$new.discovery, tmp.plat$platform)
table(tmp.plat$new.discovery, tmp.plat$platform, tmp.plat$tumor.type)
tmp.plat.450 <- subset(tmp.plat, platform == "450K")
ids.newdiscovery.450 <- c(as.character(subset(tmp.plat, platform == "450K" & new.discovery == "LGm1.Hypergroup")$longid), as.character(subset(tmp.plat, platform == "450K" & new.discovery == "LGm1.Hypogroup")$longid), as.character(subset(tmp.plat, platform == "450K" & new.discovery == "LGm2.Hypergroup")$longid))
save(ids.newdiscovery.450, file = "ids.newdiscovery.450.rda") ## sent object to Thais to get the 450K data.
load(file = "/dados/ResearchProjects/tathi/LGG.GBM_Tathi/LGG.GBM.newdiscovey.450.rda")
load(file = "/dados/ResearchProjects/tathi/LGG.GBM_Tathi/ids.newdiscovery.450.rda")
library(exactRankTests)
values <- LGG.GBM.newdiscovey.450[,c(ids.newdiscovery.450)]
values <- t(na.omit(values))
values <- data.frame(values)
length(as.character(subset(tmp.plat, platform == "450K" & new.discovery == "LGm1.Hypergroup")$longid))
length(as.character(subset(tmp.plat, platform == "450K" & new.discovery == "LGm1.Hypogroup")$longid))
w.p.values <- unlist(mclapply(values,
                              function(probe) {
                                #LGm1 hyper vs LGm1 Hypo
                                zz <- wilcox.exact(as.matrix(probe[1:34]),as.matrix(probe[35:47]), exact=TRUE, paired=FALSE) 
                                z <- zz$p.value
                                return(z)
                              }, mc.cores=8))


