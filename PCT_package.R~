
####### START getting data to test

library(stringr)
load("/home/thais/Cloud_SanDiego/TCGA/LGG.GBM/RNAseq_05-11-14/mRNA.Rda")
mrna <- merge(lgg.mrna,gbm.mrna,by="gene_id")
colnames(mrna) <- substr(colnames(mrna),1,12)
mrna <- mrna[,-which(names(mrna) %in% c("TCGA-HT-A61A"))] #remove this sample
entrezID <- str_extract(mrna$gene_id,"[^|]*$")
mrna$entrezID <- as.matrix(entrezID)
rownames(mrna) <- as.character(mrna$entrezID)
mrna <- mrna[,-c(1,669)]
load("/home/thais/Package/control_test.rda")
rm(lgg.mrna,gbm.mrna,entrezID); gc();

####### END getting data to test

####### START function

### Libraries
library(EDASeq) 
library(steFunctions) 
library(biomaRt) 
library(stringr) 


############ START TESTING DATA
load("/home/thais/Cloud_SanDiego/TCGA/LGG.GBM/NewSamples_16-09/tudao.27.450.rda")
load("/home/thais/LGG-GBM/LGG-GBM.merged.data.FB20150930v3.Rdata") #pdata from Sep 2015 (LGG-GBM project). Created by Floris Barthel
rownames(pd) <- as.character(pd$case.id)
genomeVersion <- "hg19"
meth <- LGG.GBM[,-c(1:4)]
colnames(meth) <- substr(colnames(meth),1,12)
rm(LGG.GBM); gc();
############ END TESTING DATA

meth <- meth[,intersect(colnames(meth),colnames(mrna))]
mrna <- mrna[,intersect(colnames(meth),colnames(mrna))]


####### START getting data to test

library(stringr)
load("/home/thais/Cloud_SanDiego/TCGA/LGG.GBM/RNAseq_05-11-14/mRNA.Rda")
mrna <- merge(lgg.mrna,gbm.mrna,by="gene_id")
colnames(mrna) <- substr(colnames(mrna),1,12)
mrna <- mrna[,-which(names(mrna) %in% c("TCGA-HT-A61A"))] #remove this sample
entrezID <- str_extract(mrna$gene_id,"[^|]*$")
mrna$entrezID <- as.matrix(entrezID)
rownames(mrna) <- as.character(mrna$entrezID)
mrna <- mrna[,-c(1,669)]
load("/home/thais/Package/control_test.rda")
rm(lgg.mrna,gbm.mrna,entrezID); gc();

####### END getting data to test

####### START function

### Libraries
library(EDASeq) 
library(steFunctions) 
library(biomaRt) 
library(stringr) 


############ START TESTING DATA
load("/home/thais/Cloud_SanDiego/TCGA/LGG.GBM/NewSamples_16-09/tudao.27.450.rda")
load("/home/thais/LGG-GBM/LGG-GBM.merged.data.FB20150930v3.Rdata") #pdata from Sep 2015 (LGG-GBM project). Created by Floris Barthel
rownames(pd) <- as.character(pd$case.id)
genomeVersion <- "hg19"
meth <- LGG.GBM[,-c(1:4)]
colnames(meth) <- substr(colnames(meth),1,12)
rm(LGG.GBM); gc();
############ END TESTING DATA

meth <- meth[,intersect(colnames(meth),colnames(mrna))]
mrna <- mrna[,intersect(colnames(meth),colnames(mrna))]

############ START TESTING GROUPS
groups <- as.character(pd[colnames(meth),"clustM.olddiscovery"]); rm(pd); gc();
############ END TESTING GROUPS

######### START testing data
load("/home/thais/Cloud_SanDiego/TCGA/LGG.GBM/PCT.rda")
######### END testing data

rules <- read.delim(file="/home/thais/Cloud_SanDiego/TCGA/LGG.GBM/Rules_Patient_Centric_Calls.txt",sep="\t",header=T)
rules$Class <- paste(rules$Correlation, rules$Normal.Meth, rules$Tumor.Meth, sep=".")

PCT_call <- function(meth, 
                     mrna, 
                     genomeVersion="hg38",
                     control.meth, 
                     control.mrna,
                     cores=1) {
  
  if((!genomeVersion %in% c("hg18","hg19","hg38")) | length(genomeVersion) > 1)
    stop("Genome version not supported")
  
  
  ### START select probes and genes based on genome version
  meth.info <- CpG.probes[[genomeVersion]]
  meth.info <- meth.info[(meth.info$probeID %in% rownames(meth)) & (meth.info$gene.id %in% rownames(mrna)),]
  ### END select probes and genes based on genome version
  
  ### START prepare the data and sort the genes according to CpG mapping
  meth <- meth[as.character(meth.info$probeID),intersect(colnames(meth),colnames(mrna))]
  mrna <- mrna[as.character(meth.info$gene.id),intersect(colnames(meth),colnames(mrna))]
  n.samples <- ncol(meth)
  
  cpg.gene <- cbind(meth,mrna)
  #rm(genomeVersion,meth,mrna); gc();
  ### END prepare the data and sort the genes according to CpG mapping
  
  sd <- apply(cpg.gene[,(n.samples+1):ncol(cpg.gene)],1,sd,na.rm=TRUE)
  na <- apply(cpg.gene[,1:n.samples],1,function(x) all(is.na(x)))
  cpg.gene <- cpg.gene[sd!=0 & !is.na(sd) & !na,]
  meth.info <- meth.info[rownames(cpg.gene),]
  
  if(ncol(meth)<2 | nrow(meth.info)<2) 
    stop("Not enough observations")
  
  meth.info$cor <- apply(cpg.gene,1, function(i) {
    zz <- cor(i[1:n.samples],i[(n.samples+1):ncol(cpg.gene)],method="spearman",use="complete")
    return(zz)
  })
  
  
  epiGene <- unlist(mclapply(unique(meth.info$gene.id),
                             function(i){
                               aux <- meth.info[meth.info$gene.id %in% i,]
                               aux <- as.character(aux[with(aux,order(cor))[1],"probeID"])
                               return(aux)
                             },mc.cores=cores))
  
  meth.info <- meth.info[epiGene,]						
  meth.info$correlation_call <- "NNC"
  meth.info[meth.info$cor < -0.5 & !is.na(meth.info$cor),"correlation_call"] <- "SNC"
  meth.info[meth.info$cor >= -0.5 & meth.info$cor <= -0.25 & !is.na(meth.info$cor),"correlation_call"] <- "WNC"
  
  cpg.gene <- cpg.gene[rownames(meth.info),]
  
  control.meth <- control.meth[as.character(meth.info$probeID),intersect(colnames(control.meth),colnames(control.mrna))]
  control.mrna <- control.mrna[as.character(meth.info$gene.id),intersect(colnames(control.meth),colnames(control.mrna))]
  control <- cbind(control.meth,control.mrna)
  control <- control[rownames(cpg.gene),]
  n.control <- ncol(control.meth)
  if(n.control<2 | nrow(control)<2) 
    stop("Not enough observations for control group")
  
  tum.per <- apply(cpg.gene[,1:n.samples], 1, quantile, c(.1,.5,.9), na.rm=T)
  tum.per <- t(tum.per)
  colnames(tum.per) <- c("T10","T50", "T90")
  
  con.per <- apply(control[,1:n.control], 1, quantile, c(.1,.5,.9), na.rm=T)
  con.per <- t(con.per)
  colnames(con.per) <- c("N10","N50", "N90")
  
  meth.info <- cbind(meth.info, tum.per, con.per)
  
  meth.info$N.METH <- "VMN"
  meth.info$T.METH <- "VMT"
  
  meth.info[meth.info$N90 < 0.25 & !is.na(meth.info$N90),"N.METH"] <- "CUN"
  meth.info[meth.info$N10 > 0.75 & !is.na(meth.info$N10),"N.METH"] <- "CMN"
  meth.info[meth.info$N10 > 0.25 & meth.info$N90 < 0.75 & !is.na(meth.info$N90) & !is.na(meth.info$N10),"N.METH"] <- "IMN"
  
  meth.info[meth.info$T90 < 0.25 & !is.na(meth.info$T90),"T.METH"] <- "CUT"
  meth.info[meth.info$T10 > 0.75 & !is.na(meth.info$T10),"T.METH"] <- "CMT"
  meth.info[meth.info$T10 > 0.25 & meth.info$T90 < 0.75 & !is.na(meth.info$T90) & !is.na(meth.info$T10),"T.METH"] <- "IMT"
  
  meth.info$class <- paste(meth.info$correlation_call, meth.info$N.METH, meth.info$T.METH, sep=".")
  rules.all <- merge(meth.info, rules[,c(5:11)], by.x="class",by.y="Class",all.x=TRUE)
  rownames(rules.all) <- as.character(rules.all$probeID)
  rules.all <- rules.all[rownames(meth.info),]
  patient.call <- cpg.gene[,1:n.samples]
  
  tmp.lt.25 <- patient.call < 0.25
  tmp.lt.25[is.na(tmp.lt.25)] <- FALSE
  tmp.gt.75 <- patient.call > 0.75
  tmp.gt.75[is.na(tmp.gt.75)] <- FALSE
  tmp.ge.25 <- patient.call >= 0.25 
  tmp.le.75 <- patient.call <= 0.75
  tmp <- tmp.ge.25 + tmp.le.75
  tmp.ge.25.le.75 <- tmp == 2
  tmp.ge.25.le.75[is.na(tmp.ge.25.le.75)] <- FALSE
  
  for(samp in 1:dim(patient.call)[2]){
    patient.call[tmp.lt.25[,samp],samp] <- paste(as.character(rules.all[tmp.lt.25[,samp],"Call.TM.LT.0.25"]),as.character(rules.all[tmp.lt.25[,samp],"Confidence.Score.TM.LT.0.25"]),sep=".")
    patient.call[tmp.gt.75[,samp],samp] <- paste(as.character(rules.all[tmp.gt.75[,samp],"Call.TM.GT.0.75"]),as.character(rules.all[tmp.gt.75[,samp],"Confidence.Score.TM.GT.0.75"]),sep=".")
    patient.call[tmp.ge.25.le.75[,samp],samp] <- paste(as.character(rules.all[tmp.ge.25.le.75[,samp],"Call.0.25.GE.TM.LE.0.75"]),as.character(rules.all[tmp.ge.25.le.75[,samp],"Confidence.Score.0.25.GE.TM.LE.0.75"]),sep=".")
  }
  return(patient.call)
}
################ figure

colnames(patient.score)[1:636] <- substr(colnames(patient.score)[1:636],1,12)
colnames(patient.call)[1:636] <- substr(colnames(patient.call)[1:636],1,12)
patient.call$geneName <- rownames(patient.call)
patient.score$geneName <- rownames(patient.score)
patient <- rbind(patient.call, patient.score)

patient.t <- melt(patient.call[,c(1:637)],id.vars=c("geneName"))
patient.t2 <- melt(patient.score[,c(1:637)], id.vars = c("geneName"))
patient.t3 <- merge(patient.t, patient.t2, by=c("variable","geneName"))
LGG.GBM.new.noXY.s <- LGG.GBM.250914[,5:936]
LGG.GBM.new.order <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
colnames(LGG.GBM.new.order) <- substr(colnames(LGG.GBM.new.order),1,12)
LGG.GBM.new.order <- LGG.GBM.new.order[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
a <- colnames(LGG.GBM.new.order) %in% colnames(patient.call)
teste <- LGG.GBM.new.order[,a]
patient.t3$variable <- factor(patient.t3$variable, levels = colnames(teste)) #Order the samples (like the DNA meth heatmap)

teste <- subset(patient.t3, value.x %in% c("ES") & value.y %in% c("4"))
teste <- merge(patient.t3, metadata.1129.samples.20150312[,c("id","cluster.meth")],by.x="variable",by.y="id",all.x=T)

p <- ggplot(patient.t3, aes(x=cluster.meth)) +
geom_bar(stat = "identity")
ggsave(p, filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/p.png", width=30, height=20, dpi=200, units="in")

p <- ggplot(teste[1:1000,], aes(x=value.x,fill=score)) +
#  geom_tile(aes(fill=factor(score) ) ) +
  geom_bar(stat = "bin") +
  #scale_fill_manual(values=c("blue","orange","red","purple","green","white")) +
  scale_y_log10() +
  facet_grid(. ~ cluster.meth)


g

ggsave(p, filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/teste.png", width=30, height=20, dpi=200, units="in")


ggplot(teste2, aes(y=factor(variable),x=id)) + 
  geom_tile(aes(fill=factor(value), height=as.numeric(as.character(size)))) +
  scale_fill_manual(values=c("cadetblue1","brown4","cadetblue4","red","purple","cyan","magenta","orange","green3","green3","red","black","green","red","purple","orange","yellow","blue","yellow","gray","cyan","magenta","blue","green3","red","brown1","black","cornsilk2","gray","white"), guide = guide_legend(title = "Legend")) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white"),axis.text.x = element_blank(),axis.text.y = element_text(size=8),panel.margin = unit(0, "lines")) + #, strip.background=element_rect(fill=c("green","red","purple","orange","yellow","blue")),strip.text=element_text(size=12,color=c("white"))) + 
  facet_grid(SomaticType ~ .,space="free",scales="free") +
  xlab("TCGA LGG and GBM") + 
  ylab("Gene Name") +
  guides(col = guide_legend(nrow = 8))
dev.off()

patient.m2 <- merge(patient.m,metadata.1129.samples.20150112[,1:13],by.x="variable",by.y="id",all.x=T)

p <- ggplot(meanb.naomit, aes(x=reorder(gene.name,as.numeric(value)),fill=factor(value))) + 
  geom_bar(position="fill")+ 
  #geom_bar() +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values=c("black","gray"),name = "Legend", labels = c("Mutant", "WT")) +
  #geom_bar(position = "dodge")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(y = "Percentage") + coord_flip() +
  labs(x = "Gene Name") + 
  #facet_grid(. ~ newSturmClass, scales = "free")
  facet_grid(. ~ cluster.meth2, scales = "free")


write.table(patient.score[,1:636],file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/P_Score.txt",quote=F,row.names=T,sep="\t") #de 1 a 11
write.table(patient.call[,1:636],file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/P_Call.txt",quote=F,row.names=T,sep="\t") #de 1 a 11

}
