library(minfi)

##### Pre-processing
setwd("/dados/ResearchProjects/thais/Prof_Ester/9828653150")
targets <- read.delim("/dados/ResearchProjects/thais/Prof_Ester/9828653150/SampleInfo.csv")
samples <- read.450k(paste(targets$Slide,targets$Array,sep="_"))
samples.bv <- as.data.frame(getBeta(samples))
qcReport(samples, sampNames = targets$SampleID, sampGroups = targets$Group.Info, pdf = "qcReport.pdf")

##### Column Names -> Replace File Name for IDs
colnames(samples.bv) <- as.character(targets$SampleID)

load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Hm450k.cgprobes_June19_2015.Rda")
##### 450k annotation
PE.project <- merge(Hm450k.cgprobes,samples.bv,by=0,all.y=T)
rownames(PE.project) <- as.character(PE.project$Row.names)
PE.project <- PE.project[,-1]

##### Standard Deviation
PE.project$SD <- apply(PE.project[,c(7:18)],1,sd,na.rm=T)
library(ggplot2)
ggplot(PE.project, aes(x = SD)) + stat_bin()
PE.project.SD_0.25 <- subset(PE.project, SD > 0.25) #902 probes

##### PCA Before Normalization
pca.data <- na.omit(PE.project.SD_0.25[,-c(1:6,19)])  #retirou as colunas 1:4 E omitiu NAs 

pca.data <- prcomp(t(pca.data), scales=TRUE) #calcula o principal components analysis
pca.scores <- data.frame(targets$Group.Info, pca.data$x[, 1:3]) #lgg.pca$x[, 1:3] means that estamos utilizando os 3 primeiros PCs

qplot(x=PC1, y=PC2, data=pca.scores, color=factor(targets$Group.Info)) + 
  ggtitle("Preeclampsia Project") + 
  geom_point(size = 4) +
  scale_color_manual(values=c("aquamarine4", "chocolate3"), 
                     #labels=c("Primary","Recurrent 1","Recurrent 2"),
                     name="Legend") 

##### 3D PCA (local) - Before Normalization
save(pca.scores,file="/dados/ResearchProjects/thais/Prof_Ester/Figures/pca3d_subset.rda")

##### Normalize
samples.norm <- bgcorrect.illumina(samples)
samples.norm.bv <- as.data.frame(getBeta(samples.norm))
samples.norm.bv_subset <- samples.norm.bv[rownames(PE.project.SD_0.25),]

##### PCA After Normalization
pca.data <- na.omit(samples.norm.bv_subset)  #retirou as colunas 1:4 E omitiu NAs 

pca.data <- prcomp(t(pca.data), scales=TRUE) #calcula o principal components analysis
pca.scores.norm <- data.frame(targets$Group.Info, pca.data$x[, 1:3]) #lgg.pca$x[, 1:3] means that estamos utilizando os 3 primeiros PCs

##### 3D PCA (local) - After Normalization
save(pca.scores.norm,file="/dados/ResearchProjects/thais/Prof_Ester/Figures/pca3d_norm_subset.rda")

library(rgl)
plot3d(pca.scores.norm[,2:4], col = c("chocolate3", "aquamarine4","chocolate3", "aquamarine4","chocolate3","aquamarine4","aquamarine4", "chocolate3","aquamarine4","chocolate3","chocolate3","aquamarine4"),type="s")
rgl.clear(type="lights")
rgl.light(-45, 20, ambient="black", diffuse="#dddddd", specular="white")
rgl.light(60, 30, ambient="#dddddd", diffuse="#dddddd", specular="black")
## before you animate your movie, make the window big
movie3d(spin3d(), fps=60, duration=10, dir="/home/thais")


##### Clinical Data - 17/08/2015
clin <- read.delim("/dados/ResearchProjects/thais/Prof_Ester/ClinicalData_17082015.txt")
rownames(targets) <- as.character(targets$SampleID)
rownames(clin) <- as.character(clin$ID)
clin <- clin[rownames(targets),]
clab <- clin
clab$Idade  <- jet.colors(100)[ceiling((as.numeric(clab$Idade)-min(as.numeric(clab$Idade),na.rm=T))/(max(as.numeric(clab$Idade),na.rm=T)-min(as.numeric(clab$Idade),na.rm=T))*(length(jet.colors(100))-1))+1]
clab$Idade[is.na(clab$Idade)] <- "white"
clab$Peso.RN  <- jet.colors(100)[ceiling((as.numeric(clab$Peso.RN)-min(as.numeric(clab$Peso.RN),na.rm=T))/(max(as.numeric(clab$Peso.RN),na.rm=T)-min(as.numeric(clab$Peso.RN),na.rm=T))*(length(jet.colors(100))-1))+1]
clab$Peso.RN[is.na(clab$Peso.RN)] <- "white"
levels(clab$Sexo.RN) <- c("lightpink1","lightblue1")
levels(clab$Status) <- c("aquamarine4","chocolate3")
clab$ID <- c("coral","cyan","darkorchid","darkseagreen3","firebrick1","olivedrab1","blue","magenta","darkgreen","hotpink3","indianred4","lightcyan4")

png(filename="/dados/ResearchProjects/thais/Prof_Ester/Figures/Heatmap_Norm_902p.png",res=300,width=2000,height=2000)
#pdf("/dados/ResearchProjects/thais/Prof_Ester/Figures/Heatmap_Norm_1101p.pdf")
heatmap.plus.sm(as.matrix(samples.norm.bv_subset),
                col=jet.colors(75),
                scale = "none",
                labRow = NA,
                #labCol = NA,
                #Colv = NA,
                #Rowv = NA,
                ColSideColors = as.matrix(clab)
                #RowSideColors = rlab
                
)
dev.off()


##### Subset by Standard Deviation based only on PE samples
PE.project$SD_PE <- apply(PE.project[,as.character(subset(clin,Status %in% "PE")$ID)],1,sd,na.rm=T)
PE.project.SD_PE_0.3 <- subset(PE.project, SD_PE > 0.3) #2044 probes

samples.norm.bv_subset <- samples.norm.bv[rownames(PE.project.SD_PE_0.3),]
png(filename="/dados/ResearchProjects/thais/Prof_Ester/Figures/Heatmap_Norm_600p.png",res=300,width=2000,height=2000)
#pdf("/dados/ResearchProjects/thais/Prof_Ester/Figures/Heatmap_Norm_1101p.pdf")
heatmap.plus.sm(as.matrix(samples.norm.bv_subset),
                col=jet.colors(75),
                scale = "none",
                labRow = NA,
                #labCol = NA,
                #Colv = NA,
                #Rowv = NA,
                ColSideColors = as.matrix(clab)
                #RowSideColors = rlab
                
)
dev.off()


#########2016
library(methylumi)
setwd("~/Cloud_SanDiego/Prof_Ester")
load("~/Cloud_SanDiego/TCGA/LGG.GBM/NewSamples_16-09/anno450k_CERTO.rda")

targets <- read.delim("~/Cloud_SanDiego/Prof_Ester/9828653150/SampleInfo.csv")
targets$fileName <- paste(targets$Slide,targets$Array,sep="_")

setwd("~/Cloud_SanDiego/Prof_Ester/9828653150")
files.idat <- methylumIDAT(targets$fileName)
files.proc <- stripOOB(normalizeMethyLumiSet(methylumi.bgcorr(files.idat)))
Samples.bv <- betas(files.proc)
Samples.bv <- merge(anno450k,Samples.bv,by.x="Composite.Element.REF",by.y=0)
save(Samples.bv,file="~/Cloud_SanDiego/Prof_Ester/Betas_2016.rda")


Samples.bv.NOsnp <- Samples.bv[Samples.bv$SNPs == FALSE & Samples.bv$Chromosome %in% c(1:22), ]

islands <- read.table("/home/thais/GBM WGBS/CpG_island.txt")
colnames(islands) <- c("chrom","start","end","id")
b <- Samples.bv.NOsnp
rownames(b) <- NULL
b$Chromosome <- paste0("chr",b$Chromosome)
b$CpG_type <- "CpG.openSea"
aux.GR <- makeGRangesFromDataFrame(b,start.field="Genomic_Coordinate",end.field="Genomic_Coordinate",seqnames.field="Chromosome")
aux <- findOverlaps(aux.GR,makeGRangesFromDataFrame(islands),type="within")
b[queryHits(aux),"CpG_type"] <- "CpG.island"
shore <- islands[,c(1,2,4)]
shore$startS <- shore$start - 2000
shore <- shore[,c(1,4,2,3)]
aux <- islands[,c(1,3,4)]
aux$end2 <- aux$end + 2000
aux <- aux[,c(1,2,4,3)]
colnames(shore) <- c("chrom","start","end","id")
colnames(aux) <- c("chrom","start","end","id")
shore <- rbind(shore,aux)
aux <- findOverlaps(aux.GR,makeGRangesFromDataFrame(shore),type="within")
b[queryHits(aux),"CpG_type"] <- ifelse(b[queryHits(aux),"CpG_type"] %in% "CpG.island","CpG.island","CpG.shore")
Samples.bv.NOsnp <- b

Samples.bv.NOsnp$SD <- apply(Samples.bv.NOsnp[,-c(1:5,18)],1,sd,na.rm=TRUE)
Samples.bv.NOsnp$CpG_type <- as.factor(Samples.bv.NOsnp$CpG_type)
Samples.bv.NOsnp$CpG_type <- factor(Samples.bv.NOsnp$CpG_type, levels = c("CpG.island","CpG.shore","CpG.openSea"))
library(ggplot2)
ggplot(Samples.bv.NOsnp, aes(x = SD,fill=CpG_type)) 
+ geom_histogram(binwidth = 0.0001) 
+ scale_fill_manual(values = c("darkslategrey","darkslategray4","darkslategray3")) 
ggsave("/home/thais/Cloud_SanDiego/Prof_Ester/Figures_new/hist.png",width = 18, height = 10)

Samples.bv.NOsnp.0_2SD <- Samples.bv.NOsnp[Samples.bv.NOsnp$SD > 0.2,] #3,038
pca.data <- na.omit(Samples.bv.NOsnp.0_2SD[,-c(1:5,18,19)])  #3,027

pca.data <- prcomp(t(pca.data), scales=TRUE) #calcula o principal components analysis
pca.scores <- data.frame(targets$Group.Info, pca.data$x[, 1:3]) #lgg.pca$x[, 1:3] means that estamos utilizando os 3 primeiros PCs

qplot(x=PC1, y=PC2, data=pca.scores, color=targets$Group.Info) + 
  ggtitle("Preeclampsia Project") + 
  geom_point(size = 4) +
  scale_color_manual(values=c("orchid4", "salmon"), 
                     #labels=c("Primary","Recurrent 1","Recurrent 2"),
                     name="Legend") 
ggsave("/home/thais/Cloud_SanDiego/Prof_Ester/Figures_new/pca2D_SD.png",width = 10, height = 10)
save(pca.scores,file="/home/thais/Cloud_SanDiego/Prof_Ester/Figures_new/pca_SD.rda")

library(rgl)
load("pca_SD.rda")
plot3d(pca.scores, col = c("salmon", "orchid4", "salmon", "orchid4", "salmon", "orchid4", "orchid4", "salmon","orchid4","salmon","salmon","orchid4"), type='s')
rgl.clear(type="lights")
rgl.light(-45, 20, ambient="black", diffuse="#dddddd", specular="white")
rgl.light(60, 30, ambient="#dddddd", diffuse="#dddddd", specular="black")
## before you animate your movie, make the window big
movie3d(spin3d(), fps=60, duration=10, dir="/home/thais/Cloud_SanDiego/Prof_Ester/Figures_new/",type = "gif")

#Pie chart (3D)
library(plotrix)
table(Samples.bv.NOsnp.0_2SD$CpG_type)
slices <- c(1795,  569, 674) 
lbls <- c("CpG island", "CpG shore", "Open sea")
pie3D(slices,labels=lbls,explode=0.1,
      main="Genomic Location of most variant CpG probes",col=c("darkslategrey","darkslategray4","darkslategray3"))

#Chromosome distribution
Samples.bv.NOsnp.0_2SD$Chromosome <- as.factor(Samples.bv.NOsnp.0_2SD$Chromosome)
Samples.bv.NOsnp.0_2SD$Chromosome <- factor(Samples.bv.NOsnp.0_2SD$Chromosome, levels = paste0("chr",c(1:22,"Y")))

ggplot(Samples.bv.NOsnp.0_2SD, aes(x = Chromosome,fill=CpG_type)) + 
  geom_bar() + scale_fill_manual(values = c("darkslategrey","darkslategray4","darkslategray3")) +
  labs(title="CpG distribution across chromosomes") 
ggsave("/home/thais/Cloud_SanDiego/Prof_Ester/Figures_new/chrom_dist_SD.png",width =12, height = 8)

#Hclust SD
targets <- read.table("~/Cloud_SanDiego/Prof_Ester/9828653150/SampleInfo.csv",header = TRUE,sep=",")
targets$ID <- paste(targets$Slide,targets$Array,sep="_")

aux <- as.character(subset(targets, Group.Info %in% c("Normal"))$ID)
normal <- Samples.bv.NOsnp.0_2SD[,colnames(Samples.bv.NOsnp.0_2SD) %in% aux]
aux <- hclust(dist(t(normal)))
normal <- normal[,aux$order]
aux <- as.character(subset(targets, Group.Info %in% c("PE"))$ID)
pe <- Samples.bv.NOsnp.0_2SD[,colnames(Samples.bv.NOsnp.0_2SD) %in% aux]
aux <- hclust(dist(t(pe)))
pe <- pe[,aux$order]
aux <- cbind(normal,pe)

#Heatmap SD
rownames(targets) <- as.character(targets$ID)
targets <- targets[colnames(aux),]
colnames(aux) <- as.character(targets$SampleID)
png(filename="/home/thais/Cloud_SanDiego/Prof_Ester/Figures_new/Heatmap_SD_0.2_3038p.png",width=1200,height=1200)
ha1 = HeatmapAnnotation(df=targets[,4:7])
Heatmap(as.matrix(aux), name = "DNA methylation", top_annotation = ha1,col=jet.colors(75), cluster_columns = FALSE,cluster_rows = TRUE,show_column_names = TRUE,show_row_names=FALSE)
dev.off()

##### Supervised analysis
info <- Samples.bv[,1:5]
Samples.bv <- Samples.bv[,-c(1:5)]
Samples.bv <- Samples.bv[,rownames(targets)]
rownames(Samples.bv) <- as.character(Samples.bv$Composite.Element.REF)
Samples.bv <- cbind(info, Samples.bv)
Samples.bv$p.value <- NA
for(i in 1:nrow(Samples.bv)){
  aux <- Samples.bv[i,]
  aux <- melt(aux[,6:17])
  levels(aux$variable) <- as.character(targets$Group.Info)
  if(sum(is.na(aux$value)) < 12)
    Samples.bv[i,"p.value"] <- pvalue(wilcox_test(value ~ variable,aux))
  
}
Samples.bv$p.value.adj <- p.adjust(Samples.bv$p.value,"fdr")

#Volcano
Samples.bv$meanNormal <- apply(Samples.bv[,6:11],1,mean,na.rm=T)
Samples.bv$meanPE <- apply(Samples.bv[,12:17],1,mean,na.rm=T)
Samples.bv$DiffMean <- Samples.bv$meanPE - Samples.bv$meanNormal
Samples.bv$threshold <- "1"
a <- subset(Samples.bv, p.value < 0.01 & !is.na(p.value))
b <- subset(a, DiffMean < -0.25 ) #hyper no da direita
Samples.bv[rownames(b),"threshold"] <- "2"
b <- subset(a, DiffMean > 0.25) #hyper no da esquerda
Samples.bv[rownames(b),"threshold"] <- "3"

ggplot(data=Samples.bv,aes(x=DiffMean, y=-1*log10(p.value),color=threshold)) +
  geom_point() +
  xlab("DNA methylation mean difference") + ylab("-1 * log10 of the Significance") +
  labs(title = "Volcano Plot") +
  scale_color_manual(breaks=c("1","2","3"), # color scale (for points)
                     values=c("grey", "green", "red"),
                     labels=c("Not Significant","Hypomethylated","Hypermethylated"),name="Legend")

#Heatmap
targets <- read.table("~/Cloud_SanDiego/Prof_Ester/9828653150/SampleInfo.csv",header = TRUE,sep=",")
targets$ID <- paste(targets$Slide,targets$Array,sep="_")

Samples.DMR <- Samples.bv[Samples.bv$threshold %in% c("2","3"),]

aux <- as.character(subset(targets, Group.Info %in% c("Normal"))$ID)
normal <- Samples.DMR[,colnames(Samples.DMR) %in% aux]
aux <- hclust(dist(t(normal)))
normal <- normal[,aux$order]
aux <- as.character(subset(targets, Group.Info %in% c("PE"))$ID)
pe <- Samples.DMR[,colnames(Samples.DMR) %in% aux]
aux <- hclust(dist(t(pe)))
pe <- pe[,aux$order]
order <- cbind(normal,pe)
rownames(order) <- as.character(Samples.DMR$Composite.Element.REF)

rownames(Samples.DMR) <- as.character(Samples.DMR$Composite.Element.REF)
aux <- as.character(Samples.bv[Samples.bv$threshold %in% c("2"),"Composite.Element.REF"])
hypo <- Samples.DMR[aux,6:17]
aux <- hclust(dist(t(hypo)))
hypo <- hypo[,aux$order]
aux <- as.character(Samples.bv[Samples.bv$threshold %in% c("3"),"Composite.Element.REF"])
hyper <- Samples.DMR[aux,6:17]
aux <- hclust(dist(t(hyper)))
hyper <- hyper[,aux$order]
order <- order[c(rownames(hypo),rownames(hyper)),]

rownames(targets) <- as.character(targets$ID)
targets <- targets[colnames(order),]
colnames(order) <- as.character(targets$SampleID)
png(filename="/home/thais/Cloud_SanDiego/Prof_Ester/Figures_new/Heatmap_DMR.png",width=1200,height=1200)
ha2 = rowAnnotation(df=data.frame(Probe=c(rep("Hypomethylated",nrow(hypo)),rep("Hypermethylated",nrow(hyper)))))
ha1 = HeatmapAnnotation(df=targets[,4:7])
ha <- Heatmap(as.matrix(order), 
        name = "DNA methylation", 
        top_annotation = ha1,
        col=jet.colors(75), 
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_column_names = TRUE,
        show_row_names=TRUE)
ha2+ha

islands <- read.table("/home/thais/GBM WGBS/CpG_island.txt")
colnames(islands) <- c("chrom","start","end","id")
b <- Samples.DMR
rownames(b) <- NULL
b$Chromosome <- paste0("chr",b$Chromosome)
b$CpG_type <- "CpG.openSea"
aux.GR <- makeGRangesFromDataFrame(b,start.field="Genomic_Coordinate",end.field="Genomic_Coordinate",seqnames.field="Chromosome")
aux <- findOverlaps(aux.GR,makeGRangesFromDataFrame(islands),type="within")
b[queryHits(aux),"CpG_type"] <- "CpG.island"
shore <- islands[,c(1,2,4)]
shore$startS <- shore$start - 2000
shore <- shore[,c(1,4,2,3)]
aux <- islands[,c(1,3,4)]
aux$end2 <- aux$end + 2000
aux <- aux[,c(1,2,4,3)]
colnames(shore) <- c("chrom","start","end","id")
colnames(aux) <- c("chrom","start","end","id")
shore <- rbind(shore,aux)
aux <- findOverlaps(aux.GR,makeGRangesFromDataFrame(shore),type="within")
b[queryHits(aux),"CpG_type"] <- ifelse(b[queryHits(aux),"CpG_type"] %in% "CpG.island","CpG.island","CpG.shore")
Samples.DMR <- b

pca.data <- na.omit(Samples.DMR[,-c(1:5,18:24)])  #3,027
targets <- targets[colnames(pca.data),]

pca.data <- prcomp(t(pca.data), scales=TRUE) #calcula o principal components analysis
pca.scores <- data.frame(targets$Group.Info, pca.data$x[, 1:3]) #lgg.pca$x[, 1:3] means that estamos utilizando os 3 primeiros PCs

qplot(x=PC1, y=PC2, data=pca.scores, color=targets$Group.Info) + 
  ggtitle("Preeclampsia Project") + 
  geom_point(size = 4) +
  scale_color_manual(values=c("orchid4", "salmon"), 
                     #labels=c("Primary","Recurrent 1","Recurrent 2"),
                     name="Legend") 
ggsave("/home/thais/Cloud_SanDiego/Prof_Ester/Figures_new/pca2D_DMR.png",width = 10, height = 10)
save(pca.scores,file="/home/thais/Cloud_SanDiego/Prof_Ester/Figures_new/pca_DMR.rda")

library(rgl)
load("pca_DMR.rda")
plot3d(pca.scores, col = c("orchid4", "orchid4", "orchid4", "orchid4", "orchid4", "orchid4", "salmon", "salmon","salmon","salmon","salmon","salmon"), type='s')
rgl.clear(type="lights")
rgl.light(-45, 20, ambient="black", diffuse="#dddddd", specular="white")
rgl.light(60, 30, ambient="#dddddd", diffuse="#dddddd", specular="black")
## before you animate your movie, make the window big
movie3d(spin3d(), fps=60, duration=10, dir="/home/thais/Cloud_SanDiego/Prof_Ester/Figures_new/",type = "gif")

#Pie chart (3D)
library(plotrix)
table(Samples.DMR$threshold,Samples.DMR$CpG_type)
#slices <- c(4,3,19) #hypo
slices <- c(21,9,6) #hyper
lbls <- c("CpG island", "CpG shore", "Open sea")
pie3D(slices,labels=lbls,explode=0.1,
      main="Genomic Location of most variant CpG probes",col=c("darkslategrey","darkslategray4","darkslategray3"))

#Chromosome distribution
Samples.DMR$Chromosome <- as.factor(Samples.DMR$Chromosome)
Samples.DMR$Chromosome <- factor(Samples.DMR$Chromosome, levels = paste0("chr",c(1:22,"X")))
Samples.DMR$CpG_type <- as.factor(Samples.DMR$CpG_type)
Samples.DMR$CpG_type <- factor(Samples.DMR$CpG_type, levels = c("CpG.island","CpG.shore","CpG.openSea"))
Samples.DMR$threshold <- as.factor(Samples.DMR$threshold)
levels(Samples.DMR$threshold) <- c("Hypomethylated","Hypermethylated")

ggplot(Samples.DMR, aes(x = Chromosome,fill=CpG_type)) + 
  geom_bar() + scale_fill_manual(values = c("darkslategrey","darkslategray4","darkslategray3")) +
  labs(title="CpG distribution across chromosomes") + facet_grid(. ~ threshold) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/home/thais/Cloud_SanDiego/Prof_Ester/Figures_new/chrom_dist_DMR.png",width =16, height = 8)

#nearest gene
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- genes(txdb_hg19)
symbols <- unlist(mapIds(org.Hs.eg.db, genes$gene_id, "SYMBOL", "ENTREZID", multiVals = "first"))
genes <- as.data.frame(genes)
genes$symbol <- as.matrix(symbols)

gene.GR <- makeGRangesFromDataFrame(genes,TRUE)
aux.GR <- makeGRangesFromDataFrame(Samples.DMR,start.field="Genomic_Coordinate",end.field="Genomic_Coordinate",seqnames.field="Chromosome")
aux <- distanceToNearest(aux.GR,gene.GR)
'Samples.DMR' <- cbind(Samples.DMR,genes[subjectHits(aux),c("gene_id","symbol")],data.frame(distance=mcols(aux)))
write.table(Samples.DMR[,c(1,3,4,18:27)],file="/home/thais/Cloud_SanDiego/Prof_Ester/Figures_new/genes_DMR.txt",quote=F,row.names=F,sep="\t")
