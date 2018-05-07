
############# GBM IDHwt H3K27ac
## data: GSE54047 10/19/2017
bigWigToBedGraph GSM1306369_MGH27PT.H3K27ac.bigWig GSM1306369_MGH27PT.H3K27ac.bedGraph
bigWigToBedGraph GSM1306371_MGH28PT.H3K27ac.bigWig GSM1306371_MGH28PT.H3K27ac.bedGraph
bigWigToBedGraph GSM1306373_MGH30PT.H3K27ac.bigWig GSM1306373_MGH30PT.H3K27ac.bedGraph


#cd ~/Downloads/GSE54047
a <- fread("GSM1306369_MGH27PT.H3K27ac.bedGraph")
a$id <- paste0("A",1:nrow(a))
fwrite(a,col.names=F,row.names=F,quote=F,sep="\t",file="GSM1306369_MGH27PT.H3K27ac.bedGraph")
#cut -f 1-3,5 GSM1306369_MGH27PT.H3K27ac.bedGraph | liftOver /dev/stdin /media/data1/biosoftwares/hg19ToHg38.over.chain GSM1306369_MGH27PT.H3K27ac_hg38.bedGraph
a <- fread("GSM1306369_MGH27PT.H3K27ac.bedGraph")
b <- fread("GSM1306369_MGH27PT.H3K27ac_hg38.bedGraph")
b <- merge(b,a[,4:5],by.x="V4",by.y="V5",sort=F)
fwrite(b[,c(2,3,4,1,5)],col.names=F,row.names=F,quote=F,sep="\t",file="GSM1306369_MGH27PT.H3K27ac_hg38.bedGraph")

#makeTagDirectory GSM1306369 GSM1306369_MGH27PT.H3K27ac_hg38.bedGraph -format bed -force5th -genome hg38


#cd ~/Downloads/GSE54047
a <- fread("GSM1306371_MGH28PT.H3K27ac.bedGraph")
a$id <- paste0("A",1:nrow(a))
fwrite(a,col.names=F,row.names=F,quote=F,sep="\t",file="GSM1306371_MGH28PT.H3K27ac.bedGraph")
#cut -f 1-3,5 GSM1306371_MGH28PT.H3K27ac.bedGraph | liftOver /dev/stdin /media/data1/biosoftwares/hg19ToHg38.over.chain GSM1306371_MGH28PT.H3K27ac_hg38.bedGraph unMapped
a <- fread("GSM1306371_MGH28PT.H3K27ac.bedGraph")
b <- fread("GSM1306371_MGH28PT.H3K27ac_hg38.bedGraph")
b <- merge(b,a[,4:5],by.x="V4",by.y="V5",sort=F)
fwrite(b[,c(2,3,4,1,5)],col.names=F,row.names=F,quote=F,sep="\t",file="GSM1306371_MGH28PT.H3K27ac_hg38.bedGraph")

#makeTagDirectory GSM1306371 GSM1306371_MGH28PT.H3K27ac_hg38.bedGraph -format bed -force5th -genome hg38

#cd ~/Downloads/GSE54047
a <- fread("GSM1306373_MGH30PT.H3K27ac.bedGraph")
a$id <- paste0("A",1:nrow(a))
fwrite(a,col.names=F,row.names=F,quote=F,sep="\t",file="GSM1306373_MGH30PT.H3K27ac.bedGraph")
#cut -f 1-3,5 GSM1306373_MGH30PT.H3K27ac.bedGraph | liftOver /dev/stdin /media/data1/biosoftwares/hg19ToHg38.over.chain GSM1306373_MGH30PT.H3K27ac_hg38.bedGraph unMapped
a <- fread("GSM1306373_MGH30PT.H3K27ac.bedGraph")
b <- fread("GSM1306373_MGH30PT.H3K27ac_hg38.bedGraph")
b <- merge(b,a[,4:5],by.x="V4",by.y="V5",sort=F)
fwrite(b[,c(2,3,4,1,5)],col.names=F,row.names=F,quote=F,sep="\t",file="GSM1306373_MGH30PT.H3K27ac_hg38.bedGraph")

#makeTagDirectory GSM1306373 GSM1306373_MGH30PT.H3K27ac_hg38.bedGraph -format bed -force5th -genome hg38

makeMultiWigHub.pl H3K27ac_IDHwt hg38 -webdir /media/data2/Visualization/ -url http://epigenome.fmrp.usp.br/public/gcimp/20171017/normal/ -d GSM1306369/ GSM1306371/ GSM1306373/


######### CNV track
a <- read.delim("/media/data2/CNV/GISTIC2_hg38/GCIMP.high/HIGH.all_lesions.conf_95.txt")
a <- a[,c(1,5)]

aux <- strsplit(as.character(a$Region.Limits), "\\:" )
a$Chr <- unlist(lapply(aux, "[[", 1 ))
aux <- unlist(lapply(aux, "[[", 2 ))
aux <- strsplit(aux, "\\-" )
a$Start <- unlist(lapply(aux, "[[", 1 ))
aux <- unlist(lapply(aux, "[[", 2 ))
aux <- strsplit(aux, "\\(" )
a$End <- unlist(lapply(aux, "[[", 1 ))
aux <- strsplit(as.character(a$Unique.Name), " " )
a$Event <- unlist(lapply(aux, "[[", 1 ))
a$color <- NA
a[a$Event %in% "Amplification","color"] <- "255,0,0" 
a[a$Event %in% "Deletion","color"] <- "0,0,255" 
a$score <- 0
a$strand <- "+"
a$id <- paste0("peak",1:nrow(a))
#write.table(a[,c("Chr","Start","End","Event")],quote=F,col.names=F,row.names = F,sep="\t",file="/media/data2/CNV/GISTIC2_hg38/GCIMP_high_CNV.txt")
write.table(a[1:41,c(3,4,5,10,8,9,4,5,7)],quote=F,col.names=F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/CNV/GCIMP_high_CNV.txt")

a <- read.delim("/media/data2/CNV/GISTIC2_hg38/GCIMP.low/LOW.all_lesions.conf_95.txt")
a <- a[,c(1,5)]

aux <- strsplit(as.character(a$Region.Limits), "\\:" )
a$Chr <- unlist(lapply(aux, "[[", 1 ))
aux <- unlist(lapply(aux, "[[", 2 ))
aux <- strsplit(aux, "\\-" )
a$Start <- unlist(lapply(aux, "[[", 1 ))
aux <- unlist(lapply(aux, "[[", 2 ))
aux <- strsplit(aux, "\\(" )
a$End <- unlist(lapply(aux, "[[", 1 ))
aux <- strsplit(as.character(a$Unique.Name), " " )
a$Event <- unlist(lapply(aux, "[[", 1 ))
a$color <- NA
a[a$Event %in% "Amplification","color"] <- "255,0,0" 
a[a$Event %in% "Deletion","color"] <- "0,0,255" 
a$score <- 0
a$strand <- "+"
a$id <- paste0("peak",1:nrow(a))
#write.table(a[,c("Chr","Start","End","Event")],quote=F,col.names=F,row.names = F,sep="\t",file="/media/data2/CNV/GISTIC2_hg38/GCIMP_low_CNV.txt")
write.table(a[1:12,c(3,4,5,10,8,9,4,5,7)],quote=F,col.names=F,row.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/CNV/GCIMP_low_CNV.txt")


####### 450k track
load("/media/data1/thais.projects/GCIMP-low/GBM.WGBS/low_high_450.rda")
lowXhigh_450k <- cbind(low[,-c(1:4)],high[,-c(1:4)])

my.t.wilcoxon.p.value <- function(...) {
  obj<-try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error"))
    return(NA)
  else
    return(obj$p.value)
}

values <- t(lowXhigh_450k)
values <- data.frame(values)
require(parallel)
w.p.values <- unlist(mclapply(values,
                              function(probe) {
                                z <- my.t.wilcoxon.p.value(as.matrix(probe[1:12]),as.matrix(probe[13:252]))
                                return(z)
                              }, mc.cores=12))
w.p.values.adj <- p.adjust(w.p.values, method = "fdr")
lowXhigh_450k$meanM1 <- apply(lowXhigh_450k[,1:12],1,mean,na.rm=T)
lowXhigh_450k$meanM2 <- apply(lowXhigh_450k[,13:252],1,mean,na.rm=T)
lowXhigh_450k$DiffMean <- lowXhigh_450k$meanM1 - lowXhigh_450k$meanM2
lowXhigh_450k$p.value.adj <- w.p.values.adj
lowXhigh_450k$p.value <- w.p.values

lowXhigh_450k$threshold <- "1"
a <- subset(lowXhigh_450k, p.value.adj < 0.01)
b <- subset(a, DiffMean < -0.3 ) #hyper no da direita
lowXhigh_450k[rownames(b),"threshold"] <- "2"
b <- subset(a, DiffMean > 0.3) #hyper no da esquerda
lowXhigh_450k[rownames(b),"threshold"] <- "3"

ggplot(lowXhigh_450k,aes(x=DiffMean,y=-1*log10(p.value.adj),color=threshold)) + geom_point()

clab <- cbind(c(rep("darkgreen",12),rep("firebrick",240)),c(rep("darkgreen",12),rep("firebrick",240)))

heatmap.plus(as.matrix(lowXhigh_450k[lowXhigh_450k$threshold %in% 3,1:252]),
             Colv=NA,
             col=jet.colors(75),
             cexCol = 0.5,
             labRow = NA,
             labCol = NA,
             ColSideColors = clab,
             scale = "none")

save(lowXhigh_450k,file="/media/data1/thais.projects/GCIMP-low/GBM.WGBS/lowXhigh_450k_wilcoxon.rda")
load("/media/data1/450k_hg38_GDC.rda")
hg38.450k <- hg38.450k[!(hg38.450k$Chromosome %in% "*"),]
hg38.450k$color <- "0,0,0"
hg38.450k[hg38.450k$threshold %in% "2","color"] <- "50,205,50"
hg38.450k[hg38.450k$threshold %in% "3","color"] <- "255,0,0"
hg38.450k$strand <- "+"
write.table(hg38.450k[,c(2,3,4,1,8,10,3,4,9)],quote=F,row.names=F,col.names=F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/probes450k/low_high_450k_new.bed")

#hg38.450k <- merge(hg38.450k,lowXhigh_450k[,c(255,256,258)],by.x="Composite.Element.REF",by.y=0)
#write.table(hg38.450k[hg38.450k$threshold %in% 1,c(2,3,4)],quote=F,row.names=F,col.names=F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/probes450k/low_high_450k.bed")

#write.table(hg38.450k[hg38.450k$threshold %in% 2,c(2,3,4)],quote=F,row.names=F,col.names=F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/probes450k/tmp.bed")

#write.table(hg38.450k[hg38.450k$threshold %in% 3,c(2,3,4)],quote=F,row.names=F,col.names=F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/probes450k/tmp.bed")





#### Blood H3K27ac
bigWigToBedGraph GSM1841087.bw GSM1841087.bw

a <- fread("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/blood/GSM1841087.bedGraph")
a$id <- paste0("A",1:nrow(a))
fwrite(a[,c(1,2,3,5,4)],col.names=F,row.names=F,quote=F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/blood/GSM1841087.bedGraph")

#makeTagDirectory GSM1841087 GSM1841087.bedGraph -format bed -force5th -genome hg38
#makeUCSCfile GSM1841087 -o auto

#### Pancreas H3K27ac
bigWigToBedGraph GSM1013129.bw GSM1013129.bedGraph


a <- fread("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/pancreas/files/GSM1013129.bedGraph")
a$id <- paste0("A",1:nrow(a))
fwrite(a[,c(1,2,3,5,4)],col.names=F,row.names=F,quote=F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/pancreas/files/GSM1013129.bedGraph")

#makeTagDirectory GSM1013129 GSM1013129.bedGraph -format bed -force5th -genome hg38
#makeUCSCfile GSM1013129 -o auto

makeMultiWigHub.pl H3K27ac_Control hg38 -webdir /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/control -url http://epigenome.fmrp.usp.br/public/gcimp/20171017/normal/ -d blood/files/GSM1841087/ pancreas/files/GSM1013129/

#### Normal Brain
bigWigToBedGraph GSM1112810.bw GSM1112810.bedGraph


a <- fread("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files/GSM1112810.bedGraph")
a$id <- paste0("A",1:nrow(a))
fwrite(a[,c(1,2,3,5,4)],col.names=F,row.names=F,quote=F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files/GSM1112810.bedGraph")

#makeTagDirectory GSM1112810 GSM1112810.bedGraph -format bed -force5th -genome hg38

bigWigToBedGraph GSM916035.bw GSM916035.bedGraph


a <- fread("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files/GSM916035.bedGraph")
a$id <- paste0("A",1:nrow(a))
fwrite(a[,c(1,2,3,5,4)],col.names=F,row.names=F,quote=F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/files/GSM916035.bedGraph")

#makeMultiWigHub.pl H3K27ac_NormalBrain hg38 -webdir /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/ -url http://epigenome.fmrp.usp.br/public/gcimp/20171017/normal/ -d GSM1112810/ GSM916035/


#### Gene track DEG
aux <- genes
aux$score <- 0
aux$colors <- NA
aux[aux$RNAseq.hg38.summary %in% "No.difference","color"] <- "169,169,169"
aux[aux$RNAseq.hg38.summary %in% "No.info","color"] <- "255,255,255"
aux[aux$RNAseq.hg38.summary %in% "Upregulated","color"] <- "255,0,0"
aux[aux$RNAseq.hg38.summary %in% "Downregulated","color"] <- "0,255,0"
write.table(aux[,c(2,3,4,8,29,5,3,4,31)],quote=F,row.names = F,col.names = F,sep="\t",file = "/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/genes/DEG.txt")


#### Glioblastoma Sample 1 (H3K4me3)
#Shell
wget http://dc2.cistrome.org/data5/browser/57069_treat.bw
mv 57069_treat.bw GSM1866119.bw
bigWigToBedGraph GSM1866119.bw GSM1866119.bedGraph
#end Shell

#R
library(data.table)
a <- fread("GSM1866119.bedGraph")
a$id <- paste0("A",1:nrow(a))
fwrite(a[,c(1,2,3,5,4)],col.names=F,row.names=F,quote=F,sep="\t",file="GSM1866119.bedGraph")
#end R

#Shell
makeTagDirectory GSM1866119 GSM1866119.bedGraph -format bed -force5th -genome hg38
#end Shell
#### End Glioblastoma Sample 1 (H3K4me3)

#### Glioblastoma Sample 2 (H3K4me3)
#Shell
wget http://dc2.cistrome.org/data5/browser/57077_treat.bw
mv 57077_treat.bw GSM1866127.bw
bigWigToBedGraph GSM1866127.bw GSM1866127.bedGraph
#end Shell

#R
library(data.table)
a <- fread("GSM1866127.bedGraph")
a$id <- paste0("A",1:nrow(a))
fwrite(a[,c(1,2,3,5,4)],col.names=F,row.names=F,quote=F,sep="\t",file="GSM1866127.bedGraph")
#end R

#Shell
makeTagDirectory GSM1866127 GSM1866127.bedGraph -format bed -force5th -genome hg38
#end Shell
#### End Glioblastoma Sample 2 (H3K4me3)

#### Glioblastoma Sample 3 (H3K4me3)
#Shell
wget http://dc2.cistrome.org/data5/browser/57067_treat.bw
mv 57067_treat.bw GSM1866117.bw
bigWigToBedGraph GSM1866117.bw GSM1866117.bedGraph
#end Shell

#R
library(data.table)
a <- fread("GSM1866117.bedGraph")
a$id <- paste0("A",1:nrow(a))
fwrite(a[,c(1,2,3,5,4)],col.names=F,row.names=F,quote=F,sep="\t",file="GSM1866117.bedGraph")
#end R

#Shell
makeTagDirectory GSM1866117 GSM1866117.bedGraph -format bed -force5th -genome hg38
#end Shell
#### End Glioblastoma Sample 3 (H3K4me3)


makeMultiWigHub.pl H3K4me3_GBM hg38 -webdir . -url http://epigenome.fmrp.usp.br/public/gcimp/20171017/XXXXXX/ -d GSM1866119/ GSM1866127/ GSM1866117/
  
  
##### Normal Brain Sample 1 (H3K4me3)
#Shell
wget http://dc2.cistrome.org/data5/browser/8269_treat.bw
mv 8269_treat.bw GSM772996.bw
bigWigToBedGraph GSM772996.bw GSM772996.bedGraph
#end Shell

#R
library(data.table)
a <- fread("GSM772996.bedGraph")
a$id <- paste0("A",1:nrow(a))
fwrite(a[,c(1,2,3,5,4)],col.names=F,row.names=F,quote=F,sep="\t",file="GSM772996.bedGraph")
#end R

#Shell
makeTagDirectory GSM772996 GSM772996.bedGraph -format bed -force5th -genome hg38
#end Shell
#### End Normal Brain Sample 1 (H3K4me3)

#### Normal Brain Sample 2 (H3K4me3)
#Shell
wget http://dc2.cistrome.org/data5/browser/3863_treat.bw
mv 3863_treat.bw GSM669992.bw
bigWigToBedGraph GSM669992.bw GSM669992.bedGraph
#end Shell

#R
library(data.table)
a <- fread("GSM669992.bedGraph")
a$id <- paste0("A",1:nrow(a))
fwrite(a[,c(1,2,3,5,4)],col.names=F,row.names=F,quote=F,sep="\t",file="GSM669992.bedGraph")
#end R

#Shell
makeTagDirectory GSM669992 GSM669992.bedGraph -format bed -force5th -genome hg38
#end Shell
#### End Normal Brain Sample 2 (H3K4me3)

#### Normal Brain Sample 3 (H3K4me3)
#Shell
wget http://dc2.cistrome.org/data5/browser/34339_treat.bw
mv 34339_treat.bw GSM773012.bw
bigWigToBedGraph GSM773012.bw GSM773012.bedGraph
#end Shell

#R
library(data.table)
a <- fread("GSM773012.bedGraph")
a$id <- paste0("A",1:nrow(a))
fwrite(a[,c(1,2,3,5,4)],col.names=F,row.names=F,quote=F,sep="\t",file="GSM773012.bedGraph")
#end R

#Shell
makeTagDirectory GSM773012 GSM773012.bedGraph -format bed -force5th -genome hg38
#end Shell
#### End Normal Brain Sample 3 (H3K4me3)


makeMultiWigHub.pl H3K4me3_NormalBrain hg38 -webdir /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/ -url http://epigenome.fmrp.usp.br/public/gcimp/20171017/normal/ -d GSM772996/ GSM669992/ GSM773012/
  
#### Normal Brain (H3K27ac)
#Shell
wget http://dc2.cistrome.org/data5/browser/38243_treat.bw
mv 38243_treat.bw GSM1112812.bw
bigWigToBedGraph GSM1112812.bw GSM1112812.bedGraph
#end Shell

#R
library(data.table)
a <- fread("GSM1112812.bedGraph")
a$id <- paste0("A",1:nrow(a))
fwrite(a[,c(1,2,3,5,4)],col.names=F,row.names=F,quote=F,sep="\t",file="GSM1112812.bedGraph")
#end R

#Shell
makeTagDirectory GSM1112812 GSM1112812.bedGraph -format bed -force5th -genome hg38
#end Shell
#### End Normal Brain (H3K27ac)

makeMultiWigHub.pl H3K27ac_NormalBrain hg38 -webdir /media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/normal.brain/ -url http://epigenome.fmrp.usp.br/public/gcimp/20171017/normal/ -d GSM1112810/ GSM916035/ GSM1112812/
  
### H3K27ac comparison
setwd("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer")
normal <- read.delim("outputPeaks_H3K27ac_gain_vsNormalBrain.txt")
gbm <- read.delim("outputPeaks_H3K27ac_gain_vsGBM.IDHwt.txt")
high <- read.delim("outputPeaks_histone_H3K27ac_gain.txt")
aux <- Reduce(subsetByOverlaps,list(makeGRangesFromDataFrame(high,TRUE),makeGRangesFromDataFrame(normal),makeGRangesFromDataFrame(gbm)))
a <- suppressWarnings(distanceToNearest(aux,tss,ignore.strand=T))
aux <- as.data.frame(aux)
aux <- aux[aux$bg.vs..target.adj..p.value < 0.015,]
aux$distance <- NA
aux[queryHits(a),"distance"] <- mcols(a)[,1]
aux <- aux[aux$Annotation %in% "Intergenic",]
aux <- aux[with(aux,order(bg.vs..target.adj..p.value)),]

###HiC data GBM IDHwt
setwd("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer")
load("/media/data2/HiC/HiC_commomRegions.rda")
gbm <- read.delim("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/H3K27ac/homer/outputPeaks_H3K27ac_nodiff_vsGBM.IDHwt.txt")
gbm <- gbm[gbm$bg.vs..target.adj..p.value > 0.05,]
aux <- findOverlaps(makeGRangesFromDataFrame(gbm),HiC.common)
gbm$HiC <- FALSE
gbm[queryHits(aux),"HiC"] <- TRUE
gbm <- gbm[unique(queryHits(aux)),]

normal <- read.delim("outputPeaks_H3K27ac_gain_vsNormalBrain.txt")
high <- read.delim("outputPeaks_histone_H3K27ac_gain.txt")
#aux <- findOverlaps(makeGRangesFromDataFrame(gbm)),HiC.)
aux <- Reduce(subsetByOverlaps,list(makeGRangesFromDataFrame(high,TRUE),makeGRangesFromDataFrame(gbm[gbm$HiC,])))
aux <- as.data.frame(aux)
#aux <- Reduce(subsetByOverlaps,list(makeGRangesFromDataFrame(high,TRUE),makeGRangesFromDataFrame(normal),makeGRangesFromDataFrame(gbm)))
aux <- as.data.frame(aux)


#### WGBS track
aux <- na.omit(LGG.GBM.WGBS[,c(1,2,3,17)])
aux <- aux[with(aux,order(chrom,start)),]
aux$end <- aux$start
aux$start <- as.character(aux$start)
aux$end <- as.character(aux$end)
aux$strand <- "+"
fwrite(aux[,c(1,2,3,5,4)],col.names = F,row.names = F,sep="\t",quote=F,file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/WGBS/GBM_IDHwt_3477.bedGraph")

makeTagDirectory GCIMP_low_0128 GCIMP_low_0128.bedGraph -format bed -force5th -genome hg38


#### GeneHancer track
load("/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GeneHancer/genehancer_table.rda")
a <- genehancer
a$geneHancer.score <- round(a$geneHancer.score*100)
write.table(a[,c(1,4,5,9,6)],quote=F,row.names = F,col.names = F,sep="\t",file="/media/data1/thais.projects/GCIMP-low/ActiveMotif_Data/Visualization/UCSC_track/GeneHancer/geneHancer_track.txt")
