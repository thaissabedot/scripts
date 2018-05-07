#LGG New Samples. Batches: 

## distribuicao do SD
# 3 heatmaps
# age plot
# survival
# volcano
# heatmap LGG e GBM

x = 12 #batch
lgg.patients = NULL
normal = NULL
control = NULL
i=0
end = ".10.0"
beginning <- "/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_New_Batches_10-06/jhu-usc.edu_LGG.HumanMethylation450.Level_3." #folder
pattern = 'TCGA-([0-9A-Z]{2})-([0-9A-Z]{1,})-([0-9A-Z]{3})-([0-9A-Z]{3})-([0-9A-Z]{4}-([0-9]{2}))' #barcode
while (x<17){ #from 3.12.10.0 ate 3.16.10.0
  dir = (paste(c(paste(c(paste(c(beginning,x), collapse=''),""),collapse=''),end),collapse=''))
  files = list.files(dir,pattern="txt",full.names=T) 
  z = 1 #first file
  while (!is.na(files[z])){ 
    a = strsplit(files[z],"/") 
    b = a[[1]][[9]] #file name
    c <- str_extract(b, pattern)
    patientID <- c
    if(!is.na(c)){ #do not read metadata files
      if (substr(c,14,15) > "19"){ #control samples
        if(is.null(control)){ #read the first file
          control = read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1) 
          control = subset(control,select=c(Composite.Element.REF,Gene_Symbol,Chromosome,Genomic_Coordinate,Beta_value))
          colnames(control)[5] = patientID
        }
        else{
          aux = read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1)
          colnames(aux)[2] = patientID
          control = merge(control, aux[,1:2], by.x ="Composite.Element.REF", by.y = "Composite.Element.REF") 
        }
      } #end if control
      else if(substr(c,14,15) > "9"){ #normal samples
        if(is.null(normal)){ #read the first file
          normal = read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1) 
          normal = subset(normal,select=c(Composite.Element.REF,Gene_Symbol,Chromosome,Genomic_Coordinate,Beta_value)) 
          colnames(normal)[5] = patientID
        }
        else{
          aux = read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1)
          colnames(aux)[2] = patientID
          normal = merge(normal, aux[,1:2], by.x ="Composite.Element.REF", by.y = "Composite.Element.REF")
        }
      }
      else { #tumor samples
        if(is.null(lgg.patients)){ #read the first sample
          lgg.patients = read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1) 
          lgg.patients = subset(lgg.patients,select=c(Composite.Element.REF,Gene_Symbol,Chromosome,Genomic_Coordinate,Beta_value)) 
          colnames(lgg.patients)[5] = patientID
        }
        else{
          aux = read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1)
          colnames(aux)[2] = patientID
          lgg.patients = merge(lgg.patients, aux[,1:2], by.x ="Composite.Element.REF", by.y = "Composite.Element.REF")
        }
      }
    } #end !is.na(c)
    z = z + 1 #next file
  }#end while
  print(x)
  x = x + 1 #next folder
}
#save(lgg.patients,control,file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_New_Batches_10-06/lgg_new_samples.Rda")
library(stringr)
pattern.primary <- 'TCGA-([0-9A-Z]{2})-([0-9A-Z]{1,})-([0-1A-Z]{3})-([0-9A-Z]{3})-([0-9A-Z]{4}-([0-9]{2}))'
metadata <- lgg.patients[,1:4]
primary <- str_extract(colnames(lgg.patients)[5:ncol(lgg.patients)], pattern.primary) #remove the recurrent samples
primary <- na.omit(primary) #5 recurrents
lgg.patients <-lgg.patients[,primary]
lgg.patients <- cbind(metadata,lgg.patients) #220 samples
load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG.GBM.losangeles/lgg.gbm.normals.plusAnno.rda")
load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG.GBM.losangeles/lgg_289_samples.Rda")
LGG.509 <- merge(dat,lgg.patients[,c(1,5:224)], by.x="Composite.Element.REF", by.y="Composite.Element.REF")
dimnames(LGG.509)[[1]] <- LGG.509$Composite.Element.REF #dim 509
LGG.509.noXY <- LGG.509[LGG.509$Chromosome != "X" & LGG.509$Chromosome != "Y",]


## dichotomizing the data
LGG.509.noXY.dic <- LGG.509.noXY[ ,5:513]
LGG.509.noXY.dic <- LGG.509.noXY.dic>0.3
storage.mode(LGG.509.noXY.dic) <- "numeric"
## Use probes methylated more than 10% of the tumor samples
LGG.509.noXY.dic <- LGG.509.noXY.dic[(rowSums(LGG.509.noXY[,5:513])/ncol(LGG.509.noXY[,5:513]))*100 >10,]

#ead(Fb.HM450.oecg.3kb)
FDb.HM450.oecg.3kb.s <- subset(FDb.HM450.oecg.3kb, oecg.3kb>1.5)
LGG.509.noXY.dic.oecg <- LGG.509.noXY.dic[rownames(LGG.509.noXY.dic) %in% (FDb.HM450.oecg.3kb.s$probeid),]

#TODO:
## need to select tissue specific probes (>0.2)
avg.lgg.noXY <- apply(LGG.509.noXY[,5:513], 1, mean, na.rm=TRUE) 
a <- as.data.frame(avg.lgg.noXY)
avg.norm.lgg.noXY <- merge(a, avg.normals, by = 0)
library(LSD)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_New_Batches_10-06/comparisonplot.png",res=200,width=2000,height=2000)
comparisonplot(avg.norm.lgg.noXY$avg.normals, avg.norm.lgg.noXY$avg.lgg.noXY, pimp=TRUE, xlab="Media do valor B em pacientes normais",ylab="Media do valor B em pacientes com tumor")


lgg.tumor.specific <- (subset(avg.norm.lgg.noXY, avg.normals < 0.3))

##selecting tumor specific probes for lgg.gbm
dat.lgg.noXY.dic.oecg <- LGG.509.noXY.dic.oecg[(rownames(LGG.509.noXY.dic.oecg) %in% lgg.tumor.specific$Row.names), ]
#6011 probes

library(ConsensusClusterPlus)
title="/dados/ResearchProjects/thais/TCGA/LGG.GBM/ConsensusCluster.LGG.509.samples"
setwd(title)
cc.lgg <- ConsensusClusterPlus(d=dat.lgg.noXY.dic.oecg, 
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
ccICL = calcICL(cc.lgg, plot='pdf')
cc.lgg.2 <- cbind(as.matrix(cc.lgg[[2]][[3]]), 
                                        as.matrix(cc.lgg[[3]][[3]]), 
                                        as.matrix(cc.lgg[[4]][[3]]), 
                                        as.matrix(cc.lgg[[5]][[3]]), 
                                        as.matrix(cc.lgg[[6]][[3]]),
                                        as.matrix(cc.lgg[[7]][[3]]),
                                        as.matrix(cc.lgg[[8]][[3]]),
                                        as.matrix(cc.lgg[[9]][[3]]),
                                        as.matrix(cc.lgg[[10]][[3]]))

cc.lgg.2 <- as.data.frame(cc.lgg.2)
names(cc.lgg.2) <- c("K2",  "K3",  "K4" , "K5" , "K6" , "K7",  "K8" , "K9" , "K10")

cc.lgg.2$TCGAIDnew <- substr(rownames(cc.lgg.2),1,12)
cc.lgg.2$TCGAID <- rownames(cc.lgg.2)

#save(LGG.509,dat.lgg.noXY.dic.oecg,rpmm,cc.lgg,cc.lgg.2,file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_New_Batches_10-06/lgg_tcc.Rda")

### RPMM Cluster
library("RPMM")
lgg.rpmm <- LGG.509[rownames(dat.lgg.noXY.dic.oecg),5:513]
rpmm <- blcTree(t(lgg.rpmm), verbose = 0, splitCriterion = blcSplitCriterionLevelWtdBIC) #RPMM Clusters
rpmmClass <- blcTreeLeafClasses(rpmm)
a <- lapply(levels(rpmmClass), function(k) { colnames(dat.lgg.noXY.dic.oecg)[which(rpmmClass == k)]})
names(a) <- levels(rpmmClass)
column.order <- matrix(as.character(unlist(a))) #Get the column order by RPMM Cluster
rpmm.class <- NULL
rpmm.colors <- NULL
aux <- NULL
for(i in 1:length(a)){
  if(is.null(rpmm.class))
    rpmm.class <- matrix(c(rep(names(a)[[i]], length(a[[i]]))))
  else{
    aux <- 	matrix(c(rep(names(a)[[i]], length(a[[i]]))))
    rpmm.class <- rbind(rpmm.class,aux)
  }
}
#table(a)
rownames(rpmm.class) <- c(a$rLL,a$rLR,a$rRL,a$rRRL,a$rRRR)
rpmm.class <- str_replace(rpmm.class,"rLL","1")
rpmm.class <- str_replace(rpmm.class,"rLR","2")
rpmm.class <- str_replace(rpmm.class,"rRL","3")
rpmm.class <- str_replace(rpmm.class,"rRRL","4")
rpmm.class <- str_replace(rpmm.class,"rRRR","5")

### Hieraral Clustering
hc <- hclust(dist(t(dat.lgg.noXY.dic.oecg), method = "binary"), method = "ward")
hc.4.clusters <- as.matrix(cutree(hc,k=4))
hc.5.clusters <- as.matrix(cutree(hc,k=5))
hc.6.clusters <- as.matrix(cutree(hc,k=6))

##load in LGG clinical annotations
### March 14, 2014 from synapse:
#library("synapseClient")
#synapseLogin()
#syn2369104 <- synGet(id='syn2369104')
#syn2369104 <- synGet(id='syn2369104', load=T) #updated march 14
load(file="~/.synapseCache/615/490615/LGG-Annotation.3-14-2014.RData")
recurrent.cases <- data.frame(
  c("TCGA-DH-A669",  ## not in data freeze
    "TCGA-TM-A7CF",  ## not in data freeze
    "TCGA-FG-5963",
    "TCGA-FG-5965",
    "TCGA-FG-A4MT"),
  c(rep("Recurrent", times=5))
)
colnames(recurrent.cases) <- c("Tumor", "Recurrent.status")
Clinical <- merge(Clinical, recurrent.cases, by = "Tumor", all.x = T)
clinical.cluster <- merge(Clustering, Clinical, by = "Tumor")
clinical.cluster.genetics <- merge(clinical.cluster, Genetics, by = "Tumor")
oncosign <- read.delim(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/oncosign_assignments.txt", sep="\t", header=T) ## new oncosign assignments provided by giovanni (April 5, 2014)
clinical.cluster.genetics <- merge(clinical.cluster.genetics, oncosign, by.x = "Tumor", by.y = "sample", all.x =T)

#############new
library(downloader)
download("https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/lgg/bcr/biotab/clin/nationwideldrens.org_clinical_patient_lgg.txt", "clinical_data.txt")
cc.lgg.2$TCGAID <- rownames(cc.lgg.2)
lgg.anno <- merge(cc.lgg.2, clinical.cluster.genetics[,c("Tumor","MethylationCluster","COCCluster","histological_type","neoplasm_histologic_grade", "gender","age_at_initial_pathologic","death01","daystolastordeath","IDH1","TP53","PTEN")], by.x = "TCGAIDnew", by.y = "Tumor", all.x = T)
rownames(lgg.anno) <- lgg.anno$TCGAID
clinical.data <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_New_Batches_10-06/clinical_data.txt")
rownames(clinical.data) <- clinical.data$bcr_patient_barcode
clinicaldata.all <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/pdata.516lgg.606gbm.txt")

#Clinical Information for new LGG samples
lgg.newid <- colnames(lgg.patients)[5:224]
lgg.newid <- substr(lgg.newid,1,12)
new.lgg.patients <- clinical.data[rownames(clinical.data) %in% lgg.newid,]
new.lgg.patients$vital_status <- str_replace(new.lgg.patients$vital,"Alive","0")
new.lgg.patients$vital_status <- str_replace(new.lgg.patients$vital,"Dead","1")
a <- new.lgg.patients[,c("gender","tumor_grade","histologic_diagnosis","vital_status","age_at_initial_pathologic_diagnosis","last_contact_days_to","death_days_to")]
a <- subset(a, rownames(a) != "TCGA-TM-A7CF")
#write.table(a,file="lgg_tcc.txt",quote=F,row.names=T,sep="\t")

merged <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_New_Batches_10-06/lgg_tcc.txt")
a <- a[,1:6]
a$last_contact_days_to <- merged$Merge
rownames(lgg.anno) <- lgg.anno$TCGAIDnew
lgg.anno[rownames(a) ,c("gender","neoplasm_histologic_grade","histological_type","death01","age_at_initial_pathologic")] <- a[rownames(a) %in% lgg.newid,c("gender","tumor_grade","histologic_diagnosis","vital_status","age_at_initial_pathologic_diagnosis")]
#lgg.anno[rownames(lgg.anno) %in% lgg.newid,c("gender","neoplasm_histologic_grade","histological_type","death01","age_at_initial_pathologic")] <- lgg.anno2[rownames(lgg.anno2) %in% lgg.newid,c("gender","tumor_grade","histologic_diagnosis","vital_status","age_at_initial_pathologic_diagnosis")]

clinicaldata.all <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/pdata.516lgg.606gbm.txt")

mut.new <- clinicaldata.all[rownames(clinicaldata.all) %in% lgg.newid,c("mut.TP53","mut.IDH1","mut.PTEN")]

mut.new$mut.IDH1 <- str_replace(mut.new$mut.IDH1,"Mut","Mutant")
mut.new$mut.IDH1 <- str_replace(mut.new$mut.IDH1,"WT","wt") 
mut.new$mut.PTEN <- str_replace(mut.new$mut.PTEN,"Mut","Mutant")
mut.new$mut.PTEN <- str_replace(mut.new$mut.PTEN,"WT","wt") 
mut.new$mut.TP53 <- str_replace(mut.new$mut.TP53,"Mut","Mutant")
mut.new$mut.TP53 <- str_replace(mut.new$mut.TP53,"WT","wt") 
lgg.anno[rownames(mut.new),c("IDH1","PTEN","TP53")] <- mut.new[,c("mut.IDH1","mut.PTEN","mut.TP53")]

###merge clustering methods
cl.lgg.anno <- merge(lgg.anno,rpmm.class,by.x="TCGAID",by.y=0)
#cl.lgg.anno <- merge(cl.lgg.anno, hc.4.clusters,by.x="TCGAID",by.y=0)
cl.lgg.anno <- merge(cl.lgg.anno, hc.5.clusters,by.x="TCGAID",by.y=0)
colnames(cl.lgg.anno)[23:24] <- c("RPMM","HC.5")
#cl.lgg.anno <- merge(cl.lgg.anno, hc.6.clusters,by.x="TCGAID",by.y=0)
#colnames(cl.lgg.anno)[26] <- "HC.6"
##add in age
agecat <- cut(cl.lgg.anno$age_at_initial_pathologic, c(14,32.50,41,53,87))
cl.lgg.anno$age.strat <- as.character(agecat)
rownames(cl.lgg.anno) <- cl.lgg.anno$TCGAID

library(heatmap.plus)
library(matlab)


colors.label <- function(cluster,lgg.cluster,hclust, rpmm, age,coc.c,IDH1, TP53,grade,type){
  for(i in 1:length(cluster)){
    if(is.na(cluster[i]))
      cluster[i] <- "white"
    else{ 
      if(cluster[i] == "1")
        cluster[i] <- "red"
      else if(cluster[i] == "2")
        cluster[i] <- "yellow"
      else if(cluster[i] == "3")
        cluster[i] <- "green"
      else if(cluster[i] == "4")
        cluster[i] <- "blue"
      else if(cluster[i] == "5")
        cluster[i] <- "orange"
      else if(cluster[i] == "6")
        cluster[i] <- "purple"
      else if(cluster[i] == "7")
        cluster[i] <- "pink"
      else 
        cluster[i] <- "white"
    }
  }
  
  for(i in 1:length(hclust)){
    if(is.na(hclust[i]))
      hclust[i] <- "white"
    else{ 
      if(hclust[i] == "1")
        hclust[i] <- "red"
      else if(hclust[i] == "2")
        hclust[i] <- "yellow"
      else if(hclust[i] == "3")
        hclust[i] <- "green"
      else if(hclust[i] == "4")
        hclust[i] <- "blue"
      else if(hclust[i] == "5")
        hclust[i] <- "orange"
      else 
        hclust[i] <- "white"

    }
  }
  
  for(i in 1:length(rpmm)){
    if(is.na(rpmm[i]))
      rpmm[i] <- "white"
    else{ 
      if(rpmm[i] == "1")
        rpmm[i] <- "red"
      else if(rpmm[i] == "2")
        rpmm[i] <- "yellow"
      else if(rpmm[i] == "3")
        rpmm[i] <- "green"
      else if(rpmm[i] == "4")
        rpmm[i] <- "blue"
      else if(rpmm[i] == "5")
        rpmm[i] <- "orange"
      else 
        rpmm[i] <- "white"
    }
  }

  
  for(i in 1:length(lgg.cluster)){
    if(is.na(lgg.cluster[i]))
      lgg.cluster[i] <- "white"
    else{
      if(lgg.cluster[i] == "M5")
        lgg.cluster[i] <- "darkord4"
      else if(lgg.cluster[i] == "M1")
        lgg.cluster[i] <- "black"
      else if(lgg.cluster[i] == "M2")
        lgg.cluster[i] <- "red"
      else if(lgg.cluster[i] == "M3")
        lgg.cluster[i] <- "blue"
      else if(lgg.cluster[i] == "M4")
        lgg.cluster[i] <- "forestgreen"
      else 
        lgg.cluster[i] <- "white"
    }
  }

  
  for(i in 1:length(IDH1)){
    if(is.na(IDH1[i]))
      IDH1[i] <- "white"
    else{
      if(IDH1[i] == "Mutant")
        IDH1[i] <- "black"
      else if(IDH1[i] == "wt")
        IDH1[i] <- "grey"
      else
        IDH1[i] <- "white"
      
    }
  }
  
  for(i in 1:length(TP53)){
    if(is.na(TP53[i]))
      TP53[i] <- "white"
    else{
      if(TP53[i] == "Mutant")
        TP53[i] <- "black"
      else if(TP53[i] == "wt")
        TP53[i] <- "grey"
      else
        TP53[i] <- "white"
      
    }
  }
  
  
  for(i in 1:length(age)){
    if(is.na(age[i]))
      age[i] <- "white"
    else{
      if(age[i] == "(14,32.5]")
        age[i] <- "red"
      else if(age[i] == "(32.5,41]")
        age[i] <- "yellow"
      else if(age[i] == "(41,53]")
        age[i] <- "green"
      else if(age[i] == "(53,87]")
        age[i] <- "blue"
      else 
        age[i] <- "white"
    }
  }

  
  for(i in 1:length(coc.c)){
    if(is.na(coc.c[i]))
      coc.c[i] <- "white"
    else{
      if(coc.c[i] == "coc1")
        coc.c[i] <- "darkblue"
      else if(coc.c[i] == "coc2")
        coc.c[i] <- "darkcyan"
      else if(coc.c[i] == "coc3")
        coc.c[i] <- "cyan1"
      else
        coc.c[i] <- "white"
      
    }
  }
  
  
  for(i in 1:length(grade)){
    if(is.na(grade[i]))
      grade[i] <- "white"
    else{
      if(grade[i] == "G2")
        grade[i] <- "cyan"
      else if(grade[i] == "G3")
        grade[i] <- "darkblue"
      else
        grade[i] <- "white"
      
    }
  }
  
  for(i in 1:length(type)){
    if(is.na(type[i]))
      type[i] <- "white"
    else{
      if(type[i] == "Astrocytoma")
        type[i] <- "darkblue"
      else if(type[i] == "Oligoastrocytoma")
        type[i] <- "darkcyan"
      else if(type[i] == "Oligodendroglioma")
        type[i] <- "cyan1"
      else
        type[i] <- "white"
      
    }
  }
  
  #cluster,lgg.cluster,age,coc.c,IDH1,PTEN, TP53,grade,type
  cluster <- cbind(age,IDH1,TP53,grade,type,coc.c,lgg.cluster,rpmm,hclust,cluster)
  colnames(cluster) <- c("Age", "IDH1","TP53","Tumor Grade","Histological Type", "COC Clusters [NEJM, 2014]", "LGG Clusters [NEJM, 2014]", "RPMM"," Hierarquico","Consensus Cluster")
  return(cluster)
}


#Normal samples 
library(Biobase)
load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/normals.gse41826.rda")
norms <- exprs(controls.rm)
normals.sel <- norms[rownames(lgg.order),]
normals.sel <- normals.sel[heatmap.lgg$rowInd,]

source("/dados/ResearchProjects/thais/heatmap.plus.R")

#Heatmap for RPMM Clusters
LGG.509.noXY.s <- LGG.509.noXY[rownames(dat.lgg.noXY.dic.oecg) ,5:513]
lgg.order <- LGG.509.noXY.s[,rownames(rpmm.class)]
cl.lgg.anno <- cl.lgg.anno[colnames(lgg.order),]
#function(cluster,lgg.cluster,hclust, rpmm, age,coc.c,IDH1, TP53,grade,type){
clab.rpmm <- colors.label(as.matrix(cl.lgg.anno$K7),as.matrix(cl.lgg.anno$MethylationCluster),as.matrix(cl.lgg.anno$HC.5),as.matrix(cl.lgg.anno$RPMM),as.matrix(cl.lgg.anno$age.strat),as.matrix(cl.lgg.anno$COCCluster),as.matrix(cl.lgg.anno$IDH1),as.matrix(cl.lgg.anno$TP53),as.matrix(cl.lgg.anno$neoplasm_histologic_grade),as.matrix(cl.lgg.anno$histological_type))
clab.rpmm <- clab.rpmm[c(489:509,1:244,245:413,454:488,414:453),]
lgg.order <- lgg.order[heatmap.lgg$rowInd,c(489:509,1:244,245:413,454:488,414:453)]
normals.sel <- normals.sel[rownames(lgg.order),]
lgg.order <- cbind(normals.sel,lgg.order)
norm.label <- matrix("white", nrow = 77, ncol = 10)
clab.rpmm <- rbind(norm.label,clab.rpmm)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/heatmaps_tcc/heatmap_RPMM_2normals.png",res=200,width=2000,height=2000)
heatmap.plus.sm(as.matrix(lgg.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab.rpmm
)
dev.off()


#Heatmap for 4 CC Clusters
LGG.509.noXY.s <- LGG.509.noXY[rownames(dat.lgg.noXY.dic.oecg) ,5:513]
lgg.order <- LGG.509.noXY.s[,cc.lgg[[4]][[2]][[3]]]
cl.lgg.anno <- cl.lgg.anno[colnames(lgg.order),]
clab.4 <- colors.label(as.matrix(cl.lgg.anno$K4),as.matrix(cl.lgg.anno$MethylationCluster),as.matrix(cl.lgg.anno$age.strat),as.matrix(cl.lgg.anno$COCCluster),as.matrix(cl.lgg.anno$IDH1),as.matrix(cl.lgg.anno$TP53),as.matrix(cl.lgg.anno$neoplasm_histologic_grade),as.matrix(cl.lgg.anno$histological_type))
lgg.order <- lgg.order[heatmap.lgg$rowInd,]
normals.sel <- normals.sel[rownames(lgg.order),]
lgg.order <- cbind(normals.sel,lgg.order)
norm.label <- matrix("white", nrow = 77, ncol = 8)
clab.4 <- rbind(norm.label,clab.4)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/heatmaps_tcc/heatmap_K4_normals_3.png",res=200,width=2000,height=2000)
heatmap.plus.sm(as.matrix(lgg.order),
                                   col = jet.colors(75),
                                   scale = "none",
                                   labRow = NA,
                                   labCol = NA,
                                   Colv = NA,
                                   Rowv = NA,
                                   ColSideColors = clab.4
)
dev.off()

#Heatmap for 5 CC Clusters
LGG.509.noXY.s <- LGG.509.noXY[rownames(dat.lgg.noXY.dic.oecg) ,5:513]
lgg.order <- LGG.509.noXY.s[,cc.lgg[[5]][[2]][[3]]]
cl.lgg.anno <- cl.lgg.anno[colnames(lgg.order),]
clab.5 <- colors.label(as.matrix(cl.lgg.anno$K5),as.matrix(cl.lgg.anno$MethylationCluster),as.matrix(cl.lgg.anno$age.strat),as.matrix(cl.lgg.anno$COCCluster),as.matrix(cl.lgg.anno$IDH1),as.matrix(cl.lgg.anno$TP53),as.matrix(cl.lgg.anno$neoplasm_histologic_grade),as.matrix(cl.lgg.anno$histological_type))
lgg.order <- lgg.order[heatmap.lgg$rowInd,]
normals.sel <- normals.sel[rownames(lgg.order),]
lgg.order <- cbind(normals.sel,lgg.order)
norm.label <- matrix("white", nrow = 77, ncol = 8)
clab.5 <- rbind(norm.label,clab.5)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/heatmaps_tcc/heatmap_K5_normals_3.png",res=200,width=2000,height=2000)
heatmap.plus.sm(as.matrix(lgg.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab.5
)
dev.off()

#Heatmap for 6 CC Clusters
LGG.509.noXY.s <- LGG.509.noXY[rownames(dat.lgg.noXY.dic.oecg) ,5:513]
lgg.order <- LGG.509.noXY.s[,cc.lgg[[6]][[2]][[3]]]
cl.lgg.anno <- cl.lgg.anno[colnames(lgg.order),]
clab.6 <- colors.label(as.matrix(cl.lgg.anno$K6),as.matrix(cl.lgg.anno$MethylationCluster),as.matrix(cl.lgg.anno$age.strat),as.matrix(cl.lgg.anno$COCCluster),as.matrix(cl.lgg.anno$IDH1),as.matrix(cl.lgg.anno$TP53),as.matrix(cl.lgg.anno$neoplasm_histologic_grade),as.matrix(cl.lgg.anno$histological_type))
lgg.order <- lgg.order[heatmap.lgg$rowInd,]
normals.sel <- normals.sel[rownames(lgg.order),]
lgg.order <- cbind(normals.sel,lgg.order)
norm.label <- matrix("white", nrow = 77, ncol = 8)
clab.6 <- rbind(norm.label,clab.6)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/heatmaps_tcc/heatmap_K6_normals_3.png",res=200,width=2000,height=2000)
heatmap.plus.sm(as.matrix(lgg.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab.6
)
dev.off()

#Heatmap for 7 CC Clusters
LGG.509.noXY.s <- LGG.509.noXY[rownames(dat.lgg.noXY.dic.oecg) ,5:513]
lgg.order <- LGG.509.noXY.s[,cc.lgg[[7]][[2]][[3]]]
cl.lgg.anno <- cl.lgg.anno[colnames(lgg.order),]
clab.7 <- colors.label(as.matrix(cl.lgg.anno$K7),as.matrix(cl.lgg.anno$MethylationCluster),as.matrix(cl.lgg.anno$HC.5),as.matrix(cl.lgg.anno$RPMM),as.matrix(cl.lgg.anno$age.strat),as.matrix(cl.lgg.anno$COCCluster),as.matrix(cl.lgg.anno$IDH1),as.matrix(cl.lgg.anno$TP53),as.matrix(cl.lgg.anno$neoplasm_histologic_grade),as.matrix(cl.lgg.anno$histological_type))
clab.7 <- clab.7[c(452:468,93:248,249:355,1:92,356:415,469:509,416:451),]
lgg.order <- lgg.order[heatmap.lgg$rowInd,c(452:468,93:248,249:355,1:92,356:415,469:509,416:451)]
normals.sel <- normals.sel[rownames(lgg.order),]
lgg.order <- cbind(normals.sel,lgg.order)
norm.label <- matrix("white", nrow = 77, ncol = 10)
clab.7 <- rbind(norm.label,clab.7)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/heatmaps_tcc/heatmap_K7_normals_4.png",res=200,width=2000,height=2000)
heatmap.plus.sm(as.matrix(lgg.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab.7
)
dev.off()


#### HC Heatmaps
#Heatmap for 5 HC clusters
LGG.509.noXY.s <- LGG.509.noXY[rownames(dat.lgg.noXY.dic.oecg) ,5:513]
lgg.order <- LGG.509.noXY.s[,hc$order]
cl.lgg.anno <- cl.lgg.anno[colnames(lgg.order),]
clab.hc.5 <- colors.label(as.matrix(cl.lgg.anno$K7),as.matrix(cl.lgg.anno$MethylationCluster),as.matrix(cl.lgg.anno$HC.5),as.matrix(cl.lgg.anno$RPMM),as.matrix(cl.lgg.anno$age.strat),as.matrix(cl.lgg.anno$COCCluster),as.matrix(cl.lgg.anno$IDH1),as.matrix(cl.lgg.anno$TP53),as.matrix(cl.lgg.anno$neoplasm_histologic_grade),as.matrix(cl.lgg.anno$histological_type))
clab.hc.5 <- clab.hc.5[c(1:17,95:292,354:509,293:353,18:94),]
lgg.order <- lgg.order[heatmap.lgg$rowInd,c(1:17,95:292,354:509,293:353,18:94)]
normals.sel <- normals.sel[rownames(lgg.order),]
lgg.order <- cbind(normals.sel,lgg.order)
norm.label <- matrix("white", nrow = 77, ncol = 10)
clab.hc.5 <- rbind(norm.label,clab.hc.5)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/heatmaps_tcc/heatmap_hc5_2normals.png",res=200,width=2000,height=2000)
heatmap.plus.sm(as.matrix(lgg.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab.hc.5
)
dev.off()

## plot age at diagnosis
library(ggplot2)
x <- cl.lgg.anno
x$K7 <- factor(x$K7,c("4","1","2","7","6","3","5")) #change the order
ages <- x[,c("K7", "age_at_initial_pathologic")]
ages <- na.omit(ages)
p <- ggplot(ages, aes(factor(K7), age_at_initial_pathologic))
p + geom_boxplot(aes(fill = factor(K7)), notchwidth=0.25) + 
  geom_jitter(height = 0, position = position_jitter(width = .1), size=2)  + 
  scale_fill_manual(values=c("blue","red","pink", "yellow", "purple","green","orange"
  ), labels = c("M1", "M2", "M3","M4", "M5", "M6", "M7"
  )) +
  #scale_fill_manual(values = rev(c("darkord","darkgreen","blue","red","DARKGREY"))) +
  ylab("Age At Diagnosis") +
  xlab("DNA Methylation Clusters") + 
  labs(title = "DNA Methylation Clusters by Age at Diagnosis", fill = "DNA Meth. Cluster")
ggsave(file="age_diagnosis_K7_LGG_GBM_boxplot.png")

#p.value t test

#Kaplan-Meier Ajusted

mod.cl.lgg.anno <- cl.lgg.anno
aux <- subset(mod.cl.lgg.anno, daystolastordeath > 1825)
mod.cl.lgg.anno[rownames(aux),"daystolastordeath"] <- 1825
mod.cl.lgg.anno[rownames(aux), "death01"] <- "0"
clinical.cluster.genetics <- mod.cl.lgg.anno
age.quint <- quantile(clinical.cluster.genetics$age_at_initial_pathologic,
                      probs=seq(0,1,.2), na.rm=T)
#check values and make whole number version
age.quint <- c(14,32.50,41,53,87)
agecat <- cut(clinical.cluster.genetics$age_at_initial_pathologic, age.quint)
fdata <- clinical.cluster.genetics

fdata$age <- fdata$age_at_initial_pathologic
fdata$age2 <- agecat
fdata$sex <- fdata$gender
fdata$group <- as.factor(fdata$K7) ### also do $K5 $K6
library(survival)
data(uspop2)
refpop <- uspop2[as.character(floor(min(age.quint)):ceiling(max(age.quint))), , "2000"]
temp <- as.numeric(cut(floor(min(age.quint)):ceiling(max(age.quint)), age.quint, include.lowest=T))
pi.us<- tapply(refpop, list(temp[row(refpop)], col(refpop)), sum)/sum(refpop)
tab2 <- with(fdata, table(age2, sex, group))/ nrow(fdata)
pi.us2 <- pi.us[,2:1] # this puts male/female in the right order (check data)
us.wt <- rep(pi.us2, length(levels(fdata$group)))/ tab2
index <- with(fdata, cbind(as.numeric(age2), as.numeric(sex), as.numeric(group)))
fdata$uswt <- us.wt[index]
sfit.default <- survfit(Surv(daystolastordeath.x.5years, death01.x.5years) ~ group, fdata)
sfit3.adjusted <-survfit(Surv(as.numeric(daystolastordeath), as.numeric(death01)) ~ group, data=fdata, weight=uswt)

survLGG <- sfit3.adjusted
png(filename = "/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_New_Batches_10-06/survival.png", bg="white", res=300, width=2000, height=2000)
#par(mfrow=c(1,2))
plot(survLGG, 
     lwd=4,
     col=c("blue",
           "red",
           "pink",
           "yellow",
           "purple",
           "green",
           "orange"
           #"orange",
           #"purple"
     ),
     #col=c("black","red","blue","darkgreen","darkord"),
     main="Overall Survival\nKaplan-Meier Survival Curves\nTCGA LGG SAMPLES (DNA Methylation Clusters)", 
     #main="Progression Free\nKaplan-Meier Survival Curves\nTCGA LGG SAMPLES (DNA Methylation Clusters)", 
     xlab="TIME SINCE DIAGNOSIS (YEARS)", 
     ylab="PERCENT PROBABILITY OF SURVIVAL", 
     yscale=100,
     xscale=365.25, 
     bg="black"
)
#lines(sfit1a.d, mark.time=T, col=c("black","red","blue","darkgreen","darkorchid"), lty=2, lwd=1, xscale=365.25)
box(col="black", lwd=3);
#abline(v=5*52, lwd=5, lty=2); text(4.9*52, 0.75, "5 YEARS", srt=90)
abline(h=0, lwd=4, lty=1)
legend("bottomleft", 
       legend=c(paste("M1 (n=",survLGG$n[4],")",sep=""),
                paste("M2 (n=",survLGG$n[1],")",sep=""),
                paste("M3 (n=",survLGG$n[7],")",sep=""),
                paste("M4 (n=",survLGG$n[2],")",sep=""),
                paste("M5 (n=",survLGG$n[6],")",sep=""),
                paste("M6 (n=",survLGG$n[3],")",sep=""),
                paste("M7 (n=",survLGG$n[5],")",sep="")
       ),
       col=c("blue",
             "red",
             "pink",
             "yellow",
             "purple",
             "green",
             "orange"),
       lwd=3,
       title="DNA METHYLATION Clusters",
       box.lwd=3,
       bg="white")
#text(2,0.05,paste("Log-Rank p-value:",round(1 - psq(gg$sq, length(gg$n) - 1),10)), col="black")
dev.off()

sum(is.na(a[,c(14)]) & is.na(a[,c(15)]) & is.na(a[,c(16)]) & is.na(a[,c(17)]) & is.na(a[,c(18)]) & is.na(a[,c(19)]))

#save(cl.lgg.anno,LGG.509,file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_New_Batches_10-06/volcano.Rda")
#save(volcano,file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_New_Batches_10-06/fazervolcanoplot.Rda")
#|
#Beta Value Diff M1 x M2
a <- LGG.509[,5:513]
M1 <- cl.lgg.anno[cl.lgg.anno$K7 == "3" ,]
M2 <- cl.lgg.anno[cl.lgg.anno$K7 == "5" ,]
#Selecting the Samples
M1 <- a[,M1$TCGAID]
M2 <- a[,M2$TCGAID]
volcano <- cbind(M1,M2)
volcano$meanM1 <- apply(volcano[,1:41],1,mean,na.rm=T)
volcano$meanM2 <- apply(volcano[,42:77],1,mean,na.rm=T)
volcano$DiffMean <- volcano$meanM1 - volcano$meanM2

my.wilcox.test.p.value <- function(...) {
  obj<-try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error"))
    return(NA)
  else
    return(obj$p.value)
}

volcano$p.value <- apply(volcano, 1, function(i) my.wilcox.test.p.value(i[1:41],i[42:77]) )
volcano$p.value.adj = p.adjust(volcano$p.value, method = "fdr")

#a eh d , c eh b
c <- subset(volcano, p.value.adj < 0.05)
dim(c)
d <- subset(c,DiffMean < 0) #hyper no grupo da direita
sum(c$DiffMean > 0) #hyper no grupo da esquerda
b[1:15,c("meanM1","meanM2","DiffMean")]

volcano$threshold <- "1"
a <- subset(volcano, p.value.adj < 0.05)
b <- subset(a, DiffMean < 0) #hyper no 3 e 5
volcano[rownames(b),"threshold"] <- "2"
b <- subset(a, DiffMean > 0) #hyper no 1 7 2 6 
volcano[rownames(b),"threshold"] <- "3"
p <- ggplot(data=volcano, aes(x=DiffMean, y=-1*log10(p.value), colour=threshold)) +
  geom_point() +
  #scale_color_manual(values = c("black", "red", "green")) + 
  xlim(c(-1,1)) + ylim(c(0,25)) +
  xlab("DNA Methylation Diff") + ylab("-1 * log10 of the Significance") + 
  labs(title = "Volcano Plot")

p +     scale_color_manual(breaks=c("1","2","3"), # color scale (for points) 
                           values=c("black", "green", "red"),
                           labels=c("Not Significant","Hypermethylated in M3 e M5","Hypermethylated in M1, M7, M2 and M6"),
                           name="Legend") 

t.test(cl.lgg.anno[cl.lgg.anno$K7 == 3,"age_at_initial_pathologic"],cl.lgg.anno[cl.lgg.anno$K7 == 4,"age_at_initial_pathologic"])

library(LSD)
a <- cc.lgg.2
a1 <- a[a$K7 == "4",]
a2 <- a[a$K7 == "1",]
a3 <- a[a$K7 == "7",]
a4 <- a[a$K7 == "2",]
a5 <- a[a$K7 == "6",]
a6 <- a[a$K7 == "3",]
a7 <- a[a$K7 == "5",]
c <- LGG.509[,5:513]
normal <- as.data.frame(norms[rownames(norms) %in% rownames(LGG.509),])
M1 <- c[rownames(normal),as.character(a1$TCGAID)]
M2 <- c[rownames(normal),as.character(a2$TCGAID)]
M3 <- c[rownames(normal),as.character(a3$TCGAID)]
M4 <- c[rownames(normal),as.character(a4$TCGAID)]
M5 <- c[rownames(normal),as.character(a5$TCGAID)]
M6 <- c[rownames(normal),as.character(a6$TCGAID)]
M7 <- c[rownames(normal),as.character(a7$TCGAID)]
M1$meanM1 <- apply(M1,1,mean,na.rm=T)
M2$meanM2 <- apply(M2,1,mean,na.rm=T)
M3$meanM2 <- apply(M3,1,mean,na.rm=T)
M4$meanM2 <- apply(M4,1,mean,na.rm=T)
M5$meanM2 <- apply(M5,1,mean,na.rm=T)
M6$meanM2 <- apply(M6,1,mean,na.rm=T)
M7$meanM2 <- apply(M7,1,mean,na.rm=T)
normal$mean <- apply(normal,1,mean,na.rm=T)
avg.probes.cl <- cbind(normal$mean,M1$meanM1,M2$meanM2,M3$meanM2,M4$meanM2,M5$meanM2,M6$meanM2,M7$meanM2)
rownames(avg.probes.cl) <- rownames(M1)
colnames(avg.probes.cl) <- c("normal","M1","M2","M3","M4","M5","M6","M7")
png(filename = "/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_New_Batches_10-06/heatpairs_allprobes.png", bg="white", res=300, width=2000, height=2000)
heatpairs(avg.probes.cl)