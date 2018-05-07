#LGG Samples from 3.1.12.0 to 3.16.12.0
#GBM 450 sam from 3.1.5.0 to 3.12.5.0
#GBM 27k sam from 3.1.5.0 to 3.9.5.0
x = 1 #batch
gbm27 = NULL
normal = NULL
control = NULL
i=0
end = ".5.0"
beginning <- "/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/jhu-usc.edu_GBM.HumanMethylation27.Level_3." #folder
pattern = 'TCGA-([0-9A-Z]{2})-([0-9A-Z]{1,})-([0-9A-Z]{3})-([0-9A-Z]{3})-([0-9A-Z]{4}-([0-9]{2}))' #barcode
pattern.primary = 'TCGA-([0-9A-Z]{2})-([0-9A-Z]{1,})-([0-1A-Z]{3})-([0-9A-Z]{3})-([0-9A-Z]{4}-([0-9]{2}))' #barcode
while (x<10){ #from 3.1.5.0 ate 3.12.5.0
  dir = (paste(c(paste(c(paste(c(beginning,x), collapse=''),""),collapse=''),end),collapse=''))
  files = list.files(dir,pattern="txt",full.names=T) 
  z = 1 #first file
  while (!is.na(files[z])) {
                            a = strsplit(files[z],"/") 
                            b = a[[1]][[9]] #file name
                            c <- str_extract(b, pattern.primary)
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
                              else if(substr(c,14,15) > "09"){ #normal samples
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
                                if(is.null(gbm27)){ #read the first sample
                                  gbm27 = read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1) 
                                  gbm27 = subset(gbm27,select=c(Composite.Element.REF,Gene_Symbol,Chromosome,Genomic_Coordinate,Beta_value)) 
                                  colnames(gbm27)[5] = patientID
                                }
                                else{
                                  aux = read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1)
                                  colnames(aux)[2] = patientID
                                  gbm27 = merge(gbm27, aux[,1:2], by.x ="Composite.Element.REF", by.y = "Composite.Element.REF")
                                }
                              }
                            } #end !is.na(c)
                            z = z + 1 #next file
  }#end while
  print(x)
  x = x + 1 #next folder
}

LGG.GBM.new <- merge(lgg.450k,gbm.450k[,-c(2,3,4)],by="Composite.Element.REF")
LGG.GBM.new <- merge(LGG.GBM.new,gbm27[,-c(2,3,4)],by="Composite.Element.REF")
LGG.GBM.new <- na.omit(LGG.GBM.new)
load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG.GBM.losangeles/lgg.gbm.normals.plusAnno.rda")
rownames(LGG.GBM.new) <- LGG.GBM.new$Composite.Element.REF
LGG.GBM.new.noXY <- LGG.GBM.new[LGG.GBM.new$Chromosome != "X" & LGG.GBM.new$Chromosome != "Y",]
LGG.GBM.new.noXY <- LGG.GBM.250914[LGG.GBM.250914$Chromosome != "X" & LGG.GBM.250914$Chromosome != "Y",]


## dichotomizing the data
LGG.GBM.new.noXY.dic <- LGG.GBM.new.noXY[ ,5:936]
LGG.GBM.new.noXY.dic <- LGG.GBM.new.noXY.dic>0.3
storage.mode(LGG.GBM.new.noXY.dic) <- "numeric"
## Use probes methylated more than 10% of the tumor samples
LGG.GBM.new.noXY.dic <- LGG.GBM.new.noXY.dic[(rowSums(LGG.GBM.new.noXY.dic)/ncol(LGG.GBM.new.noXY.dic))*100 >10,]
#LGG.GBM.new.noXY.dic <- LGG.GBM.new.noXY.dic[(rowSums(LGG.GBM.new.noXY[,5:936])/ncol(LGG.GBM.new.noXY[,5:936]))*100 >10,]

#head(FDb.HM450.oecg.3kb)
FDb.HM450.oecg.3kb.s <- subset(FDb.HM450.oecg.3kb, oecg.3kb>1.25)
LGG.GBM.new.noXY.dic.oecg <- LGG.GBM.new.noXY.dic[rownames(LGG.GBM.new.noXY.dic) %in% (FDb.HM450.oecg.3kb.s$probeid),]

#TODO:
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
title="/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/ConsensusCluster.LGG.GBM.new_932_samples"
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


#load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/metadata.Rda")
#load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/LGG_GBM_new.Rda")

metadata.new$new.cluster <- NA
aux <- subset(metadata.new, cluster.mRNA == "LGr4" | cluster.mRNA == "LGr5")
aux <- subset(aux, cluster.meth2 == "LGm1" | cluster.meth2 == "LGm2")
metadata.new[rownames(aux),"new.cluster"] <- "IDHmut"
aux <- subset(metadata.new, cluster.mRNA == "LGr1" | cluster.mRNA == "LGr2")
aux <- subset(aux, cluster.meth2 == "LGm3")
metadata.new[rownames(aux),"new.cluster"] <- "IDHmut.codel"
aux <- subset(metadata.new, cluster.mRNA == "LGr6" | cluster.mRNA == "LGr7")
aux <- subset(aux, cluster.meth2 == "LGm4" | cluster.meth2 == "LGm5" | cluster.meth2 == "LGm6")
metadata.new[rownames(aux),"new.cluster"] <- "IDHwt"

metadata.new$cluster.meth2 <- NA
metadata.new[metadata.new$K6.new %in% "4","cluster.meth2"] <- "LGm1"
metadata.new[metadata.new$K6.new %in% "2","cluster.meth2"] <- "LGm2"
metadata.new[metadata.new$K6.new %in% "3","cluster.meth2"] <- "LGm3"
metadata.new[metadata.new$K6.new %in% "6","cluster.meth2"] <- "LGm4"
metadata.new[metadata.new$K6.new %in% "1","cluster.meth2"] <- "LGm5"
metadata.new[metadata.new$K6.new %in% "5","cluster.meth2"] <- "LGm6"

library(heatmap.plus)
library(matlab)


colors.label <- function(lgg.cluster.new,gbm.cluster.new,p19q,age,cluster,tumor,IDH,rna.dna.c,rna.cluster,coc.c,leukocyte,estimates,sturm,anno,cnv,tp53,egfr,braf,histology,grade,subtype){
  
  for(i in 1:length(cluster)){
    if(is.na(cluster[i]))
      cluster[i] <- "white"
    else{ 
      if(cluster[i] == "LGm1")
        cluster[i] <- "green"
      else if(cluster[i] == "LGm2")
        cluster[i] <- "red"
      else if(cluster[i] == "LGm3")
        cluster[i] <- "purple"
      else if(cluster[i] == "LGm4")
        cluster[i] <- "orange"
      else if(cluster[i] == "LGm5")
        cluster[i] <- "yellow"
      else if(cluster[i] == "LGm6")
        cluster[i] <- "blue"
      else 
        cluster[i] <- "white"
    }
  }
  
  for(i in 1:length(tp53)){
    if(is.na(tp53[i]))
      tp53[i] <- "white"
    else{
      if(tp53[i] == "Mut")
        tp53[i] <- "black"
      else if(tp53[i] == "WT")
        tp53[i] <- "grey"
      else
        tp53[i] <- "white"
      
    }
  }
  
  
  for(i in 1:length(egfr)){
    if(is.na(egfr[i]))
      egfr[i] <- "white"
    else{
      if(egfr[i] == "Hemi-del")
        egfr[i] <- "blue"
      # else if(egfr[i] == "Homo-del")
      #  egfr[i] <- "cadetblue4"
      else if(egfr[i] == "Low-amp")
        egfr[i] <- "pink"
      else if(egfr[i] == "High-amp")
        egfr[i] <- "red"
      else if(egfr[i] == "Neutral")
        egfr[i] <- "gray"
      else
        egfr[i] <- "white"
      
    }
  }
  
  for(i in 1:length(braf)){
    if(is.na(braf[i]))
      braf[i] <- "white"
    else{
      if(braf[i] == "Mut")
        braf[i] <- "black"
      else if(braf[i] == "WT")
        braf[i] <- "grey"
      else
        braf[i] <- "white"
      
    }
  }
  
  
  for(i in 1:length(histology)){
    if(is.na(histology[i]))
      histology[i] <- "white"
    else{
      if(histology[i] == "astrocytoma")
        histology[i] <- "red"
      else if(histology[i] == "glioblastoma")
        histology[i] <- "purple"
      else if(histology[i] == "oligoastrocytoma")
        histology[i] <- "cyan"
      else if(histology[i] == "oligodendroglioma")
        histology[i] <- "green3"
      #else if(histology[i] == "pilocytic astrocytoma")
      # histology[i] <- "deeppink3"
      else
        histology[i] <- "white"
      
    }
  }
  
  for(i in 1:length(grade)){
    if(is.na(grade[i]))
      grade[i] <- "white"
    else{
      if(grade[i] == "G2")
        grade[i] <- "green3"
      else if(grade[i] == "G3")
        grade[i] <- "red"
      else if(grade[i] == "G4")
        grade[i] <- "blue"
      else
        grade[i] <- "white"
      
    }
  }
  
  # for(i in 1:length(gbm.cluster)){
  #  if(is.na(gbm.cluster[i]))
  #   gbm.cluster[i] <- "white"
  #else{
  #      if(gbm.cluster[i] == "G-CIMP")
  #       gbm.cluster[i] <- "yellow"
  #    else if(gbm.cluster[i] == "M1")
  #     gbm.cluster[i] <- "green"
  #  else if(gbm.cluster[i] == "M2")
  #        gbm.cluster[i] <- "black"
  #     else if(gbm.cluster[i] == "M3")
  #      gbm.cluster[i] <- "blue"
  #   else if(gbm.cluster[i] == "M4")
  #    gbm.cluster[i] <- "red"
  #      else if(gbm.cluster[i] == "M6")
  #       gbm.cluster[i] <- "darksalmon"
  #    else 
  #     gbm.cluster[i] <- "white"
  #}
  #}
  
  for(i in 1:length(cnv)){
    if(is.na(cnv[i]))
      cnv[i] <- "white"
    else{ 
      if(cnv[i] == "LGc1")
        cnv[i] <- "red"
      else if(cnv[i] == "LGc2")
        cnv[i] <- "purple"
      else if(cnv[i] == "LGc3")
        cnv[i] <- "cyan"
      else if(cnv[i] == "LGc4")
        cnv[i] <- "magenta"
      else if(cnv[i] == "LGc5")
        cnv[i] <- "orange"
      else if(cluster[i] == "LGc6")
        cnv[i] <- "green3"
      else 
        cnv[i] <- "white"
    }
  }
  
  
  
  for(i in 1:length(gbm.cluster.new)){
    if(is.na(gbm.cluster.new[i]))
      gbm.cluster.new[i] <- "white"
    else{
      if(gbm.cluster.new[i] == "Classical")
        gbm.cluster.new[i] <- "red"
      else if(gbm.cluster.new[i] == "G-CIMP")
        gbm.cluster.new[i] <- "yellow"
      else if(gbm.cluster.new[i] == "Mesenchymal")
        gbm.cluster.new[i] <- "green"
      else if(gbm.cluster.new[i] == "Neural")
        gbm.cluster.new[i] <- "blue"
      else if(gbm.cluster.new[i] == "Proneural")
        gbm.cluster.new[i] <- "orange"
      else 
        gbm.cluster.new[i] <- "white"
    }
  }
  
  #  for(i in 1:length(lgg.cluster)){
  #   if(is.na(lgg.cluster[i]))
  #    lgg.cluster[i] <- "white"
  #  else{
  #   if(lgg.cluster[i] == "M5")
  #    lgg.cluster[i] <- "darkorchid4"
  #  else if(lgg.cluster[i] == "M1")
  #   lgg.cluster[i] <- "black"
  #  else if(lgg.cluster[i] == "M2")
  #      lgg.cluster[i] <- "red"
  #   else if(lgg.cluster[i] == "M3")
  #   lgg.cluster[i] <- "blue"
  #else if(lgg.cluster[i] == "M4")
  #        lgg.cluster[i] <- "forestgreen"
  #     else 
  #      lgg.cluster[i] <- "white"
  # }
  #}
  
  for(i in 1:length(subtype)){
    if(is.na(subtype[i]))
      subtype[i] <- "white"
    else{ 
      if(subtype[i] == "Classical")
        subtype[i] <- "red"
      else if(subtype[i] == "G-CIMP")
        subtype[i] <- "black"
      else if(subtype[i] == "IDHmut-codel")
        subtype[i] <- "cyan"
      else if(subtype[i] == "IDHmut-non-codel")
        subtype[i] <- "tomato"
      else if(subtype[i] == "IDHwt")
        subtype[i] <- "gold"
      else if(subtype[i] == "Mesenchymal")
        subtype[i] <- "green"
      else if(subtype[i] == "Neural")
        subtype[i] <- "blue"
      else if(subtype[i] == "Proneural")
        subtype[i] <- "orange"
      else 
        subtype[i] <- "white"
    }
  }
  
  
  for(i in 1:length(lgg.cluster.new)){
    if(is.na(lgg.cluster.new[i]))
      lgg.cluster.new[i] <- "white"
    else{
      if(lgg.cluster.new[i] == "M5")
        lgg.cluster.new[i] <- "darkorchid4"
      else if(lgg.cluster.new[i] == "M1")
        lgg.cluster.new[i] <- "black"
      else if(lgg.cluster.new[i] == "M2")
        lgg.cluster.new[i] <- "red"
      else if(lgg.cluster.new[i] == "M3")
        lgg.cluster.new[i] <- "blue"
      else if(lgg.cluster.new[i] == "M4")
        lgg.cluster.new[i] <- "forestgreen"
      else 
        lgg.cluster.new[i] <- "white"
    }
  }
  
  
  for(i in 1:length(sturm)){
    if(is.na(sturm[i]))
      sturm[i] <- "white"
    else{
      if(sturm[i] == "IDH")
        sturm[i] <- "red3"
      else if(sturm[i] == "Mesenchymal")
        sturm[i] <- "deeppink"
      else if(sturm[i] == "RTK I 'PDGFRA'")
        sturm[i] <- "tan3"
      else if(sturm[i] == "RTK II 'Classic'")
        sturm[i] <- "darkslateblue"
      else if(sturm[i] == "G34")
        sturm[i] <- "green"
      else if(sturm[i] == "K27")
        sturm[i] <- "cyan"
      else 
        sturm[i] <- "white"
    }
  }
  
  for(i in 1:length(tumor)){
    if(is.na(tumor[i]))
      tumor[i] <- "white"
    else{ 
      if(tumor[i] == "LGG")
        tumor[i] <- "black"
      else if(tumor[i] == "GBM")
        tumor[i] <- "darkgrey"
      else
        tumor[i] <- "khaki3"
    }
    
  }
  
  for(i in 1:length(IDH)){
    if(is.na(IDH[i]))
      IDH[i] <- "white"
    else{
      if(IDH[i] == "Mutant")
        IDH[i] <- "black"
      else if(IDH[i] == "WT")
        IDH[i] <- "grey"
      else
        IDH[i] <- "white"
      
    }
  }
  
  for(i in 1:length(p19q)){
    if(is.na(p19q[i]))
      p19q[i] <- "white"
    else{
      if(p19q[i] == "codel")
        p19q[i] <- "darkblue"
      else if(p19q[i] == "non-codel")
        p19q[i] <- "cyan"
      else
        p19q[i] <- "white"
      
    }
  }
  
  #  for(i in 1:length(idh2)){
  #   if(is.na(idh2[i]))
  #    idh2[i] <- "white"
  # else{ 
  #      if(idh2[i] == "R172G")
  #       idh2[i] <- "black"
  #    else if(idh2[i] == "R172K")
  #     idh2[i] <- "black"
  #      else if(idh2[i] == "R172M")
  #       idh2[i] <- "black"
  #    else if(idh2[i] == "R172S")
  #     idh2[i] <- "black"
  #  else if(idh2[i] == "R172W")
  #   idh2[i] <- "black"
  #      else if(idh2[i] == "WT")
  #       idh2[i] <- "gray"
  #    else 
  #     idh2[i] <- "white"
  #  }
  #}
  
  # for(i in 1:length(idh1)){
  #  if(is.na(idh1[i]))
  #   idh1[i] <- "white"
  # else{ 
  #    if(idh1[i] == "R132C")
  #     idh1[i] <- "black"
  #  else if(idh1[i] == "R132G")
  #        idh1[i] <- "black"
  #     else if(idh1[i] == "R132H")
  #      idh1[i] <- "black"
  #   else if(idh1[i] == "R132S")
  #        idh1[i] <- "black"
  #     else if(idh1[i] == "WT")
  #      idh1[i] <- "gray"
  #      else 
  #       idh1[i] <- "white"
  #  }
  #}
  
  estimates <- jet.colors(100)[ceiling((estimates-min(estimates,na.rm=T))/(max(estimates,na.rm=T)-min(estimates,na.rm=T))*(length(jet.colors(100))-1))+1]
  estimates[is.na(estimates)] <- "white"
  leukocyte <- jet.colors(100)[ceiling((leukocyte-min(leukocyte,na.rm=T))/(max(leukocyte,na.rm=T)-min(leukocyte,na.rm=T))*(length(jet.colors(100))-1))+1]
  leukocyte[is.na(leukocyte)] <- "white"
  age <- jet.colors(100)[ceiling((age-min(age,na.rm=T))/(max(age,na.rm=T)-min(age,na.rm=T))*(length(jet.colors(100))-1))+1]
  age[is.na(age)] <- "white"
  
  for(i in 1:length(rna.dna.c)){
    if(is.na(rna.dna.c[i]))
      rna.dna.c[i] <- "white"
    else{
      if(rna.dna.c[i] == "Classical")
        rna.dna.c[i] <- "red"
      else if(rna.dna.c[i] == "G-CIMP")
        rna.dna.c[i] <- "black"
      else if(rna.dna.c[i] == "Mesenchymal")
        rna.dna.c[i] <- "green"
      else if(rna.dna.c[i] == "Neural")
        rna.dna.c[i] <- "blue"
      else if(rna.dna.c[i] == "Proneural")
        rna.dna.c[i] <- "orange"
      else 
        rna.dna.c[i] <- "white"
    }
  }
  
  
  for(i in 1:length(rna.cluster)){
    if(is.na(rna.cluster[i]))
      rna.cluster[i] <- "white"
    else{
      if(rna.cluster[i] == "LGr1")
        rna.cluster[i] <- "yellow"
      else if(rna.cluster[i] == "LGr2")
        rna.cluster[i] <- "cyan"
      else if(rna.cluster[i] == "LGr3")
        rna.cluster[i] <- "blue"
      else if(rna.cluster[i] == "LGr4")
        rna.cluster[i] <- "red"
      else if(rna.cluster[i] == "unclassified")
        rna.cluster[i] <- "brown"
      else 
        rna.cluster[i] <- "white"
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
  
  
  
  for(i in 1:length(anno)){
    if(is.na(anno[i]))
      anno[i] <- "white"
    else{
      if(anno[i] == "Case submitted is found to be a recurrence after submission")
        anno[i] <- "yellow"
      else if(anno[i] == "Case submitted is found to be a recurrence after submission/Item does not meet study protocol")
        anno[i] <- "yellow"
      else if(anno[i] == "History of unacceptable prior treatment related to a prior/other malignancy")
        anno[i] <- "yellow"
      else if(anno[i] == "Item in special subset")
        anno[i] <- "gray"
      else if(anno[i] == "Molecular analysis outside specification")
        anno[i] <- "gray"
      else if(anno[i] == "Neoadjuvant therapy")
        anno[i] <- "cyan"
      else if(anno[i] == "Neoadjuvant therapy/Case submitted is found to be a recurrence after submission")
        anno[i] <- "yellow"
      else if(anno[i] == "Normal tissue origin incorrect")
        anno[i] <- "gray"
      else if(anno[i] == "Item may not meet study protocol")
        anno[i] <- "gray"
      else if(anno[i] == "Pathology outside specification")
        anno[i] <- "gray"
      else if(anno[i] == "Prior malignancy")
        anno[i] <- "yellow"
      else if(anno[i] == "Prior malignancy/History of acceptable prior treatment related to a prior/other malignancy")
        anno[i] <- "yellow"
      else if(anno[i] == "Qualified in error")
        anno[i] <- "yellow"
      else if(anno[i] == "Qualified in error/Item in special subset")
        anno[i] <- "gray"
      else if(anno[i] == "Subject withdrew consent")
        anno[i] <- "magenta"
      
      else 
        anno[i] <- "white"
    }
  }
  
  
  
  #cluster <- cbind(age,anno,p19q,idh1,idh2,IDH,sturm,rna.dna.c,gbm.cluster.new,gbm.cluster,coc.c,lgg.cluster.new,lgg.cluster,rna.cluster,cluster,new.cl,estimates,leukocyte,tumor)
  cluster <- cbind(age,anno,p19q,IDH,sturm,rna.dna.c,gbm.cluster.new,coc.c,lgg.cluster.new,rna.cluster,cluster,estimates,leukocyte,tumor,cnv,tp53,egfr,braf,histology,grade,subtype)
  colnames(cluster) <- c("Age","Annotation","1p19q codel","IDH","Sturm et al","GBM Clusters [Brennan, 2013]","GBM DNA Meth Predicted Cluster","LGG COC [NEJM, unpublished]","LGG Predicted Clusters", "RNA Clusters","Cluster meth","RNA Leukocyte", "DNA Meth Leukocyte", "Sample Type","Cluster CNV","TP53","EGFR CNV","BRAF","Histology","Grade","Subtype")
  #colnames(cluster) <- c("Age","Annotation","1p19q codel","IDH1","IDH2", "IDH","Sturm et al","GBM Clusters [Brennan, 2013]","GBM DNA Meth Predicted Cluster","GBM DNA Meth Cl [Brennan, 2013]","COC Clusters [NEJM, unpublished]","LGG Predicted Clusters", "LGG Clusters [NEJM, unpublished]", "RNA Clusters","Cluster meth","New Master Cluster","RNA Leukocyte", "DNA Meth Leukocyte", "Sample Type")
  
  return(cluster)
}

source("/dados/ResearchProjects/thais/heatmap.plus.R")

### CpG island (downloaded from Encode Aug, 28th)
cpg <- read.table("/dados/ResearchProjects/thais/TCGA/CpGisland.txt")
library(GenomicRanges)
cpg.is <- GRanges(seqnames = cpg$V1,ranges = IRanges(start = cpg$V2, end=cpg$V3), ID = cpg$V4)
aux <- LGG.GBM.250914[rownames(dat.lgg.gbm.new.noXY.dic.oecg),]
probe <- GRanges(seqnames = paste0("chr",aux$Chromosome),ranges = IRanges(start = aux$Genomic_Coordinate, end=aux$Genomic_Coordinate), ID = aux$Composite.Element.REF)
overlap <- as.data.frame(findOverlaps(probe,cpg.is))
rlab <- matrix("white", nrow = nrow(aux), ncol = 2)
rlab[overlap$queryHits,] <- "green"

#normal.age <- read.table("/dados/ResearchProjects/thais/TCGA/LGG.GBM/normal_age.txt",header=T)

####### Pilocytic Astrocytoma GSE44684
pilocytic.astro.GSE44684 <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Pilocytic_astrocytoma_GSE44684/GSE44684_beta.txt")
#remove the p value columns
index <- grep("Pval",colnames(pilocytic.astro.GSE44684))
pilocytic.astro.GSE44684 <- pilocytic.astro.GSE44684[,-index]
rownames(pilocytic.astro.GSE44684) <- pilocytic.astro.GSE44684[,1]
pilocytic.astro.GSE44684 <- pilocytic.astro.GSE44684[,c(1:38,40:68,39)] #coloca os 6 controles pro final
pilocytic.astro.GSE44684 <- pilocytic.astro.GSE44684[,c(1:62,65:67,63:64,68)] #manter a ordem que esta no GEO para coincidir com os dados de idade
PA.age <- t(read.table("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Pilocytic_astrocytoma_GSE44684/PA_age.txt"))
colnames(pilocytic.astro.GSE44684)[2:68] <- as.character(RF.PA.on.GBM$PA.sampleID)

####### mRNA for GSE15745 (control)
normalBrain_dnameth <- read.delim("/dados/ResearchProjects/thais/TCGA/Normals/mRNA_Brain/07-03-14/GSE15745-GPL8490_series_matrix.txt",skip=96)
normalBrain_dnameth_id <- scan("/dados/ResearchProjects/thais/TCGA/Normals/mRNA_Brain/07-03-14/id_dnameth.txt",sep="\t") 
colnames(normalBrain_dnameth)[2:507] <- as.character(normalBrain_dnameth_id)
rownames(normalBrain_dnameth) <- normalBrain_dnameth$ID_REF

#RF.LGG.new.2 <- RF.LGG.new[!duplicated(RF.LGG.new$TCGAIDnew), ]
#a <- merge(metadata.new,RF.LGG.new,by.x="TCGAID",by.y="ID",all.x=T)

colnames(metadata.new)[265] <- "GBM.Expres.Clusters"
colnames(metadata.new)[270] <- "Master.Clusters"
metadata.new <- metadata.new[, !(colnames(metadata.new) %in% c("leukocyte.strat"))]
metadata.new <- metadata.new[, !(colnames(metadata.new) %in% c("age.strat"))]
a <- metadata.new
a$newLGGcluster <- a$LGG.DNA.Methyl.Clusters
a[as.character(RF.LGG.new$ID),"newLGGcluster"] <- as.character(RF.LGG.new$newlgg.labels)
a$newGBMcluster <- a$GBM.DNA.Methyl.Clusters
a[as.character(RF.GBM.new$ID),"newGBMcluster"] <- as.character(RF.GBM.new$RFtype.all)



#a <- merge(a,RF.GBM.new,by.x="TCGAID",by.y="ID",all.x=T)
#a <- merge(a,RF.LGG.new,by.x="TCGAIDnew",by.y="TCGAIDnew",all.x=T)
#c <- metadata.new
#metadata.new <- a
#rownames(metadata.new) <- metadata.new$TCGAID
#Normal samples 

########## Ordem das amostras Houtan
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/LGG_GBM_Sturm_PA_dnaMethyl-ORDERED.Rda")
#obj full.heatmap.orderd 1129 samples = 136+61+932

PA.age <- t(read.table("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Pilocytic_astrocytoma_GSE44684/PA_age.txt"))

a <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/pdata.516lgg.606gbm.txt")

## novo Random Forest file
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/pa.sturm.RF.s1.rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/all.lgggbm.tested.on.sturm.rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/RF.full.LGG.GBM.PA.Sturm.classification.rda")

pa.sturm.RF.s1.onlyTumor <- subset(pa.sturm.RF.s1, Tumor.Type %in% c("Pilocytic astrocytoma tumor","Primary brain tumor tissue"))
pa.sturm.RF.s1.onlyTumor$mut.BRAF <- as.character(pa.sturm.RF.s1.onlyTumor$PA.mutation.status)
aux <- grepl("BRAF",pa.sturm.RF.s1.onlyTumor$mut.BRAF)
pa.sturm.RF.s1.onlyTumor$mut.BRAF[aux] <- "Mut"
pa.sturm.RF.s1.onlyTumor$mut.BRAF[!aux & !is.na(pa.sturm.RF.s1.onlyTumor$mut.BRAF)] <- "WT"
levels(pa.sturm.RF.s1.onlyTumor$Tumor.Type) <- c("","PA","","","GBM")
levels(pa.sturm.RF.s1.onlyTumor$Gender) <- c("FEMALE","MALE","FEMALE","MALE","FEMALE")
pa.sturm.RF.s1.onlyTumor$RF.GBM.class <- str_replace(pa.sturm.RF.s1.onlyTumor$RF.GBM.class,"GCIMP","G-CIMP")
levels(pa.sturm.RF.s1.onlyTumor$Sturm.IDH.1.Mutation.Status) <- c("Mutant","WT")
levels(pa.sturm.RF.s1.onlyTumor$Sturm.TP53.Mutation.Status) <- c("Mut","WT")
pa.sturm.RF.s1.onlyTumor$Sturm.CDKN2A.Deletion <- str_replace(pa.sturm.RF.s1.onlyTumor$Sturm.CDKN2A.Deletion,"0","Non-Del")
pa.sturm.RF.s1.onlyTumor$Sturm.CDKN2A.Deletion <- str_replace(pa.sturm.RF.s1.onlyTumor$Sturm.CDKN2A.Deletion,"1","Del")
pa.sturm.RF.s1.onlyTumor$Sturm.EGFR.Amplification <- str_replace(pa.sturm.RF.s1.onlyTumor$Sturm.EGFR.Amplification,"0","Non-Amp")
pa.sturm.RF.s1.onlyTumor$Sturm.EGFR.Amplification <- str_replace(pa.sturm.RF.s1.onlyTumor$Sturm.EGFR.Amplification,"1","Amp")
pa.sturm.RF.s1.onlyTumor$Sturm.PDGFRA.Amplification <- str_replace(pa.sturm.RF.s1.onlyTumor$Sturm.PDGFRA.Amplification,"0","Non-Amp")
pa.sturm.RF.s1.onlyTumor$Sturm.PDGFRA.Amplification <- str_replace(pa.sturm.RF.s1.onlyTumor$Sturm.PDGFRA.Amplification,"1","Amp")
pa.sturm.RF.s1.onlyTumor$Sturm.Death <- str_replace(pa.sturm.RF.s1.onlyTumor$Sturm.Death,"0","alive")
pa.sturm.RF.s1.onlyTumor$Sturm.Death <- str_replace(pa.sturm.RF.s1.onlyTumor$Sturm.Death,"1","dead")
levels(a$IDH1.status) <- c("Mutant","Mutant","Mutant","Mutant","WT")
pa.sturm.RF.s1.onlyTumor$Tumor.Type <- as.character(pa.sturm.RF.s1.onlyTumor$Tumor.Type)
a$tumor.type <- as.character(a$tumor.type)
a$gender <- as.character(a$gender)
aux <- matrix(NA, nrow = nrow(pa.sturm.RF.s1.onlyTumor), ncol = ncol(a))
colnames(aux) <- colnames(a)
a <- rbind(a,aux)
a$Study <- "TCGA"
a[1123:1319,"id"] <- as.character(pa.sturm.RF.s1.onlyTumor$nonTCGA.ID)
a[1123:1319,"Study"] <- as.character(pa.sturm.RF.s1.onlyTumor$Study)
a[1123:1319,"tumor.type"] <- as.character(pa.sturm.RF.s1.onlyTumor$Tumor.Type)
a[1123:1319,"gender"] <- as.character(pa.sturm.RF.s1.onlyTumor$Gender)
a[1123:1319,"age"] <- as.character(pa.sturm.RF.s1.onlyTumor$Age.at.diagnosis.years)
a[1123:1319,"cluster.meth.lgg"] <- as.character(pa.sturm.RF.s1.onlyTumor$RF.LGG.class)
a[1123:1319,"cluster.meth.gbm"] <- as.character(pa.sturm.RF.s1.onlyTumor$RF.GBM.class)
a[1123:1319,"IDH1.status"] <- as.character(pa.sturm.RF.s1.onlyTumor$Sturm.IDH.1.Mutation.Status)
a[1123:1319,"mut.TP53"] <- as.character(pa.sturm.RF.s1.onlyTumor$Sturm.TP53.Mutation.Status)
#a[933:1129,"del.CDKN2A"] <- as.character(pa.sturm.RF.s1.onlyTumor$Sturm.CDKN2A.Deletion)
#a[933:1129,"amp.EGFR"] <- as.character(pa.sturm.RF.s1.onlyTumor$Sturm.EGFR.Amplification)
#a[933:1129,"amp.PDGFRA"] <- as.character(pa.sturm.RF.s1.onlyTumor$Sturm.PDGFRA.Amplification)
a[1123:1319,"mut.BRAF"] <- as.character(pa.sturm.RF.s1.onlyTumor$mut.BRAF)
a$PA.BRAF.fusion <- NA
pa.sturm.RF.s1.onlyTumor$Study.Class.original <- as.factor(pa.sturm.RF.s1.onlyTumor$Study.Class.original)
y <- levels(pa.sturm.RF.s1.onlyTumor$Study.Class.original)
levels(pa.sturm.RF.s1.onlyTumor$Study.Class.original) <- c("1-PA", "2-PA", y[3:8])
a$Sturm.PA.merged.clusters <- NA
a[1123:1319,"Sturm.PA.merged.clusters"] <- as.character(pa.sturm.RF.s1.onlyTumor$Study.Class.original)
a[1123:1319,"PA.BRAF.fusion"] <- as.character(pa.sturm.RF.s1.onlyTumor$PA.mutation.status)
a[1123:1319,"vital"] <- as.character(pa.sturm.RF.s1.onlyTumor$Sturm.Death)
a[1123:1319,"os.months"] <- as.character(pa.sturm.RF.s1.onlyTumor$Sturm.OS.months)
#a[,"GBM.Cluster.Sturm.et.al"] <- as.character(a[,"GBM.Cluster.Sturm.et.al"])
pa.sturm.RF.s1.onlyTumor$Study.Class.original <- as.character(pa.sturm.RF.s1.onlyTumor$Study.Class.original)
#a[933:1068,"GBM.Cluster.Sturm.et.al"] <- as.character(pa.sturm.RF.s1.onlyTumor$Study.Class.original[1:136])
a[1123:1258,"histology"] <- "glioblastoma"
a[1123:1258,"grade"] <- "G4"
a[1123:1319,"Methylation.Platform"] <- "450k"
a$histology <- as.character(a$histology)
a$grade <- as.character(a$grade)
a[1259:1319,"histology"] <- "pilocytic astrocytoma"
a[1259:1319,"grade"] <- "G1"
a$PA.original.clusters <- NA
a[1259:1319,"PA.original.clusters"] <- as.character(pa.sturm.RF.s1.onlyTumor$Study.Class.original[137:197])
rownames(a) <- a$id
a$cluster.meth.all <- NA
a[as.character(RF.class$TCGAID),"cluster.meth.all"] <- as.character(RF.class$cluster.meth)
rownames(a) <- a$id

#all.lgggbm.tested.on.sturm <- data.frame(all.lgggbm.tested.on.sturm)
#all.lgggbm.tested.on.sturm$ID <- rownames(all.lgggbm.tested.on.sturm)
#all.lgggbm.tested.on.sturm$all.lgggbm.tested.on.sturm <- str_replace(all.lgggbm.tested.on.sturm$all.lgggbm.tested.on.sturm,"Classical","RTK II 'Classic'")
#all.lgggbm.tested.on.sturm$all.lgggbm.tested.on.sturm <- str_replace(all.lgggbm.tested.on.sturm$all.lgggbm.tested.on.sturm,"PDGFRA","RTK I 'PDGFRA'")
#a$newSturmClass <- a$GBM.Cluster.Sturm.et.al
#a[as.character(all.lgggbm.tested.on.sturm$ID),"newSturmClass"] <- as.character(all.lgggbm.tested.on.sturm$all.lgggbm.tested.on.sturm)
#rownames(a) <- a$TCGAID
## FIM novo Random Forest file

metadata.1129.samples.20141016 <- a
metadata.1129.samples.20141016.annotation <- read.table(file = "/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG-GBM.Annotation.Table-20141016.csv", na.strings="NA", stringsAsFactors=FALSE, header = T, sep = ",")

metadata.1129.samples.20141016[is.na(metadata.1129.samples.20141016$newSturmClass) & metadata.1129.samples.20141016$Study == "TCGA", "newSturmClass"] <- "Unclassified-TCGA-27K"
metadata.1129.samples.20141016$Sturm.PA.merged.clusters[1:932] <- "TCGA"
metadata.1129.samples.20141016$age <- as.numeric(metadata.1129.samples.20141016$age)
metadata.1129.samples.20141016$os.months <- as.numeric(metadata.1129.samples.20141016$os.months)

save(metadata.1129.samples.20141016,metadata.1129.samples.20141016.annotation,file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/metadata_LGG_GBM_PA_Sturm-20141016.Rda")

b <- subset(metadata.1129.samples.20141016, Study == "TCGA")
b <- subset(b, !is.na(newSturmClass))
c <- cbind(b$TCGAID, paste(b$newSturmClass, b$cluster.meth2, sep=" - "))
write.table(metadata.1129.samples.20141016[1:932,c("TCGAID","newLGGcluster","newGBMcluster")],file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_GBM_classified.as.Sturm.membership.txt",quote=F,row.names=F,sep="\t")

d <- subset(metadata.1129.samples.20141016, Study == "TCGA" & newSturmClass=="Mesenchymal" & cluster.meth.all == "LGm6")
e <- subset(metadata.1129.samples.20141016, Study == "TCGA" & newSturmClass=="Mesenchymal" & cluster.meth.all == "LGm5")


load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/gbm450.Rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/lgg.Rda")
sturm.probes <- read.table("/dados/ResearchProjects/thais/TCGA/LGG.GBM/hccpgs_order.txt")

c <- subset(metadata.1129.samples.20141016, Methylation.Platform == "450k")
c <- subset(c, Study != "PA")
lgg.gbm.450 <- merge(lgg.patients[,c(1,5:520)],gbm.patients[,c(1,5:133)],by="Composite.Element.REF")
lgg.gbm.sturm.450 <- merge(lgg.gbm.450,sturm.samples.nocontrol,by.x="Composite.Element.REF",by.y="ID_REF")
rownames(lgg.gbm.sturm.450) <- lgg.gbm.sturm.450$Composite.Element.REF
lgg.gbm.sturm.450.sel <- lgg.gbm.sturm.450[as.character(sturm.probes$V1),]

sturm.samples.nocontrol <- sturm.samples[,-c(77:82)]
rownames(sturm.samples.nocontrol) <- sturm.samples.nocontrol$ID_REF
colnames(sturm.samples.nocontrol)[2:137] <- as.character(rownames(RF.Sturm.on.GBM)[c(1:75,82:142)])

aux <- subset(c, newSturmClass == "K27")
k27 <- lgg.gbm.sturm.450.sel[,as.character(aux$TCGAID)]
hc.k27 <- hclust(dist(t(k27)))
aux <- subset(c, newSturmClass == "G34")
g34 <- lgg.gbm.sturm.450.sel[,as.character(aux$TCGAID)]
hc.g34 <- hclust(dist(t(g34)))
aux <- subset(c, newSturmClass == "Mesenchymal")
mesenchymal <- lgg.gbm.sturm.450.sel[,as.character(aux$TCGAID)]
hc.mesenchymal <- hclust(dist(t(mesenchymal)))
aux <- subset(c, newSturmClass == "RTK I 'PDGFRA'")
pdgfra <- lgg.gbm.sturm.450.sel[,as.character(aux$TCGAID)]
hc.pdgfra <- hclust(dist(t(pdgfra)))
aux <- subset(c, newSturmClass == "RTK II 'Classic'")
classic <- lgg.gbm.sturm.450.sel[,as.character(aux$TCGAID)]
hc.classic <- hclust(dist(t(classic)))
aux <- subset(c, newSturmClass == "IDH")
idh <- lgg.gbm.sturm.450.sel[,as.character(aux$TCGAID)]
hc.idh <- hclust(dist(t(idh)))

lgg.gbm.sturm.order <- cbind(idh[,hc.idh$order],k27[,hc.k27$order],g34[,hc.g34$order],pdgfra[,hc.pdgfra$order],mesenchymal[,hc.mesenchymal$order],classic[,hc.classic$order])

cores <- c[colnames(lgg.gbm.sturm.order),c("newSturmClass","cluster.meth.all","Study")]
cores$newSturmClass <- as.factor(cores$newSturmClass)
cores$cluster.meth.all <- as.factor(cores$cluster.meth.all)
cores$Study <- as.factor(cores$Study)
levels(cores$newSturmClass) <- c("blue","red3","green3","darkgoldenrod3","darkorange3","brown4")
levels(cores$cluster.meth.all) <- c("green","red","purple","orange","yellow","blue")
levels(cores$Study) <- c("#fa7f71","tan")


png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/testeSturm.png",res=300,width=2000,height=2000)


heatmap.plus.sm(as.matrix(lgg.gbm.sturm.order[rev(1:8000),]),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = as.matrix(cores)
                #RowSideColors = rlab
                
)
dev.off()

######################### OLD

RF.Sturm.on.GBM$ID <- rownames(RF.Sturm.on.GBM)
RF.Sturm.on.GBM.onlyTumor <- RF.Sturm.on.GBM[c(1:75,82:142),] #Tumor Samples
RF.Sturm.on.LGG.onlyTumor <- RF.Sturm.on.LGG[c(1:75,82:142),] #Tumor Samples
RF.Sturm.on.GBM.onlyTumor$rftype <- str_replace(RF.Sturm.on.GBM.onlyTumor$rftype,"GCIMP","G-CIMP")
aux <- matrix(NA, nrow = nrow(RF.Sturm.on.GBM.onlyTumor), ncol = ncol(a))
colnames(aux) <- colnames(a)
a <- rbind(a,aux)
a[933:1068,"TCGAID"] <- as.character(RF.Sturm.on.GBM.onlyTumor$ID)
a[933:1068,"newGBMcluster"] <- as.character(RF.Sturm.on.GBM.onlyTumor$rftype)
a[933:1068,"newLGGcluster"] <- as.character(RF.Sturm.on.LGG.onlyTumor$Sturm.RFtype)
a[933:1068,"tumor.type"] <- "GBM"
levels(RF.Sturm.on.LGG.onlyTumor$Sturm.characteristics_ch1.3) <- c("","FEMALE","MALE")
a[933:1068,"age"] <- as.character(RF.Sturm.on.LGG.onlyTumor$Sturm.Age)
a[933:1068,"gender"] <- as.character(RF.Sturm.on.LGG.onlyTumor$Sturm.characteristics_ch1.3)
RF.PA.on.GBM.onlyTumor <- RF.PA.on.GBM[1:61,]
RF.PA.on.LGG.onlyTumor <- RF.PA.on.GBM[1:61,]
aux <- matrix(NA, nrow = nrow(RF.PA.on.GBM.onlyTumor), ncol = ncol(a))
colnames(aux) <- colnames(a)
a <- rbind(a,aux)
a[1069:1129,"TCGAID"] <- as.character(RF.PA.on.GBM.onlyTumor$PA.sampleID)

a[1069:1129,"newGBMcluster"] <- as.character(RF.PA.on.GBM.onlyTumor$PA.cluster)
levels(RF.PA.on.GBM.onlyTumor$PA.characteristics_ch1.3) <- c("FEMALE","MALE")
a[1069:1129,"gender"] <- as.character(RF.PA.on.GBM.onlyTumor$PA.characteristics_ch1.3)
a[1069:1129,"newLGGcluster"] <- as.character(RF.PA.on.LGG.onlyTumor$PA.cluster)
a$tumor.type <- as.character(a$tumor.type)
a[1069:1129,"tumor.type"] <- "PA"
a[1069:1129,"age"] <- PA.age[1:61]
sturm.et.al <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Sturm_Cluster.txt")
a$GBM.Cluster.Sturm.et.al <- as.character(a$GBM.Cluster.Sturm.et.al)
a[933:1068,"GBM.Cluster.Sturm.et.al"] <- as.character(sturm.et.al$Cluster)[1:136]
rownames(a) <- a$TCGAID
#LGm clusters para todas as amostras (PA e Sturm)
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/RF.full.LGG.GBM.PA.Sturm.classification.rda")
a$cluster.meth.all <- NA
a[as.character(RF.class$TCGAID),"cluster.meth.all"] <- as.character(RF.class$cluster.meth)
rownames(a) <- a$TCGAID

a <- a[colnames(full.heatmap.orderd),]

load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/all.lgggbm.tested.on.sturm.rda")
all.lgggbm.tested.on.sturm <- data.frame(all.lgggbm.tested.on.sturm)
all.lgggbm.tested.on.sturm$ID <- rownames(all.lgggbm.tested.on.sturm)
all.lgggbm.tested.on.sturm$all.lgggbm.tested.on.sturm <- str_replace(all.lgggbm.tested.on.sturm$all.lgggbm.tested.on.sturm,"Classical","RTK II 'Classic'")
all.lgggbm.tested.on.sturm$all.lgggbm.tested.on.sturm <- str_replace(all.lgggbm.tested.on.sturm$all.lgggbm.tested.on.sturm,"PDGFRA","RTK I 'PDGFRA'")
a$newSturmClass <- a$GBM.Cluster.Sturm.et.al
a[as.character(all.lgggbm.tested.on.sturm$ID),"newSturmClass"] <- as.character(all.lgggbm.tested.on.sturm$all.lgggbm.tested.on.sturm)
metadata.1129.samples <- a
#save(metadata.1129.samples,file="metadata_LGG_GBM_PA_Sturm.Rda")

############################ FIM OLD

## adicionar age, gender, mimimi

clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth.all),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.DNA.Methyl.Clusters),as.matrix(a$cluster.mRNA),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$Master.Clusters),as.matrix(a$Annotation.Category))


normals.sel <- norms[rownames(full.heatmap.orderd),]
pilocytic.astro.GSE44684.order <- pilocytic.astro.GSE44684[rownames(full.heatmap.orderd),]
normalBrain_dnameth.order <- normalBrain_dnameth[rownames(full.heatmap.orderd),]
colnames(normalBrain_dnameth.order)[507] <- "UMARY.1648.TCTX"

normalBrain_dnameth.order <- normalBrain_dnameth.order[,-1]
crbl.mean <- apply(as.matrix(normalBrain_dnameth.order[,1:121]),1,mean,na.rm=T)
frctx.mean <- apply(as.matrix(normalBrain_dnameth.order[,122:254]),1,mean,na.rm=T)
pons.mean <- apply(as.matrix(normalBrain_dnameth.order[,255:380]),1,mean,na.rm=T)
tpctx.mean <- apply(as.matrix(normalBrain_dnameth.order[,381:506]),1,mean,na.rm=T)
regions.mean <- cbind(crbl.mean,frctx.mean,pons.mean,tpctx.mean)

norm.label <- matrix("white", nrow = 77, ncol = 15)
norm.label[,15] <- "darkblue"
region.label <- matrix("white", nrow = 4, ncol = 15)
region.label[1,15] <- "aquamarine" #cerebelo
region.label[2,15] <- "pink" #frontal cortex
region.label[3,15] <- "deeppink" #pons 
region.label[4,15] <- "brown" #temporal contex

clab.6 <- rbind(norm.label,clab.6)
clab.6 <- rbind(region.label,clab.6)

full.heatmap.orderd.normal <- cbind(normals.sel,full.heatmap.orderd)
full.heatmap.orderd.normal <- cbind(regions.mean,full.heatmap.orderd.normal)

png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/teste8.png",res=200,width=2000,height=2000)

a <- metadata.1129.samples.20141016[1:932,]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth.all),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.DNA.Methyl.Clusters),as.matrix(a$cluster.mRNA),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$Master.Clusters),as.matrix(a$Annotation.Category))

LGG.GBM.new.noXY.s <- LGG.GBM.new.noXY[rownames(dat.lgg.gbm.new.noXY.dic.oecg) ,5:936]
LGG.GBM.new.order <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
a <- metadata.1129.samples.20141016[colnames(LGG.GBM.new.order),]
LGG.GBM.new.order <- LGG.GBM.new.order[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/teste.png",res=200,width=2000,height=2000)
heatmap.plus.sm(as.matrix(LGG.GBM.new.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                #Rowv = NA,
                ColSideColors = clab.6
                #RowSideColors = rlab
                
)
dev.off()

library(Biobase)
load(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/normals.gse41826.rda")
norms <- exprs(controls.rm)

LGG.GBM.new.noXY.s <- LGG.GBM.new.noXY[rownames(dat.lgg.gbm.new.noXY.dic.oecg) ,5:936]
LGG.GBM.new.order <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
metadata.new <- metadata.new[colnames(LGG.GBM.new.order),]
clab.6 <- colors.label(as.matrix(metadata.new$newLGGcluster),as.matrix(metadata.new$newGBMcluster),as.matrix(metadata.new$codel1p19q),as.matrix(metadata.new$IDH1.status),as.matrix(metadata.new$IDH2.status),as.matrix(metadata.new$age),as.matrix(metadata.new$cluster.meth2),as.matrix(metadata.new$GBM.DNA.Methyl.Clusters),as.matrix(metadata.new$LGG.DNA.Methyl.Clusters),as.matrix(metadata.new$tumor.type),as.matrix(metadata.new$IDH.status),as.matrix(metadata.new$GBM.Clusters),as.matrix(metadata.new$cluster.mRNA),as.matrix(metadata.new$LGG.COC.Clusters),as.matrix(metadata.new$leukocyte.DNA.methyl),as.matrix(metadata.new$leukocyte.RNA),as.matrix(metadata.new$GBM.Cluster.Sturm.et.al),as.matrix(metadata.new$new.cluster),as.matrix(metadata.new$Annotation.Category))

clab.6 <- clab.6[c(1:63,64:327,405:531,532:682,683:932,328:404),]
#metadata.new <- metadata.new[c(1:63,64:327,405:531,532:682,683:932,328:404),]
LGG.GBM.new.order <- LGG.GBM.new.order[,c(1:63,64:327,405:531,532:682,683:932,328:404)]

sturm.order <- sturm.samples.nocontrol[rownames(LGG.GBM.new.order),2:137]
normals.sel <- norms[rownames(LGG.GBM.new.order),]
pilocytic.astro.GSE44684.order <- pilocytic.astro.GSE44684[rownames(LGG.GBM.new.order),]
normalBrain_dnameth.order <- normalBrain_dnameth[rownames(LGG.GBM.new.order),]
colnames(normalBrain_dnameth.order)[507] <- "UMARY.1648.TCTX"

normalBrain_dnameth.order <- normalBrain_dnameth.order[,-1]
crbl.mean <- apply(as.matrix(normalBrain_dnameth.order[,1:121]),1,mean,na.rm=T)
frctx.mean <- apply(as.matrix(normalBrain_dnameth.order[,122:254]),1,mean,na.rm=T)
pons.mean <- apply(as.matrix(normalBrain_dnameth.order[,255:380]),1,mean,na.rm=T)
tpctx.mean <- apply(as.matrix(normalBrain_dnameth.order[,381:506]),1,mean,na.rm=T)
regions.mean <- cbind(crbl.mean,frctx.mean,pons.mean,tpctx.mean)

LGG.GBM.new.order <- cbind(normals.sel,LGG.GBM.new.order)
LGG.GBM.new.order <- cbind(LGG.GBM.new.order,sturm.order)
LGG.GBM.new.order <- cbind(LGG.GBM.new.order,pilocytic.astro.GSE44684.order[,-c(1,63:68)]) #remove the control samples from PA
LGG.GBM.new.order <- cbind(regions.mean,LGG.GBM.new.order)

sturm.label <- matrix("white", nrow = 136, ncol = 19)
sturm.label[,19] <- "darkgray"

norm.label <- matrix("white", nrow = 77, ncol = 19)
#norm.label[,1] <- jet.colors(100)[ceiling((as.numeric(normal.age[,colnames(normals.sel)])-min(as.numeric(normal.age[,colnames(normals.sel)]),na.rm=T))/(max(as.numeric(normal.age[,colnames(normals.sel)]),na.rm=T)-min(as.numeric(normal.age[,colnames(normals.sel)]),na.rm=T))*(length(jet.colors(100))-1))+1]
norm.label[,19] <- "darkblue"
region.label <- matrix("white", nrow = 4, ncol = 19)
region.label[1,19] <- "aquamarine" #cerebelo
region.label[2,19] <- "pink" #frontal cortex
region.label[3,19] <- "deeppink" #pons 
region.label[4,19] <- "brown" #temporal contex

clab.6 <- rbind(norm.label,clab.6)
pa.label <- matrix("white", nrow = 61, ncol = 19)
pa.label[,19] <- "khaki3" #tumor type
#pa.label[,1] <- jet.colors(100)[ceiling((PA.age[1:61]-min(PA.age[1:61],na.rm=T))/(max(PA.age[1:61],na.rm=T)-min(PA.age[1:61],na.rm=T))*(length(jet.colors(100))-1))+1]
clab.6 <- rbind(clab.6,sturm.label)
clab.6 <- rbind(clab.6,pa.label)
clab.6 <- rbind(region.label,clab.6)
#clab.6 <- cbind(clab.6[,1:8],"white",clab.6[,9:19])
clab.6[1150:1210,9] <- "green" #RF.PA.on.GBM$PA.cluster
#clab.6 <- cbind(clab.6[,1:9],"white",clab.6[,10:20])
aux <- RF.Sturm.on.GBM
levels(aux$rftype) <- c("yellow","green","black","blue","red","darksalmon") #M1 e M4
clab.6[1014:1149,9] <-  as.character(aux$rftype[c(1:75,82:142)])


#clab.6 <- cbind(clab.6[,1:13],"white",clab.6[,14:21])
aux <- RF.PA.on.LGG
levels(aux$PA.cluster) <- c("black","forestgreen") #M1 e M4
clab.6[1150:1210,12] <- as.character(aux$PA.cluster[1:61])

#clab.6 <- cbind(clab.6[,1:14],"white",clab.6[,15:22])
aux <- RF.Sturm.on.LGG
levels(aux$Sturm.RFtype) <- c("black","red","blue","forestgreen","purple") #M1 e M4
clab.6[1014:1149,12] <-  as.character(aux$Sturm.RFtype[c(1:75,82:142)])
#colnames(clab.6) <- 1:13
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/teste5.png",res=200,width=2000,height=2000)


heatmap.plus.sm(as.matrix(LGG.GBM.new.order[order$rowInd,]),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab.6,
                RowSideColors = rlab
                
)
dev.off()

## plot age at diagnosis
library(ggplot2)
x <- subset(metadata.1129.samples.20141030, Study == "TCGA")
#x$K7 <- factor(x$K7,c("3","1","5","7","6","2","4")) #change the order
x$K6 <- factor(x$cluster.meth) #change the order
ages <- x
#ages <- na.omit(ages)
p <- ggplot(ages, aes(factor(K6), age))
png("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/age2.png",res=200,width=4000,height=2000)
p + geom_boxplot(aes(fill = factor(K6)), notchwidth=0.25) + 
  geom_jitter(height = 0, position = position_jitter(width = .1), size=2)  + 
  scale_fill_manual(values=c("LGm1"="green","red","purple", "orange", "yellow","blue"
  ), labels = c("LGm1", "LGm2", "LGm3","LGm4", "LGm5", "LGm6"
  )) +
  #scale_fill_manual(values = rev(c("darkorchid","darkgreen","blue","red","DARKGREY"))) +
  ylab("Age At Diagnosis") +
  xlab("DNA Methylation Clusters") + 
  facet_grid(. ~ tumor.type) +
  labs(title = "DNA Methylation Clusters by Age At Diagnosis", fill = "DNA Meth. Cluster")# +
  #theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white"), panel.border = element_rect(colour = "white"), legend.position="none")
dev.off()

t.test(metadata.1129.samples.20141030[metadata.1129.samples.20141030$cluster.meth %in% "LGm6" & metadata.1129.samples.20141030$tumor.type %in% "LGG","age"],metadata.1129.samples.20141030[metadata.1129.samples.20141030$cluster.meth %in% "LGm6" & metadata.1129.samples.20141030$tumor.type %in% "GBM","age"])


### Gene expression: RNA-seq + microarrray
mergedRNA <- read.table("/dados/ResearchProjects/thais/TCGA/LGG.GBM/MergedData-allGenes.csv")

a <- metadata.1129.samples.20141016[,c("Annotation.Category","cluster.meth2","tumor.type")]
a <- na.omit(a)
a$cluster.meth2 <- factor(a$cluster.meth2)
a$Annotation.Category <- factor(a$Annotation.Category)
a$tumor.type <- factor(a$tumor.type)
library(scales)
library(reshape)
tcga <- subset(metadata.1129.samples.20141030, Study == "TCGA")
MeanTCGA <- melt(tcga[,c(49:151,269)],id="cluster.meth2")


b <- subset(metadata.1129.samples.20141030, Study == "TCGA")
b <- subset(b, !is.na(newSturmClass))
b$ID <- paste(b$newSturmClass, b$cluster.meth, sep=" - ")
#b <- subset(b, cluster.meth2 == "LGm6")
Meanb <- melt(b[,c(93:148,8,241,242)],id=c("cluster.meth","newSturmClass","ID"))
Meanb$gene.name <- substr(Meanb$variable, 5, length(Meanb$variable))
#png(filename = "Mut_LGm5.png", res=200, width=5000, height=3000)
meanb.naomit <- na.omit(Meanb)
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
ggsave(p, filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Mut_AllOrder.png", width=30, height=20, dpi=200, units="in")


LGG.GBM.new.noXY.s <- LGG.GBM.new.noXY[rownames(dat.lgg.gbm.new.noXY.dic.oecg) ,5:936]
LGG.GBM.new.order <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
LGG.GBM.new.order <- LGG.GBM.new.order[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
heatmap.mut <- metadata.1129.samples.20141030[,names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"High-amp"], decreasing=F))]
heatmap.mut <- heatmap.mut[substr(colnames(LGG.GBM.new.order),1,12),]
heatmap.mut <- as.matrix(sapply(heatmap.mut, as.numeric))
rownames(heatmap.mut) <- substr(colnames(LGG.GBM.new.order),1,12)
heatmap.mut <- na.omit(heatmap.mut)
a <- metadata.1129.samples.20141030[rownames(heatmap.mut),]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$Master.cluster),as.matrix(a$flag.categories))
clab.6 <- clab.6[,-c(4,13,14)]

#### order by Hclust within each LGm
mut.naomit <- na.omit(metadata.1129.samples.20141030[,c(names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"High-amp"], decreasing=F)),"cluster.meth")])
aux <- subset(mut.naomit, cluster.meth == "LGm1")
lgm1 <- mut.naomit[as.character(rownames(aux)),names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"High-amp"], decreasing=F))]
lgm1 <- as.matrix(sapply(lgm1, as.numeric))
lgm1 <- lgm1 == 1
storage.mode(lgm1) <- "numeric"
rownames(lgm1) <- rownames(aux)
hc.lgm1 <- hclust(dist(lgm1,"binary"),"ward")
aux <- subset(mut.naomit, cluster.meth == "LGm2")
lgm2 <- mut.naomit[as.character(rownames(aux)),names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"High-amp"], decreasing=F))]
lgm2 <- as.matrix(sapply(lgm2, as.numeric))
lgm2 <- lgm2 == 1
storage.mode(lgm2) <- "numeric"
rownames(lgm2) <- rownames(aux)
hc.lgm2 <- hclust(dist(lgm2,method="binary"),method="ward")
aux <- subset(mut.naomit, cluster.meth == "LGm3")
lgm3 <- mut.naomit[as.character(rownames(aux)),names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"High-amp"], decreasing=F))]
lgm3 <- as.matrix(sapply(lgm3, as.numeric))
lgm3 <- lgm3 == 1
storage.mode(lgm3) <- "numeric"
rownames(lgm3) <- rownames(aux)
hc.lgm3 <- hclust(dist(lgm3,method="binary"),method="ward")
aux <- subset(mut.naomit, cluster.meth == "LGm4")
lgm4 <- mut.naomit[as.character(rownames(aux)),names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"High-amp"], decreasing=F))]
lgm4 <- as.matrix(sapply(lgm4, as.numeric))
lgm4 <- lgm4 == 1
storage.mode(lgm4) <- "numeric"
rownames(lgm4) <- rownames(aux)
hc.lgm4 <- hclust(dist(lgm4,method="binary"),method="ward")
aux <- subset(mut.naomit, cluster.meth == "LGm5")
lgm5 <- mut.naomit[as.character(rownames(aux)),names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"High-amp"], decreasing=F))]
lgm5 <- as.matrix(sapply(lgm5, as.numeric))
lgm5 <- lgm5 == 1
storage.mode(lgm5) <- "numeric"
rownames(lgm5) <- rownames(aux)
hc.lgm5 <- hclust(dist(lgm5,method="binary"),method="ward")
aux <- subset(mut.naomit, cluster.meth == "LGm6")
lgm6 <- mut.naomit[as.character(rownames(aux)),names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"High-amp"], decreasing=F))]
lgm6 <- as.matrix(sapply(lgm6, as.numeric))
lgm6 <- lgm6 == 1
storage.mode(lgm6) <- "numeric"
rownames(lgm6) <- rownames(aux)
hc.lgm6 <- hclust(dist(lgm6,method="binary"),method="ward")


heatmap.mut <- rbind(lgm1[hc.lgm1$order,],lgm2[hc.lgm2$order,],lgm3[hc.lgm3$order,],lgm4[hc.lgm4$order,],lgm5[hc.lgm5$order,],lgm6[hc.lgm6$order,])
a <- metadata.1129.samples.20141030[rownames(heatmap.mut),]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$Master.cluster),as.matrix(a$flag.categories))
clab.6 <- clab.6[,-c(4,13,14)]

#### order the gene frequency by cluster
b <- subset(metadata.1129.samples.20141030, Study == "TCGA")
b <- subset(b, !is.na(newSturmClass) & !is.na(cluster.meth))
b$ID <- paste(b$newSturmClass, b$cluster.meth, sep=" - ")
#b <- subset(b, cluster.meth2 == "LGm6")
Meanb <- melt(b[,c(53:92,8,241,242)],id=c("cluster.meth","newSturmClass","ID"))
Meanb$gene.name <- substr(Meanb$variable, 5, length(Meanb$variable))
#png(filename = "Mut_LGm5.png", res=200, width=5000, height=3000)
meanb.naomit <- na.omit(Meanb)
mut.naomit <- na.omit(metadata.1129.samples.20141030[,c(names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"Mut"], decreasing=F)),"cluster.meth")])
aux <- subset(b, cluster.meth == "LGm1")
c <- subset(meanb.naomit,cluster.meth == "LGm1")
#c <- meanb.naomit
heatmap.mut <- rbind(lgm1[hc.lgm1$order,names(sort(table(c$variable, c$value)[,"High-amp"], decreasing=F))],lgm2[hc.lgm2$order,names(sort(table(c$variable, c$value)[,"High-amp"], decreasing=F))],lgm3[hc.lgm3$order,names(sort(table(c$variable, c$value)[,"High-amp"], decreasing=F))],lgm4[hc.lgm4$order,names(sort(table(c$variable, c$value)[,"High-amp"], decreasing=F))],lgm5[hc.lgm5$order,names(sort(table(c$variable, c$value)[,"High-amp"], decreasing=F))],lgm6[hc.lgm6$order,names(sort(table(c$variable, c$value)[,"High-amp"], decreasing=F))])
a <- metadata.1129.samples.20141030[rownames(heatmap.mut),]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$Master.cluster),as.matrix(a$flag.categories))
clab.6 <- clab.6[,-c(4,13,14)]

### lgm6
heatmap.mut <- lgm6[hc.lgm6$order,names(sort(table(c$variable, c$value)[,"Mut"], decreasing=F))]
a <- metadata.1129.samples.20141030[rownames(heatmap.mut),]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$Master.cluster),as.matrix(a$flag.categories))
clab.6 <- clab.6[,-c(4,13,14)]
#meanb.naomit$variable <- factor(meanb.naomit$variable, levels=names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"Mut"], decreasing=T)))


### sem ordenar 
heatmap.mut <- rbind(lgm1[,names(sort(table(c$variable, c$value)[,"High-amp"], decreasing=F))],lgm2[,names(sort(table(c$variable, c$value)[,"High-amp"], decreasing=F))],lgm3[,names(sort(table(c$variable, c$value)[,"High-amp"], decreasing=F))],lgm4[,names(sort(table(c$variable, c$value)[,"High-amp"], decreasing=F))],lgm5[,names(sort(table(c$variable, c$value)[,"High-amp"], decreasing=F))],lgm6[,names(sort(table(c$variable, c$value)[,"High-amp"], decreasing=F))])
a <- metadata.1129.samples.20141030[rownames(heatmap.mut),]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$Master.cluster),as.matrix(a$flag.categories))
clab.6 <- clab.6[,-c(4,13,14)]

## ordenar caracter
library(cluster)
LGG.GBM.new.noXY.s <- LGG.GBM.new.noXY[rownames(dat.lgg.gbm.new.noXY.dic.oecg) ,5:936]
LGG.GBM.new.order <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
LGG.GBM.new.order <- LGG.GBM.new.order[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
mut.naomit <- metadata.1129.samples.20141030[substr(colnames(LGG.GBM.new.order),1,12),93:148]
mut.naomit <- na.omit(mut.naomit)
a <- as.hclust(agnes(daisy(data.frame(t(mut.naomit)), metric="gower"), method="ward"))
heatmap.mut <- mut.naomit[,a$order]
heatmap.mut <- as.matrix(sapply(heatmap.mut, as.numeric))
rownames(heatmap.mut) <- rownames(mut.naomit)

a <- metadata.1129.samples.20141030[rownames(heatmap.mut),]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$Master.cluster),as.matrix(a$flag.categories))
clab.6 <- clab.6[,-c(4,13,14)]


png(filename = "Heatmap_CNV.png", res=300, width=5000, height=3000)
heatmap.plus.sm(as.matrix(t(heatmap.mut)),
                col = c("darkblue","red4","deepskyblue","red","lightgray"), 
                scale = "none",
                labRow = substr(colnames(heatmap.mut),5,length(colnames(heatmap.mut))),
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab.6)
dev.off()

##### Melt Mutation Status + CNV
coc <- read.csv("CoCassigments845Samples3Cl_v2.csv")
rownames(coc) <- coc$Sample
b <- merge(metadata.1129.samples.20141030,coc[,2:3],by.x="id",by.y="Sample",all.x=T)
m <- subset(metadata.1129.samples.20141030, Study == "TCGA")
m <- m[,c(1,53:92)]
c <- subset(metadata.1129.samples.20141030, Study == "TCGA")
c <- c[,c(1,93:148)]
labels <- subset(b, Study == "TCGA")
labels <- labels[,c("id","CoCluster","cluster.meth","cluster.expr","cluster.cnv")]
aux <- matrix(NA, nrow=932, ncol = ncol(labels))
colnames(aux) <- colnames(labels)
aux[,1] <- as.character(labels$id)
labels <- rbind(labels,aux)
colnames(m)[2:41] <- substr(colnames(m)[2:41],5,length(colnames(m)[2:41]))
colnames(c)[2:57] <- substr(colnames(c)[2:57],5,length(colnames(c)[2:57]))
g.both <- intersect(colnames(m)[2:41], colnames(c)[2:57])
both <- rbind(m[,c("id",g.both)],c[,c("id",g.both)])
m.s <- m[,c("id",setdiff(colnames(m)[2:41], colnames(c)[2:57]))]
c.s <- c[,c("id",setdiff(colnames(c)[2:57], colnames(m)[2:41]))]
aux <- matrix(NA, nrow=932, ncol = ncol(m.s))
colnames(aux) <- colnames(m.s)
aux[,1] <- as.character(m.s$id)
m.s <- rbind(m.s, aux)
aux <- matrix(NA, nrow=932, ncol = ncol(c.s))
colnames(aux) <- colnames(c.s)
aux[,1] <- as.character(c.s$id)
c.s <- rbind(c.s, aux)

both.na <- both
both.na <- as.data.frame(lapply(both.na,function(z){z <- as.character(z); z[is.na(z)] <- "NA"; return(z)}))
hc.both <- as.hclust(agnes(daisy(data.frame(t(both.na[,-1])), metric="gower"), method="ward"))
both <- both[,c(1,hc.both$order+1)]

m.s.na <- m.s
m.s.na <- as.data.frame(lapply(m.s.na,function(z){z <- as.character(z); z[is.na(z)] <- "NA"; return(z)}))
hc.m.s <- as.hclust(agnes(daisy(data.frame(t(m.s.na[,-1])), metric="gower"), method="ward"))
m.s <- m.s[,c(1,hc.m.s$order+1)]

c.s.na <- c.s
c.s.na <- as.data.frame(lapply(c.s.na,function(z){z <- as.character(z); z[is.na(z)] <- "NA"; return(z)}))
hc.c.s <- as.hclust(agnes(daisy(data.frame(t(c.s.na[,-1])), metric="gower"), method="ward"))
c.s <- c.s[,c(1,hc.c.s$order+1)]

colnames(labels) <- c("id","LGcc","LGm","LGr","LGc")
labels <- labels[,c(1,2,5,4,3)]
all <- cbind(labels,both,m.s[,-1],c.s[,-1])

teste <- melt(all,id="id")
teste1 <- merge(na.omit(teste),metadata.1129.samples.20141030,by.x="id",by.y="id")
#teste2 <- subset(teste1, cluster.meth == "LGm5")
#teste2 <- teste2[teste2$variable %in% colnames(both)[2:9],]
#teste1 <- na.omit(teste1)
teste2 <- teste1
teste2$SomaticType <- NA
teste2[teste2$variable %in% names(labels)[-1],"SomaticType"] <- "Clust"
teste2[teste2$variable %in% names(both)[-1],"SomaticType"] <- "CNV & Mut"
teste2[teste2$variable %in% names(m.s)[-1],"SomaticType"] <- "Mut"
teste2[teste2$variable %in% names(c.s)[-1],"SomaticType"] <- "CNV"
teste2$SomaticType <- factor(teste2$SomaticType)
teste2$size <- factor(teste2$value)
#levels(teste2$size) <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1, 1, .75, .75, .50, .50, .50)
levels(teste2$size) <- c(.75, .50, .75, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,.50,1, .75, 1,1 )
teste3 <- subset(teste2, cluster.meth == "LGm1")
teste2$id <- factor(teste2$id)
#LGG.GBM.new.order <- LGG.GBM.new.order[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
ID_order <- as.character(substr(colnames(LGG.GBM.new.order),1,12))
#LGG.GBM.new.order <- LGG.GBM.new.order[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
teste2$id <- factor(teste2$id, levels = ID_order)
ID_order <- c("Clust","CNV","CNV & Mut","Mut")
teste2$SomaticType <- factor(teste2$SomaticType, levels = ID_order)
png(filename = "tileOrderedByCluster2.png", bg="white", res=200, width=4000, height=2000)
#pdf("tileOrderedByCluster2.pdf",width=30, height=15)
ggplot(teste2, aes(y=factor(variable),x=id)) + 
  geom_tile(aes(fill=factor(value), height=as.numeric(as.character(size)))) +
  scale_fill_manual(values=c("cadetblue1","brown4","cadetblue4","red","purple","cyan","magenta","orange","green3","green3","red","black","green","red","purple","orange","yellow","blue","yellow","gray","cyan","magenta","blue","green3","red","brown1","black","cornsilk2","gray","white"),guide = guide_legend(title = "Legend")) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white"),axis.text.x = element_blank(),axis.text.y = element_text(size=8),panel.margin = unit(0, "lines")) + #, strip.background=element_rect(fill=c("green","red","purple","orange","yellow","blue")),strip.text=element_text(size=12,color=c("white"))) + 
  facet_grid(SomaticType ~ .,space="free",scales="free") +
  xlab("TCGA LGG and GBM") + 
  ylab("Gene Name") +
  guides(col = guide_legend(nrow = 8))
dev.off()


mut.naomit <- na.omit(metadata.1129.samples.20141030[,c(names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"High-amp"], decreasing=F)),"cluster.meth")])
aux <- subset(mut.naomit, cluster.meth == "LGm1")
lgm1 <- mut.naomit[as.character(rownames(aux)),names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"High-amp"], decreasing=F))]
rownames(lgm1) <- rownames(aux)
hc.lgm1 <- as.hclust(agnes(daisy(data.frame(t(lgm1)), metric="gower"), method="ward"))
aux <- subset(mut.naomit, cluster.meth == "LGm2")
lgm2 <- mut.naomit[as.character(rownames(aux)),names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"High-amp"], decreasing=F))]
rownames(lgm2) <- rownames(aux)
hc.lgm2 <- as.hclust(agnes(daisy(data.frame(t(lgm2)), metric="gower"), method="ward"))
aux <- subset(mut.naomit, cluster.meth == "LGm3")
lgm3 <- mut.naomit[as.character(rownames(aux)),names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"High-amp"], decreasing=F))]
rownames(lgm3) <- rownames(aux)
hc.lgm3 <- as.hclust(agnes(daisy(data.frame(t(lgm3)), metric="gower"), method="ward"))
aux <- subset(mut.naomit, cluster.meth == "LGm4")
lgm4 <- mut.naomit[as.character(rownames(aux)),names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"High-amp"], decreasing=F))]
rownames(lgm4) <- rownames(aux)
hc.lgm4 <- as.hclust(agnes(daisy(data.frame(t(lgm4)), metric="gower"), method="ward"))
aux <- subset(mut.naomit, cluster.meth == "LGm5")
lgm5 <- mut.naomit[as.character(rownames(aux)),names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"High-amp"], decreasing=F))]
rownames(lgm5) <- rownames(aux)
hc.lgm5 <- as.hclust(agnes(daisy(data.frame(t(lgm5)), metric="gower"), method="ward"))
aux <- subset(mut.naomit, cluster.meth == "LGm6")
lgm6 <- mut.naomit[as.character(rownames(aux)),names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"High-amp"], decreasing=F))]
rownames(lgm6) <- rownames(aux)
hc.lgm6 <- as.hclust(agnes(daisy(data.frame(t(lgm6)), metric="gower"), method="ward"))

heatmap.mut <- rbind(lgm1[,hc.lgm1$order],lgm2[,hc.lgm2$order],lgm3[,hc.lgm3$order],lgm4[,hc.lgm4$order],lgm5[,hc.lgm5$order],lgm6[,hc.lgm6$order])
a <- metadata.1129.samples.20141030[rownames(heatmap.mut),]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$Master.cluster),as.matrix(a$flag.categories))
clab.6 <- clab.6[,-c(4,13,14)]
heatmap.mut <- as.matrix(sapply(heatmap.mut, as.numeric))

png(filename = "Heatmap_CNV_LGm1_2.png", res=300, width=5000, height=3000)
heatmap.plus.sm(as.matrix(t(heatmap.mut)),
             col = c("deepskyblue","red4","darkblue","red","lightgray"), 
             scale = "none",
             labRow = substr(colnames(heatmap.mut),5,length(colnames(heatmap.mut))),
             labCol = NA,
             Colv = NA,
             Rowv = NA,
             ColSideColors = clab.6)
dev.off()

aux <- substr(metadata.1129.samples.20141030$id,6,7)
a <- metadata.1129.samples.20141030
a$center <- as.character(aux)
b <- subset(a, center == "TQ")

png(filename = "teste.png", bg="white", res=200, width=4000, height=2000)
ggplot(a, aes(x=factor(cluster.meth2),fill=factor(Annotation.Category))) + 
  #geom_bar(position="fill")+ 
  geom_bar(aes(y = (..count..)/sum(..count..))) +
  scale_y_continuous(labels = percent) +
  #geom_bar(position = "dodge")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  #facet_grid(tumor.type ~ ., scales = "free") + 
  facet_grid(tumor.type ~ .) + 
  labs(y = "Total Count") + 
  labs(x = "'Annotation.Classification'")


### Replicates (patients with more than one sample available)
#teste <- metadata.new[duplicated(metadata.new$TCGAIDnew),] 
LGG.GBM.250914 <- LGG.GBM.new[,!colnames(LGG.GBM.new) %in% c("TCGA-06-0145-01A-02D-0218-05","TCGA-06-0145-01A-03D-0218-05","TCGA-06-0145-01A-04D-0218-05","TCGA-06-0145-01A-05D-0218-05","TCGA-06-0145-01A-06D-0218-05","TCGA-06-0137-01A-02D-0218-05","TCGA-06-0137-01A-03D-0218-05","TCGA-06-0137-01B-02D-0218-05","TCGA-02-0002-01A-01D-0186-05")]
#sel <- LGG.GBM.new[,c("TCGA-06-0145-01A-05D-0218-05","TCGA-06-0145-01A-01D-0218-05","TCGA-06-0145-01A-03D-0218-05","TCGA-06-0145-01A-04D-0218-05","TCGA-06-0145-01A-06D-0218-05","TCGA-06-0145-01A-02D-0218-05")]
#sel$mean <- apply(sel,1,mean,na.rm=T)
#sel2 <- LGG.GBM.new[,c("TCGA-06-0137-01B-02D-0218-05","TCGA-06-0137-01A-03D-0218-05","TCGA-06-0137-01A-01D-0218-05","TCGA-06-0137-01A-02D-0218-05")]
#se2l$mean <- apply(sel2,1,mean,na.rm=T)
#aux <- cbind(aux,sel$mean)
#aux <- cbind(aux,sel2$mean)
#colnames(aux)[934:935] <- c("TCGA-06-0145-01A-mean","TCGA-06-0137-01A-mean")

#status[status$TCGAIDnew %in% "TCGA-02-0002",1:10] #Genotype mismatch
#status[status$TCGAIDnew %in% "TCGA-32-2498",1:10] #Withdrew consent


### Info about sample status
status.lgg <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/tcga_annotations_LGG_24Sept2014.txt")
status.gbm <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/tcga_annotations_GBM_24Sept2014.txt")
a <- subset(status.gbm, Annotation.Classification != "CenterNotification")
a <- subset(teste, Annotation.Classification != "CenterNotification")
ggplot(d, aes(x=factor(Annotation.Classification), fill=factor(Annotation.Category))) + 
  geom_bar(position = "fill")+ 
  #geom_bar(position = "dodge")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(Disease ~ ., scales = "free") + 
  labs(y = "Total Percent") + 
  labs(x = "'Annotation.Classification'")

status <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/tcga_annotations_LGGxGBM_25Sept2014.txt")
status$TCGAIDnew <- substr(status$Item.Barcode,1,12)
mutation.clinical.data.synapse <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/new_mut_info.txt") #table from synapse
teste <- merge(status,mutation.clinical.data.synapse,by.x="TCGAIDnew",by.y="id")
b <- subset(teste, Annotation.Classification != "CenterNotification")
c <- b[substr(b$Item.Barcode,14,15) < "10",] #only tumor samples
d <- subset(c, Item.Type == "Patient" | Item.Type == "Sample" & Annotation.Notes != "No Matching Normal")
mutation.clinical.data.synapse$Annotation.Category <- NA
mutation.clinical.data.synapse$Annotation.Notes <- NA
mutation.clinical.data.synapse[mutation.clinical.data.synapse$id %in% as.character(d$TCGAIDnew),"Annotation.Category"] <- as.character(d$Annotation.Category)
mutation.clinical.data.synapse[mutation.clinical.data.synapse$id %in% as.character(d$TCGAIDnew),"Annotation.Notes"] <- as.character(d$Annotation.Notes)
#save(mutation.clinical.data.synapse,file="mutation.clinical.data.synapse.250914.Rda")

length(unique(c$TCGAIDnew))
d <- subset(c,tumor.type == "LGG")
length(unique(d$TCGAIDnew))

ggplot(b, aes(x=factor(Annotation.Classification), fill=factor(Annotation.Category))) + geom_bar(position="fill")+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(. ~ tumor.type) + scale_colour_discrete(drop = FALSE) #+ scale_fill_brewer(palette="Set1") 

clinical.tgca.lgg <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/lgg_clinical_290914.txt")
clinical.tgca.gbm <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/gbm_clinical_290914.txt")

a <- metadata.new
metadata.new$gender <- NA
rownames(metadata.new) <- metadata.new$TCGAIDnew
clinical.tgca.lgg <- clinical.tgca.lgg[clinical.tgca.lgg$bcr_patient_barcode %in% metadata.new$TCGAIDnew,]
clinical.tgca.gbm <- clinical.tgca.gbm[clinical.tgca.gbm$bcr_patient_barcode %in% metadata.new$TCGAIDnew,]
metadata.new[as.character(clinical.tgca.lgg$bcr_patient_barcode),"gender"] <- as.character(clinical.tgca.lgg$gender)
metadata.new[as.character(clinical.tgca.gbm$bcr_patient_barcode),"gender"] <- as.character(clinical.tgca.gbm$gender)
rownames(metadata.new) <- metadata.new$TCGAID

## KM plot
mod.cl.lgg.anno <- subset(metadata.1129.samples.20150312, Study == "TCGA") #objeto com as informacoes dos pacientes
#mod.cl.lgg.anno <- subset(mod.cl.lgg.anno, cluster.meth == "LGm6")
aux <- subset(mod.cl.lgg.anno, os.months > 60)
mod.cl.lgg.anno[rownames(aux),"os.months"] <- 60 #colocar a coluna com a data do ultimo contato com o paciente
mod.cl.lgg.anno[rownames(aux), "vital"] <- "0" #colocar a coluna com informacao sobre o estado do paciente 
clinical.cluster.genetics <- mod.cl.lgg.anno  
age.quint <- quantile(clinical.cluster.genetics$age,
                      probs=seq(0,1,.25), na.rm=T) #colocar a coluna com a idade no momento do diagnostico
#check values and make whole number version
age.quint <- c(10,38,51,62,89) #verificar quais valores estao em age.quint
agecat <- cut(clinical.cluster.genetics$age, age.quint)
fdata <- clinical.cluster.genetics 

fdata$age <- fdata$age #colocar a coluna com a idade no momento do diagnostico
fdata$age2 <- agecat
fdata$sex <- fdata$gender #colocar a coluna com o sexo do paciente
fdata$group <- as.factor(fdata$cluster.meth) #colocar a coluna com o resultado do Consensus Cluster
#library(survival)
data(uspop2)
refpop <- uspop2[as.character(floor(min(age.quint)):ceiling(max(age.quint))), , "2000"]
temp <- as.numeric(cut(floor(min(age.quint)):ceiling(max(age.quint)), age.quint, include.lowest=T))
pi.us<- tapply(refpop, list(temp[row(refpop)], col(refpop)), sum)/sum(refpop)
tab2 <- with(fdata, table(age2, sex, group))/ nrow(fdata)
pi.us2 <- pi.us[,2:1] # this puts male/female in the right order (check data)
us.wt <- rep(pi.us2, length(levels(fdata$group)))/ tab2
fdata$sex <- ifelse(fdata$sex == "female", 1, 2) #convert female to 1 and male to 2
index <- with(fdata, cbind(as.numeric(age2), as.numeric(sex), as.numeric(group)))
fdata$uswt <- us.wt[index]
sfit.default <- survfit(Surv(as.numeric(os.days), as.numeric(vital)) ~ group, fdata) #colocar a coluna com a data do ultimo contato com o paciente e a coluna com informacao sobre o estado do paciente. Utilizar este objeto caso nao deseje fazer o ajuste para idade e sexo do paciente
sfit3.adjusted <-survfit(Surv(as.numeric(os.days), as.numeric(vital)) ~ group, data=fdata, weight=uswt) #colocar a coluna com a data do ultimo contato com o paciente e a coluna com informacao sobre o estado do paciente. Utilizar este objeto caso deseje fazer o ajuste para idade e sexo do paciente
aux2 <- subset(fdata, group %in% c("LGm1","LGm3"))
gg.s<-coxph(Surv(as.numeric(os.days), as.numeric(vital)) ~ group, data = fdata)

fdata$s <- fdata$vital == "1"
f.m = formula(Surv(os.months,s) ~ cluster.meth)
fit.lgm = survfit(f.m, data=fdata)

survLGG <- fit.lgm#ou sfit.default
png(filename = "/dados/ResearchProjects/thais/TCGA/LGG.GBM/survival_all.png", bg="white", res=300, width=2000, height=2000)
plot(survLGG, 
     lwd=4,
     col=c("green",
       "red", #colocar a cor de acordo com a ordem de 1 a 6
           "purple",
           "orange",
           "yellow",
           "blue"

     ),
     main="Overall Survival\nKaplan-Meier Survival Curves\nTCGA GLIOMA SAMPLES ", 
     xlab="TIME SINCE DIAGNOSIS (YEARS)", 
     ylab="PERCENT PROBABILITY OF SURVIVAL", 
     yscale=100,
    # xscale=365.25, 
     bg="black"
)
box(col="black", lwd=3);
abline(h=0, lwd=4, lty=1)
legend("bottomleft", 
       legend=c(
         paste("LGm1 (n=",survLGG$n[1],")",sep=""),
                paste("LGm2 (n=",survLGG$n[2],")",sep=""),
         paste("LGm3 (n=",survLGG$n[3],")",sep=""),
         paste("LGm4 (n=",survLGG$n[4],")",sep=""),
         paste("LGm5 (n=",survLGG$n[5],")",sep=""),
         paste("LGm6 (n=",survLGG$n[6],")",sep="")

       ),
       col=c("green",
             "red", #colocar a cor de acordo com a ordem de 1 a 6
             "purple",
             "orange",
             "yellow",
             "blue" #colocar a cor de acordo com a ordem final

       ),
       lwd=3,
       title="DNA Meth Clusters",
       box.lwd=3,
       bg="white")
dev.off()

aux2 <- subset(metadata.1129.samples.20141016, LGG.DNA.Methyl.Clusters == "M1")

platform <- read.csv("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/LGG-GBM_sampleByPlatform_dataAvailability_09222014.csv")
metadata.new <- merge(metadata.new,platform[,c("id","Methylation.Platform")],by.x="TCGAIDnew",by.y="id")


#### Sturm samples (74 TCGA samples, we already have .. 59 pediatric and 77 adult GBMs in this file) Consensus Cluster + Heatmap
sturm.samples <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Sturm_GBM_GSE36278.txt",skip=139)
sturm.samples <- na.omit(sturm.samples)
sturm.samples.nocontrol <- sturm.samples[,-c(77:82)]
rownames(sturm.samples.nocontrol) <- sturm.samples.nocontrol$ID_REF
colnames(sturm.samples.nocontrol)[2:137] <- as.character(rownames(RF.Sturm.on.GBM)[c(1:75,82:142)])
#sturm.probes <- read.table("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/sturm2012_top8000probes.txt") # 204 probes no merged e 103 probes no 1300
#inters.sturm.ours.pb <- intersect(rownames(LGG.GBM.250914),sturm.probes$V1)
#sturm.samples.nocontrol.sel <- sturm.samples.nocontrol[inters.sturm.ours.pb,]
LGG.GBM.sturm <- merge(LGG.GBM.250914,sturm.samples.nocontrol,by.x="Composite.Element.REF",by.y="ID_REF")

#### Selecting probes for LGG + GBM + Sturm samples
LGG.GBM.sturm.noXY <- LGG.GBM.sturm[LGG.GBM.sturm$Chromosome != "X" & LGG.GBM.sturm$Chromosome != "Y",]
rownames(LGG.GBM.sturm.noXY) <- LGG.GBM.sturm.noXY$Composite.Element.REF

## dichotomizing the data
LGG.GBM.sturm.noXY.dic <- LGG.GBM.sturm.noXY[ ,5:1072]
LGG.GBM.sturm.noXY.dic <-LGG.GBM.sturm.noXY.dic>0.3
storage.mode(LGG.GBM.sturm.noXY.dic) <- "numeric"
## Use probes methylated more than 10% of the tumor samples
LGG.GBM.sturm.noXY.dic <- LGG.GBM.sturm.noXY.dic[(rowSums(LGG.GBM.sturm.noXY.dic)/ncol(LGG.GBM.sturm.noXY.dic))*100 >10,]
#LGG.GBM.new.noXY.dic <- LGG.GBM.new.noXY.dic[(rowSums(LGG.GBM.new.noXY[,5:936])/ncol(LGG.GBM.new.noXY[,5:936]))*100 >10,]

#head(FDb.HM450.oecg.3kb)
FDb.HM450.oecg.3kb.s <- subset(FDb.HM450.oecg.3kb, oecg.3kb>1.25)
LGG.GBM.sturm.noXY.dic.oecg <- LGG.GBM.sturm.noXY.dic[rownames(LGG.GBM.sturm.noXY.dic) %in% (FDb.HM450.oecg.3kb.s$probeid),]

#TODO:
## need to select tissue specific probes (>0.3)
avg.lgg.gbm.new.noXY <- apply(LGG.GBM.sturm.noXY[,5:1072], 1, mean, na.rm=TRUE)
avg.norm.lgg.gbm.new.noXY <- merge(as.data.frame(avg.lgg.gbm.new.noXY), avg.normals, by = 0, all.x = T)
#library(LSD)
#comparisonplot(avg.norm.lgg.gbm.27.450.noXY$avg.normals, avg.norm.lgg.gbm.27.450.noXY$avg.lgg.gbm.27.450.noXY, pimp=TRUE)
lgg.gbm.new.tumor.specific <- (subset(avg.norm.lgg.gbm.new.noXY, avg.normals < 0.3))

##selecting tumor specific probes for lgg.gbm
dat.lgg.gbm.sturm.noXY.dic.oecg <- LGG.GBM.sturm.noXY.dic.oecg[(rownames(LGG.GBM.sturm.noXY.dic.oecg) %in% lgg.gbm.new.tumor.specific$Row.names), ] ##1193 probes

library(ConsensusClusterPlus)
title="/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/ConsensusCluster.LGG.GBM.Sturm.Samples"
setwd(title)
cc.lgg.gbm.sturm <- ConsensusClusterPlus(d=dat.lgg.gbm.sturm.noXY.dic.oecg, 
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
ccICL = calcICL(cc.lgg.gbm.sturm, plot='pdf')

cc.lgg.gbm.sturm.2 <- cbind(as.matrix(cc.lgg.gbm.sturm[[2]][[3]]), 
                                     as.matrix(cc.lgg.gbm.sturm[[3]][[3]]), 
                                     as.matrix(cc.lgg.gbm.sturm[[4]][[3]]), 
                                     as.matrix(cc.lgg.gbm.sturm[[5]][[3]]), 
                                     as.matrix(cc.lgg.gbm.sturm[[6]][[3]]),
                                     as.matrix(cc.lgg.gbm.sturm[[7]][[3]]),
                                     as.matrix(cc.lgg.gbm.sturm[[8]][[3]]),
                                     as.matrix(cc.lgg.gbm.sturm[[9]][[3]]),
                                     as.matrix(cc.lgg.gbm.sturm[[10]][[3]]))

cc.lgg.gbm.sturm.2 <- as.data.frame(cc.lgg.gbm.sturm.2)
names(cc.lgg.gbm.sturm.2) <- c("K2.s",  "K3.s",  "K4.s" , "K5.s" , "K6.s" , "K7.s",  "K8.s" , "K9.s" , "K10.s")
cc.lgg.gbm.sturm.2$Sturm.IDnew <- substr(rownames(cc.lgg.gbm.sturm.2),1,12)
cc.lgg.gbm.sturm.2$Sturm.ID <- rownames(cc.lgg.gbm.sturm.2)
cc.lgg.gbm.sturm.2$TCGAID <- rownames(cc.lgg.gbm.sturm.2)
a <- merge(cc.lgg.gbm.sturm.2,metadata.new,by.x="TCGAID",by.y="TCGAID",all.x=T)
sturm.et.al <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Sturm_Cluster.txt")
c <- merge(a, sturm.et.al,by.x="Sturm.IDnew",by.y="ID",all.x=T)
rownames(c) <- c$Sturm.ID

LGG.GBM.sturm.noXY.s <- LGG.GBM.sturm.noXY[rownames(dat.lgg.gbm.sturm.noXY.dic.oecg) ,5:1072]
LGG.GBM.sturm.noXY.order <- LGG.GBM.sturm.noXY.s[,cc.lgg.gbm.sturm[[6]][[2]][[3]]]
c <-c[colnames(LGG.GBM.sturm.noXY.order),]
# Alterar a funcao para colocar o K6.s
clab.6 <- colors.label(as.matrix(c$K6.s),as.matrix(c$cluster.meth2),as.matrix(c$GBM.DNA.Methyl.Clusters),as.matrix(c$LGG.DNA.Methyl.Clusters),as.matrix(c$tumor.type),as.matrix(c$IDH.status),as.matrix(c$GBM.Clusters),as.matrix(c$cluster.mRNA),as.matrix(c$LGG.COC.Clusters),as.matrix(c$leukocyte.DNA.methyl),as.matrix(c$leukocyte.RNA),as.matrix(c$Cluster),as.matrix(c$new.cluster),as.matrix(c$Annotation.Category))

normals.sel <- norms[rownames(LGG.GBM.sturm.noXY.order),]
pilocytic.astro.GSE44684.order <- pilocytic.astro.GSE44684[rownames(LGG.GBM.sturm.noXY.order),]
normalBrain_dnameth.order <- normalBrain_dnameth[rownames(LGG.GBM.sturm.noXY.order),]
colnames(normalBrain_dnameth.order)[507] <- "UMARY.1648.TCTX"

normalBrain_dnameth.order <- normalBrain_dnameth.order[,-1]
crbl.mean <- apply(as.matrix(normalBrain_dnameth.order[,1:121]),1,mean,na.rm=T)
frctx.mean <- apply(as.matrix(normalBrain_dnameth.order[,122:254]),1,mean,na.rm=T)
pons.mean <- apply(as.matrix(normalBrain_dnameth.order[,255:380]),1,mean,na.rm=T)
tpctx.mean <- apply(as.matrix(normalBrain_dnameth.order[,381:506]),1,mean,na.rm=T)
regions.mean <- cbind(crbl.mean,frctx.mean,pons.mean,tpctx.mean)

LGG.GBM.sturm.noXY.order <- cbind(normals.sel,LGG.GBM.sturm.noXY.order)
LGG.GBM.sturm.noXY.order <- cbind(LGG.GBM.sturm.noXY.order,pilocytic.astro.GSE44684.order[,-c(1,63:68)]) #remove the control samples from PA
LGG.GBM.sturm.noXY.order <- cbind(regions.mean,LGG.GBM.sturm.noXY.order)


norm.label <- matrix("white", nrow = 77, ncol = 14)
#norm.label[,1] <- jet.colors(100)[ceiling((as.numeric(normal.age[,colnames(normals.sel)])-min(as.numeric(normal.age[,colnames(normals.sel)]),na.rm=T))/(max(as.numeric(normal.age[,colnames(normals.sel)]),na.rm=T)-min(as.numeric(normal.age[,colnames(normals.sel)]),na.rm=T))*(length(jet.colors(100))-1))+1]
norm.label[,14] <- "darkblue"
region.label <- matrix("white", nrow = 4, ncol = 14)
region.label[1,14] <- "aquamarine" #cerebelo
region.label[2,14] <- "pink" #frontal cortex
region.label[3,14] <- "deeppink" #pons 
region.label[4,14] <- "brown" #temporal contex

clab.6 <- rbind(norm.label,clab.6)
pa.label <- matrix("white", nrow = 61, ncol = 14)
pa.label[,14] <- "khaki3" #tumor type
#pa.label[,1] <- jet.colors(100)[ceiling((PA.age[1:61]-min(PA.age[1:61],na.rm=T))/(max(PA.age[1:61],na.rm=T)-min(PA.age[1:61],na.rm=T))*(length(jet.colors(100))-1))+1]
clab.6 <- rbind(clab.6,pa.label)
clab.6 <- rbind(region.label,clab.6)
#colnames(clab.6) <- 1:13
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/heatmap_sturm_ccProbes.png",res=200,width=2000,height=2000)

heatmap.plus.sm(as.matrix(LGG.GBM.sturm.noXY.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                #Rowv = NA,
                ColSideColors = clab.6,
                RowSideColors = rlab
                
)
dev.off()




## Heatpairs
K1 <- metadata.new[metadata.new$cluster.meth2 == "LGm1",]
K2 <- metadata.new[metadata.new$cluster.meth2 == "LGm2",]
K3 <- metadata.new[metadata.new$cluster.meth2 == "LGm3",]
K4 <- metadata.new[metadata.new$cluster.meth2 == "LGm4",]
K5 <- metadata.new[metadata.new$cluster.meth2 == "LGm5",]
K6 <- metadata.new[metadata.new$cluster.meth2 == "LGm6",]
b <- LGG.GBM.250914[rownames(dat.lgg.gbm.new.noXY.dic.oecg),5:936]
normal <- as.data.frame(norms[rownames(norms) %in% rownames(b),])
M1 <- b[rownames(normal),as.character(K1$TCGAID)]
M2 <- b[rownames(normal),as.character(K2$TCGAID)]
M3 <- b[rownames(normal),as.character(K3$TCGAID)]
M4 <- b[rownames(normal),as.character(K4$TCGAID)]
M5 <- b[rownames(normal),as.character(K5$TCGAID)]
M6 <- b[rownames(normal),as.character(K6$TCGAID)]
M1$meanM1 <- apply(M1,1,mean,na.rm=T)
M2$meanM2 <- apply(M2,1,mean,na.rm=T)
M3$meanM3 <- apply(M3,1,mean,na.rm=T)
M4$meanM4 <- apply(M4,1,mean,na.rm=T)
M5$meanM5 <- apply(M5,1,mean,na.rm=T)
M6$meanM6 <- apply(M6,1,mean,na.rm=T)
normal$mean <- apply(normal,1,mean,na.rm=T)
avg.probes.cl <- cbind(normal$mean,M1$meanM1,M2$meanM2,M3$meanM3,M4$meanM4,M5$meanM5,M6$meanM6)
rownames(avg.probes.cl) <- rownames(M1)
colnames(avg.probes.cl) <- c("normal","LGm1","LGm2","LGm3","LGm4","LGm5","LGm6")
png(filename = "/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/heatpairs_spec_probes.png", bg="white", res=300, width=2000, height=2000)
heatpairs(avg.probes.cl)
dev.off()

## Key
z <- seq(min(LGG.GBM.250914[,5:936],na.rm=TRUE), max(LGG.GBM.250914[,5:936],na.rm=TRUE), length = 200) #min(lgg.gbm.27.450.order) and #max(lgg.gbm.27.450.order)
n <- max(LGG.GBM.250914[,5:936],na.rm=TRUE)
png("/dados/ResearchProjects/thais/TCGA/LGG.GBM/key.png",res=600,width=5000,height=2000)
image(matrix(z, ncol = 1), col = jet.colors(75),
      xaxt = "n", yaxt = "n", main = "Overlap count Key")
box()
par(usr = c(0, n, 0, n))
axis(1, at = c(min(z),max(z)))


write.table(metadata.new[,c("TCGAID","TCGAIDnew","cluster.meth2")],file="DNA_methylation_new_clusters.txt",quote=F,row.names=F,sep="\t")

#### Only 450k samples (GBM + LGG) and Sturm GBM samples
sturm.et.al <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Sturm_Cluster.txt") #136 samples
a <- metadata.new
rownames(a) <- a$TCGAID
teste <- matrix(data="NA",nrow=136,ncol=ncol(a))
colnames(teste) <- colnames(a)
rownames(teste) <- sturm.et.al$ID[1:136]
teste[,"GBM.Cluster.Sturm.et.al"] <- as.character(sturm.et.al$Cluster)[1:136]
a <- rbind(a,teste)
lgg.gbm.450 <- merge(lgg.patients,gbm.patients[,-c(2,3,4)],by="Composite.Element.REF") #516 LGG 129 GBM
lgg.gbm.sturm <- merge(lgg.gbm.450,sturm.samples,by.x="Composite.Element.REF",by.y="ID_REF")
rownames(lgg.gbm.sturm) <- as.character(lgg.gbm.sturm$Composite.Element.REF)
probes.sturm <- read.table("/dados/ResearchProjects/thais/TCGA/LGG.GBM/hccpgs_order.txt")
lgg.gbm.sturm.sel <- lgg.gbm.sturm[as.character(probes.sturm$V1),]
b <- a[colnames(lgg.gbm.sturm.sel)[5:785],]

### CpG island (downloaded from Encode Aug, 28th)
cpg <- read.table("/dados/ResearchProjects/thais/TCGA/CpGisland.txt")
library(GenomicRanges)
cpg.is <- GRanges(seqnames = cpg$V1,ranges = IRanges(start = cpg$V2, end=cpg$V3), ID = cpg$V4)
aux <- lgg.gbm.sturm.sel
probe <- GRanges(seqnames = paste0("chr",aux$Chromosome),ranges = IRanges(start = aux$Genomic_Coordinate, end=aux$Genomic_Coordinate), ID = aux$Composite.Element.REF)
overlap <- as.data.frame(findOverlaps(probe,cpg.is))
rlab <- matrix("white", nrow = nrow(aux), ncol = 2)
rlab[overlap$queryHits,] <- "green"


clab.6 <- colors.label(as.matrix(b$cluster.meth2),as.matrix(b$GBM.DNA.Methyl.Clusters),as.matrix(b$LGG.DNA.Methyl.Clusters),as.matrix(b$tumor.type),as.matrix(b$IDH.status),as.matrix(b$GBM.Clusters),as.matrix(b$cluster.mRNA),as.matrix(b$LGG.COC.Clusters),as.numeric(b$leukocyte.DNA.methyl),as.numeric(b$leukocyte.RNA),as.matrix(b$GBM.Cluster.Sturm.et.al),as.matrix(b$new.cluster),as.matrix(b$Annotation.Category))

png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/heatmap_lgg_gbm_sturm_8000.png",res=200,width=2000,height=2000)

heatmap.plus.sm(as.matrix(lgg.gbm.sturm.sel[,5:785]),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                #Colv = NA,
                #Rowv = NA,
                ColSideColors = clab.6,
                RowSideColors = rlab
                
)
dev.off()

RF.PA.on.LGG <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/RF.PA.clusters.on.LGG")
RF.PA.on.GBM <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/RF.PA.clusters.on.GBM")
RF.Sturm.on.LGG <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/RF.STURM.clusters.on.LGG")
RF.Sturm.on.GBM <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/RF.STURM.clusters.on.GBM")
RF.GBM.new <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/RF.GBM.clusters.20")
RF.GBM.new$RFtype.all <- str_replace(RF.GBM.new$RFtype.all,"GCIMP","G-CIMP")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/RF.newLGG.lgg.rda")
RF.LGG.new <- as.data.frame(newlgg.labels)
RF.LGG.new$ID <- rownames(RF.LGG.new)


lgg.gbm.sturm.sel.dic <- lgg.gbm.sturm.sel[ ,5:785]
lgg.gbm.sturm.sel.dic <- lgg.gbm.sturm.sel.dic>0.3
storage.mode(lgg.gbm.sturm.sel.dic) <- "numeric"


library(ConsensusClusterPlus)
title="/dados/ResearchProjects/thais/TCGA/LGG.GBM/ConsensusCluster.450.LGG.GBM.Sturm"
setwd(title)
cc.450.lgg.gbm.sturm <- ConsensusClusterPlus(d=lgg.gbm.sturm.sel.dic, 
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
ccICL = calcICL(cc.450.lgg.gbm.sturm, plot='pdf')

cc.450.lgg.gbm.sturm.2 <- cbind(as.matrix(cc.450.lgg.gbm.sturm[[2]][[3]]), 
                            as.matrix(cc.450.lgg.gbm.sturm[[3]][[3]]), 
                            as.matrix(cc.450.lgg.gbm.sturm[[4]][[3]]), 
                            as.matrix(cc.450.lgg.gbm.sturm[[5]][[3]]), 
                            as.matrix(cc.450.lgg.gbm.sturm[[6]][[3]]),
                            as.matrix(cc.450.lgg.gbm.sturm[[7]][[3]]),
                            as.matrix(cc.450.lgg.gbm.sturm[[8]][[3]]),
                            as.matrix(cc.450.lgg.gbm.sturm[[9]][[3]]),
                            as.matrix(cc.450.lgg.gbm.sturm[[10]][[3]]))

cc.lgg.gbm.sturm.2 <- as.data.frame(cc.lgg.gbm.sturm.2)
names(cc.lgg.gbm.sturm.2) <- c("K2.s",  "K3.s",  "K4.s" , "K5.s" , "K6.s" , "K7.s",  "K8.s" , "K9.s" , "K10.s")
cc.lgg.gbm.sturm.2$Sturm.IDnew <- substr(rownames(cc.lgg.gbm.sturm.2),1,12)
cc.lgg.gbm.sturm.2$Sturm.ID <- rownames(cc.lgg.gbm.sturm.2)
cc.lgg.gbm.sturm.2$TCGAID <- rownames(cc.lgg.gbm.sturm.2)
a <- merge(cc.lgg.gbm.sturm.2,metadata.new,by.x="TCGAID",by.y="TCGAID",all.x=T)
sturm.et.al <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Sturm_Cluster.txt")
c <- merge(a, sturm.et.al,by.x="Sturm.IDnew",by.y="ID",all.x=T)
rownames(c) <- c$Sturm.ID

### Boxplot (mean beta value) all probes
LGG.GBM.sturm.PA <- merge(LGG.GBM.250914,sturm.samples.nocontrol,by.x = "Composite.Element.REF",by.y="ID_REF")
LGG.GBM.sturm.PA <- merge(LGG.GBM.sturm.PA, pilocytic.astro.GSE44684[,1:62],by.x="Composite.Element.REF",by.y="ID_REF")
lgm1 <- subset(metadata.new, cluster.meth2 == "LGm1")
lgm2 <- subset(metadata.new, cluster.meth2 == "LGm2")
lgm3 <- subset(metadata.new, cluster.meth2 == "LGm3")
lgm4 <- subset(metadata.new, cluster.meth2 == "LGm4")
lgm5 <- subset(metadata.new, cluster.meth2 == "LGm5")
lgm6 <- subset(metadata.new, cluster.meth2 == "LGm6")
lgm1.s <- LGG.GBM.sturm.PA[,as.character(lgm1$TCGAID)]
lgm2.s <- LGG.GBM.sturm.PA[,as.character(lgm2$TCGAID)]
lgm3.s <- LGG.GBM.sturm.PA[,as.character(lgm3$TCGAID)]
lgm4.s <- LGG.GBM.sturm.PA[,as.character(lgm4$TCGAID)]
lgm5.s <- LGG.GBM.sturm.PA[,as.character(lgm5$TCGAID)]
lgm6.s <- LGG.GBM.sturm.PA[,as.character(lgm6$TCGAID)]
lgm1.s$mean <- apply(lgm1.s,1,mean,na.rm=T)
lgm2.s$mean <- apply(lgm2.s,1,mean,na.rm=T)
lgm3.s$mean <- apply(lgm3.s,1,mean,na.rm=T)
lgm4.s$mean <- apply(lgm4.s,1,mean,na.rm=T)
lgm5.s$mean <- apply(lgm5.s,1,mean,na.rm=T)
lgm6.s$mean <- apply(lgm6.s,1,mean,na.rm=T)
sturm.s <- LGG.GBM.sturm.PA[,937:1072]
sturm.s$mean <- apply(sturm.s,1,mean,na.rm=T)
pa.s <- LGG.GBM.sturm.PA[,1073:1133]
pa.s$mean <- apply(pa.s,1,mean,na.rm=T)

x <- data.frame(mean = c(lgm1.s$mean,lgm2.s$mean,lgm3.s$mean,lgm4.s$mean,lgm5.s$mean,lgm6.s$mean,sturm.s$mean,pa.s$mean),
                type = c(rep("LGm1",20306),rep("LGm2",20306),rep("LGm3",20306),rep("LGm4",20306),rep("LGm5",20306),rep("LGm6",20306),rep("Sturm Samples",20306),rep("PA",20306)))
png("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/meanMethyGroups.png",res=300,width=4000,height=2000)
p <- ggplot(x, aes(factor(type), mean))
p + geom_boxplot(aes(fill = factor(type)), notchwidth=0.25) + 
  #geom_jitter(height = 0, position = position_jitter(width = .1), size=2)  + 
  scale_fill_manual(values=c("green","red","purple", "orange", "yellow","blue","brown","khaki3")) +
  #), labels = c("LGm1", "LGm2", "LGm3","LGm4", "LGm5", "LGm6"
  #)) +
  #scale_fill_manual(values = rev(c("darkorchid","darkgreen","blue","red","DARKGREY"))) +
  ylab("Mean Methylation Values") +
  xlab("Groups") + 
  labs(title = "Mean Methylation Values By Groups", fill = "Groups") +
  theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white"), panel.border = element_rect(colour = "white"), legend.position="none")
dev.off()


### new metadata file - 30/10/14
met <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/pdata.516lgg.606gbm.txt")


library(xlsx)
gbm.molecular.subtype <- read.xlsx("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Molecular_subtype_classification.xlsx",1)

met2 <- merge(met,gbm.molecular.subtype[,1:2], by.x = "id", by.y = "sample.id", all.x = T)

#Information about GBM clusters from Brennan, 2013
brennan.class <- read.delim(file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Brennan2013_GBM27K-450K/DNA.methylation.k6.txt", sep="\t", header=T)
#cc.out.purity.lgg.gbm.27.450.2$TCGAIDnew <- substr(rownames(cc.out.purity.lgg.gbm.27.450.2),1,12)
brennan.class$TCGAIDnew <- substr(brennan.class$TCGAID, 1, 12)
brennan.class.s <- brennan.class[brennan.class$DNA.Methylation.Clusters != "UNKNOWN", ]
met3 <- merge(met2, brennan.class.s[,c("DNA.Methylation.Clusters","TCGAIDnew")],by.x="id", by.y = "TCGAIDnew", all.x =T)

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

met3 <- merge(met3, clinical.cluster.genetics[,c("Tumor","MethylationCluster","COCCluster")], by.x = "id", by.y = "Tumor", all.x = T)

leukocyte <- read.csv("/dados/ResearchProjects/thais/TCGA/LGG.GBM/leukocyte_for_houtan2.csv")
leukocyte <- leukocyte[,-1] #remove column with row number

met4 <- merge(met3,leukocyte, by.x="id", by.y= "x",all.x=T)

estimate.affy <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/estimate.gbm529.affy.u133a.txt",skip=2)
estimate.unc <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/estimate.gbm154.lgg463.unc.rnaseqv2.txt",skip=2)
library(Hmisc)
uniques.affy <- estimate.affy[,colnames(estimate.affy) %nin% colnames(estimate.unc)]
uniques.unc <- estimate.unc[,colnames(estimate.unc) %nin% colnames(estimate.affy)]
duplicates <- estimate.unc[,colnames(estimate.unc) %in% colnames(estimate.affy)[3:531]]
uniques.affy <- uniques.affy[3,]
uniques.unc <- uniques.unc[3,]
duplicates <- duplicates[3,]
#duplicates <- uniques.unc[colnames(estimate.unc) %nin% colnames(uniques.unc)]
rna.estimates <- cbind(uniques.affy,uniques.unc,duplicates)
rna.estimates <- as.matrix(t(rna.estimates))
rownames(rna.estimates) <- str_replace_all(rownames(rna.estimates),"[.]","-")

met4 <- merge(met4,rna.estimates,by.x="id",by.y=0,all.x=T)

sturm.et.al <- read.xlsx("/dados/ResearchProjects/thais/TCGA/LGG.GBM/mmc2.xlsx",1)
met4 <- merge(met4, sturm.et.al[,c(1,2)],by.x="id",by.y="Sample.ID",all.x=T)
rownames(met4) <- met4$id
colnames(met4)[329:335] <- c("GBM.Expres.Clusters","GBM.DNA.Methyl.Clusters","LGG.DNA.Methyl.Clusters","LGG.COC.Clusters","leukocyte.DNA.methyl","leukocyte.RNA","GBM.Cluster.Sturm.et.al")

platform <- read.csv("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/LGG-GBM_sampleByPlatform_dataAvailability_09222014.csv")
met4 <- merge(met4,platform[,c("id","Methylation.Platform")],by.x="id",by.y="id")
met4$vital <- as.factor(met4$vital)
levels(met4$vital) <- c("0","1") #TRUE dead 1, FALSE alive 0


### Random Forest
RF.PA.on.LGG <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/RF.PA.clusters.on.LGG")
RF.PA.on.GBM <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/RF.PA.clusters.on.GBM")
RF.Sturm.on.LGG <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/RF.STURM.clusters.on.LGG")
RF.Sturm.on.GBM <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/RF.STURM.clusters.on.GBM")
RF.GBM.new <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/RF.GBM.clusters.20")
RF.GBM.new$RFtype.all <- str_replace(RF.GBM.new$RFtype.all,"GCIMP","G-CIMP")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/RF.newLGG.lgg.rda")
RF.LGG.new <- as.data.frame(newlgg.labels)
RF.LGG.new$ID <- rownames(RF.LGG.new)
a <- met4
rownames(a) <- a$id
a$newLGGcluster <- as.character(a$LGG.DNA.Methyl.Clusters)
RF.LGG.new$IDnew <- substr(RF.LGG.new$ID,1,12)
RF.GBM.new$IDnew <- substr(RF.GBM.new$ID,1,12)
a[as.character(RF.LGG.new$IDnew),"newLGGcluster"] <- as.character(RF.LGG.new$newlgg.labels)
a$newGBMcluster <- as.character(a$GBM.DNA.Methyl.Clusters)
a[as.character(RF.GBM.new$IDnew),"newGBMcluster"] <- as.character(RF.GBM.new$RFtype.all)

load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/pa.sturm.RF.s1.rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/all.lgggbm.tested.on.sturm.rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/RF.full.LGG.GBM.PA.Sturm.classification.rda")

pa.sturm.RF.s1.onlyTumor <- subset(pa.sturm.RF.s1, Tumor.Type %in% c("Pilocytic astrocytoma tumor","Primary brain tumor tissue"))
pa.sturm.RF.s1.onlyTumor$mut.BRAF <- as.character(pa.sturm.RF.s1.onlyTumor$PA.mutation.status)
aux <- grepl("BRAF",pa.sturm.RF.s1.onlyTumor$mut.BRAF)
pa.sturm.RF.s1.onlyTumor$mut.BRAF[aux] <- "Mut"
pa.sturm.RF.s1.onlyTumor$mut.BRAF[!aux & !is.na(pa.sturm.RF.s1.onlyTumor$mut.BRAF)] <- "WT"
levels(pa.sturm.RF.s1.onlyTumor$Tumor.Type) <- c("","PA","","","GBM")
levels(pa.sturm.RF.s1.onlyTumor$Gender) <- c("FEMALE","MALE","FEMALE","MALE","FEMALE")
pa.sturm.RF.s1.onlyTumor$RF.GBM.class <- str_replace(pa.sturm.RF.s1.onlyTumor$RF.GBM.class,"GCIMP","G-CIMP")
levels(pa.sturm.RF.s1.onlyTumor$Sturm.IDH.1.Mutation.Status) <- c("Mutant","WT")
levels(pa.sturm.RF.s1.onlyTumor$Sturm.TP53.Mutation.Status) <- c("Mut","WT")
pa.sturm.RF.s1.onlyTumor$Sturm.CDKN2A.Deletion <- str_replace(pa.sturm.RF.s1.onlyTumor$Sturm.CDKN2A.Deletion,"0","Non-Del")
pa.sturm.RF.s1.onlyTumor$Sturm.CDKN2A.Deletion <- str_replace(pa.sturm.RF.s1.onlyTumor$Sturm.CDKN2A.Deletion,"1","Del")
pa.sturm.RF.s1.onlyTumor$Sturm.EGFR.Amplification <- str_replace(pa.sturm.RF.s1.onlyTumor$Sturm.EGFR.Amplification,"0","Non-Amp")
pa.sturm.RF.s1.onlyTumor$Sturm.EGFR.Amplification <- str_replace(pa.sturm.RF.s1.onlyTumor$Sturm.EGFR.Amplification,"1","Amp")
pa.sturm.RF.s1.onlyTumor$Sturm.PDGFRA.Amplification <- str_replace(pa.sturm.RF.s1.onlyTumor$Sturm.PDGFRA.Amplification,"0","Non-Amp")
pa.sturm.RF.s1.onlyTumor$Sturm.PDGFRA.Amplification <- str_replace(pa.sturm.RF.s1.onlyTumor$Sturm.PDGFRA.Amplification,"1","Amp")
#pa.sturm.RF.s1.onlyTumor$Sturm.Death <- str_replace(pa.sturm.RF.s1.onlyTumor$Sturm.Death,"0","alive")
#pa.sturm.RF.s1.onlyTumor$Sturm.Death <- str_replace(pa.sturm.RF.s1.onlyTumor$Sturm.Death,"1","dead")
levels(a$IDH1.status) <- c("Mutant","Mutant","Mutant","Mutant","WT")
pa.sturm.RF.s1.onlyTumor$Tumor.Type <- as.character(pa.sturm.RF.s1.onlyTumor$Tumor.Type)
a$tumor.type <- as.character(a$tumor.type)
a$gender <- as.character(a$gender)
aux <- matrix(NA, nrow = nrow(pa.sturm.RF.s1.onlyTumor), ncol = ncol(a))
colnames(aux) <- colnames(a)
a <- rbind(a,aux)
a$Study <- "TCGA"
a$id <- as.character(a$id)
a[1123:1319,"id"] <- as.character(pa.sturm.RF.s1.onlyTumor$nonTCGA.ID)
a[1123:1319,"Study"] <- as.character(pa.sturm.RF.s1.onlyTumor$Study)
a[1123:1319,"tumor.type"] <- as.character(pa.sturm.RF.s1.onlyTumor$Tumor.Type)
a[1123:1319,"gender"] <- as.character(pa.sturm.RF.s1.onlyTumor$Gender)
a[1123:1319,"age"] <- as.character(pa.sturm.RF.s1.onlyTumor$Age.at.diagnosis.years)
a[1123:1319,"newLGGcluster"] <- as.character(pa.sturm.RF.s1.onlyTumor$RF.LGG.class)
a[1123:1319,"newGBMcluster"] <- as.character(pa.sturm.RF.s1.onlyTumor$RF.GBM.class)
a[1123:1319,"IDH1.status"] <- as.character(pa.sturm.RF.s1.onlyTumor$Sturm.IDH.1.Mutation.Status)
a[1123:1319,"mut.TP53"] <- as.character(pa.sturm.RF.s1.onlyTumor$Sturm.TP53.Mutation.Status)
#a[1123:1319,"del.CDKN2A"] <- as.character(pa.sturm.RF.s1.onlyTumor$Sturm.CDKN2A.Deletion)
#a[1123:1319,"amp.EGFR"] <- as.character(pa.sturm.RF.s1.onlyTumor$Sturm.EGFR.Amplification)
#a[1123:1319,"amp.PDGFRA"] <- as.character(pa.sturm.RF.s1.onlyTumor$Sturm.PDGFRA.Amplification)
a[1123:1319,"mut.BRAF"] <- as.character(pa.sturm.RF.s1.onlyTumor$mut.BRAF)
a$PA.BRAF.fusion <- NA
pa.sturm.RF.s1.onlyTumor$Study.Class.original <- as.factor(pa.sturm.RF.s1.onlyTumor$Study.Class.original)
y <- levels(pa.sturm.RF.s1.onlyTumor$Study.Class.original)
levels(pa.sturm.RF.s1.onlyTumor$Study.Class.original) <- c("1-PA", "2-PA", y[3:9])
a$Sturm.PA.merged.clusters <- NA
a[1123:1319,"Sturm.PA.merged.clusters"] <- as.character(pa.sturm.RF.s1.onlyTumor$Study.Class.original)
a[1123:1319,"PA.BRAF.fusion"] <- as.character(pa.sturm.RF.s1.onlyTumor$PA.mutation.status)
a[1123:1319,"vital"] <- as.character(pa.sturm.RF.s1.onlyTumor$Sturm.Death)
a[1123:1319,"os.months"] <- as.character(pa.sturm.RF.s1.onlyTumor$Sturm.OS.months)
a[,"GBM.Cluster.Sturm.et.al"] <- as.character(a[,"GBM.Cluster.Sturm.et.al"])
pa.sturm.RF.s1.onlyTumor$Study.Class.original <- as.character(pa.sturm.RF.s1.onlyTumor$Study.Class.original)
a[1123:1258,"GBM.Cluster.Sturm.et.al"] <- as.character(pa.sturm.RF.s1.onlyTumor$Study.Class.original[1:136])
a[1123:1258,"histology"] <- "glioblastoma"
a[1123:1258,"grade"] <- "G4"
a[1123:1258,"Methylation.Platform"] <- "450k"
a$histology <- as.character(a$histology)
a$grade <- as.character(a$grade)
a[1259:1319,"histology"] <- "pilocytic astrocytoma"
a[1259:1319,"grade"] <- "G1"
a$PA.original.clusters <- NA
a[1259:1319,"PA.original.clusters"] <- as.character(pa.sturm.RF.s1.onlyTumor$Study.Class.original[137:197])
rownames(a) <- a$id
a$cluster.meth.all <- NA
RF.class$TCGAIDnew <- RF.class$TCGAID
RF.class$TCGAIDnew[1:932] <- substr(RF.class$TCGAIDnew[1:932],1,12)
a[as.character(RF.class$TCGAIDnew),"cluster.meth.all"] <- as.character(RF.class$cluster.meth2)
rownames(a) <- a$id

all.lgggbm.tested.on.sturm <- data.frame(all.lgggbm.tested.on.sturm)
all.lgggbm.tested.on.sturm$ID <- substr(rownames(all.lgggbm.tested.on.sturm),1,12)
all.lgggbm.tested.on.sturm$all.lgggbm.tested.on.sturm <- str_replace(all.lgggbm.tested.on.sturm$all.lgggbm.tested.on.sturm,"Classical","RTK II 'Classic'")
all.lgggbm.tested.on.sturm$all.lgggbm.tested.on.sturm <- str_replace(all.lgggbm.tested.on.sturm$all.lgggbm.tested.on.sturm,"PDGFRA","RTK I 'PDGFRA'")
a$newSturmClass <- a$GBM.Cluster.Sturm.et.al
a[as.character(all.lgggbm.tested.on.sturm$ID),"newSturmClass"] <- as.character(all.lgggbm.tested.on.sturm$all.lgggbm.tested.on.sturm)
rownames(a) <- a$id
## FIM  Random Forest 

metadata.1129.samples.20150312 <- a
#metadata.1129.samples.20141016.annotation <- read.table(file = "/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG-GBM.Annotation.Table-20141016.csv", na.strings="NA", stringsAsFactors=FALSE, header = T, sep = ",")

metadata.1129.samples.20150312 <- subset(metadata.1129.samples.20150312, !is.na(cluster.meth.all))

metadata.1129.samples.20150312[is.na(metadata.1129.samples.20150312$newSturmClass) & metadata.1129.samples.20150312$Study == "TCGA", "newSturmClass"] <- "Unclassified-TCGA-27K"
metadata.1129.samples.20150312$Sturm.PA.merged.clusters[1:932] <- as.character(metadata.1129.samples.20150312$cluster.meth[1:932])
metadata.1129.samples.20150312$age <- as.numeric(metadata.1129.samples.20150312$age)
metadata.1129.samples.20150312$os.months <- as.numeric(metadata.1129.samples.20150312$os.months)

#metadata.all <- metadata.1129.samples.20141030


save(metadata.1129.samples.20150312,file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/metadata_LGG_GBM_PA_Sturm-20150312.Rda")


auxGBM <- subset(metadata.1129.samples.20141030, cluster.meth == "LGm4")
auxGBM <- subset(auxGBM, tumor.type == "GBM")
aux <- subset(mut.naomit, cluster.meth == "LGm6")
lgm6 <- mut.naomit[as.character(rownames(auxGBM)),names(sort(table(meanb.naomit$variable, meanb.naomit$value)[,"Mut"], decreasing=F))]
lgm6 <- as.matrix(sapply(lgm6, as.numeric))
lgm6 <- lgm6 == 1
storage.mode(lgm6) <- "numeric"
rownames(lgm6) <- rownames(auxGBM)
hc.lgm6.GBM <- hclust(dist(lgm6,method="binary"),method="ward")


auxGBM <- subset(metadata.1129.samples.20141030, cluster.meth %in% "LGm6")
aux2.m1 <- subset(auxGBM, newGBMcluster != "G-CIMP" | is.na(newGBMcluster))
aux2.m2 <- subset(auxGBM, newGBMcluster != "G-CIMP" | is.na(newGBMcluster))
aux2.m3 <- subset(auxGBM, newGBMcluster != "G-CIMP" | is.na(newGBMcluster))
aux2.m4 <- subset(auxGBM, newGBMcluster == "G-CIMP")
aux2.m5 <- subset(auxGBM, newGBMcluster == "G-CIMP")
aux2.m6 <- subset(auxGBM, newGBMcluster == "G-CIMP")

aux2.ms <- rbind(aux2.m1,aux2.m2,aux2.m3)

aux <- LGG.GBM.250914
colnames(aux) <- substr(colnames(aux),1,12)
aux3 <- aux[,as.character(auxGBM$id)]
aux3.m1 <- apply(aux3,1,mean,na.rm=TRUE)
aux3.m2 <- apply(aux3,1,mean,na.rm=TRUE)
aux3.m3 <- apply(aux3,1,mean,na.rm=TRUE)
aux3.m4 <- apply(aux3,1,mean,na.rm=TRUE)
aux3.m5 <- apply(aux3,1,mean,na.rm=TRUE)
aux3.m6 <- apply(aux3,1,mean,na.rm=TRUE)

aux3.m1 <- as.data.frame(aux3.m1)
aux3.m2 <- as.data.frame(aux3.m2)
aux3.m3 <- as.data.frame(aux3.m3)
aux3.m4 <- as.data.frame(aux3.m4)
aux3.m5 <- as.data.frame(aux3.m5)
aux3.m6 <- as.data.frame(aux3.m6)
aux3.m1$type2 <- "non.gcimp"
aux3.m2$type2 <- "non.gcimp"
aux3.m3$type2 <- "non.gcimp"
aux3.m4$type2 <- "gcimp"
aux3.m5$type2 <- "gcimp"
aux3.m6$type2 <- "gcimp"
aux3.m1$type <- "lgm1"
aux3.m2$type <- "lgm2"
aux3.m3$type <- "lgm3"
aux3.m4$type <- "lgm4"
aux3.m5$type <- "lgm5"
aux3.m6$type <- "lgm6"
names(aux3.m1)[1] <- "lgm1"
names(aux3.m2)[1] <- "lgm2"
names(aux3.m3)[1] <- "lgm3"
names(aux3.m4)[1] <- "lgm4"
names(aux3.m5)[1] <- "lgm5"
names(aux3.m6)[1] <- "lgm6"
aux3.ms <- rbind(aux3.m1,aux3.m2,aux3.m3,aux3.m4,aux3.m5,aux3.m6)

png(filename = "meanMethGroup.png", res=200, width=2000, height=2000)
ggplot(aux3.ms, aes(factor(type), aux3.m2)) + geom_boxplot(aes(fill = factor(type))) + ylim(0,1) #+ facet_wrap( ~ type)
dev.off()

library(LSD)
png(filename = "heatpairsMeanProbe.png", res=200, width=2000, height=2000)
heatpairs(cbind(aux3.m1,aux3.m2,aux3.m3,aux3.m4,aux3.m5,aux3.m6))


### Heatmap

a <- subset(metadata.1129.samples.20141030, Study=="TCGA")


LGG.GBM.new.noXY.s <- LGG.GBM.new.noXY[rownames(full.heatmap.orderd) ,5:936]
LGG.GBM.new.order <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
a <- metadata.1129.samples.20141030[substr(colnames(LGG.GBM.new.order),1,12),]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$Master.cluster),as.matrix(a$flag.categories),as.matrix(a$cluster.cnv))
clab.6 <- clab.6[,c(3,4,6,8,16,10,15,11)]
LGG.GBM.new.order <- LGG.GBM.new.order[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
clab.6 <- clab.6[c(1:63,64:327,405:531,532:682,683:932,328:404),]
normals.sel <- norms[rownames(LGG.GBM.new.order),]
norm.label <- matrix("white", nrow = 77, ncol = 8)
norm.label[,7] <- "darkblue"
clab.6 <- rbind(norm.label,clab.6)
LGG.GBM.new.order <- cbind(normals.sel,LGG.GBM.new.order)
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap.png",res=800,width=4000,height=4000)
heatmap.plus.sm(as.matrix(LGG.GBM.new.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab.6
                #RowSideColors = rlab
                
)
dev.off()


## heatmap GCIMP
GCIMP <- subset(metadata.1129.samples.20141030, cluster.meth.gbm == "G-CIMP" )
LGG.GBM.new.noXY.s <- LGG.GBM.250914[rownames(full.heatmap.orderd),5:936]
LGG.GBM.new.order <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
colnames(LGG.GBM.new.order) <- substr(colnames(LGG.GBM.new.order),1,12)
LGG.GBM.new.order <- LGG.GBM.new.order[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
b <- colnames(LGG.GBM.new.order) %in% as.character(GCIMP$id)
teste <- LGG.GBM.new.order[,b]
a <- metadata.1129.samples.20141030[substr(colnames(teste),1,12),]
clab.6 <- colors.label(as.matrix(a$newLGGcluster),as.matrix(a$newGBMcluster),as.matrix(a$codel1p19q),as.matrix(a$age),as.matrix(a$cluster.meth),as.matrix(a$tumor.type),as.matrix(a$IDH.status),as.matrix(a$GBM.Expres.Clusters),as.matrix(a$cluster.expr),as.matrix(a$LGG.COC.Clusters),as.matrix(a$leukocyte.DNA.methyl),as.matrix(a$leukocyte.RNA),as.matrix(a$newSturmClass),as.matrix(a$Master.cluster),as.matrix(a$flag.categories),as.matrix(a$cluster.cnv))
clab.6 <- clab.6[,c(3,4,7,16,10,11,15)]

png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/Heatmap_GCIMP.png",res=800,width=4000,height=4000)
heatmap.plus.sm(as.matrix(teste),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
                Rowv = NA,
                ColSideColors = clab.6
                #RowSideColors = rlab
                
)
dev.off()

library(survival)
pd <- subset(metadata.1129.samples.20150312,Study == "TCGA" & cluster.meth %in% c("LGm1","LGm2","LGm3") & tumor.type == "LGG")   #8 180
pd$OS.months <- as.numeric(as.character(pd$os.months))
aux <- subset(pd, OS.months > 60) ##quem viveu mais que 5 anos
pd[rownames(aux),"OS.months"] <- 60 #colocar a coluna com a data do ultimo contato com o paciente
pd[rownames(aux), "vital"] <- "0" #colocar a coluna com informacao sobre o estado do paciente 
pd$s <- pd$vital == "1"
pd$os <- pd$OS.months
f = formula(Surv(os,s) ~ cluster.meth)
fit = survfit(f, data=pd)
pdf(file = "/dados/ResearchProjects/thais/TCGA/LGG.GBM/survival_LGGinLGm1-2-3.pdf")
plot(fit, 
     lwd=4,  
     col=c("green", "red","purple"),
     main="Overall Survival\nKaplan-Meier Survival Curves\nTCGA LGG Samples within LGm1-2-3 ", 
     xlab="TIME SINCE DIAGNOSIS (MONTHS)", 
     ylab="PERCENT PROBABILITY OF SURVIVAL", 
     yscale=100,
     xscale=1, 
     bg="black")
box(col="black", lwd=3);
#abline(h=0, lwd=4, lty=1)
legend("bottomleft", cex=1, 
       legend=c(
         paste("LGm1 LGG (n=",fit$n[1],")",sep=""),
         paste("LGm2 LGG (n=",fit$n[2],")",sep=""),
         paste("LGm3 LGG (n=",fit$n[3],")",sep="")),
       col=c("green", "red","purple"),
       lwd=3,
       title="DNA Meth Clusters",
       box.lwd=	3,
       bg="white")
dev.off()


##### LGm6