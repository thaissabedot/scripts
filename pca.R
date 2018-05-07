####pca
setwd("/dados/ResearchProjects/thais/TCGA/LGG.GBM")
load("LGGs_IDHwt.rda")
load("Turcan.DNA_Meth.Rda")
load("PA.DNA_Meth.Rda")
load("Mur_Sturm_clinical.rda")
load("Sturm_DNA_meth.rda")
load("mur.DNA_Meth.Rda")
load("LGG_GBM_new_order.rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/pa.sturm.RF.s1.rda")

pa.sturm.RF.s1.onlyTumor <- subset(pa.sturm.RF.s1, Tumor.Type %in% c("Pilocytic astrocytoma tumor","Primary brain tumor tissue"))

LGG.22.IDHwt.pdata$IDH.status <- NA
LGG.22.IDHwt.pdata[LGG.22.IDHwt.pdata$characteristics_ch1.3 == levels(LGG.22.IDHwt.pdata$characteristics_ch1.3)[1],"IDH.status"] <- "WT"
LGG.22.IDHwt.pdata[LGG.22.IDHwt.pdata$characteristics_ch1.3 == levels(LGG.22.IDHwt.pdata$characteristics_ch1.3)[2],"IDH.status"] <- "WT"
LGG.22.IDHwt.pdata[LGG.22.IDHwt.pdata$characteristics_ch1.3 == levels(LGG.22.IDHwt.pdata$characteristics_ch1.3)[3],"IDH.status"] <- "Mutant"
LGG.22.IDHwt.pdata[LGG.22.IDHwt.pdata$characteristics_ch1.3 == levels(LGG.22.IDHwt.pdata$characteristics_ch1.3)[4],"IDH.status"] <- "Mutant"

LGG.44.IDHwt.450k.pdata$IDH.status <- NA
LGG.44.IDHwt.450k.pdata[LGG.44.IDHwt.450k.pdata$characteristics_ch1.6 == levels(LGG.44.IDHwt.450k.pdata$characteristics_ch1.6)[2],"IDH.status"] <- "Mutant"
LGG.44.IDHwt.450k.pdata[LGG.44.IDHwt.450k.pdata$characteristics_ch1.6 == levels(LGG.44.IDHwt.450k.pdata$characteristics_ch1.6)[3],"IDH.status"] <- "WT"

LGG.44.IDHwt.27k.pdata$IDH.status <- NA
LGG.44.IDHwt.27k.pdata[LGG.44.IDHwt.27k.pdata$characteristics_ch1.6 == levels(LGG.44.IDHwt.27k.pdata$characteristics_ch1.6)[2],"IDH.status"] <- "Mutant"
LGG.44.IDHwt.27k.pdata[LGG.44.IDHwt.27k.pdata$characteristics_ch1.6 == levels(LGG.44.IDHwt.27k.pdata$characteristics_ch1.6)[3],"IDH.status"] <- "WT"


mur.meth <- gset.mur.m[,as.character(gset.mur.p$geo_accession)]


LGG.22.IDHwt.pdata$Study <- "Zhang"
LGG.44.IDHwt.27k.pdata$Study <- "van_den_Bent"
LGG.44.IDHwt.450k.pdata$Study <- "van_den_Bent"
pa.pdata$Study <- "Lambert" #tem os 6 controles aqui
sturm.ids$Study <- "Sturm"

###Heatmap com LGm1-hypo e LGm6-LGG-IDHwt probes do validation set com informação de tumor type, study e IDH status
sturm.ids$IDH.status <- NA
sturm.ids[sturm.ids$characteristics_ch1.5 == levels(sturm.ids$characteristics_ch1.5)[2],"IDH.status"] <- "Mutant"
sturm.ids[sturm.ids$characteristics_ch1.5 == levels(sturm.ids$characteristics_ch1.5)[3],"IDH.status"] <- "WT"

turcan.sv$IDH.status <- NA
turcan.sv[turcan.sv$CIMP.phenotype == levels(turcan.sv$CIMP.phenotype)[1],"IDH.status"] <- "Mutant"
turcan.sv[turcan.sv$CIMP.phenotype == levels(turcan.sv$CIMP.phenotype)[2],"IDH.status"] <- "WT"

pa.pdata$IDH.status <- NA
pa.pdata <- subset(pa.pdata, source_name_ch1 == "Frozen pilocytic astrocytoma tumor tissue")

levels(mur_clinical_46$IDH) <- c("WT","Mutant")

a <- read.delim("Clinical_LGG22_Nature.txt")
LGG.22.IDHwt.pdata2 <- merge(LGG.22.IDHwt.pdata,a,by.x="title",by.y="Sample_title",all.x=TRUE)

LGG.44.IDHwt.450k.pdata$age <- str_extract(as.character(LGG.44.IDHwt.450k.pdata$characteristics_ch1.13),"[0-9]{1,}")
LGG.44.IDHwt.27k.pdata$age <- str_extract(as.character(LGG.44.IDHwt.27k.pdata$characteristics_ch1.13),"[0-9]{1,}")
LGG.22.IDHwt.pdata2$age <- str_extract(as.character(LGG.22.IDHwt.pdata2$characteristics_ch1.1),"[0-9]{1,}.*")
pa.pdata <- merge(pa.pdata,pa.sturm.RF.s1.onlyTumor[,c(1,7)],by.x="geo_accession",by.y="nonTCGA.ID")
sturm.ids <- merge(sturm.ids,pa.sturm.RF.s1.onlyTumor[,c(1,7,22,23)],by.x="geo_accession",by.y="nonTCGA.ID")
LGG.44.IDHwt.450k.pdata$OS.months <- str_extract(as.character(LGG.44.IDHwt.450k.pdata$characteristics_ch1),"[0-9]{1,}.*")
LGG.44.IDHwt.450k.pdata$OS.months <- as.numeric(LGG.44.IDHwt.450k.pdata$OS.months)*12
LGG.44.IDHwt.27k.pdata$OS.months <- str_extract(as.character(LGG.44.IDHwt.27k.pdata$characteristics_ch1),"[0-9]{1,}.*")
LGG.44.IDHwt.27k.pdata$OS.months <- as.numeric(LGG.44.IDHwt.27k.pdata$OS.months)*12


ids <- data.frame(id = c(as.character(pa.pdata$geo_accession),as.character(mur_clinical_46$geo_accession),as.character(sturm.ids$geo_accession),as.character(turcan.sv$geo_accession),as.character(LGG.22.IDHwt.pdata2$geo_accession),as.character(LGG.44.IDHwt.27k.pdata$geo_accession),as.character(LGG.44.IDHwt.450k.pdata$geo_accession)),Study = c(as.character(pa.pdata$Study),as.character(mur_clinical_46$Study),as.character(sturm.ids$Study),as.character(turcan.sv$Study),as.character(LGG.22.IDHwt.pdata2$Study),as.character(LGG.44.IDHwt.27k.pdata$Study),as.character(LGG.44.IDHwt.450k.pdata$Study)),IDH.status=c(as.character(pa.pdata$IDH.status),as.character(mur_clinical_46$IDH),as.character(sturm.ids$IDH.status),as.character(turcan.sv$IDH.status),as.character(LGG.22.IDHwt.pdata2$IDH.status),as.character(LGG.44.IDHwt.27k.pdata$IDH.status),as.character(LGG.44.IDHwt.450k.pdata$IDH.status)),Histology=c(as.character(pa.pdata$characteristics_ch1),as.character(mur_clinical_46$Diagnosis),as.character(sturm.ids$characteristics_ch1),as.character(turcan.sv$Pathology),as.character(LGG.22.IDHwt.pdata2$Histopathology),as.character(LGG.44.IDHwt.27k.pdata$characteristics_ch1.11),as.character(LGG.44.IDHwt.450k.pdata$characteristics_ch1.11)),age = c(as.character(pa.pdata$Age.at.diagnosis.years),as.character(mur_clinical_46$Age),as.character(sturm.ids$Age.at.diagnosis.years),as.character(turcan.sv$Age.at.Diagnosis),as.character(LGG.22.IDHwt.pdata2$age),as.character(LGG.44.IDHwt.27k.pdata$age),as.character(LGG.44.IDHwt.450k.pdata$age)),OS.months = c(rep(NA,61),as.character(mur_clinical_46$OS.months),as.character(sturm.ids$Sturm.OS.months),as.character(turcan.sv$Follow.Up..Months.),as.character(LGG.22.IDHwt.pdata2$survival.months),as.character(LGG.44.IDHwt.27k.pdata$OS.months),as.character(LGG.44.IDHwt.450k.pdata$OS.months)),status = c(rep(NA,61),as.character(mur_clinical_46$Status),as.character(sturm.ids$Sturm.Death),as.character(turcan.sv$Deceased),as.character(LGG.22.IDHwt.pdata2$status.y),as.character(LGG.44.IDHwt.27k.pdata$characteristics_ch1.1),as.character(LGG.44.IDHwt.450k.pdata$characteristics_ch1.1)))
rownames(ids) <- as.character(ids$id)

levels(ids$Histology) <- c("Astrocytoma", "Astrocytoma","other/non-glioma", 
                           "other/non-glioma", "Glioblastoma", "Oligoastrocytoma",
                           "Oligoastrocytoma", "Oligodendroglioma","Oligoastrocytoma",
                           "Oligoastrocytoma", "Oligodendroglioma", "Oligodendroglioma",
                           "other/non-glioma", "other/non-glioma","Not Available", 
                           "Glioblastoma","Oligoastrocytoma","Oligodendroglioma",
                           "LGG", "Pilocytic astrocytoma", "Glioblastoma")
levels(ids$status) <- c(0,1,0,1,NA,0,NA,0,1,1)
#levels(ids$IDH.status) <- c("Mutant","WT","WT","Mutant")
levels(ids$age)[levels(ids$age) == "ND "] <- NA
levels(ids$OS.months)[levels(ids$OS.months) == "ND "] <- NA
ids$age <- as.numeric(as.character(ids$age))
ids$OS.months <- as.numeric(as.character(ids$OS.months))
ids$Histology <- as.character(ids$Histology)



all <- merge(LGG.GBM.250914,norms,by=0)
rownames(all) <- as.character(all$Row.names)
all <- all[,-1]
all <- merge(all,sturm.meth.tumor,by=0)
rownames(all) <- as.character(all$Row.names)
all <- all[,-1]
all <- merge(all,pa.meth.nocontrol,by=0)
rownames(all) <- as.character(all$Row.names)
all <- all[,-1]
all <- merge(all,turcan.tumor.bv,by=0)
rownames(all) <- as.character(all$Row.names)
all <- all[,-1]
all <- merge(all,mur.meth,by=0)
rownames(all) <- as.character(all$Row.names)
all <- all[,-1]

all <- all[,5:1337]
a <- substr(colnames(all)[1:932],1,12)
a <- as.character(metadata.1129.samples.20150312[as.character(a),"cluster.meth"])
b <- substr(colnames(all)[1:932],1,12)
b <- as.character(metadata.1129.samples.20150312[as.character(b),"tumor.type"])
c <- ids[colnames(all)[1010:1333],]
c$Histology <- as.factor(c$Histology)
levels(c$Histology) <- c("LGG","GBM","LGG","LGG","other/non-glioma","PA")
pca.meta <- data.frame(Study = c(a,rep("non-tumor",77),rep("Sturm et al., 2012",136),rep("Lambert et al., 2013",61),rep("Turcan et al., 2012",81),rep("Mur et al., 2013",46)), Tumor.Type = c(b,rep("non-tumor brain",77),as.character(c$Histology)))



pca.data <- na.omit(all)  #17,492 1,333

pca.data <- prcomp(t(pca.data), scales=TRUE) #calcula o principal components analysis
scores <- data.frame(pca.meta, pca.data$x[, 1:3]) #lgg.pca$x[, 1:3] means that estamos utilizando os 3 primeiros PCs
#para checar lgg.pca:
#summary(lgg.pca)
#plot(lgg.pca) (Scree Plot)

library(ggplot2)
#png((filename = "/dados/ResearchProjects/tathi/TCGA_LGG_Recurrent/Images/PCAnew", res=300, width=2000, height=2000))
qplot(x=PC1, y=PC3, data=scores, color=factor(pca.meta$Study)) + 
  ggtitle("PCA Plot\nGliomas") + 
  geom_point(size = 4,aes(shape=factor(pca.meta$Tumor.Type))) +
  #scale_color_discrete(name="recurrentVSprimmary") + 
  #scale_color_discrete(name="primary-recurrent") + 
  scale_color_manual(values=c("green", "red", "purple","orange","yellow","blue","lightblue","salmon","lightpink4","tomato4","magenta"), 
                    # breaks=c("LGm1","LGm2","LGm3","LGm4","LGm5","LGm6","non-tumor","Sturm et al., 2012","Lambert et al., 2013","Turcan et al., 2012","Mur et al., 2013"), # color scale (for points) 
                     labels=c("LGm1","LGm2","LGm3","LGm4","LGm5","LGm6","non-tumor","Sturm et al., 2012","Lambert et al., 2013","Turcan et al., 2012","Mur et al., 2013"),
                     name="Legend") 
ggsave(file="PCA_tumorType.pdf")
dev.off()

tomato4 = Turcan
salmon = Sturm
Magenta = Mur
Turquoise = TCGA
Lightpink4 = Lambert
levels(info$Study) <- c("magenta","lightpink4","salmon","turquoise","tomato4")

