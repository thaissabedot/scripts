setwd("/dados/ResearchProjects/thais/TCGA/LGG.GBM")
load("Turcan.DNA_Meth.Rda")
load("PA.DNA_Meth.Rda")
load("Mur_Sturm_clinical.rda")
load("Sturm_DNA_meth.rda")
load("mur.DNA_Meth.Rda")
load("LGG_GBM_new_order.rda")


load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/pa.sturm.RF.s1.rda")
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

levels(mur_clinical_46$IDH) <- c("WT","Mutant")

load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/all_clinical_data.Rda")

pa.pdata <- merge(pa.pdata, pa.sturm.RF.s1[143:209,c("nonTCGA.ID","Age.at.diagnosis.years","Study.Class.original","RF.LGG_GBM.class")],by.x="geo_accession",by.y="nonTCGA.ID")
pa.pdata <- subset(pa.pdata, characteristics_ch1 %in% "tissue: Pilocytic astrocytoma tumor")

sturm.ids <- merge(sturm.ids, pa.sturm.RF.s1[1:142,c("nonTCGA.ID","Age.at.diagnosis.years","Study.Class.original","RF.LGG_GBM.class","Sturm.Death","Sturm.OS.months")],by.x="geo_accession",by.y="nonTCGA.ID")

ids <- data.frame(id = c(as.character(pa.pdata$geo_accession),as.character(mur.RF$geo_accession),as.character(sturm.ids$geo_accession),as.character(turcan.sv$geo_accession)),Study = c(as.character(pa.pdata$Study),as.character(mur.RF$Study),as.character(sturm.ids$Study),as.character(turcan.sv$Study)),IDH.status=c(rep(NA,61),as.character(mur.RF$IDH),as.character(sturm.ids$IDH.status),as.character(turcan.sv$IDH.status)),Histology=c(as.character(pa.pdata$characteristics_ch1),as.character(mur.RF$Diagnosis),as.character(sturm.ids$characteristics_ch1),as.character(turcan.sv$Pathology)),age = c(as.character(pa.pdata$Age.at.diagnosis.years),as.character(mur.RF$Age),as.character(sturm.ids$Age.at.diagnosis.years),as.character(turcan.sv$Age.at.Diagnosis)),OS.months = c(rep(NA,61),as.character(mur.RF$OS.months),as.character(sturm.ids$Sturm.OS.months),as.character(turcan.sv$Follow.Up..Months.)),status = c(rep(NA,61),as.character(mur.RF$Status),as.character(sturm.ids$Sturm.Death),as.character(turcan.sv$Deceased)),grade=c(rep("G1",61),as.character(mur.RF$grade),rep("G4",136),as.character(turcan.sv$Grade)),original.cluster = c(as.character(pa.pdata$Study.Class.original),as.character(mur.RF$CIMP),as.character(sturm.ids$Study.Class.original),as.character(turcan.sv$CIMP.phenotype)),LGm = c(as.character(pa.pdata$RF.LGG_GBM.class),as.character(mur.RF$RFtype),as.character(sturm.ids$RF.LGG_GBM.class),as.character(turcan.sv$cluster.meth.all)))
#all.test.350.t$id <- rownames(all.test.350.t)
#ids <- merge(ids, all.test.350.t[,351:352],by ="id")
##  falta colocar age, pediatric glioma
rownames(ids) <- as.character(ids$id)
ids <- ids[rownames(all.test.212.t),]
#identical(rownames(ids),rownames(all.test.212.t)) #ok
ids <- cbind(ids,all.test.212.t$RFtype)
colnames(ids)[8] <- "RFtype.IDHwt"

levels(ids$Histology) <- c("Astrocytoma", "other/non-glioma","other/non-glioma", 
                           "Glioblastoma", "Oligoastrocytoma", "Oligoastrocytoma",
                           "Oligodendroglioma","Oligoastrocytoma","Oligodendroglioma", 
                           "other/non-glioma", "other/non-glioma", "Pilocytic astrocytoma", "Glioblastoma")
levels(ids$status) <- c(0,1,0,1,NA,0,1)
levels(ids$IDH.status) <- c("Mutant","WT","WT","Mutant")
levels(ids$age)[levels(ids$age) == "ND "] <- NA
levels(ids$OS.months)[levels(ids$OS.months) == "ND "] <- NA
levels(ids$grade) <- c("G2","G3","G4","G1","G4","G2","G3")
ids$age <- as.numeric(as.character(ids$age))
ids$OS.months <- as.numeric(as.character(ids$OS.months))
ids$Histology <- as.character(ids$Histology)


b <- subset(ids, age < 22 & Histology != "Pilocytic astrocytoma")
ids[rownames(b), "Histology"] <- "Pediatric.Glioma"

