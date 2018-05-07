library(xlsx)
a <- read.xlsx2("/media/data1/thais.projects/Glioma.Rec/HFH_Female_Deceased_020518.xls",1,startRow=2)
#a[duplicated(as.character(a$Name)),] #6 patients
a[as.character(a$Name) %in% "Jenovai, Kathleen ",]
a <- a[-522,]
sum(duplicated(as.character(a$MRN))) #0
b <- grep(".*[0-9].*",as.character(a$HF.))
a <- a[b,]
others <- levels(a$Histology)[-c(1,2,3,6,7,9,10,11,19,20,29,30,31,32)]
a <- a[a$Histology %in% levels(a$Histology)[c(1,2,3,6,7,8,9,10,11,18,19,20,21,29,30,31,32)],]
a$Histology <- factor(a$Histology)
levels(a$Histology)
aux <- strsplit(as.character(a$HF.),"\n")
a$MRN <- as.character(a$MRN)
matches <- regmatches(a$MRN, gregexpr("[[:digit:]]+", a$MRN))
a$MRN <- matches
names(aux) <- (matches)


list.2 <- unlist(lapply(aux, function(x) {length(x) > 1}))
b <- aux[list.2]
length(b) #87 patients
c <- unlist(b)
length(c) #210 samples
duplicated(c) 
c[duplicated(c)]

a <- a[a$MRN %in% names(b),]
aux <- strsplit(as.character(a$Tissue.Amount),"\n")
aux2 <- strsplit(as.character(a$HF.),"\n")
aux <- unlist(aux)
names(aux2) <- as.character(a$HF)
aux2 <- melt(aux2)
aux2$value <- sub("^0","",aux2$value)
tissue <- data.frame(HF.ID=aux2$value,Tissue.amount=aux,Status="Recurrent",Primary=aux2$L1,stringsAsFactors = FALSE)
tissue[tissue$HF.ID %in% tissue$Primary,"Status"] <- "Primary"

write.xlsx(a[,c(2,21,22,23)], "/media/data1/thais.projects/Glioma.Rec/female_deceased_87p_210s.xlsx", sheetName="Sheet1",col.names=TRUE, row.names=F, append=FALSE)

# 101 patients - 243 samples (one patient has the same HF ID twice (MRN_57295282))
# GLIOMAS 87 patients - 210 samples



library(xlsx)
a <- read.xlsx2("/media/data1/thais.projects/Glioma.Rec/HFH_Female_Alive_020518.xls",1,startRow=2)
#a[duplicated(as.character(a$Name)),] #40 patients
sum(duplicated(as.character(a$MRN))) #0
b <- grep(".*[0-9].*",as.character(a$HF.))
a <- a[b,]
others <- c(others,levels(a$Histology)[-c(1,2,3,6,7,10,11,32,34,46,47,48,49)])
a <- a[a$Histology %in% levels(a$Histology)[c(1,2,3,6,7,9,10,11,13,30,32,33,34,35,46,47,48,49)],]
a$Histology <- factor(a$Histology)
levels(a$Histology)
aux <- strsplit(as.character(a$HF.),"\n")
a$MRN <- as.character(a$MRN)
matches <- regmatches(a$MRN, gregexpr("[[:digit:]]+", a$MRN))
a$MRN <- matches
names(aux) <- (matches)


list.2 <- unlist(lapply(aux, function(x) {length(x) > 1}))
b <- aux[list.2]
length(b) #42 patients
c <- unlist(b)
length(c) #93 samples but #one patient has the same HF ID twice (MRN_53945380)
duplicated(c) 
c[duplicated(c)]

a <- a[a$MRN %in% names(b),]
aux <- strsplit(as.character(a$Tissue.Amount),"\n")
aux2 <- strsplit(as.character(a$HF.),"\n")
aux <- unlist(aux)
names(aux2) <- as.character(a$HF)
aux2 <- melt(aux2)
aux2$value <- sub("^0","",aux2$value)
aux <- data.frame(HF.ID=aux2$value,Tissue.amount=aux,Status="Recurrent",Primary=aux2$L1,stringsAsFactors = FALSE)
aux[aux$HF.ID %in% aux$Primary,"Status"] <- "Primary"
tissue <- rbind(tissue,aux)

write.xlsx(a[,c(2,21,22,23)], "/media/data1/thais.projects/Glioma.Rec/female_alive_42p_93s.xlsx", sheetName="Sheet1",col.names=TRUE, row.names=F, append=FALSE)

# 77 patients - 168 samples
# GLIOMAS 38 patients - 84 samples #one duplicate

library(xlsx)
a <- read.xlsx2("/media/data1/thais.projects/Glioma.Rec/HFH_Male_Deceased_020518.xls",1,startRow=2)
#a[duplicated(as.character(a$Name)),] #11 patients
sum(duplicated(as.character(a$MRN))) #0
b <- grep(".*[0-9].*",as.character(a$HF.))
a <- a[b,]
others <- c(others,levels(a$Histology)[-c(1,3,4,7,8,9,10,11,17,21,23,34,35,36,37)])
a <- a[a$Histology %in% levels(a$Histology)[c(1,2,3,4,7,8,9,10,11,17,20,21,22,23,24,34,35,36,37)],]
a$Histology <- factor(a$Histology)
levels(a$Histology)
aux <- strsplit(as.character(a$HF.),"\n")
a$MRN <- as.character(a$MRN)
matches <- regmatches(a$MRN, gregexpr("[[:digit:]]+", a$MRN))
a$MRN <- matches
names(aux) <- (matches)


list.2 <- unlist(lapply(aux, function(x) {length(x) > 1}))
b <- aux[list.2]
length(b) #184 patients
c <- unlist(b)
length(c) #417 samples
duplicated(c) 
c[duplicated(c)] #0

a <- a[a$MRN %in% names(b),]
aux <- strsplit(as.character(a$Tissue.Amount),"\n")
aux2 <- strsplit(as.character(a$HF.),"\n")
aux <- unlist(aux)
names(aux2) <- as.character(a$HF)
aux2 <- melt(aux2)
aux2$value <- sub("^0","",aux2$value)
aux <- data.frame(HF.ID=aux2$value,Tissue.amount=aux,Status="Recurrent",Primary=aux2$L1,stringsAsFactors = FALSE)
aux[aux$HF.ID %in% aux$Primary,"Status"] <- "Primary"
tissue <- rbind(tissue,aux)

write.xlsx(a[,c(2,21,22,23)], "/media/data1/thais.projects/Glioma.Rec/male_deceased_184p_417s.xlsx", sheetName="Sheet1",col.names=TRUE, row.names=F, append=FALSE)

# 201 patients - 458 samples
# GLIOMAS 177 patients - 402 samples

library(xlsx)
a <- read.xlsx2("/media/data1/thais.projects/Glioma.Rec/HFH_Male_Alive_020518.xls",1,startRow=2)
#a[duplicated(as.character(a$Name)),] #52 patients
sum(duplicated(as.character(a$MRN))) #0
b <- grep(".*[0-9].*",as.character(a$HF.))
a <- a[b,]
others <- c(others,levels(a$Histology)[-c(1,3,4,7,8,12,13,14,27,34,35,48,49,50,51)])
a <- a[a$Histology %in% levels(a$Histology)[c(1,2,3,4,7,8,9,11,12,13,14,16,27,32,34,35,36,48,49,50,51)],]
a$Histology <- factor(a$Histology)
levels(a$Histology)
aux <- strsplit(as.character(a$HF.),"\n")
a$MRN <- as.character(a$MRN)
matches <- regmatches(a$MRN, gregexpr("[[:digit:]]+", a$MRN))
a$MRN <- matches
names(aux) <- (matches)



list.2 <- unlist(lapply(aux, function(x) {length(x) > 1}))
b <- aux[list.2]
length(b) #56 patients
c <- unlist(b)
length(c) #138 samples
duplicated(c) 
c[duplicated(c)] #0

a <- a[a$MRN %in% names(b),]
aux <- strsplit(as.character(a$Tissue.Amount),"\n")
aux2 <- strsplit(as.character(a$HF.),"\n")
aux <- unlist(aux)
names(aux2) <- as.character(a$HF)
aux2 <- melt(aux2)
aux2$value <- sub("^0","",aux2$value)
aux <- data.frame(HF.ID=aux2$value,Tissue.amount=aux,Status="Recurrent",Primary=aux2$L1,stringsAsFactors = FALSE)
aux[aux$HF.ID %in% aux$Primary,"Status"] <- "Primary"
tissue <- rbind(tissue,aux)

write.xlsx(a[,c(2,21,22,23)], "/media/data1/thais.projects/Glioma.Rec/male_alive_56p_138s.xlsx", sheetName="Sheet1",col.names=TRUE, row.names=F, append=FALSE)


# 85 patients - 200 samples
# GLIOMAS 55 patients - 136 samples

others <- unique(others)
write.table(others,quote=F,row.names=F,col.names=F,sep="\t",file="Levels_Histology.txt")

tissue[!duplicated(tissue$Primary),]
tissue[!duplicated(tissue$Primary),"Status"] <- "Primary"
tissue[tissue$Primary %in% c(2194,2211,2435),]
tissue[c(158,489,816),"Status"] <- "Recurrent"

#tissue$Primary <- as.factor(tissue$Primary)
for(i in levels(factor(tissue$Primary))){
  aux <- tissue[tissue$Primary == i & !is.na(tissue$Primary),]
  p <- sum(as.numeric(aux[aux$Status %in% "Primary","Tissue.amount"]) > 5,na.rm=T)
  r <- sum(as.numeric(aux[aux$Status %in% "Recurrent","Tissue.amount"]) > 5,na.rm=T)
  if(p > 0 & r > 0) print(i)
}
ids <- c("1295","1325","1691","1789","250","2666","2742","285","3381","408","442","936")
aux <- tissue[tissue$Primary %in% ids,]
write.table(aux,quote=F,row.names=F,col.names=T,sep="\t",file="/media/data1/thais.projects/Glioma.Rec/Sample_morethan5.txt")