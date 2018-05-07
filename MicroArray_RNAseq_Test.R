setwd("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Normalization_Test/RNASeqV2//UNC__IlluminaHiSeq_RNASeqV2/Level_3/")

id.mrna <- read.delim("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Normalization_Test/file_manifest.txt")
rownames(id.mrna) <- id.mrna$File.Name
gbm.norm = NULL
pattern.primary = 'TCGA-([0-9A-Z]{2})-([0-9A-Z]{1,})-([0-1A-Z]{3})-([0-9A-Z]{3})-([0-9A-Z]{4}-([0-9]{2}))' #barcode
files = list.files("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Normalization_Test/RNASeqV2//UNC__IlluminaHiSeq_RNASeqV2/Level_3/") 
z = 1 #first file
while (!is.na(files[z])) {
  if(nchar(as.character(id.mrna[as.character(files[z]),"Barcode"]))<29) {
    aux <- id.mrna[as.character(files[z]),"Barcode"]
    c <- str_extract(as.character(aux), pattern.primary)
    if(!is.na(c)){ #do not read metadata files and non-tumor samples
      if(is.null(gbm.norm)){ #read the first sample
        gbm.norm = read.table(files[z],sep="\t",header=TRUE) 
        colnames(gbm.norm)[2] <- as.character(c)
        gbm.norm <- gbm.norm[,c(1,2)]
      }
      else{
        aux = read.table(files[z], header=TRUE,sep="\t")
        colnames(aux)[2] = c
        gbm.norm = merge(gbm.norm, aux[1:2], by="gene_id")
      }
    }
  }
  z = z + 1 #next file
}#end while
# gbm -> 140 tumor samples


setwd("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Normalization_Test/Expression-Genes//UNC__AgilentG4502A_07_1/Level_3/")

gbm.micro1 = NULL
files = list.files("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Normalization_Test/Expression-Genes//UNC__AgilentG4502A_07_1/Level_3/") 
z = 1 #first file
while (!is.na(files[z])) {
  if(nchar(as.character(id.mrna[as.character(files[z]),"Barcode"]))<29) {
    aux <- id.mrna[as.character(files[z]),"Barcode"]
    c <- str_extract(as.character(aux), pattern.primary)
    if(!is.na(c)){ #do not read metadata files and non-tumor samples
      if(is.null(gbm.micro1)){ #read the first sample
        gbm.micro1 = read.table(as.character(files[z]),skip=2,sep="\t") 
        colnames(gbm.micro1)[2] <- as.character(c)
        colnames(gbm.micro1)[1] <- "Hybridization.REF"
        
      }
      else{
        aux = read.delim(as.character(files[z]))
        colnames(aux)[2] = c
        colnames(aux)[1] <- "Hybridization.REF"
        gbm.micro1 = merge(gbm.micro1, aux, by="Hybridization.REF")
      }
    }
  }
  z = z + 1 #next file
}#end while
# gbm -> 140 tumor samples


setwd("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Normalization_Test/Expression-Genes//UNC__AgilentG4502A_07_2/Level_3/")

gbm.micro2 = NULL
files = list.files("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Normalization_Test/Expression-Genes//UNC__AgilentG4502A_07_2/Level_3/") 
z = 1 #first file
while (!is.na(files[z])) {
  if(nchar(as.character(id.mrna[as.character(files[z]),"Barcode"]))<29) {
    aux <- id.mrna[as.character(files[z]),"Barcode"]
    c <- str_extract(as.character(aux), pattern.primary)
    if(!is.na(c)){ #do not read metadata files and non-tumor samples
      if(is.null(gbm.micro2)){ #read the first sample
        gbm.micro2 = read.table(as.character(files[z]),skip=2) 
        colnames(gbm.micro2)[2] <- as.character(c)
        colnames(gbm.micro2)[1] <- "Hybridization.REF"
        
      }
      else{
        aux = read.delim(as.character(files[z]),,skip=2)
        colnames(aux)[2] = c
        colnames(aux)[1] <- "Hybridization.REF"
        gbm.micro2 = merge(gbm.micro2, aux, by="Hybridization.REF")
      }
    }
  }
  z = z + 1 #next file
}#end while
# gbm -> 140 tumor samples

gbm.micro <- merge(gbm.micro1,gbm.micro2,by="Hybridization.REF")
colnames(gbm.micro) <- substr(colnames(gbm.micro),1,12)
colnames(gbm.norm) <- substr(colnames(gbm.norm),1,12)
gbm.norm$gene_symbol <- str_extract(gbm.norm$gene_id,"*[^|]*")
a <- intersect(colnames(gbm.micro),colnames(gbm.norm))
info <- gbm.micro$Hybridizatio
gbm.micro <- gbm.micro[,as.character(a)]
gbm.micro$Hybridization.REF <- as.character(info)
info <- gbm.norm$gene_symbol
gbm.norm <- gbm.norm[,as.character(a)]
gbm.norm$gene_symbol <- as.character(info)

a <- intersect(as.character(gbm.micro$Hybridization.REF),as.character(gbm.norm$gene_symbol))
gbm.norm <- subset(gbm.norm, gene_symbol %in% a)
gbm.norm <- gbm.norm[!duplicated(gbm.norm$gene_symbol),]
gbm.micro <- subset(gbm.micro, Hybridization.REF %in% a)


#normalize the expression data
gbm.norm.t  <- t(apply(gbm.norm[,-138], 1, function(x) { (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) }))
gbm.norm.t[is.nan(gbm.norm.t)] <- 0 #colunas que soh tem zero ficam com NaN por causa da divisao por zero na normalizacao
gbm.norm.t <- as.data.frame(gbm.norm.t)
gbm.norm.t$gene_symbol <- as.character(gbm.norm$gene_symbol)


gbm.micro[gbm.micro == "null"] <- NA
gbm.micro.t <- sapply(gbm.micro[,-138], function(x) as.numeric(as.character(x)))
gbm.micro.t  <- t(apply(gbm.micro.t, 1, function(x) { (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) }))
gbm.micro.t[is.nan(gbm.micro.t)] <- 0 #colunas que soh tem zero ficam com NaN por causa da divisao por zero na normalizacao
gbm.micro.t <- as.data.frame(gbm.micro.t)
gbm.micro.t$Hybridization.REF <- as.character(gbm.micro$Hybridization.REF)


#### teste
x <- as.numeric(gbm.micro[gbm.micro$Hybridization.REF == "REST",1:137])
y <- as.numeric(gbm.norm[gbm.norm$gene_symbol == "REST",1:137])
x <- (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
y <- (y-min(y,na.rm=T))/(max(y,na.rm=T)-min(y,na.rm=T))
teste2 <- cbind(as.numeric(x),as.numeric(gbm.micro[gbm.micro$Hybridization.REF == "REST",1:137]),as.numeric(gbm.norm[gbm.norm$gene_symbol == "REST",1:137]),as.numeric(y))
colnames(teste2) <- c("Micro Norm", "Micro Raw", "RNA seq Raw", "RNA seq Norm")
teste <- cbind(as.numeric(gbm.micro.t[gbm.micro.t$Hybridization.REF == "A2M",1:137]),as.numeric(gbm.micro[gbm.micro$Hybridization.REF == "A2M",1:137]),as.numeric(gbm.norm[gbm.norm$gene_symbol == "A2M",1:137]),as.numeric(gbm.norm.t[gbm.norm.t$gene_symbol == "A2M",1:137]))
colnames(teste) <- c("Micro Norm", "Micro Raw", "RNA seq Raw", "RNA seq Norm")
png(filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/Normalization_Test/teste.png",res=300,width=2000,height=2000)
heatpairs(teste2)
dev.off()