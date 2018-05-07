load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/MicroArray_Affy/glio20140927FRMAData.RData")
require(hgu133a2.db)
geneName<-unlist(mget(colnames(geData), hgu133a2SYMBOL))
#geneID <- unlist(mget(colnames(geData), hgu133a2ENTREZID))
geneName <- as.data.frame(geneName)
geData <- t(geData)
geData.map <- merge(geData, as.data.frame(geneName), by=0) #dados de expressao genica com anotação das probes
rownames(geData.map) <- as.character(geData.map$Row.names)
colnames(geData.map)[1] <- "ProbeID" #529 patients 2 metadata columns



