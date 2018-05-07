##Normalized data using Methylumi (Thais did this part)
print(load("/dados/camila/Radiation_MS/Validation_Set_1/Costello_Data/Mazor_Cancer_Cell_2015/Dataframe/Normalized_Results/Costello_Validation_71samples_Methylumi.rda"))  ##Costello.Validation.bv  ##data.frame  ##485577 rows and 76 columns
rownames(Costello.Validation.bv) <- Costello.Validation.bv$Composite.Element.REF

##Link Patient_ID with Sample_ID
setwd("/dados/camila/Radiation_MS/Validation_Set_1/Costello_Data/Mazor_Cancer_Cell_2015/Molecular_Clinical_Features")
Clinical_Data <- read.delim("Costello_EGA_methylation450k_CancerCell_2015.csv")  ##data.frame  ##71 rows and 10 columns
Clinical_Data$Colnames <- paste("ega-box-68_",substr(Clinical_Data$files,1,17),sep = "")
Clinical_Data_s <- subset(Clinical_Data,SampleAttrPatID %in% c("Patient01","Patient02","Patient03","Patient04","Patient07","Patient08","Patient10","Patient11","Patient12","Patient13","Patient14","Patient16","Patient17","Patient18","Patient22","Patient36","Patient38","Patient68","Patient90","Control"))
Costello_Validation <- Costello.Validation.bv[,as.character(Clinical_Data_s$Colnames)]  ##485577 rows and 65 columns
Costello_Validation_s <- Costello_Validation  ##485577 rows and 65 columns

identical(colnames(Costello_Validation_s),Clinical_Data_s$Colnames)

a <- as.matrix(c("Patient01_InitialA","Patient01_InitialB","Patient01_InitialC","Patient01_InitialD","Patient01_Recurrence1A","Patient01_Recurrence1B","Patient01_Recurrence1C","Patient02_Initial","Patient02_Recurrence1","Patient03_Initial","Patient03_Recurrence1","Patient04_InitialA","Patient04_InitialB","Patient04_InitialC","Patient04_InitialD","Patient04_InitialE","Patient04_InitialF","Patient04_Recurrence1","Patient04_Recurrence2","Patient04_Recurrence3","Patient07_Initial","Patient07_Recurrence1","Patient08_Initial","Patient08_Recurrence1","Patient10_Initial","Patient10_Recurrence1","Patient11_Initial","Patient11_Recurrence1","Patient12_Initial","Patient12_Recurrence1","Patient13_Initial","Patient13_Recurrence1","Patient14_Initial","Patient14_Recurrence1","Patient16_Initial","Patient16_Recurrence1","Patient17_InitialA","Patient17_InitialB","Patient17_InitialC","Patient17_Recurrence1A","Patient17_Recurrence1B","Patient17_Recurrence1C","Patient17_Recurrence1D","Patient18_InitialA","Patient18_InitialB","Patient18_InitialC","Patient18_InitialD","Patient18_Recurrence1","Patient22_Initial","Patient22_Recurrence1","Patient36_Initial","Patient36_Recurrence1","Patient38_Initial","Patient38_Recurrence1","Patient68_Initial","Patient68_Recurrence1","Patient90_InitialA","Patient90_InitialB","Patient90_InitialC","Patient90_InitialD","Patient90_InitialE","Patient90_InitialF","Patient90_Recurrence1A","Patient90_Recurrence1B","A141_Adult_Insula"))
cbind(Clinical_Data_s[,c("Colnames","SampleTitle")],a)

colnames(Costello_Validation_s) <- as.character(Clinical_Data_s$SampleTitle)

Costello_Validation_Rad <- Costello_Validation_s[rownames(order.row),]
Costello_Validation_Rad <- Costello_Validation_Rad[,c(65,21,22,49,50,27,28,1:20,23:26,29:48,51:64)]

##Heatmap
setwd("/dados/camila/Radiation_MS/Validation_Set_1/Costello_Data/Mazor_Cancer_Cell_2015/Check_Costello_data")
pdf("teste2.pdf")
heatmap.plus.sm(as.matrix(Costello_Validation_Rad),
                col = jet.colors(75),
                scale = "none",
                # labCol = colnames(Costello_Validation_Rad),
                Colv = NA,
                cexCol = 0.3,
                Rowv = NA
                #ColSideColors = as.matrix(Costello_Labels[colnames(Costello_Validation_Rad),c(2:8,15,16)]),
                #RowSideColors = rlab.rad.rec.cg)
)
dev.off()