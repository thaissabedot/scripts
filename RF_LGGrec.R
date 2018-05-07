

# Load required data for training
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_GBM_obj_20150127.Rda") # the methylation matrix
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGGRecurrent_DNAmeth.Rda") # LGGrecurrent - DNA methylation obj
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/gbm.lgg.recurrent.meta.rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/df.gbm.rec.Rda") #df.gbm.rec - DNA methylation obj
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG-GBM.merged.data.FB20150708.Rdata") # pdata

ids <- subset(gbm.lgg.recurrent.meta, cluster.meth.new == "UNK")$Patient.ID

all <- merge(LGGrecurrent, df.gbm.rec[,c(4:27)],by=0) # 32 LGGs + 24 GBMs = 56 samples
rownames(all) <- as.character(all$Row.names)
all <- all[,-1]

# find common probes 
Rec.s <- all[,substr(colnames(all),1,12) %in% as.character(ids) | substr(colnames(all),14,15) == "02"]
LGGrec.1300 <- Rec.s[rownames(dat.lgg.gbm.new.noXY.dic.oecg),]
Rec.1298 <- na.omit(LGGrec.1300)
LGG.GBM.1298 <- LGG.GBM.250914[rownames(Rec.1298),5:936]
# do this however you see fit. Basically choose a probe set, and then make
# sure that you na.omit both training and testing, and then find the common probeset

# the lgg.gbm object below represents this common probeset
# (in this example between sturm and lgg.gbm)


# Transpose such that probe ids are colnames and samples are rownames
trainingdata <- t(LGG.GBM.1298)
# add an additional column of data that you wish to train on
# (for classifying the new data). In this case the new column is called "papercluster"
# and is using the "cluster.meth2" field from the metadata table.
trainingdata <- merge(trainingdata, pd[,c("case.id","cluster.meth")], by.x=0,by.y="case.id")
rownames(trainingdata) <- as.character(trainingdata$Row.names)
trainingdata <- trainingdata[,-1] #932 samples, 1281 probes, 1 cluster.meth


save(trainingdata, file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/trainingdata.rda")
# This must be a factor, and must include ONLY the factor levels that are in
# the subset of data that you wish to use. For example, if your subset of data 
# that is defined by rownames(trainingdata) does not include LGm4, make sure that the
# result of levels(trainingdata$papercluster) does not include LGm4
#trainingdata <- data.frame(trainingdata,
#papercluster=factor(lgg.gbm.meta[rownames(trainingdata), "cluster.meth2"]))


################### START HERE
library(caret)
library(randomForest)
library(doMC)
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/trainingdata.rda")
#teste <- trainingdata[1:31,c(1:10,1282)]
# register cores for doMC
registerDoMC(cores = 8)
# set up k-fold cross validation
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)
# you may additionally, if you wish use a different method, for validating your
# model parameters, such as oob (Out of Bag).  oob is faster.

# Set your seed so your work is repeatable
set.seed(42)
# Create a subset of your data to train your model on.  This makes sure you have
# equal representation of the 'papercluster' groups in your training set
inTraining <- createDataPartition(trainingdata$cluster.meth, p=0.8, list=FALSE, times=1)
# Training Set
myTrain <- trainingdata[inTraining, ]
# Testing Set
myTest <- trainingdata[-inTraining, ]
# Confirm seed is set
set.seed(210)
# set values for mtry
# mtry is the "Number of variables randomly sampled as candidates at each split"
# traditionally for classification you use the sqrt of the number of variables
# but here we try a range of mtry values to find the best parameters for our model
mtryVals <- floor(c(seq(100, 2000, by=100),
                    sqrt(ncol(trainingdata))))
mtryGrid <- data.frame(.mtry=mtryVals)
# Confirm seed again
set.seed(420)
# Set number of cores
registerDoMC(cores = 8)
# Run Training
lgg_gbm.RF.to.Recurrent <- train(cluster.meth ~ ., # variable to be trained on
                                 data = trainingdata, # Data we are using
                                 method = "rf", # Method we are using
                                 trControl = fitControl, # How we validate
                                 # We created this object above
                                 ntree = 5000, # number of trees
                                 # is dependent on training data size
                                 importance = TRUE, # calculate varible importance
                                 # can be omitted to speed up calc
                                 tuneGrid = mtryGrid, # set mtrys
                                 subset = inTraining # define training set
)


############# PAREI AQUI

save(lgg_gbm.RF.to.Recurrent,file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/RF_LGGrec.rda")






# make sure you have your object that you wish to reclassify
# this object must have only the features that you have trained on, and must have
# all the features that are trained on.

# predict the probability with which each sample of the test set is classified
lgg_gbm.RF.to.mur.pred <- predict(lgg_gbm.RF.to.mur, myTest, type="prob")
# record the best classification in the test set
myTest$RFtype <- predict(lgg_gbm.RF.to.mur, myTest)
# show the confusion matrix
confusionMatrix(data = myTest$RFtype, reference = myTest$cluster.meth)

# predict the classification of the new data set (sturm probes in this case)
# this data contains the same common set of probes
mur.46s.1281.t <- t(mur.46s.1281)
lgg_gbm.RF.to.mur.prob <- predict(lgg_gbm.RF.to.mur, mur.46s.1281.t, type="prob")
# record the best classification
mur.46s.t <- t(mur.46s.1281)
mur.46s.t <- as.data.frame(mur.46s.t)
mur.46s.t$RFtype <- predict(lgg_gbm.RF.to.mur, mur.46s.t)
# show how the new data falls into your classification groups
table(sturm.136$RFtype)

turcan.1281 <- turcan.tumor[rownames(LGG.GBM.1281),]
turcan.1281.t <- t(turcan.1281)
turcan.1281.t <- as.data.frame(turcan.1281.t)
#transformar para beta value
turcan.1281.t.bv <- apply(turcan.1281, c(1,2), function(x){y <- exp(x)/(1+exp(x)); return(y)})
turcan.1281.t.bv <- t(turcan.1281.t.bv)

lgg_gbm.RF.to.turcan.prob <- predict(lgg_gbm.RF.to.mur, turcan.1281.t.bv,type="prob") 

# save everything
save(list=ls(), file="sturm.probes.lgg_gbm.alldata.Rda")