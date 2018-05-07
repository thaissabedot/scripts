recurrents = NULL
z = 1 #first file
pattern = 'TCGA-([0-9A-Z]{2})-([0-9A-Z]{1,})-([0-9A-Z]{3})-([0-9A-Z]{3})-([0-9A-Z]{4}-([0-9]{2}))' #barcode
setwd("/dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09/Recurrents")
files <- list.files()
while (!is.na(files[z])) {
  c <- str_extract(as.character(files[z]), pattern)
  if(!is.na(c)){ #do not read metadata files and non-tumor samples
    if(is.null(recurrents)){ #read the first sample
      recurrents =  read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1)
      colnames(recurrents)[2] <- as.character(c)
      recurrents <- recurrents[,c(1,3,4,5,2)]
    }
    else{
      aux =  read.table(files[z], header=TRUE,quote="",comment="",sep="\t", skip=1)
      colnames(aux)[2] = c
      recurrents = merge(recurrents, aux[1:2],by ="Composite.Element.REF")
    }
  }
  z = z + 1 #next file
}#end while

load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_Recurrent_RNAseqV2/gbm.lgg.recurrent.meta.rda")


rownames(recurrents) <- as.character(recurrents$Composite.Element.REF)
#find commom probes
recurrents.1300 <- recurrents[rownames(dat.lgg.gbm.new.noXY.dic.oecg),-c(1:4)]
LGG.GBM.1300 <- LGG.GBM.250914[rownames(recurrents.1300),5:936]

load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG-GBM.merged.data.FB20150918.Rdata")

trainingdata <- t(LGG.GBM.1300)
# add an additional column of data that you wish to train on
# (for classifying the new data). In this case the new column is called "papercluster"
# and is using the "cluster.meth2" field from the metadata table.
trainingdata <- merge(trainingdata, pd[,c("case.id","cluster.meth")], by.x=0,by.y="case.id")
rownames(trainingdata) <- as.character(trainingdata$Row.names)
trainingdata <- trainingdata[,-1] #932 samples, 1281 probes, 1 cluster.meth
trainingdata$cluster.meth <- factor(trainingdata$cluster.meth)

save(trainingdata,file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/RF_1300p_932LGGnGBM.rda")

library(caret)
library(randomForest)
library(doMC)
load("RF_1300p_932LGGnGBM.rda")
#teste <- trainingdata[1:31,c(1:10,1282)]
# register cores for doMC
registerDoMC(cores = 12)
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
registerDoMC(cores = 12)
# Run Training
lgg_gbm.RF.to.Rec.p <- train(cluster.meth ~ ., # variable to be trained on
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

# make sure you have your object that you wish to reclassify
# this object must have only the features that you have trained on, and must have
# all the features that are trained on.

# predict the probability with which each sample of the test set is classified
lgg_gbm.RF.to.Rec.p.pred <- predict(lgg_gbm.RF.to.Rec.p, myTest, type="prob")
# record the best classification in the test set
myTest$RFtype <- predict(lgg_gbm.RF.to.Rec.p, myTest)
# show the confusion matrix
confusionMatrix(data = myTest$RFtype, reference = myTest$cluster.meth)