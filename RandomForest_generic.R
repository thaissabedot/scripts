
# Load in required R packages
library(caret)
library(randomForest)
library(doMC)


# Load required data for training
load("LGG_GBM_Cell_paper.Rda") #obj: LGG.GBM

# load data you wish to apply the model
load("new_gliomas.Rda") #obj: new.gliomas


# find common probes. Can't have NAs
identical(nrow(na.omit(LGG.GBM)),nrow(na.omit(new.gliomas))) #must be TRUE. If not, subset the data accordingly

identical(rownames(LGG.GBM),rownames(new.gliomas)) #They must also be in the same order

# Transpose such that probe ids are colnames and samples are rownames
trainingdata <- t(LGG.GBM)

# add an additional column of data that you wish to train on
# (for classifying the new data). In this case the new column is called "cartoon"
load("metadata.rda")
trainingdata <- merge(trainingdata, pd[,c("case.id","cartoon")], by.x=0,by.y="case.id")
rownames(trainingdata) <- as.character(trainingdata$Row.names)
trainingdata <- trainingdata[,-1] #X samples, Y probes, 1 metadata

# The metadata column must be a factor, and must include ONLY the factor levels that are in
# the subset of data that you wish to use. For example, if your subset of data 
# that is defined by rownames(trainingdata) does not include LGm4, make sure that the
# result of levels(trainingdata$papercluster) does not include LGm4
trainingdata$cartoon <- factor(trainingdata$cartoon)

save(trainingdata,file="RF_training.rda")


# register cores for doMC
registerDoMC(cores = 8) #Number of CPUs

# set up k-fold cross validation
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

# Set your seed so your work is repeatable
set.seed(42)

# Create a subset of your data to train your model on.  This makes sure you have
# equal representation of the 'cartoon' groups in your training set
inTraining <- createDataPartition(trainingdata$cartoon, p=0.8, list=FALSE, times=1) #using 80% of data to train and the remaining 20% for testing

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
registerDoMC(cores = 8) #Number of CPUs

# Run Training
lgg_gbm.RF <- train(cartoon ~ ., # variable to be trained on
                           data = trainingdata, # Data we are using
                           method = "rf", # Method we are using
                           trControl = fitControl, # How we validate
                           ntree = 5000, # number of trees is dependent on training data size
                           importance = TRUE, # calculate varible importance
                           tuneGrid = mtryGrid, # set mtrys
                           subset = inTraining # define training set
)

# make sure you have your object that you wish to reclassify
# this object must have only the features that you have trained on, and must have
# all the features that are trained on.

# predict the probability with which each sample of the test set is classified
lgg_gbm.RF.pred <- predict(lgg_gbm.RF, myTest, type="prob")

# record the best classification in the test set
myTest$RF.classification <- predict(lgg_gbm.RF, myTest)

# show the confusion matrix
confusionMatrix(data = myTest$RF.classification, reference = myTest$cartoon)

# predict the classification of the new data set (new gliomas in this case)
# this data contains the same common set of probes
new.gliomas <- t(new.gliomas)
lgg_gbm.RF.pred <- predict(lgg_gbm.RF, new.gliomas, type="prob")

# record the best classification
new.gliomas$RF.classification <- predict(lgg_gbm.RF, new.gliomas)

# show how the new data falls into your classification groups
table(new.gliomas$RF.classification)


# save everything
save(list=ls(), file="RF_LGG.GBM.Rda")