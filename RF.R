Random Forest in R
---
  There are a few resources for doing random forests in R.
The implementation we will be using is [`randomForest`](http://www.stat.berkeley.edu/~breiman/RandomForests/%22randomForest%22) from CRAN. In order to allow it to make it run in parallel we use the package [`caret`](http://topepo.github.io/caret/index.html%22caret%22) which implements `doMC` or `doMPI` for parallel backends.  `caret` also provides a consistent interface across many machine learning algorithms, and allows for easy training and tuning of parameters before finalizing the model.

Some example code follows, but one important thing to realize is that the features that are going to be in your test set (the set you wish to reclassify, have to be fully and exclusively in the training set - aka you have to have exactly the same features in both training and testing data sets).

```R
# Load in required packages
library(caret)
library(randomForest)
library(doMC)


# Load required data for training
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG_GBM_obj_20150127.Rda") # the methylation matrix

# load sturm methylation matrix
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/mur.DNA_Meth.Rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Turcan.DNA_Meth.Rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/PA.DNA_Meth.Rda")
load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/Sturm_DNA_meth.rda")

load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/LGG-GBM.merged.data.FB20150918.Rdata")

IDHwt.914p <- read.table("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/p914_IDHwt.txt")
IDHmut.1308p <- read.table("/dados/ResearchProjects/thais/TCGA/LGG.GBM/RandomForest.calls/p1308_IDHmut.txt")

load("/dados/ResearchProjects/thais/TCGA/LGG.GBM/clinical_validation.rda")
IDHwt.s <- subset(ids, IDH.status %in% "WT" | is.na(IDH.status))
IDHmut.s <- subset(ids, IDH.status %in% "Mutant")
mur.meth <- gset.mur.m[,as.character(gset.mur.p$geo_accession)]
# find common probes

all <- merge(mur.meth,pa.meth.nocontrol,by=0)
rownames(all) <- as.character(all$Row.names)
all <- all[,-1]
all <- merge(all, sturm.meth.tumor,by=0)
rownames(all) <- as.character(all$Row.names)
all <- all[,-1]
all <- merge(all, turcan.tumor.bv,by=0)
rownames(all) <- as.character(all$Row.names)
all <- all[,-1]

all.IDHwt <- all[,as.character(IDHwt.s$id)]
all.914 <- all.IDHwt[as.character(IDHwt.914p$V1),]
all.820 <- na.omit(all.914)

aux <- subset(pd, cluster.meth2 %in% c("IDHwt-K1","IDHwt-K2","IDHwt-K3"))
LGG.GBM.1280 <- LGG.GBM.250914[rownames(all.1280),colnames(LGG.GBM.250914) %in% as.character(aux$case.id)]
# do this however you see fit. Basically choose a probe set, and then make
# sure that you na.omit both training and testing, and then find the common probeset

# the lgg.gbm object below represents this common probeset
# (in this example between sturm and lgg.gbm)

aux <- as.character(subset(ids, RF.IDH %in% c("IDHmut-hypo"))$id
a1 <- all.1308[,aux]
a1.s <- hclust(dist(t(a1)))
aux <- as.character(subset(ids, RF.IDH %in% c("IDHmut-nonCodels"))$id)
a2 <- all.1308[,aux]
a2.s <- hclust(dist(t(a2)))
aux <- as.character(subset(ids, RF.IDH %in% c("IDHmut-codels"))$id)
a3 <- all.1308[,aux]
a3.s <- hclust(dist(t(a3)))

order <- cbind(a1[,a1.s$order],a2[,a2$order],a3[,a3$order])


#order <- rbind(ESg1[ESg1.p$order,],ESg2[ESg2.p$order,],ESg3[ESg3.p$order,],ESg4[ESg4.p$order,],ESg5[ESg5.p$order,])
#a <- metadata.1129.samples.20150312[colnames(order),]
order <- rbind(ESg5[ESg5.p$order,],ESg4[ESg4.p$order,],ESg3[ESg3.p$order,],ESg2[ESg2.p$order,],ESg1[ESg1.p$order,])

# Transpose such that probe ids are colnames and samples are rownames
trainingdata <- t(LGG.GBM.820)
# add an additional column of data that you wish to train on
# (for classifying the new data). In this case the new column is called "papercluster"
# and is using the "cluster.meth2" field from the metadata table.
trainingdata <- merge(trainingdata, pd[,c("case.id","cluster.meth2")], by.x=0,by.y="case.id")
rownames(trainingdata) <- as.character(trainingdata$Row.names)
trainingdata <- trainingdata[,-1] #932 samples, 1281 probes, 1 cluster.meth
trainingdata$cluster.meth2 <- factor(trainingdata$cluster.meth2)

save(trainingdata,file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/RF_training_IDHwt_914.rda")

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
load("RF_training_IDHmut_1308.rda")
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
inTraining <- createDataPartition(trainingdata$cluster.meth2, p=0.8, list=FALSE, times=1)
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
lgg_gbm.RF.to.IDHwt.p <- train(cluster.meth2 ~ ., # variable to be trained on
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
lgg_gbm.RF.to.IDHwt.p.pred <- predict(lgg_gbm.RF.to.IDHwt.p, myTest, type="prob")
# record the best classification in the test set
myTest$RFtype <- predict(lgg_gbm.RF.to.IDHwt.p, myTest)
# show the confusion matrix
confusionMatrix(data = myTest$RFtype, reference = myTest$cluster.meth2)

# predict the classification of the new data set (sturm probes in this case)
# this data contains the same common set of probes
all.820.t <- t(all.820)
lgg_gbm.RF.to.IDHwt.p <- predict(lgg_gbm.RF.to.IDHwt.p, all.1280.t, type="prob")
# record the best classification
all.820.t <- as.data.frame(all.820.t)
all.820.t$RFtype <- predict(lgg_gbm.RF.to.IDHwt.p, all.820.t)

all.1280.t$RFtype <- predict(lgg_gbm.RF.to.IDHmut.p, all.1280.t)
# show how the new data falls into your classification groups
table(all.820.t$RFtype)

ids$RF.IDH <- NA
ids[rownames(all.820.t),"RF.IDH"] <- as.character(all.820.t$RFtype)
ids[rownames(all.1280.t),"RF.IDH"] <- as.character(all.1280.t$RFtype)
ids$RF.IDH <- as.factor(ids$RF.IDH)

turcan.1281 <- turcan.tumor[rownames(LGG.GBM.1281),]
turcan.1281.t <- t(turcan.1281)
turcan.1281.t <- as.data.frame(turcan.1281.t)
#transformar para beta value
turcan.1281.t.bv <- apply(turcan.1281, c(1,2), function(x){y <- exp(x)/(1+exp(x)); return(y)})
turcan.1281.t.bv <- t(turcan.1281.t.bv)

lgg_gbm.RF.to.turcan.prob <- predict(lgg_gbm.RF.to.mur, turcan.1281.t.bv,type="prob") 

# save everything
save(list=ls(), file="sturm.probes.lgg_gbm.alldata.Rda")
```

Running this script
---
  In order to run this script you need save it as (filename).R and then make a bash type script to run it.
This is what the bash script looks like.
```sh
#!/bin/bash

#PBS -q long
#PBS -l nodes=1:ppn=40

R CMD BATCH --no-restore --max-ppsize=500000 /vol3/fmrp/updated_lgg/LGG_NORM/LGG_11977probes/sturm.R /vol3/fmrp/updated_lgg/LGG_NORM/LGG_11977probes/sturm.out
```
The form of the R command is `R CMD BATCH` to run it non interactively `--no-restore` so that it doesn't load any old data.  `--max-ppsize=500000` so that it can use maximum memory and use maximum number of pointers. `/vol3/fmrp/updated_lgg/LGG_NORM/LGG_11977probes/sturm.R` is the R file above. `/vol3/fmrp/updated_lgg/LGG_NORM/LGG_11977probes/sturm.out` is the output of any messages created, basically what you would see if you ran the script interactively.
The two `#PBS` comments are not normal comments, but are read by torque to tell it how many cores and how much time to devote to the job.
Save this as `job.sh` (or any descriptive name)  then submit using the `qsub` command.
```sh
qsub job.sh
```
To see what is happening with your jobs you can do this:
```sh
# show submitted jobs
qstat -l
# see detailed information about a particular job
qstat -f 697.bm-fmrp-1 # example job
```
Make certain all paths are referenced absolutely. both in the R script and in the torque script. This is essential since the command may execute on a different machine.
