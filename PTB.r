ssh cloudsandiego
tmux
cd /dados/ResearchProjects/thais/TCGA/LGG.GBM/NewSamples_16-09
R
load("lgg.Rda")

######## NAO USAR ESSE
teste <- c()
for(i in 1:nrow(lgg.450k)){
    obj <- lgg.450k[i,5:520]
    if(sum(is.na(obj)) == ncol(obj))
        teste <- c(teste,TRUE)
    else
        teste <- c(teste,FALSE)
    print(i)
   

}
####### USAR ESSE
numNAs <- apply(lgg.450k[,5:520], 1, function(z) sum(is.na(z)))
teste <- numNAs == 516

tmux a #sandiego
anno450k <- cbind(lgg.450k[,1:4],data.frame(SNPs=teste))
save(anno450k,file="anno450k.rda")

#baixar o arquivo anno450k.rda para o meu computador e transferir para seattle
rsync -av -e ssh --progress anno450k.rda tsabedot@gateway.systemsbiology.net:/home/ISB/tsabedot
rsync -av -e ssh --progress anno450k.rda tsabedot@romulus:/users/tsabedot

ssh gateway
ssh romulus
load("/titan/ITMI1/projects/gamcop/data/methylation/products/data_METH_2014-11-25_raw")
clinical <- read.table("/titan/ITMI1/workspaces/users/selasady/FM_PW/df5/2015_02_25_df5/2015_02_25_clinical_admx_hilevel_allBlood_output/2015_02_25_clinical_admx_hilevel_allBlood_hilevel.fm",header=T)
rownames(methsub) <- as.character(methsub$.)
methsub <- methsub[,-1]
load("anno450k.rda")

anno450k.noSNPs <- subset(anno450k,SNPs == FALSE)
PBT.noSNPs <- methsub[rownames(methsub) %in% as.character(anno450k.noSNPs$Composite.Element.REF),]
PBT.noSNPs.noNA <- na.omit(PBT.noSNPs)
PBT.noSNPs.noNA$sd <- apply(PBT.noSNPs.noNA,1,sd,na.rm=TRUE)
png(filename="hist_SD2.png")
hist(PBT.noSNPs.noNA$sd)
dev.off()

rsync -av -e ssh --progress metadata.rda tsabedot@gateway.systemsbiology.net:/home/ISB/tsabedot
scp gateway:/home/ISB/tsabedot/metadata.rda .

PTB.128p <- subset(PBT.noSNPs.noNA,sd > 0.25) #antes 1.5
PTB.128p <- PTB.128p[,-786] #remove a coluna de SD
save(PTB.128p,file="PBT_128p.rda")

library(ConsensusClusterPlus)
title="/users/tsabedot/cc_128p"
setwd(title)
cc.out.ptb.128 <- ConsensusClusterPlus(		     d=as.matrix(PTB.128p), 
                                                     maxK=10, 
                                                     reps=1000, 
                                                     clusterAlg="km",
                                                     distance="euclidean",
                                                     seed=1022,
                                                     plot="pdf",
						     title=title,
                                                     verbose = TRUE)
ccICL = calcICL(cc.out.ptb.128, plot='pdf')
save(cc.out.ptb.128,file="cc.128.rda")

cc.out.ptb.128.2 <- cbind(as.matrix(cc.out.ptb.128[[2]][[3]]), 
                                        as.matrix(cc.out.ptb.128[[3]][[3]]), 
                                        as.matrix(cc.out.ptb.128[[4]][[3]]), 
                                        as.matrix(cc.out.ptb.128[[5]][[3]]), 
                                        as.matrix(cc.out.ptb.128[[6]][[3]]),
                                        as.matrix(cc.out.ptb.128[[7]][[3]]),
                                        as.matrix(cc.out.ptb.128[[8]][[3]]),
                                        as.matrix(cc.out.ptb.128[[9]][[3]]),
                                        as.matrix(cc.out.ptb.128[[10]][[3]]))

cc.out.ptb.128.2 <- as.data.frame(cc.out.ptb.128.2)
names(cc.out.ptb.128.2) <- c("K2",  "K3",  "K4" , "K5" , "K6" , "K7",  "K8" , "K9" , "K10")
#cc.out.ptb.128.2$ID <- rownames(cc.out.ptb.128.2)
rownames(cc.out.ptb.128.2) <- substr(rownames(cc.out.ptb.128.2),1,8)

rownames(clinical) <- as.character(clinical[,1])
clinical <- clinical[,-1]
clinical <- t(clinical)
clinical <- as.data.frame(clinical)
clinical$ID <- rownames(clinical)
rownames(clinical) <- substr(rownames(clinical),1,8)
clinical.s <- clinical[rownames(clinical) %in% substr(colnames(PTB.128p),1,8),]

metadata <- merge(cc.out.ptb.128.2, clinical.s, by=0)
rownames(metadata) <- as.character(metadata$Row.names)
metadata <- metadata[,-1]
save(metadata,file="metadata")

PTB.order <- PTB.128p[,cc.out.ptb.128[[4]][[2]][[3]]]
info <- metadata[substr(colnames(PTB.order),1,8),]
clab <- info
clab$K2 <- as.factor(clab$K2)
clab$K3 <- as.factor(clab$K3)
clab$K4 <- as.factor(clab$K4)
clab$K5 <- as.factor(clab$K5)
levels(clab$K2) <- c("red","blue")
levels(clab$K3) <- c("red","blue","green")
levels(clab$K4)<- c("red","blue","green","orange")
levels(clab$K5) <- c("red","green","darkblue","orange","violet") #
source("heatmap.plus.R")
library(matlab)
a <- cbind(as.character(clab$K4),as.character(clab$K4))
png(filename="Heatmap_4k.png",res=300,width=1500,height=1500)
heatmap.plus.sm(as.matrix(PTB.order),
                col = jet.colors(75),
                scale = "none",
                labRow = NA,
                labCol = NA,
                Colv = NA,
               # Rowv = NA,
                ColSideColors = a
                #RowSideColors = rlab
                
)
dev.off()


qplot(factor(K4), data=metadata, geom="bar", fill=factor(orig.m))
library(scales)
ggplot(data=metadata,aes(x=factor(K2),fill=factor(orig.m))) + geom_bar(position="fill") + scale_y_continuous(labels = percent) + labs(x="Clusters",y="Percentage",fill = "Ancestry")


load("metadata.rda")
colnames(metadata)[26] <- "ancestry_m"
metadata$clusters.5 <- as.factor(metadata$K5)  
levels(metadata$clusters.5) <- c("K2","K5","K4","K1","K3")
metadata$ID_mother <- paste(substr(rownames(metadata),1,8),".M",sep="")
save(metadata,file="metadata.rda")


## Random Forest (K5) - Classify the "mixed" samples
load("PBT_128p.rda")
load("metadata.rda")

aux <- subset(metadata,ancestry_m != "mixed")$ID_mother
knw.clas.PTB <- PTB.128p[,colnames(PTB.128p) %in% as.character(aux)]
aux <- subset(metadata,ancestry_m == "mixed")$ID_mother
unknw.clas.PTB <- PTB.128p[,colnames(PTB.128p) %in% as.character(aux)]

library(caret)
library(randomForest)
library(doMC)

trainingdata <- t(knw.clas.PTB)
trainingdata <- merge(trainingdata, metadata[,c("ID_mother","ancestry_m")], by.x=0,by.y="ID_mother")
rownames(trainingdata) <- as.character(trainingdata$Row.names)
trainingdata <- trainingdata[,-1] # 617 samples, 128 probes, 1 ancestry_m
trainingdata$ancestry_m <- factor(trainingdata$ancestry_m)
save(trainingdata, file="RF_trainingdata.rda")

load("RF_trainingdata.rda")
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
inTraining <- createDataPartition(trainingdata$ancestry_m, p=0.8, list=FALSE, times=1)
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
RF.ancestry.mother <- train(ancestry_m ~ ., # variable to be trained on
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
save(RF.ancestry.mother,file="RF.ancestry.mother.rda")
RF.ancestry.mother.pred <- predict(RF.ancestry.mother, myTest, type="prob")
# record the best classification in the test set
myTest$RFtype <- predict(RF.ancestry.mother, myTest)
# show the confusion matrix
confusionMatrix(data = myTest$RFtype, reference = myTest$ancestry_m)

# predict the classification of the new data set (sturm probes in this case)
# this data contains the same common set of probes
unknw.clas.PTB.t <- t(unknw.clas.PTB)
RF.ancestry.mother.prob <- predict(RF.ancestry.mother, unknw.clas.PTB.t, type="prob")
# record the best classification
unknw.clas.PTB.t <- as.data.frame(unknw.clas.PTB.t)
unknw.clas.PTB.t$RFtype <- predict(RF.ancestry.mother, unknw.clas.PTB.t)
# show how the new data falls into your classification groups
table(unknw.clas.PTB.t$RFtype)
RF.mixed <- data.frame(ID=rownames(unknw.clas.PTB.t),RF=unknw.clas.PTB.t$RFtype)
metadata <- merge(metadata, RF.mixed,by.x="ID_mother",by.y="ID",all.x=T)
rownames(metadata) <- as.character(metadata$ID_mother)
save(metadata,file="metadata.rda")
