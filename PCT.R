load("PCT.rda") #Objects: cpg.gene, cpg.gene.normal, patient.score, patient.call

#cpg.gene 8:643 -> DNA methylation 
#         644:1279 -> expression
cor.cpg.gene <- cpg.gene[,c(2,6)]
cor.cpg.gene$cor <- apply(cpg.gene[,8:1279], 1, function(i) {
  #zz <- cor(i[1:636],i[637:1272])
  zz <- cor(i[637:1272],i[1:636])
  return(zz)
})
cor.cpg.gene <- na.omit(cor.cpg.gene)
## first, categorize each gene and probe to three classes
## SNC (strong negative correlation), WNC (weak negative correlation), NNC (no negative correlation).
cor.cpg.gene$correlation_call <- c("UNK")
for(i in 1:dim(cor.cpg.gene)[1]){
  if(!is.na(cor.cpg.gene[i,"cor"])){
    if(cor.cpg.gene[i,"cor"] < c(-0.4)){
      cor.cpg.gene[i,"correlation_call"]<-c("SNC")
    } else {
      if(cor.cpg.gene[i,"cor"] >= c(-0.4) && cor.cpg.gene[i,"cor"] <= c(-0.2) ){
        cor.cpg.gene[i,"correlation_call"]<-c("WNC")
      } else {
        cor.cpg.gene[i,"correlation_call"]<-c("NNC")
      }
      
      
    }
  }
}

cor.cpg.gene[,"correlation_call"] <- as.factor(cor.cpg.gene[,"correlation_call"])

tum.per <- apply(cpg.gene[rownames(cor.cpg.gene),8:643], 1, quantile, c(.1,.5,.9), na.rm=T)
tum.per <- t(tum.per)
dimnames(tum.per)[[2]] <- c("T10","T50", "T90")
#cpg.gene.normal 8:117 -> DNA methylation 
#                118:227 -> expression
nor.per <- apply(cpg.gene.normal[rownames(cor.cpg.gene),8:117], 1, quantile, c(.1,.5,.9), na.rm=T)
nor.per <- t(nor.per)
dimnames(nor.per)[[2]] <- c("N10","N50", "N90")


genenames <- cbind(cor.cpg.gene, tum.per, nor.per)


### now, make conditions in the following manner

genenames$N.METH <- c("UNK")
genenames$T.METH <- c("UNK")
## for normals:


for(n in 1:dim(genenames)[1]){
  if(is.na(genenames[n, "N10"] < 0.25 && genenames[n, "N50"] < 0.25 && genenames[n, "N90"] < 0.25)){
    genenames[n, "N.METH"] <- c("UNK")
  } else{
    if(genenames[n, "N90"] < 0.25){
      genenames[n, "N.METH"] <- c("CUN")
    } else {
      if(genenames[n, "N10"] > 0.75){
        genenames[n, "N.METH"] <- c("CMN")
      } else {
        if(genenames[n, "N10"] > 0.25 && genenames[n, "N90"] < 0.75){
          genenames[n, "N.METH"] <- c("IMN")
        } else {
          genenames[n, "N.METH"] <- c("VMN")
          
        }
        
        
      }	
      
      
    }
  }
}#endfor

genenames[,"N.METH"] <- as.factor(genenames[,"N.METH"])



## for tumors
for(n in 1:dim(genenames)[1]){
  if(is.na(genenames[n, "T10"] < 0.25 && genenames[n, "T50"] < 0.25 && genenames[n, "T90"] < 0.25)){
    genenames[n, "T.METH"] <- c("UNK")
  } else{
    if(genenames[n, "T90"] < 0.25){
      genenames[n, "T.METH"] <- c("CUT")	
    } else {
      if(genenames[n, "T10"] > 0.75){
        genenames[n, "T.METH"] <- c("CMT")
      } else {
        if(genenames[n, "T10"] > 0.25 && genenames[n, "T90"] < 0.75){
          genenames[n, "T.METH"] <- c("IMT")
        } else {
          genenames[n, "T.METH"] <- c("VMT")
          
        }	
        
      }
    }
  }
}#endfor

genenames[,"T.METH"] <- as.factor(genenames[,"T.METH"])


rules <- read.delim(file="Rules_Patient_Centric_Calls.txt",sep="\t",header=T)
rules[,"Label"] <- paste(rules[,"Correlation"], rules[,"Normal.Meth"], rules[,"Tumor.Meth"], sep=".")

genenames$probeID <- rownames(genenames)
genenames$Label <- paste(genenames[,"correlation_call"], genenames[,"N.METH"], genenames[,"T.METH"], sep=".")
cpg.gene.m <- cbind(cpg.gene[rownames(genenames),],genenames[,"Label"])
rules.all <- merge(genenames[,c("probeID","Label")], rules[,c(1,5:10)],by.x="Label", by.y="Label",all.x=T)
dimnames(rules.all)[[1]] <- rules.all[,"probeID"]
rules.all <- rules.all[,-2]
rules.all <- rules.all[dimnames(genenames)[[1]],]
patient.call <- cpg.gene[rownames(genenames),8:643]
patient.score <- cpg.gene[rownames(genenames),8:643]

tmp.lt.25 <- patient.call < 0.25
tmp.gt.75 <- patient.call > 0.75
tmp.ge.25 <- patient.call >= 0.25 
tmp.le.75 <- patient.call <= 0.75
tmp <- tmp.ge.25 + tmp.le.75
tmp.ge.25.le.75 <- tmp == 2

for(samp in 1:dim(patient.call)[2]){
  
  ## convert NAs to FALSE
  tmp.lt.25[is.na(tmp.lt.25[,samp]),samp] <- FALSE
  tmp.gt.75[is.na(tmp.gt.75[,samp]),samp] <- FALSE
  tmp.ge.25.le.75[is.na(tmp.ge.25.le.75[,samp]),samp] <- FALSE
  
  
  patient.call[tmp.lt.25[,samp],samp] <- as.character(rules.all[tmp.lt.25[,samp],"Call.TM.LT.0.25"])
  patient.call[tmp.gt.75[,samp],samp] <- as.character(rules.all[tmp.gt.75[,samp],"Call.TM.GT.0.75"])
  patient.call[tmp.ge.25.le.75[,samp],samp] <- as.character(rules.all[tmp.ge.25.le.75[,samp],"Call.0.25.GE.TM.LE.0.75"])
  
  patient.score[tmp.lt.25[,samp],samp] <- as.character(rules.all[tmp.lt.25[,samp],"Confidence.Score.TM.LT.0.25"])
  patient.score[tmp.gt.75[,samp],samp] <- as.character(rules.all[tmp.gt.75[,samp],"Confidence.Score.TM.GT.0.75"])
  patient.score[tmp.ge.25.le.75[,samp],samp] <- as.character(rules.all[tmp.ge.25.le.75[,samp],"Confidence.Score.0.25.GE.TM.LE.0.75"])
  
  patient.call[,samp] <- as.factor(patient.call[,samp])
  patient.score[,samp] <- as.factor(patient.score[,samp])
  cat(date(), " :: Finished Sample : ",samp,"out of ", dim(patient.call)[2],"\n")
  
}

################ figure

colnames(patient.score)[1:636] <- substr(colnames(patient.score)[1:636],1,12)
colnames(patient.call)[1:636] <- substr(colnames(patient.call)[1:636],1,12)
patient.call$geneName <- rownames(patient.call)
patient.score$geneName <- rownames(patient.score)
patient <- rbind(patient.call, patient.score)

patient.t <- melt(patient.call[,c(1:637)],id.vars=c("geneName"))
patient.t2 <- melt(patient.score[,c(1:637)], id.vars = c("geneName"))
patient.t3 <- merge(patient.t, patient.t2, by=c("variable","geneName"))
LGG.GBM.new.noXY.s <- LGG.GBM.250914[,5:936]
LGG.GBM.new.order <- LGG.GBM.new.noXY.s[,cc.out.purity.lgg.gbm.new[[6]][[2]][[3]]]
colnames(LGG.GBM.new.order) <- substr(colnames(LGG.GBM.new.order),1,12)
LGG.GBM.new.order <- LGG.GBM.new.order[,c(1:63,64:327,405:531,532:682,683:932,328:404)]
a <- colnames(LGG.GBM.new.order) %in% colnames(patient.call)
teste <- LGG.GBM.new.order[,a]
patient.t3$variable <- factor(patient.t3$variable, levels = colnames(teste)) #Order the samples (like the DNA meth heatmap)

teste <- subset(patient.t3, value.x %in% c("ES") & value.y %in% c("4"))
teste <- merge(patient.t3, metadata.1129.samples.20150312[,c("id","cluster.meth")],by.x="variable",by.y="id",all.x=T)

p <- ggplot(patient.t3, aes(x=cluster.meth)) +
geom_bar(stat = "identity")
ggsave(p, filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/p.png", width=30, height=20, dpi=200, units="in")

p <- ggplot(teste[1:1000,], aes(x=value.x,fill=score)) +
#  geom_tile(aes(fill=factor(score) ) ) +
  geom_bar(stat = "bin") +
  #scale_fill_manual(values=c("blue","orange","red","purple","green","white")) +
  scale_y_log10() +
  facet_grid(. ~ cluster.meth)


g

ggsave(p, filename="/dados/ResearchProjects/thais/TCGA/LGG.GBM/teste.png", width=30, height=20, dpi=200, units="in")


ggplot(teste2, aes(y=factor(variable),x=id)) + 
  geom_tile(aes(fill=factor(value), height=as.numeric(as.character(size)))) +
  scale_fill_manual(values=c("cadetblue1","brown4","cadetblue4","red","purple","cyan","magenta","orange","green3","green3","red","black","green","red","purple","orange","yellow","blue","yellow","gray","cyan","magenta","blue","green3","red","brown1","black","cornsilk2","gray","white"), guide = guide_legend(title = "Legend")) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white"),axis.text.x = element_blank(),axis.text.y = element_text(size=8),panel.margin = unit(0, "lines")) + #, strip.background=element_rect(fill=c("green","red","purple","orange","yellow","blue")),strip.text=element_text(size=12,color=c("white"))) + 
  facet_grid(SomaticType ~ .,space="free",scales="free") +
  xlab("TCGA LGG and GBM") + 
  ylab("Gene Name") +
  guides(col = guide_legend(nrow = 8))
dev.off()

patient.m2 <- merge(patient.m,metadata.1129.samples.20150112[,1:13],by.x="variable",by.y="id",all.x=T)

p <- ggplot(meanb.naomit, aes(x=reorder(gene.name,as.numeric(value)),fill=factor(value))) + 
  geom_bar(position="fill")+ 
  #geom_bar() +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values=c("black","gray"),name = "Legend", labels = c("Mutant", "WT")) +
  #geom_bar(position = "dodge")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(y = "Percentage") + coord_flip() +
  labs(x = "Gene Name") + 
  #facet_grid(. ~ newSturmClass, scales = "free")
  facet_grid(. ~ cluster.meth2, scales = "free")


write.table(patient.score[,1:636],file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/P_Score.txt",quote=F,row.names=T,sep="\t") #de 1 a 11
write.table(patient.call[,1:636],file="/dados/ResearchProjects/thais/TCGA/LGG.GBM/P_Call.txt",quote=F,row.names=T,sep="\t") #de 1 a 11
