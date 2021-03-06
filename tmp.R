load("/home/thais/Cloud_Sandiego")
TCGA.06.0128 <- TCGA.06.0128[-1,]
TCGA.16.1460 <- TCGA.16.1460[-1,]
TCGA.19.1788 <- TCGA.19.1788[-1,]

colnames(TCGA.06.0128) <- c("chrom","start","end","name","score","strand","percentMeth","numCTreads")
colnames(TCGA.16.1460) <- c("chrom","start","end","name","score","strand","percentMeth","numCTreads")
colnames(TCGA.19.1788) <- c("chrom","start","end","name","score","strand","percentMeth","numCTreads")
TCGA.06.0128$numCTreads <- as.numeric(as.character(TCGA.06.0128$numCTreads))
TCGA.16.1460$numCTreads <- as.numeric(as.character(TCGA.16.1460$numCTreads))
TCGA.19.1788$numCTreads <- as.numeric(as.character(TCGA.19.1788$numCTreads))

TCGA.06.0128$percentMeth <- as.numeric(as.character(TCGA.06.0128$percentMeth))
TCGA.16.1460$percentMeth <- as.numeric(as.character(TCGA.16.1460$percentMeth))
TCGA.19.1788$percentMeth <- as.numeric(as.character(TCGA.19.1788$percentMeth))

TCGA.06.0128 <- subset(TCGA.06.0128, numCTreads > 4)
TCGA.16.1460 <- subset(TCGA.16.1460, numCTreads > 4)
TCGA.19.1788 <- subset(TCGA.19.1788, numCTreads > 4)

rownames(TCGA.06.0128) <- paste(TCGA.06.0128$chrom,TCGA.06.0128$start,sep=".")
rownames(TCGA.16.1460) <- paste(TCGA.16.1460$chrom,TCGA.16.1460$start,sep=".")
rownames(TCGA.19.1788) <- paste(TCGA.19.1788$chrom,TCGA.19.1788$start,sep=".")

tudao <- subset(tudao, chrom %in% paste0("chr",c(1:22,"X","Y")))
tudao <- droplevels(tudao)
tudao$start <- as.numeric(as.character(tudao$start))

win.5CpG.0128 <- list()
for(j in c("chr4","chrX")){
  aux1 <- c()
  aux <- subset(tudao,chrom %in% j)
  aux <- aux[with(aux,order(start)),]
  if(nrow(aux) > 5){
    for(i in 1:(nrow(aux)-4)){
      aux1 <- c(aux1,mean(aux[c(i:(i+4)),c("percentMeth.0128")],na.rm = T))
    }
  }
  else{
    aux1 <- c(aux1,mean(aux[,c("percentMeth.0128")],na.rm = T))
  }
  win.5CpG.0128 <- append(win.5CpG.0128,list(aux1))
  print(j)
}


win.5CpG.high <- list()
for(j in c("chr4","chrX")){
  aux1 <- c()
  aux <- subset(tudao,chrom %in% j)
  aux <- aux[with(aux,order(start)),]
  if(nrow(aux) > 5){
    for(i in 1:(nrow(aux)-4)){
      aux1 <- c(aux1,mean(aux[c(i:(i+4)),c("mean_GCIMP.high")],na.rm = T))
    }
  }
  else{
    aux1 <- c(aux1,mean(aux[,c("mean_GCIMP.high")],na.rm = T))
  }
  win.5CpG.high <- append(win.5CpG.high,list(aux1))
  print(j)
}


ids <- c(as.character(pd[pd$clustM.supervised2 %in% "G-CIMP-low","case.id"]),as.character(pd[pd$clustM.supervised2 %in% "G-CIMP-high","case.id"]))


all$percentMeth.0128 <- as.numeric(as.character(all$percentMeth.0128))
all$percentMeth.1460 <- as.numeric(as.character(all$percentMeth.1460))
all$percentMeth.1788 <- as.numeric(as.character(all$percentMeth.1788))

ids <- read.table("/home/thais/Downloads/OncoPrintPatients_all.txt",header=T)
ids$group <- c(rep("G1",26),rep("G2",248))
aux <- merge(ids,pd[,c("case.id","tumor.type","os.months","vital","clustM.supervised2")],by.x="Order",by.y="case.id",all.x=T,sort=F)

f.m = formula(Surv(os.months,vital) ~ group)
fit.lgm = survfit(f.m, data=aux)

plot(fit.lgm, 
     lwd=4,
     col=c("blue",
       "purple" #colocar a cor de acordo com a ordem de 1 a 6

     ),
     main="Overall Survival\nKaplan-Meier Survival Curves\nTCGA GLIOMA SAMPLES ", 
     xlab="TIME SINCE DIAGNOSIS (MONTHS)", 
     ylab="PERCENT PROBABILITY OF SURVIVAL", 
     yscale=100,
    # xscale=365.25, 
     bg="black"
)
box(col="black", lwd=3);
abline(h=0, lwd=4, lty=1)
legend("topright", 
       legend=c(
         paste("Amp. group (n=",fit.lgm$n[1],")",sep=""),
         paste("Intact group (n=",fit.lgm$n[2],")",sep="")

       ),
       col=c("blue",
             "purple", #colocar a cor de acordo com a ordem de 1 a 6
             "purple"

       ),
       lwd=3,
       title="Legend",
       box.lwd=3,
       bg="white")



aux$clustM.supervised2 <- factor(aux$clustM.supervised2,levels=c("G-CIMP-low","G-CIMP-high"))
f.m = formula(Surv(os.months,vital) ~ clustM.supervised2)
fit.lgm = survfit(f.m, data=aux)

plot(fit.lgm, 
     lwd=4,
     col=c("green3",
           "firebrick" #colocar a cor de acordo com a ordem de 1 a 6
           
     ),
     main="Overall Survival\nKaplan-Meier Survival Curves\nTCGA GLIOMA SAMPLES ", 
     xlab="TIME SINCE DIAGNOSIS (MONTHS)", 
     ylab="PERCENT PROBABILITY OF SURVIVAL", 
     yscale=100,
     # xscale=365.25, 
     bg="black"
)
box(col="black", lwd=3);
abline(h=0, lwd=4, lty=1)
legend("topright", 
       legend=c(
           paste("G-CIMP-low (n=",fit.lgm$n[1],")",sep=""),
           paste("G-CIMP-high (n=",fit.lgm$n[2],")",sep="")
           
       ),
       col=c("green3",
             "firebrick" #colocar a cor de acordo com a ordem de 1 a 6
             
       ),
       lwd=3,
       title="Legend",
       box.lwd=3,
       bg="white")
