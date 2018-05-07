library(ConsensusClusterPlus)
title="/pasta"
setwd(title)
cc.out <- ConsensusClusterPlus(d=dados, 
                               maxK=10, 
                               reps=1000, 
                               pItem=0.8, 
                               pFeature=1,
                               title=title, 
                               clusterAlg="hc",   #ou: km, etc.
                               distance="binary", #ou: pearson, euclidean, spearman, etc.
                               seed=1022,
                               plot="pdf",
                               verbose = TRUE)
ccICL = calcICL(cc.out, plot='pdf')



### survival
library(survival)
f.m = formula(Surv(OS.months,Status) ~ type)
fit.m = survfit(f.m, data=metadata)

pdf(file = "survival.pdf")
plot(fit.m, 
     lwd=4,
     col=c("green",
       "firebrick4", #colocar a cor de acordo com a ordem dos grupos
           "orange3",
           "blue"
           
           
     ),
     main="Kaplan-Meier Overall Survival Curves", 
     xlab="TIME SINCE DIAGNOSIS (MONTHS)", 
     ylab="PERCENT PROBABILITY OF SURVIVAL", 
     yscale=100,
     # xscale=365.25, 
     bg="black"
)
box(col="black", lwd=3);
legend("bottomleft", 
       legend=c(
         paste("Group1 (n=",fit.m$n[1],")",sep=""),
         paste("Group2 (n=",fit.m$n[2],")",sep=""),
         paste("Group3 (n=",fit.m$n[3],")",sep=""),
         paste("Group4 (n=",fit.m$n[4],")",sep="")
         
       ),
       col=c("green",
             "firebrick4",
             "orange3",   #colocar a cor na mesma ordem que está acima
             "blue"
             
       ),
       lwd=3,
       title="Legend",
       box.lwd=3,
       bg="white")
dev.off()

