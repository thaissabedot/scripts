## Heatpairs
scores$Study <- as.character(scores$Study)
scores[scores$Study == "LGm6" & scores$Tumor.Type == "LGG","Study"] <- "LGm6.LGG"
scores$id <- rownames(scores)
K1 <- scores[scores$Study == "LGm1","id"]
K2 <- scores[scores$Study == "LGm2","id"]
K3 <- scores[scores$Study == "LGm3","id"]
K4 <- scores[scores$Study == "LGm4","id"]
K5 <- scores[scores$Study == "LGm5","id"]
K6 <- scores[scores$Study == "LGm6","id"]
K7 <- scores[scores$Study == "Lambert et al., 2013","id"]
K8 <- scores[scores$Study == "LGm6.LGG","id"]
K9 <- scores[scores$Study == "non-tumor","id"]
K10 <- scores[scores$Study == "Sturm et al., 2012","id"]
K11 <- scores[scores$Study == "Turcan et al., 2012","id"]
K12 <- scores[scores$Study == "Mur et al., 2013","id"]
M1 <- all[,as.character(K1)]
M2 <- all[,as.character(K2)]
M3 <- all[,as.character(K3)]
M4 <- all[,as.character(K4)]
M5 <- all[,as.character(K5)]
M6 <- all[,as.character(K6)]
M7 <- all[,as.character(K7)]
M8 <- all[,as.character(K8)]
M9 <- all[,as.character(K9)]
M10 <- all[,as.character(K10)]
M11 <- all[,as.character(K11)]
M12 <- all[,as.character(K12)]
M1$meanM1 <- apply(M1,1,mean,na.rm=T)
M2$meanM2 <- apply(M2,1,mean,na.rm=T)
M3$meanM3 <- apply(M3,1,mean,na.rm=T)
M4$meanM4 <- apply(M4,1,mean,na.rm=T)
M5$meanM5 <- apply(M5,1,mean,na.rm=T)
M6$meanM6 <- apply(M6,1,mean,na.rm=T)
M7$meanM6 <- apply(M7,1,mean,na.rm=T)
M8$meanM6 <- apply(M8,1,mean,na.rm=T)
M9$meanM6 <- apply(M9,1,mean,na.rm=T)
M10$meanM6 <- apply(M10,1,mean,na.rm=T)
M11$meanM6 <- apply(M11,1,mean,na.rm=T)
M12$meanM6 <- apply(M12,1,mean,na.rm=T)
avg.probes.cl <- cbind(M1$meanM1,M2$meanM2,M3$meanM3,M4$meanM4,M5$meanM5,M6$meanM6,M7$meanM6,M8$meanM6,M9$meanM6,M10$meanM6,M11$meanM6,M12$meanM6)
rownames(avg.probes.cl) <- rownames(M1)
colnames(avg.probes.cl) <- c("LGm1","LGm2","LGm3","LGm4","LGm5","LGm6","Lambert et al., 2013","LGm6.LGG","non-tumor","Sturm et al., 2012","Turcan et al., 2012","Mur et al., 2013")
avg.probes.cl <- avg.probes.cl[,c(1:6,8,9,7,10:12)]
png(filename = "/dados/ResearchProjects/thais/TCGA/LGG.GBM/Heatmaps_932/heatpairs_all.png", bg="white", res=300, width=2000, height=2000)
heatpairs(avg.probes.cl)
dev.off()




