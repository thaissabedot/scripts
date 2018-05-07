m <- subset(metadata.1129.samples.20141030, Study == "TCGA")
m <- m[,c(1,53:92)] #colunas de ID do paciente e de CNV
labels <- subset(b, Study == "TCGA")
labels <- labels[,c("id","CoCluster","cluster.meth","cluster.expr","cluster.cnv")] #Colunas que voce quer colocar nas tracks do heatmap

m.s.na <- m
m.s.na <- as.data.frame(lapply(m.s.na,function(z){z <- as.character(z); z[is.na(z)] <- "NA"; return(z)})) #Transforma NA em "NA"
hc.m.s <- as.hclust(agnes(daisy(data.frame(t(m.s.na[,-1])), metric="gower"), method="ward")) #Fazer a clusterização
m.s <- m.s[,c(1,hc.m.s$order+1)] #Ordena as colunas com o resultado da clusterização
colnames(labels) <- c("id","LGcc","LGm","LGr","LGc") #Nome que voce quer que apareça ao lado da track no heatmap
labels <- labels[,c(1,2,5,4,3)]
all <- cbind(labels,m.s[,-1]) #label + dados

teste <- melt(all,id="id")
teste1 <- merge(na.omit(teste),metadata.1129.samples.20141030,by.x="id",by.y="id")
teste1$id <- factor(teste1$id)
ID_order <- as.character(substr(colnames(LGG.GBM.new.order),1,12)) #Aqui você coloca a ordem das amostras que voce quer no heatmap
teste1$id <- factor(teste1$id, levels = ID_order)
png(filename = "Mutation_Heatmap.png", bg="white", res=200, width=4000, height=2000)
ggplot(teste1, aes(y=factor(variable),x=id)) + 
  geom_tile(aes(fill=factor(value))) +
  scale_fill_manual(values=c("cadetblue1","brown4","cadetblue4","red","purple","cyan","magenta","orange","green3","green3","red","black","green","red","purple","orange","yellow","blue","yellow","gray","cyan","magenta","blue","green3","red","brown1","black","cornsilk2","gray","white"),guide = guide_legend(title = "Legend")) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white"),axis.text.x = element_blank(),axis.text.y = element_text(size=8),panel.margin = unit(0, "lines")) + #, strip.background=element_rect(fill=c("green","red","purple","orange","yellow","blue")),strip.text=element_text(size=12,color=c("white"))) + 
  facet_grid(SomaticType ~ .,space="free",scales="free") +
  xlab("TCGA LGG and GBM") + 
  ylab("Gene Name") +
  guides(col = guide_legend(nrow = 8))
dev.off()

