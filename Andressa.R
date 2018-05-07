#### Tests between different treatment for the same condition
setwd("/home/thais/Andressa")
files <- list.files()[grep("*CS.txt",list.files())]
a <- NULL
for(i in 1:length(files)){
	if(is.null(a)){
		aux <- read.delim(files[i])
		aux <- aux[grep("*S",aux$Sample.Name),]
		a <- unlist(strsplit(as.character(aux$Sample.Name), " "))
		aux$groups <- a[which(1:length(a) %% 2 == 0)]
		a <- aux[,c("Target.Name","groups","Delta.Delta.Ct")]
		colnames(a)[3] <- unlist(strsplit(files[i], " "))[1]
		
	}
	else{
		aux <- read.delim(files[i])
		aux <- aux[grep("*S",aux$Sample.Name),]
		b <- unlist(strsplit(as.character(aux$Sample.Name), " "))
		aux$groups <- b[which(1:length(b) %% 2 == 0)]
		a <- cbind(a,aux[,"Delta.Delta.Ct"])
		colnames(a)[ncol(a)] <- unlist(strsplit(files[i], " "))[1]
	}

}
a$Target.Name <- factor(a$Target.Name) #aqui são os nomes das probes
anova_groups <- data.frame(miRNA=levels(a$Target.Name)) #aqui vai ser o seu data frame com os p-value, ai criei uma coluna com o nome de cada probe pra você saber a ordem
anova_groups$pvalue <- NA #coluna de p-values
for(i in 1:nlevels(a$Target.Name)){ #a variável i vai de 1 até o número total de probes diferentes que você tem
	b <- subset(a, Target.Name %in% levels(a$Target.Name)[i]) #seleciona somente as instâncias dessa probe no seu objeto
	b <- melt(b) #faz um melt só para essa probe. Em uma coluna você vai ter os valores (beta, por exemplo) e na outra a qual grupo a amostra pertecen
	if(nrow(na.omit(b)) > 1){ #se você só tiver NA no seu objeto ele não realiza o teste
		test <- summary(aov(value ~ groups, data=b)) #faz o teste anova. A função summary ajuda a extrair o resultado do teste.
		if(ncol(test[[1]]) > 3) 
			anova_groups[i,2] <- test[[1]][[5]][[1]] #dimensão onde o p-value está armazenado
	}

}
a[a$Target.Name %in% anova_groups[anova_groups$pvalue < 0.05 & !is.na(anova_groups$pvalue),"miRNA"],] #ajuste do p-value


b <- melt(a)
b$Target.Name <- factor(b$Target.Name, levels=c(as.character(anova_groups[anova_groups$pvalue < 0.05 & !is.na(anova_groups$pvalue),"miRNA"]),as.character(anova_groups[anova_groups$pvalue > 0.05 & !is.na(anova_groups$pvalue),"miRNA"]),as.character(anova_groups[is.na(anova_groups$pvalue),"miRNA"])))

b$type <- paste(b$groups,b$variable,sep=".")
b <- dcast(b, Target.Name ~ type)
rownames(b) <- as.character(b$Target.Name)
b <- b[(c(as.character(anova_groups[anova_groups$pvalue < 0.05 & !is.na(anova_groups$pvalue),"miRNA"]),as.character(anova_groups[anova_groups$pvalue > 0.05 & !is.na(anova_groups$pvalue),"miRNA"]),as.character(anova_groups[is.na(anova_groups$pvalue),"miRNA"]))),]


c <- b
c[is.na(c)] <- min(c[,-1],na.rm=TRUE)+100
b <- b[hclust(dist(c))$order,]
aux <- colnames(c)[2:11]
order <- aux[hclust(dist(t(c[,2:11])))$order]
aux <- colnames(c)[12:21]
order <- c(order,aux[hclust(dist(t(c[,12:21])))$order])
aux <- colnames(c)[22:31]
order <- c(order,aux[hclust(dist(t(c[,22:31])))$order])
aux <- colnames(c)[32:41]
order <- c(order,aux[hclust(dist(t(c[,32:41])))$order])
b <- b[,order]
order <- c(rownames(b)[10:48],rownames(b)[1:9])
b <- b[order,]
ha1 = HeatmapAnnotation(df = data.frame(type1 = c(rep("CS",10),rep("RS",10),rep("RTS",10),rep("TS",10))), 
    col = list(type1 = c("CS" =  "#F8766D", "RS" = "#7CAE00","RTS"="#00BFC4","TS"="#C77CFF")))

pdf("/home/thais/Andressa/Sus_heatmap.pdf",width=10, height=10)
Heatmap(as.matrix(b), name = "Treatment", col=greenred(75),cluster_rows = FALSE, cluster_columns = FALSE,top_annotation = ha1)
dev.off()

b <- melt(a)
b$Target.Name <- factor(b$Target.Name, levels=order)
ggplot(b, aes(groups, value,fill=groups)) + 
  geom_boxplot() + 
  #geom_jitter()  + 
  facet_wrap(~ Target.Name, scales = "free") 
ggsave("/home/thais/Andressa/Sus_boxplot.pdf",width = 20, height = 11)

anova_groups$adj_pvalue <- p.adjust(anova_groups$pvalue, method = "fdr")
write.table(anova_groups,quote=F,row.names=F,sep="\t",file="/home/thais/Andressa/Sus_ANOVA.txt")

 
#### Tests between different conditions for the same treatment
setwd("/home/thais/Andressa")
files <- list.files()[grep("*CA.txt",list.files())]
c <- NULL
for(i in 1:length(files)){
	if(is.null(c)){
		aux <- read.delim(files[i])
		aux <- aux[grep("*A",aux$Sample.Name),]
		c <- unlist(strsplit(as.character(aux$Sample.Name), " "))
		aux$groups <- c[which(1:length(c) %% 2 == 0)]
		c <- aux[,c("Target.Name","groups","Delta.Delta.Ct")]
		colnames(c)[3] <- unlist(strsplit(files[i], " "))[1]
		
	}
	else{
		aux <- read.delim(files[i])
		aux <- aux[grep("*A",aux$Sample.Name),]
		b <- unlist(strsplit(as.character(aux$Sample.Name), " "))
		aux$groups <- b[which(1:length(b) %% 2 == 0)]
		c <- cbind(c,aux[,"Delta.Delta.Ct"])
		colnames(c)[ncol(c)] <- unlist(strsplit(files[i], " "))[1]
	}

}
c$Target.Name <- factor(c$Target.Name)


ttest_groups <- data.frame(miRNA=levels(a$Target.Name))
ttest_groups$pvalue <- NA
sus <- NULL
adh <- NULL
for(i in 1:nlevels(a$Target.Name)){
	aux <- a[a$Target.Name %in% levels(a$Target.Name)[i] & a$groups %in% "TS",3:12]
	aux2 <- c[c$Target.Name %in% levels(a$Target.Name)[i] & c$groups %in% "TA",3:12]
	if(is.null(sus)){
		sus <- cbind(miRNA=levels(a$Target.Name)[i],a[a$Target.Name %in% levels(a$Target.Name)[i] & a$groups %in% "TS",3:12])
		adh <- cbind(miRNA=levels(a$Target.Name)[i],c[c$Target.Name %in% levels(a$Target.Name)[i] & c$groups %in% "TA",3:12])
	}
	else{
		sus <- rbind(sus,cbind(miRNA=levels(a$Target.Name)[i],a[a$Target.Name %in% levels(a$Target.Name)[i] & a$groups %in% "TS",3:12]))
		adh <- rbind(adh,cbind(miRNA=levels(a$Target.Name)[i],c[c$Target.Name %in% levels(a$Target.Name)[i] & c$groups %in% "TA",3:12]))
	}
	if(sum(is.na(aux)) < 9 & sum(is.na(aux2)) < 9)
		ttest_groups[i,2] <- t.test(aux,aux2)$p.value
}



ttest_groups$adj_pvalue <- p.adjust(ttest_groups$pvalue, method = "fdr")
write.table(ttest_groups,quote=F,row.names=F,sep="\t",file="/home/thais/Andressa/T_Ttest.txt")


all <- cbind(adh[,-1],sus[,-1])
rownames(all) <- adh[,1]
colnames(all)[1:10] <- paste("adh",colnames(all)[1:10],sep=".")
colnames(all)[11:20] <- paste("sus",colnames(all)[11:20],sep=".")

c <- all
c[is.na(c)] <- min(c[,-1],na.rm=TRUE)+100
all <- all[hclust(dist(c))$order,]
aux <- colnames(c)[1:10]
order <- aux[hclust(dist(t(c[,1:10])))$order]
aux <- colnames(c)[11:20]
order <- c(order,aux[hclust(dist(t(c[,11:20])))$order])
all <- all[,order]
order <- c(rownames(all)[15:48],rownames(all)[1:14]) #R c(rownames(all)[1:21],rownames(all)[48:22]) #T c(rownames(all)[15:48],rownames(all)[1:14])
all <- all[order,]
ha1 = HeatmapAnnotation(df = data.frame(type1 = c(rep("Adhesion",10),rep("Suspension",10))), 
    col = list(type1 = c("Adhesion" =  "#F8766D", "Suspension" = "#00BFC4")))

pdf("/home/thais/Andressa/T_heatmap.pdf",width=10, height=10)
Heatmap(as.matrix(all), name = "Expression values", col=greenred(75),cluster_rows = FALSE, cluster_columns = FALSE,top_annotation = ha1)
dev.off()


aux <- melt(sus)
aux$Condition <- "Suspension"
aux2 <- melt(adh)
aux2$Condition <- "Adhesion"
aux <- rbind(aux,aux2)

aux$miRNA <- factor(aux$miRNA, levels=order)
ggplot(aux, aes(Condition, value,fill=Condition)) + 
  geom_boxplot() + 
  #geom_jitter()  + 
  facet_wrap(~ miRNA, scales = "free") 
ggsave("/home/thais/Andressa/T_boxplot.pdf",width = 20, height = 11)


heatmap.plus(as.matrix(aux),
	     Colv=NA,
	     Rowv=NA,
	     col=greenred(75),
	     scale = "row")




b <- stack(a[,-c(1:2)])
b$subject = rep(a$Target.Name, 10) 
colnames(b) = c("price", "store", "subject")
with(b, pairwise.t.test(x=price, g=store, p.adjust.method="none", paired=T))

aov.out = aov(price ~ store + Error(subject/store), data=b)

for(i in 1:length(files)){
	if(is.null(a)){
		aux <- read.delim(files[i])
		aux <- aux[grep("*A",aux$Sample.Name),]
		a <- unlist(strsplit(as.character(aux$Sample.Name), " "))
		aux$groups <- a[which(1:length(a) %% 2 == 0)]
		a <- aux[,c("groups","Delta.Delta.Ct")]
		colnames(a)[2] <- unlist(strsplit(files[i], " "))[1]
		
	}
	else{
		aux <- read.delim(files[i])
		aux <- aux[grep("*A",aux$Sample.Name),]
		b <- unlist(strsplit(as.character(aux$Sample.Name), " "))
		aux$groups <- b[which(1:length(b) %% 2 == 0)]
		a <- cbind(a,aux[,"Delta.Delta.Ct"])
		colnames(a)[ncol(a)] <- unlist(strsplit(files[i], " "))[1]
	}

}


#T-test instead of ANOVA
setwd("/home/thais/Andressa")
files <- list.files()[grep("*CS.txt",list.files())]
c <- NULL
for(i in 1:length(files)){
	if(is.null(c)){
		aux <- read.delim(files[i])
		aux <- aux[grep("*S",aux$Sample.Name),]
		c <- unlist(strsplit(as.character(aux$Sample.Name), " "))
		aux$groups <- c[which(1:length(c) %% 2 == 0)]
		c <- aux[,c("Target.Name","groups","Delta.Delta.Ct")]
		colnames(c)[3] <- unlist(strsplit(files[i], " "))[1]
		
	}
	else{
		aux <- read.delim(files[i])
		aux <- aux[grep("*S",aux$Sample.Name),]
		b <- unlist(strsplit(as.character(aux$Sample.Name), " "))
		aux$groups <- b[which(1:length(b) %% 2 == 0)]
		c <- cbind(c,aux[,"Delta.Delta.Ct"])
		colnames(c)[ncol(c)] <- unlist(strsplit(files[i], " "))[1]
	}

}
c$Target.Name <- factor(c$Target.Name)


ttest_groups <- data.frame(miRNA=levels(c$Target.Name))
ttest_groups$pvalue <- NA
gr1 <- NULL
gr2 <- NULL
for(i in 1:nlevels(c$Target.Name)){
	aux <- c[c$Target.Name %in% levels(c$Target.Name)[i] & c$groups %in% "RTS",3:12]
	aux2 <- c[c$Target.Name %in% levels(c$Target.Name)[i] & c$groups %in% "TS",3:12]
	if(sum(is.na(aux)) < 9 & sum(is.na(aux2)) < 9)
		ttest_groups[i,2] <- t.test(aux,aux2)$p.value
}

ttest_groups$adj_pvalue <- p.adjust(ttest_groups$pvalue, method = "fdr")
write.table(ttest_groups,quote=F,row.names=F,sep="\t",file="/home/thais/Andressa/Sus_Ttest_RT_vs_T.txt")

