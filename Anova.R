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

