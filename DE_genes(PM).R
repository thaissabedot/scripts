norm.edger <- DGEList(counts= "counts", genes= "seus genes")

egSYMBOL <- toTable(org.Hs.egSYMBOL)
idfound <- norm.edger$genes$genes %in% egSYMBOL$gene_id
norm.edger.clean <- norm.edger[idfound,]

norm.edger.n <- calcNormFactors(norm.edger.clean, method = "TMM")
celltype <- factor(targets.new$CellType)
levels(celltype)[2:3] <- c("NFTE", "NOSE")
design.new <- model.matrix( ~ celltype)
v.norm <- voom(norm.edger.n, design.new, plot = F, normalize.method= "quantile")
v.norm.stdev <- sort(apply(v.norm$E, 1 , sd, na.rm=T), decreasing = T) ### a partir daqui com seu  "v.norm$E ja da pra plotar heatmaps"

# Comparison
contrast.matrix <- makeContrasts(celltypeNOSE-celltypeNFTE, levels = design.new[,2:3])
fit.norm <- lmFit(v.norm, design.new[,2:3])
fit.norm$genes <- data.frame(fit.norm$genes, select(x=org.Hs.eg.db,
                                             keys=as.character(fit.norm$genes$genes),
                                             columns=c("SYMBOL"),
                                             keytype=c("ENTREZID")))



fit.norm <- contrasts.fit(fit.norm, contrast.matrix)
fit.norm.e <- eBayes(fit.norm, robust = T)          

#p/ volcano
cutoff.norm <- topTable(fit.norm.e, number = "Inf", coef = "celltypeNOSE - celltypeNFTE") 
cutoff.norm$STATUS <- "NOT SIG"
cutoff.norm[cutoff.norm$logFC < -3 & cutoff.norm$adj.P.Val < 0.000001, ]$STATUS <- "DOWN"
cutoff.norm[cutoff.norm$logFC > 3 & cutoff.norm$adj.P.Val < 0.000001, ]$STATUS <- "UP"                                   


#p/ heatmap
nose <- v.norm$E[as.character(diff.genes), c(1:51,76,79,81:100)]
nfte <- v.norm$E[as.character(diff.genes), c(55:75,77:78,80)]
nose.o <- hclust(dist(t(nose)))
nfte.o <- hclust(dist(t(nfte)))
rnaseq.order.all <- cbind(nose[,nose.o$order], nfte[,nfte.o$order])
collabell <- as.data.frame(colnames(rnaseq.order.all));names(collabell) <- "colors"
collabell$colors <- c(rep("#377EB8",73),rep("red",24))
collabell$ID <- colnames(rnaseq.order.all)
rowlabel <- matrix(c(rep("red",399), rep("green", 902)))
