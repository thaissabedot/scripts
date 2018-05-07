
setwd("/dados/normalized.geneExp/rename")
patients <- list.files(pattern="TCGA*")

gene.id <- read.delim(patients[1], header=TRUE, sep="\t")[,1]

require(plyr)
gene.express <- sapply(patients, function(x) {
  gene.expres <- read.delim(x, header=TRUE, sep="\t")[,2]
})

require(stringr)
gene.id <- str_split(as.character(gene.id), pattern="\\|")
gene.id <- do.call(rbind.data.frame, gene.id)
colnames(gene.id)  <- c("SYM", "GENEID")
gene.matrix <- cbind(gene.id, gene.express)
colnames(gene.matrix)[-c(1:2)] <- substr(colnames(gene.matrix)[-c(1:2)], 1, 12)
crap.genes.map <- read.delim("crapgenesmap.txt", header=T, sep="\t") ## inputted GeneID into pubmed and extracted this table.
x <- gene.id[1:29, ]
xx <- merge(x, crap.genes.map, by.x="GENEID", by.y="GeneID")[, c("Symbol", "GENEID")]
gene.matrix[1:29, c("SYM", "GENEID")]
gene.matrix$SYM <- as.character(gene.matrix$SYM)
gene.matrix$SYM[1:29] <- as.character(xx[, 1]); rm(x); rm(xx)

load("/dados/methylation.LGG.289.rda")
#load(file="~houtan/.synapseCache/313/488313//LGG-Annotation.3-10-2014.RData")

cpg.matrix <- dat
rm(dat)
colnames(cpg.matrix)[-c(1:4)] <- substr(colnames(cpg.matrix)[-c(1:4)], 1, 12)

rsid <- grepl("cg", cpg.matrix$Composite.Element.REF)
cpg.matrix.s <- cpg.matrix[rsid,]
cpg.matrix.sm <- na.omit(cpg.matrix.s)
cpg.matrix.sm$Gene_Symbol[cpg.matrix.sm$Gene_Symbol==""] <- NA


cpg.matrix.sm <- cpg.matrix.sm[order(cpg.matrix.sm$Gene_Symbol, cpg.matrix.sm$Composite.Element.REF), ]
gene.matrix.sm <- gene.matrix[order(gene.matrix$SYM), ]

cpg.sample <- colnames(cpg.matrix.sm[-c(1:4)])
gene.sample <- colnames(gene.matrix.sm[-c(1:2)])
both.sample <- intersect(cpg.sample, gene.sample)
both.sample <- both.sample[order(both.sample)]

gene.matrix.samples <- gene.matrix.sm[, both.sample]
cpg.matrix.samples <- cpg.matrix.sm[, both.sample]

gene.matrix.sm <- cbind(gene.matrix.sm[, c(1:2)], gene.matrix.samples)
cpg.matrix.sm <- cbind(cpg.matrix.sm[, c(1:4)], cpg.matrix.samples)

load("cpg.gene.loc.Rda") ## Adds cpg.gene.loc.df
gene.table <- read.delim('~scoetzee/Galaxy1-[Homo_sapiens_genes_(GRCh37.p13)].tabular', header=T, sep="\t")
gene.table <- gene.table[, c("EntrezGene.ID", "Associated.Gene.Name")]
dups <- duplicated(gene.table$EntrezGene.ID)
uniques <- gene.table[!dups, ]
cpg.gene.loc.df <- merge(cpg.gene.loc.df, uniques, by.x="feature", by.y="EntrezGene.ID", all.x=T)
cpg.gene.loc.df$names <- cpg.gene.loc.df$Associated.Gene.Name
cpg.gene.loc.df$Associated.Gene.Name <- NULL
cpg.matrix.sm.b <- merge(cpg.gene.loc.df, cpg.matrix.sm, by.x = "peak", by.y = "Composite.Element.REF")
cpg.matrix.sm.b <- cpg.matrix.sm.b[, -c(3:6,8:11,13:17)]
x <- colnames(cpg.matrix.sm.b)
x[2] <- "GeneID"
colnames(cpg.matrix.sm.b) <- x
cpg.matrix.sm <- cpg.matrix.sm.b
rm(x); rm(gene.matrix.samples); rm(cpg.matrix.samples); rm(cpg.matrix.sm.b)

save(list=c("cpg.matrix.sm", "gene.matrix.sm"), file="cpg.gene.matricies.Rda")

#### CUSTOM BIT FOR A SINGLE GENE
# exon.matrix.CDKN2A.1a <- read.delim(file="../exon/rename/CDKN2A//CDKN2A.1a", header=F,sep="\t")
# colnames(exon.matrix.CDKN2A.1a) <- c("SampleID", "Position", "RawCount", "LengthNorm", "RPKM")
# cpgID <- "cg13601799"
# geneSymbol <- "CDKN2A"
# cpg.for.gene <- cpg.matrix[cpg.matrix$Gene_Symbol == geneSymbol, -c(2:4)]
# cpg.for.gene$Composite.Element.REF <- as.character(cpg.for.gene$Composite.Element.REF)
# cpg.for.gene <- cpg.for.gene[-1]
# expression.for.exon <- t(exon.matrix.CDKN2A.1a[, "RPKM"])
# names(expression.for.exon) <- substr(exon.matrix.CDKN2A.1a[, "SampleID"],1,12)
# expression.for.exon <- expression.for.exon[intersect(names(cpg.for.gene), names(expression.for.exon))]
# expression.for.exon <- as.matrix(expression.for.exon)
# z <- cbind(t(cpg.for.gene), expression.for.exon)
# colnames(z) <- c("meth", "expression")
# x <- read.delim("CDKN2A_GISTIC_ThresholdedByGene_withSubtypes.txt", header=T, sep="\t")
# zz <- merge(z, x, by.x=0, by.y = "Tumor", all.x=T)
# zz <- as.data.frame(zz)
# zz$methcut <- NA
# zz[zz$meth > 0.3, "methcut"] <- TRUE
# zz[zz$meth <= 0.3, "methcut"] <- FALSE
# 
# unmethmean <- mean(zz[zz[, "meth"] <= 0.3, "expression"], na.rm = T)
# ggplot(zz, aes(x=meth, y=expression)) +
#   geom_point(size=3, aes(colour = factor(CDKN2A_GISTIC), shape = factor(Subtype))) + 
#   geom_hline(aes(yintercept=unmethmean), linetype="dashed") + 
#   geom_vline(aes(xintercept=0.3), linetype="dashed") + scale_x_continuous(limits=c(0,1), breaks=c(0, 0.3, 0.5, 0.7, 1)) +
# #  scale_y_log10() +
#   ylab("RPKM") +
#   xlab("DNA Methylation Beta Values") + 
#   labs(title = paste("DNA Methylation vs Expression\n for cg13601799|CDKN2A (exon 1a)"), colour = "CDKN2A GISTIC", shape = "IDH + 1p19qStatus")

# load(file="~houtan/.synapseCache/313/488313//LGG-Annotation.3-10-2014.RData") #updated march 10
# ## pull in the 7 miRNA cluster
# #syn2362918 <- synGet(id='syn2362918', load=T)
# miR7 <- read.delim("~houtan/.synapseCache/713/484713//LGG.miRNAseq.NMF-7groups.cluster-calls.BCGSC.20140115.txt")
# miR7$new <- substr(miR7$ID, 1, 12)
# names(miR7)[4] <- "miRNACluster7"
# Clustering <- merge(Clustering, miR7[,c(7,4)], by.x="Tumor", by.y = "new")
# 
# ##recurrent status
# recurrent.cases <- data.frame(
#   c("TCGA-DH-A669",  ## not in data freeze
#     "TCGA-TM-A7CF",  ## not in data freeze
#     "TCGA-FG-5963", 
#     "TCGA-FG-5965", 
#     "TCGA-FG-A4MT"), 
#   c(rep("Recurrent", times=5))
# )
# colnames(recurrent.cases) <- c("Tumor", "Recurrent.status")
# Clinical <- merge(Clinical, recurrent.cases, by = "Tumor", all.x = T)
# clinical.cluster <- merge(Clustering, Clinical, by = "Tumor")
# clinical.cluster.genetics <- merge(clinical.cluster, Genetics, by = "Tumor")
agecat <- cut(clinical.cluster.genetics$age_at_initial_pathologic, c(13,31.25,40,54,75))
clinical.cluster.genetics$age.strat <- as.factor(agecat)

save(list=c("cpg.matrix.sm", "gene.matrix.sm", "clinical.cluster.genetics"), file="cpg.gene.matricies.Rda")



load("cpg.gene.matricies.Rda")

InterrogateGene <- function(geneID = NULL, geneSymbol = NULL, gene.matrix = gene.matrix.sm, cpg.matrix = cpg.matrix.sm,
                            only.pass=FALSE, enrich.cluster="K4") {
  if(is.null(geneID)) {
    cpg.for.gene <- cpg.matrix[cpg.matrix$names == geneSymbol, -c(2:4)]
    cpg.for.gene$peak <- as.character(cpg.for.gene$peak)
    expression.for.gene <- gene.matrix.sm[gene.matrix.sm$SYM == geneSymbol, -c(1:2)]
    row.names(expression.for.gene) <- geneSymbol
    if(only.pass==FALSE) {
      x <- lapply(split(cpg.for.gene, cpg.for.gene$peak), function(m, e = expression.for.gene) {
        m <- as.numeric(m[-1])
        e <- as.numeric(e)
        cor.pvalue <- cor.test(e, m, method="spearman",alternative="two.sided")
        return(c(rho=cor.pvalue$estimate, p.value=cor.pvalue$p.value))
      })
    } else {
      x <- lapply(split(cpg.for.gene, cpg.for.gene$peak), function(m, e = expression.for.gene) {
        m <- as.numeric(m[-1])
        e <- as.numeric(e)
        z <- cbind(m, e)
        unmeth <- z[,1] <= 0.3
        meth <- z[,1] > 0.3
        if(sum(meth) > 1 & sum(unmeth) > 1) {
          if(mean(z[meth, "e"], na.rm=TRUE) < 1.28*sd(z[unmeth, "e"], na.rm=TRUE)) {
            unmeth.mean.expression <- mean(z[unmeth, "e"], na.omit=TRUE)
            meth.samples.count <- sum(meth)
            meth.silenced.count <- sum(meth & z[, "e"] < unmeth.mean.expression)
            if(meth.silenced.count > 0.8 * meth.samples.count) {
              cor.pvalue <- cor.test(e, m, method="spearman",alternative="two.sided")
            } else {
              cor.pvalue <- NULL
            }
          } else {
            cor.pvalue <- NULL
          } 
        } else {
          cor.pvalue <- NULL
        }
        return(c(rho=cor.pvalue$estimate, p.value=cor.pvalue$p.value))
      })
    }
    x.df <- do.call(rbind.data.frame, x)
    if(length(x.df) >= 2) {
      colnames(x.df) <- c("rho", "p.value")
      return(x.df[order(x.df$rho), ])
    } else {
      return(x.df)
    }
  } else {
    cpg.for.gene <- cpg.matrix[cpg.matrix$GeneID == geneID, -c(2:4)]
    cpg.for.gene$peak <- as.character(cpg.for.gene$peak)
    expression.for.gene <- gene.matrix.sm[gene.matrix.sm$GENEID == geneID, -c(1:2)]
    row.names(expression.for.gene) <- geneID
    if(only.pass==FALSE) {
      x <- lapply(split(cpg.for.gene, cpg.for.gene$peak), function(m, e = expression.for.gene) {
        m <- as.numeric(m[-1])
        e <- as.numeric(e)
        cor.pvalue <- cor.test(e, m, method="spearman",alternative="two.sided")
        return(c(rho=cor.pvalue$estimate, p.value=cor.pvalue$p.value))
      })
    } else {
      x <- lapply(split(cpg.for.gene, cpg.for.gene$peak), function(m, e = expression.for.gene) {
        m <- as.numeric(m[-1])
        e <- as.numeric(e)
        z <- cbind(m, e)
        unmeth <- z[,1] <= 0.3
        meth <- z[,1] > 0.3
        if(sum(meth) > 1 & sum(unmeth) > 1) {
          if(mean(z[meth, "e"], na.rm=TRUE) < 1.28*sd(z[unmeth, "e"], na.rm=TRUE)) {
            unmeth.mean.expression <- mean(z[unmeth, "e"], na.omit=TRUE)
            meth.samples.count <- sum(meth)
            meth.silenced.count <- sum(meth & z[, "e"] < unmeth.mean.expression)
            if(meth.silenced.count > 0.8 * meth.samples.count) {
              cor.pvalue <- cor.test(e, m, method="spearman",alternative="two.sided")
            } else {
              cor.pvalue <- NULL
            }
          } else {
            cor.pvalue <- NULL
          } 
        } else {
          cor.pvalue <- NULL
        }
        return(c(rho=cor.pvalue$estimate, p.value=cor.pvalue$p.value))
      })
    }
    x.df <- do.call(rbind.data.frame, x)
    if(length(x.df) >= 2) {
      colnames(x.df) <- c("rho", "p.value")
      return(x.df[order(x.df$rho), ])
    } else {
      return(x.df)
    }
  }
}

PlotGeneCpG <- function(geneSymbol = NULL, cpgID = NULL, gene.matrix = gene.matrix.sm, cpg.matrix = cpg.matrix.sm, 
                        annotation = clinical.cluster.genetics, annotate.color="MethylationCluster", annotate.shape="IDH.1p19q.Subtype.x") {
  require(ggplot2)
  cpg.for.gene <- subset(cpg.matrix, peak == cpgID & names == geneSymbol)[, -c(1:4)]
  expression.for.gene <- subset(gene.matrix, SYM == geneSymbol)[, -c(1:2)]
  expression.for.gene <- (expression.for.gene - min(gene.matrix[, -c(1,2)], na.rm=T))/(max(gene.matrix[, -c(1,2)], na.rm=T)- min(gene.matrix[, -c(1,2)], na.rm=T))
  x <- cbind(t(cpg.for.gene), t(expression.for.gene))
  colnames(x) <- c("meth", "expression")
  x <- as.data.frame(x)
  x <- merge(x, annotation, by.x=0, by.y="Tumor", all.x=TRUE)
  message("calculating mean")
  x.df <<- x
  unmethmean <<- mean(x[x[, "meth"] <= 0.3, "expression"], na.rm=T)
  message(paste0("mean is ", unmethmean, " calculated ", mean(x[x[, "meth"] <= 0.3, "expression"], na.rm=T)))
  x[, annotate.shape] <- factor(x[, annotate.shape])
  x[, annotate.color] <- factor(x[, annotate.color])
  ggplot(x, aes(x=meth, y=expression)) +
    geom_point(size=3, aes_string(colour = annotate.color, shape = annotate.shape)) +
    #scale_colour_manual(values = c("K1" = "black","K2" = "red","K3" = "blue", "K4" = "darkgreen", "K5" = "darkorchid")) + 
    scale_y_log10() +
    ylab("Normalized RNA seq expression (log10)") +
    xlab("DNA Methylation Values") + 
    geom_hline(aes(yintercept=unmethmean), linetype="dashed") + 
    geom_vline(aes(xintercept=0.3), linetype="dashed") + scale_x_continuous(limits=c(0,1), breaks=c(0, 0.3, 0.5, 0.7, 1)) +
    labs(title = paste0("DNA Methylation vs Expression\n for ", cpgID, "|", geneSymbol))
}

GroupAnnotation <- function(geneSymbol = NULL, cpgID = NULL, gene.matrix = gene.matrix.sm, cpg.matrix = cpg.matrix.sm, 
         annotation = clinical.cluster.genetics, annotate.color="MethylationCluster", annotate.shape="IDH.1p19q.Subtype.x") {
cpg.for.gene <- cpg.matrix[cpg.matrix$names == geneSymbol & cpg.matrix$peak == cpgID, -c(1:4)]
expression.for.gene <- gene.matrix.sm[gene.matrix.sm$SYM == geneSymbol, -c(1:2)]
x <- cbind(t(cpg.for.gene), t(expression.for.gene))
colnames(x) <- c("meth", "expression")
x <- as.data.frame(x)
x <- merge(x, annotation, by.x=0, by.y="Tumor", all.x=TRUE)
x.df <<- x
x.df$meanstatus <- NA
x.df[x.df$MethylationCluster=="K1", "meanstatus"] <- log2(x.df[x.df$MethylationCluster == "K1", "expression"]) < log2(mean(x.df[x.df$MethylationCluster != "K1", "expression"]))
x.df.tf <- is.na(x.df$meanstatus)
x.df <- x.df[!x.df.tf, ]
#x.df$meanstatus <- x.df$meanstatus * 1
x.df <- x.df[order(x.df$meanstatus), ]
return(x.df)
}
KHDRBS2 cg16587616
CallSilencedSample <- function(geneSymbol = NULL, cpgID = NULL, gene.matrix = gene.matrix.sm,
                               cpg.matrix = cpg.matrix.sm, 
                               annotation = clinical.cluster.genetics) {
  cpg.for.gene <- as.numeric(as.character(cpg.matrix[cpg.matrix$names == geneSymbol & cpg.matrix$peak == cpgID, -c(1:4)]))
  expression.for.gene <- gene.matrix.sm[gene.matrix.sm$SYM == geneSymbol, -c(1:2)]
  x <- as.data.frame(cbind(cpg.for.gene, t(expression.for.gene)), stringsAsFactors=F)
  colnames(x) <- c("meth", "expression")
  x$es <- FALSE
  x[(x$meth > 0.3) & (x$expression < mean(x[x$meth < 0.3, "expression"], na.rm=T)), "es"] <- TRUE
  z <- x$es
  names(z) <- row.names(x)
  return(z)
}

enrichmentList.geneNames2$uid <- paste(enrichmentList.geneNames2$cpgID, enrichmentList.geneNames2$GeneName, enrichmentList.geneNames2$cluster, sep="|")
x <- lapply(split(enrichmentList.geneNames2, enrichmentList.geneNames2$uid), function(x, gm=gene.matrix.sm, cm=cpg.matrix.sm, a=clinical.cluster.genetics) {
  SYM <<- as.character(x[["GeneName"]])
  CPG <<- as.character(x[["cpgID"]])
  CallSilencedSample(SYM, CPG, gm, cm, a)
})

x <- do.call(rbind, x)

# unmethmean <- mean(zz[zz[, "meth"] <= 0.3, "expression"], na.rm = T)
# ggplot(zz, aes(x=meth, y=expression)) +
#   geom_point(size=3, aes(colour = factor(CDKN2A_GISTIC), shape = factor(Subtype))) + 
#   geom_hline(aes(yintercept=unmethmean), linetype="dashed") + 
#   geom_vline(aes(xintercept=0.3), linetype="dashed") + scale_x_continuous(limits=c(0,1), breaks=c(0, 0.3, 0.5, 0.7, 1)) +
#   scale_y_log10() +
#   ylab("RNA seq expression (log2)") +
#   xlab("DNA Methylation Values") + 
#   labs(title = paste("DNA Methylation vs Expression\n for cg13601799|CDKN2A"), colour = "CDKN2A GISTIC", shape = "IDH + 1p19qStatus")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# gene.choice <- cpg.matrix.sm[cpg.matrix.sm$Gene_Symbol == "CDKN2A", -c(1:4)]
# 
# x <- cbind(t(cpg.matrix.sm[cpg.matrix.sm$Gene_Symbol == "CDKN2A", -c(1:4)]), t(gene.matrix.sm[gene.matrix.sm$SYM == "CDKN2A", -c(1:2)]))
# colnames(x) <- c("cg13601799", "CDKN2A")
# x <- as.data.frame(x)
# x$methstatus <- "NONE"
# x[x$cg13601799 > 0.3, "methstatus"] <- "methylated"
# x[x$cg13601799 <= 0.3, "methstatus"] <- "unmethylated"
# unmethmean <- mean(x[x$methstatus == "unmethylated", "CDKN2A"])
# 
# x$methstatus <- as.factor(x$methstatus)
# 
# 
# require(ggplot2)
# ggplot(x, aes(x=cg13601799, y=CDKN2A, color=methstatus)) +
#   geom_smooth(method=lm, se=FALSE) +
#   geom_point() + geom_hline(aes(yintercept=unmethmean), linetype="dashed") + 
#   geom_vline(aes(xintercept=0.3), linetype="dashed")
# 
