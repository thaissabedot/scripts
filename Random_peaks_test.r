library(GenomicRanges)


######################## 
buildTable <- function(rd_peaks.GR, p450k.GR) {
  aux <- suppressWarnings(findOverlaps(rd_peaks.GR,p450k.GR))
  aux <- length(unique(queryHits(aux)))/length(rd_peaks.GR)
  return(aux)
}

##################### Random peaks

buildRandomPeaks <- function(p5hmC.info) {
  
  chrom.size = read.delim("/media/data1/hg38.chrom.sizes.txt", header = F)
  random.peaks = data.frame(seqnames= character(0), start = numeric(0),  end = numeric(0), name = character(0))
  
  n.peaks = dim(p5hmC.info)[1]
  limits <- summary(p5hmC.info$end - p5hmC.info$start)
  levels(random.peaks$seqnames) <- levels(p5hmC.info$seqnames)
  
  ## get the proportion of peaks per chr
  chr.freq <- data.frame(table(p5hmC.info$seqnames))
  
  random.chr = NULL
  random.pos = NULL
  
  for (i in 1: dim(chr.freq)[1]) {
    chr = rep(as.character(chr.freq$Var1[i]), times=chr.freq$Freq[i])
    random.chr = c(random.chr, chr)
    
    chr_size = chrom.size$V2[which(chrom.size$V1 == as.character(chr.freq$Var1[i]))]
    random.pos = c(random.pos, sample(1:chr_size, chr.freq$Freq[i]))
  }
  
  hundredths <- seq(from=limits[2], to=limits[5], by=1)
  lens <- sample(hundredths, size=n.peaks, replace=TRUE)
  
  random.peaks = data.frame(seqnames = random.chr, start = random.pos, end = random.pos + lens, name = paste0("rd_peak_", 1:n.peaks))
  
  rd_peaks.GR <- makeGRangesFromDataFrame(random.peaks, ignore.strand=TRUE)
  return(rd_peaks.GR)
}


ninter = 1000
#################

load("cat_5hmc_v4.rda")
head(cat.5hmC)
dim(cat.5hmC)

####
prop1 <- NULL
for (inter in 1: ninter) {
  print(paste0("Iteraction: ", inter))
  rd.peaks = buildRandomPeaks(cat.5hmC)
  temp_rd = buildTable(rd.peaks, p450k.GR) #p450k.GR could be any set of probes you wish to test the overlap with 5hmC peaks. Just make sure the reference genome is hg38.
  prop1 = c(prop1,temp_rd)
  
}

n <- nrow(cat.5hmC[cat.5hmC$status %in% "loss",])
p <- suppressWarnings(findOverlaps(makeGRangesFromDataFrame(cat.5hmC[cat.5hmC$status %in% "loss",]),p450k.GR))
p <- length(unique(queryHits(p)))/n

#### test
#Null hypothesis: the proportion of random peaks that overlap 450k is equal to the proportion of differentially bound peaks that overlap 450k probes
alpha <- 0.05
p.hat <- mean(prop1)
z.calc <- function(alpha,p.hat,p,n){
	z <- (p.hat - p)/sqrt((p*(1-p))/n)
	p.value <- 2*pnorm(-abs(z)) #pnorm(-abs(z)) se o teste for unilateral. Para 0.05 unilateral, z = -1.645. Para 0.05 bilateral, z = -1.96.
	print(paste0("p-value: ",p.value))
	if(p.value <= alpha)
		print("Reject the null hypothesis") #which means, the proportion is different
	else
		print("Fail to reject the null hypothesis") #which means, the proportion does not change
}
z.calc(alpha,p,p.hat,n)
p/p.hat
