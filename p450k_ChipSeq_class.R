load("/media/data1/Recurrent_hg38.rda")


tss <- promoters(makeGRangesFromDataFrame(transcripts,TRUE), upstream = 0, downstream = 1)
a <-  distanceToNearest(makeGRangesFromDataFrame(Recurrent_CIMP),tss)
a <- as.data.frame(a)
a <- a[a$distance < 1000,]
a$entrezID <- mcols(tss[a$subjectHits,])[,3]
Recurrent_CIMP$on.TSS <- FALSE
Recurrent_CIMP[a$queryHits,"on.TSS"] <- as.character(a$entrezID)

h3k4me3.r <- dba.report(h3k4me3.a)
a <-  distanceToNearest(makeGRangesFromDataFrame(Recurrent_CIMP),h3k4me3.r)
a <- as.data.frame(a)
a <- a[a$distance == 0,]
a$Fold <- mcols(h3k4me3.r[a$subjectHits,])[,4]


Recurrent_CIMP$on.H3K4me3.gained.in.GCIMPlow <- FALSE
Recurrent_CIMP[a[a$Fold < 0, "queryHits"],"on.H3K4me3.gained.in.GCIMPlow"] <- TRUE

Recurrent_CIMP$on.H3K4me3.lost.in.GCIMPlow <- FALSE
Recurrent_CIMP[a[a$Fold > 0, "queryHits"],"on.H3K4me3.lost.in.GCIMPlow"] <- TRUE


h3k27ac.r <- dba.report(h3k27ac.a)
a <-  distanceToNearest(makeGRangesFromDataFrame(Recurrent_CIMP),h3k27ac.r)
a <- as.data.frame(a)
a <- a[a$distance == 0,]
a$Fold <- mcols(h3k27ac.r[a$subjectHits,])[,4]


Recurrent_CIMP$on.H3K27ac.gained.in.GCIMPlow <- FALSE
Recurrent_CIMP[a[a$Fold < 0, "queryHits"],"on.H3K27ac.gained.in.GCIMPlow"] <- TRUE

Recurrent_CIMP$on.H3K27ac.lost.in.GCIMPlow <- FALSE
Recurrent_CIMP[a[a$Fold > 0, "queryHits"],"on.H3K27ac.lost.in.GCIMPlow"] <- TRUE

ctcf.r <- dba.report(ctcf.a,method = DBA_EDGER)
a <-  distanceToNearest(makeGRangesFromDataFrame(Recurrent_CIMP),ctcf.r)
a <- as.data.frame(a)
a <- a[a$distance == 0,]
a$Fold <- mcols(ctcf.r[a$subjectHits,])[,4]


Recurrent_CIMP$on.CTCF.gained.in.GCIMPlow <- FALSE
Recurrent_CIMP[a[a$Fold < 0, "queryHits"],"on.CTCF.gained.in.GCIMPlow"] <- TRUE

Recurrent_CIMP$on.CTCF.lost.in.GCIMPlow <- FALSE
Recurrent_CIMP[a[a$Fold > 0, "queryHits"],"on.CTCF.lost.in.GCIMPlow"] <- TRUE

hmec.r <- dba.report(hmec.a)
a <-  distanceToNearest(makeGRangesFromDataFrame(Recurrent_CIMP),hmec.r)
a <- as.data.frame(a)
a <- a[a$distance == 0,]
a$Fold <- mcols(hmec.r[a$subjectHits,])[,4]


Recurrent_CIMP$on.5hmC.gained.in.GCIMPlow <- FALSE
Recurrent_CIMP[a[a$Fold < 0, "queryHits"],"on.5hmC.gained.in.GCIMPlow"] <- TRUE

Recurrent_CIMP$on.5hmC.lost.in.GCIMPlow <- FALSE
Recurrent_CIMP[a[a$Fold > 0, "queryHits"],"on.5hmC.lost.in.GCIMPlow"] <- TRUE

save(Recurrent_CIMP,file="/media/data1/Recurrent_CIMP_MACSsummit.rda")