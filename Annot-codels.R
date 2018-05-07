## Annotate codels

## Case list
#cases = read.delim("data/freeze/caseid_by_tss.txt", header=T)

## Case list
cases = read.delim("data/freeze/LGG-GBM.final.cases.txt", header=T)

## Load LGG subtypes. Updated by Hailei Zhang on 6/26/2014.
suplgg = openxlsx::read.xlsx("data/clin/IDHsubtypes.v20141219.xlsx")

## Load GBM codels. Annotated by Hailei Zhang june 2014. 
codelgbm = read.delim("data/cn/GBM_1p_19q_subtypes.txt", header=T, as.is=F, row.names=1)

## Transform GBM codeletion table
rownames(codelgbm) <- gsub("GBM\\-(\\w{2})\\-(\\w{4})", "TCGA-\\1-\\2", rownames(codelgbm), perl=T)
levels(codelgbm$codel) = c("N", "Y")
codelgbm$codel[which(rownames(codelgbm) == "TCGA-19-1390")] = "N" # AWG agreed this sample is not a codel

## Match LGG subtypes and cases
idx = match(cases$case, suplgg$id)
cases$subtype = suplgg$subtype[idx]

################
## 1p 19q codel
################

cases$codel1p19q = factor(NA, levels=c("non-codel", "codel"))

idx = which(cases$subtype == "IDHmut-codel")
cases$codel1p19q[idx] = "codel"
idx = which(cases$subtype %in% c("IDHmut-non-codel", "IDHwt"))
cases$codel1p19q[idx] = "non-codel"

## Manually add three non-codels
idx = which(cases$case %in% c("TCGA-DU-7014", "TCGA-DU-A7TI", "TCGA-HW-7493"))
cases$codel1p19q[idx] = "non-codel"

codelgbm = subset(codelgbm, complete.cases(codelgbm$codel))
idx = match(rownames(codelgbm), cases$case)
cases$codel1p19q[idx] = factor(codelgbm$codel, levels=c("N", "Y"), labels=c("non-codel", "codel"))

## save
codel = data.frame(case.id=cases$case, Codel.1p.19q=cases$codel1p19q)
write.table(codel, "results/cn/LGG-GBM.codel.1p.19q.status.txt", sep="\t", quote=F, row.names=F, col.names=T)
save(codel, file="results/cn/LGG-GBM.codel.1p.19q.status.Rdata")
