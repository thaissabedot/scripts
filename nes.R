samples <- rownames(annot.order.cancer[annot.order.cancer$disease_type %in% "Glioblastoma Multiforme",])

norm.order.normal <- sort(norm.order.normal[1,])
b <- norm.order.normal
b[2,] <- norm.order.normal

norm.order.cancer <- sort(norm.order.cancer[1,])
a <- norm.order.cancer
a[2,] <- norm.order.cancer

summary(t(a[1,]))

metadata <- annot.order.cancer[substr(colnames(a),1,15),]
identical(rownames(metadata),substr(colnames(a),1,15)) #TRUE
metadata$cut.NES.all <- NA

aux <- as.data.frame(t(a))
metadata[substr(rownames(aux[aux$NES < quantile(aux$NES,prob=5/100),]),1,15),"cut.NES.all"] <- "Top 5% with lower expression" #menor valor de expressao
metadata[substr(rownames(aux[aux$NES > quantile(aux$NES,prob=1-5/100),]),1,15),"cut.NES.all"] <- "Top 5% with higher expression" #maior valor de expressao
metadata[substr(rownames(aux[aux$NES > quantile(aux$NES,prob=5/100) & aux$NES < quantile(aux$NES,prob=1-5/100),]),1,15),"cut.NES.all"] <- "Intermediate expression" #resto

metadata <- metadata[metadata$cut.NES.all %in% c("Top 5% with lower expression","Top 5% with higher expression"),]


TCGAanalyze_survival(metadata,clusterCol = "cut.NES.all", filename = "all_gliomas_NES_LowHigh.png")

ha = HeatmapAnnotation(metadata[,c(2:7, 20, 29:35,37:42, 48:54, 63,
                                             grep("gistic2",colnames(annot.order.cancer)),
                                             grep("mut",colnames(annot.order.cancer)))],
                       show_annotation_name = TRUE,
                       annotation_name_side = "left",
                       annotation_name_gp = gpar(fontsize = 6),
                       #annotation_height = unit(5, "mm"),
                       annotation_legend_param = list(nrow = 1, labels_gp = gpar(fontsize = 4), title_gp = gpar(fontsize = 4),grid_height = unit(2, "mm"),
                               grid_width = unit(2, "mm")),
                       col = c(cnv.col,
                               mut.col,
                               list("definition"=c("Primary solid Tumor"="blue","Recurrent Solid Tumor"= "green","Solid Tissue Normal"="black")))
)
aux <- norm.order.cancer[,substr(colnames(norm.order.cancer),1,15) %in% rownames(metadata)]
identical(rownames(metadata),substr(colnames(aux),1,15)) #TRUE
#heatmap <- Heatmap(norm.order.cancer, col = greenred(255),
heatmap <- Heatmap(aux, col = greenred(255),
                   show_column_names = F,
                   top_annotation = ha,row_names_side = "left",
                   show_row_names = TRUE, cluster_columns = F , cluster_rows = T,
                   heatmap_legend_param = list(color_bar = "continuous"),row_names_gp =  gpar(fontsize = 6)) +
    Heatmap(norm.order.normal,
            col = greenred(255),
            column_title_gp = gpar(fontsize = 10, fontface = "bold"),
            heatmap_legend_param = list(color_bar = "continuous"),
            #top_annotation = ha.normal,
            #column_title = "Normal samples",
            show_column_names = F,show_row_names = F, cluster_columns = F,
            show_heatmap_legend = F)

pdf("heatmap_all_gliomas_LowHigh.pdf", width = 15, height = 25)
draw(heatmap, annotation_legend_side = "left")
dev.off()

########## GBM ONLY

table(annot.order.cancer$disease_type)

metadata <- annot.order.cancer[substr(colnames(a),1,15),]
identical(rownames(metadata),substr(colnames(a),1,15)) #TRUE

samples <- rownames(metadata[metadata$disease_type %in% "Glioblastoma Multiforme",])
gbm <- a[,substr(colnames(a),1,15) %in% samples]

gbm.meta <- metadata[samples,]
identical(rownames(gbm.meta),substr(colnames(gbm),1,15)) #TRUE

gbm.meta$cut.NES.all <- NA

aux <- as.data.frame(t(gbm))
gbm.meta[substr(rownames(aux[aux$NES < quantile(aux$NES,prob=5/100),]),1,15),"cut.NES.all"] <- "Top 5% with lower expression" #menor valor de expressao
gbm.meta[substr(rownames(aux[aux$NES > quantile(aux$NES,prob=1-5/100),]),1,15),"cut.NES.all"] <- "Top 5% with higher expression" #maior valor de expressao
gbm.meta[substr(rownames(aux[aux$NES > quantile(aux$NES,prob=5/100) & aux$NES < quantile(aux$NES,prob=1-5/100),]),1,15),"cut.NES.all"] <- "Intermediate expression" #resto

gbm.meta <- gbm.meta[gbm.meta$cut.NES.all %in% c("Top 5% with lower expression","Top 5% with higher expression"),]

TCGAanalyze_survival(gbm.meta,clusterCol = "cut.NES.all", filename = "gbm_NES_5_LowHigh.png")

ha = HeatmapAnnotation(gbm.meta[,c(2:7, 20, 29:35,37:42, 48:54, 63,
                                             grep("gistic2",colnames(annot.order.cancer)),
                                             grep("mut",colnames(annot.order.cancer)))],
                       show_annotation_name = TRUE,
                       annotation_name_side = "left",
                       annotation_name_gp = gpar(fontsize = 6),
                       #annotation_height = unit(5, "mm"),
                       annotation_legend_param = list(nrow = 1, labels_gp = gpar(fontsize = 4), title_gp = gpar(fontsize = 4),grid_height = unit(2, "mm"),
                               grid_width = unit(2, "mm")),
                       col = c(cnv.col,
                               mut.col,
                               list("definition"=c("Primary solid Tumor"="blue","Recurrent Solid Tumor"= "green","Solid Tissue Normal"="black")))
)
aux <- gbm[,substr(colnames(gbm),1,15) %in% rownames(gbm.meta)]
identical(rownames(gbm.meta),substr(colnames(aux),1,15)) #TRUE
#heatmap <- Heatmap(gbm[1,], col = greenred(255),
heatmap <- Heatmap(aux[1,], col = greenred(255),
                   show_column_names = F,
                   top_annotation = ha,row_names_side = "left",
                   show_row_names = TRUE, cluster_columns = F , cluster_rows = T,
                   heatmap_legend_param = list(color_bar = "continuous"),row_names_gp =  gpar(fontsize = 6)) +
    Heatmap(norm.order.normal,
            col = greenred(255),
            column_title_gp = gpar(fontsize = 10, fontface = "bold"),
            heatmap_legend_param = list(color_bar = "continuous"),
            #top_annotation = ha.normal,
            #column_title = "Normal samples",
            show_column_names = F,show_row_names = F, cluster_columns = F,
            show_heatmap_legend = F)

pdf("heatmap_clustered_in_gbm_LowHigh.pdf", width = 15, height = 25)
draw(heatmap, annotation_legend_side = "left")
dev.off()


########## IDHwt ONLY

metadata <- annot.order.cancer[substr(colnames(a),1,15),]
identical(rownames(metadata),substr(colnames(a),1,15)) #TRUE


samples <- rownames(metadata[metadata$subtype_IDH.status %in% "WT",])
wt <- a[,substr(colnames(a),1,15) %in% samples]

wt.meta <- metadata[samples,]
identical(rownames(wt.meta),substr(colnames(wt),1,15)) #TRUE

wt.meta$cut.NES.all <- NA

aux <- as.data.frame(t(wt))
wt.meta[substr(rownames(aux[aux$NES < quantile(aux$NES,prob=5/100),]),1,15),"cut.NES.all"] <- "Top 5% with lower expression" #menor valor de expressao
wt.meta[substr(rownames(aux[aux$NES > quantile(aux$NES,prob=1-5/100),]),1,15),"cut.NES.all"] <- "Top 5% with higher expression" #maior valor de expressao
wt.meta[substr(rownames(aux[aux$NES > quantile(aux$NES,prob=5/100) & aux$NES < quantile(aux$NES,prob=1-5/100),]),1,15),"cut.NES.all"] <- "Intermediate expression" #resto

wt.meta <- wt.meta[wt.meta$cut.NES.all %in% c("Top 5% with lower expression","Top 5% with higher expression"),]

TCGAanalyze_survival(wt.meta,clusterCol = "cut.NES.all", filename = "IDHwt_NES_5_LowHigh.png")

ha = HeatmapAnnotation(wt.meta[,c(2:7, 20, 29:35,37:42, 48:54, 63,
                                             grep("gistic2",colnames(annot.order.cancer)),
                                             grep("mut",colnames(annot.order.cancer)))],
                       show_annotation_name = TRUE,
                       annotation_name_side = "left",
                       annotation_name_gp = gpar(fontsize = 6),
                       #annotation_height = unit(5, "mm"),
                       annotation_legend_param = list(nrow = 1, labels_gp = gpar(fontsize = 4), title_gp = gpar(fontsize = 4),grid_height = unit(2, "mm"),
                               grid_width = unit(2, "mm")),
                       col = c(cnv.col,
                               mut.col,
                               list("definition"=c("Primary solid Tumor"="blue","Recurrent Solid Tumor"= "green","Solid Tissue Normal"="black")))
)
aux <- wt[,substr(colnames(wt),1,15) %in% rownames(wt.meta)]
identical(rownames(wt.meta),substr(colnames(aux),1,15)) #TRUE
#heatmap <- Heatmap(wt[1,], col = greenred(255),
heatmap <- Heatmap(aux[1,], col = greenred(255),
                   show_column_names = F,
                   top_annotation = ha,row_names_side = "left",
                   show_row_names = TRUE, cluster_columns = F , cluster_rows = T,
                   heatmap_legend_param = list(color_bar = "continuous"),row_names_gp =  gpar(fontsize = 6)) +
    Heatmap(norm.order.normal,
            col = greenred(255),
            column_title_gp = gpar(fontsize = 10, fontface = "bold"),
            heatmap_legend_param = list(color_bar = "continuous"),
            #top_annotation = ha.normal,
            #column_title = "Normal samples",
            show_column_names = F,show_row_names = F, cluster_columns = F,
            show_heatmap_legend = F)

pdf("heatmap_clustered_in_IDHwt_LowHigh.pdf", width = 15, height = 25)
draw(heatmap, annotation_legend_side = "left")
dev.off()
