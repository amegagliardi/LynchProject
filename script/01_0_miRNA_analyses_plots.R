library("tidyverse")
library("ggpubr")
require("gplots")
library("ggfortify")
library("circlize")
library("ComplexHeatmap")

path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

###### Figure 1A - PCA ----
count <- "0_Count_table_miRNAs_subjects_no_novel.txt_DESeq2_Normalized_adjusted.tsv"
pdata <- "0_Sample_data_subjects.txt"

DE <- "0_DE_miRNA_subjects_all_IT.txt"

countData <- read.delim(count, row.names=1, check.names = F)
phenoData <- read.delim(pdata, check.names=FALSE)
DE <- read.delim(DE, check.names=FALSE)

phenoData <- dplyr::filter(phenoData, id_mirna != "NO")
phenoData <- dplyr::filter(phenoData, to_remove %in% c("No", "WMS"))
phenoData <- column_to_rownames(phenoData, "id_mirna")
phenoData <- dplyr::filter(phenoData, time %in% c("T0"))
phenoData <- dplyr::filter(phenoData, cohort %in% c("IIGM-INT"))
phenoData <- dplyr::filter(phenoData, diagnosis %in% c("LS"))
phenoData <- dplyr::filter(phenoData, class_sample %in% c("Adenoma", "Cancer", "Negative", "Negative_Hist+"))

countData <- countData[, row.names(phenoData)]
countData <- countData[rowSums(countData) > 0, ]
countData <- countData[DE$ID, ]

cdatat <- data.frame(t(countData), check.names = F)

cdata.pca<-prcomp(log(cdatat+0.1,2), center=T, scale=T)

dat <- data.frame(cbind(phenoData, cdatat, ID=row.names(phenoData)), check.names=F)

dat$diagnosis <- factor(dat$diagnosis, levels=c("Control", "LS"))
dat$cohort <- factor(dat$cohort, levels=c("IIGM-INT", "RUN", "RUN2"))

colorc<-rev(colorRampPalette(brewer.pal(11,"RdBu"))(12))[c(1,3,4,12,10,9)]
colorc<-rev(colorRampPalette(brewer.pal(11,"RdBu"))(12))[c(1,3,12,10)]

# Mutation
colorc<-rev(colorRampPalette(brewer.pal(11,"RdBu"))(12))[c(2,12,10,9,8)]

# Lesion
colorc<-rev(colorRampPalette(brewer.pal(11,"RdBu"))(12))[c(2,12,10,8)]

# Class sample
colorc<-rev(colorRampPalette(brewer.pal(11,"RdBu"))(12))[c(4,2,10,12)]

colorcohort<-colorRampPalette(brewer.pal(8,"Set2"))(8)

dat <- phenoData %>%
  mutate(
    mutation2 = case_when(
      mutation %in% c("MLH1", "MSH2") ~ "HR",
      mutation %in% c("MSH6", "PMS2") ~ "LR",
      mutation %in% c("Control") ~ "Control",
    )
  )

dat$class_sample <- factor(dat$class_sample, levels=c("Negative", "Negative_Hist+", "Adenoma", "Cancer"))

autoplot(cdata.pca, data = dat, colour = "class_sample", size=1.5, 
         alpha=0.75, label=F, label.vjust=2, label.size=1) +
  theme_bw()+
  scale_color_manual(values=colorc)+
  scale_fill_manual(values=colorc)+
  stat_ellipse(geom="polygon", level=0.95, alpha=0.2,aes(fill=class_sample, color=class_sample)) +
  theme(legend.position="right")

### 1B - Scatter plot DE miRNAs ----
a <- read.delim("0_DE_miRNA_subjects_LS_vs_controls_coherent_log2FC.txt", check.names=F)[,c(1:3)]

ggplot(a, aes(x=`INT-IIGM`, y=RUN))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()+
  stat_cor(method="spearman", label.y = 2, label.x = -2.2,
           cor.coef.name = c("rho"), digits = 3, p.digits = 3) +
  labs(x="log2FC(LS-CRL-H- / Control) - IT", 
       y="log2FC(LS-CRL-H- / Control) - NL")+
  geom_hline(yintercept = 0, linetype="dashed", color="red")+
  geom_vline(xintercept = 0, linetype="dashed", color="red")

### Figure 1C - Heatmap DE miRNAs log2FC ------
library("RColorBrewer")
library("ComplexHeatmap")

DE <- read.delim("0_DE_miRNA_subjects_LS_vs_controls_coherent_log2FC.txt", row.names=1, check.names=F)

colorh<-rev(colorRampPalette(brewer.pal(11,"RdBu"))(256))

log2FC <- DE[, c(1:8)]
pval <- DE[, c(9:16)]
  
h1 <- Heatmap(log2FC,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 6),
              cluster_columns = F,
              cluster_rows = T,
              column_split = c(rep("1.CRL-", 2), rep("2.CRL-LP", 2), rep("3.CRL+", 2), rep("4.All LS", 2)),
              name = "log2FC",
              col=colorRamp2(breaks=c(-4,-2,0,2,4), colors=c(colorh[1], colorh[64], "white", colorh[192], colorh[256])),
              border = "black",
              column_names_rot = 45,
              heatmap_legend_param = list(border = "black"),
              width = 90,
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2",
              cell_fun = function(j,i,x,y,w,h,fill) {
                if(pval[i, j] < 0.001) {
                  grid.text("*", x, y, vjust = 0.75)
                } else if(pval[i, j] < 0.01) {
                  grid.text("*", x, y, vjust = 0.75)
                }else if(pval[i, j] < 0.05) {
                  grid.text("*", x, y, vjust = 0.75)
                }},
        clustering_distance_rows = "euclidean")

countData <- read.delim("0_Count_table_miRNAs_subjects_no_novel.txt_DESeq2_Normalized_adjusted.tsv", row.names=1)

phenoData <- read.delim("0_Sample_data_subjects.txt")

phenoData <- dplyr::filter(phenoData, id_mirna != "NO")
phenoData <- dplyr::filter(phenoData, to_remove %in% c("No", "WMS"))
phenoData <- dplyr::filter(phenoData, time %in% c("T0", "Control"))
phenoData <- dplyr::filter(phenoData, is.na(age) == F)
                           
phenoData <- column_to_rownames(phenoData, "id_mirna")

phenoData <- filter(phenoData, cohort %in% c("IIGM-INT", "RUN"))
phenoData <- filter(phenoData, phenoData$class_sample %in% c(c("Control","NEG","LGD","HGD","ADK")))

phenoData2 <-filter(phenoData, diagnosis %in% c("LS"))

count <- countData[row.names(log2FC), row.names(phenoData)]
count2 <- countData[row.names(log2FC), row.names(phenoData2)]

count <- scale(t(log10(count+1)))
count2 <- scale(t(log10(count2+1)))

colorc<-rev(colorRampPalette(brewer.pal(11,"RdBu"))(12))

colorm<-rev(colorRampPalette(brewer.pal(11,"Spectral"))(11))

colorcohort<-colorRampPalette(brewer.pal(8,"Set2"))(8)

colclass <- c("Control_INT"=colorc[2], "LS_INT_T0"=colorc[11], "Control_RUN"=colorc[2], "LS_RUN_T0"=colorc[5], "Control_RUN-V"=colorc[3], "LS_RUN-V"=colorc[6])

colclass <- c("Control"=colorc[2], "LS"=colorc[11])

colrisk <- c("Control"=colorc[2], "High"=colorc[11], "Low"=colorc[9])

colorcohort<-colorRampPalette(brewer.pal(8,"Set2"))(8)

colcohort <- c("IIGM-INT"=colorcohort[1], "RUN"=colorcohort[2], "RUN-V"=colorcohort[3])

collesion.c<-colorRampPalette(brewer.pal(5,"Greys"))(10)

collesion <- c("Control"="white", "NEG"=collesion.c[3], "LGD"=collesion.c[6], "HGD"=collesion.c[8], "ADK"=collesion.c[10], "Other_cancer"="black")

phenoData$mutation <- ifelse(phenoData$mutation=="Control", "NO", phenoData$mutation)

phenoData$mutation <- factor(phenoData$mutation, levels=c("MSH2", "MSH6", "MLH1", "PMS2", "NO"))

colmutation <- c("NO"=colorm[1], "PMS2"=colorm[2], "MSH6"=colorm[4], "MLH1"=colorm[8], "MSH2"=colorm[10], "EPCAM/MSH2"=colorm[6])

phenoData$mutation <- factor(phenoData$mutation, levels=c("MSH2", "MLH1", "EPCAM/MSH2", "MSH6", "PMS2", "NO"))

column_ha = HeatmapAnnotation(Class=phenoData$diagnosis,
                              Lesion=phenoData$class_sample,
                              Mutation=phenoData$mutation,
                              Cohort=phenoData$cohort,
                              col = list(Class = colclass,
                                         Lesion = collesion,
                                         Cohort = colcohort,
                                         Mutation = colmutation),
                              border = TRUE,
                               simple_anno_size = unit(0.3, "cm"),
                              annotation_name_gp = gpar(fontsize =7),
                              annotation_legend_param = list(border="black"))

h2 <- Heatmap(t(count),
              top_annotation = column_ha,
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 1),
              cluster_columns = T,
              cluster_rows = F,
              name = "Z-score",
              col=colorRamp2(breaks=c(-2,-1,0,1,2), colors=c(colorh[1], colorh[64], "white", colorh[192], colorh[256])),
              border = "black",
              column_names_rot = 45,
              heatmap_legend_param = list(border = "black"),
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2", 
              clustering_distance_columns = "euclidean")
h1 + h2

### Figure 1E - Dotplot Enrichment ------
library("RColorBrewer")
library("tidyverse")

a <- read.delim("~/Desktop/a", check.names=F)

colcor <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(256))
a$GS <- factor(a$GS, levels=rev(unique(a$GS)))
a$Class <- factor(a$Class, levels=c("Up", "Down"))

ggplot(a, aes(x=-log10(adj.p.val), y=GS, fill=coef, size=(gene.tested))) +
  geom_point(pch=21) +
  scale_fill_gradientn(colcor, colors=colcor, limits=c(-0.3,0.3), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  facet_wrap(~Class, scale="free", nrow=1) +
  theme_bw() +
  labs(y="Gene set", x="-log10(Adj. p-value)", size="N°target genes", fill="Coefficient", col="Coefficient") +
  theme(text = element_text(size=9),
        axis.text.y = element_text(size=6))

### DE miRNA boxplot
library("ggfortify")
count <- "0_Count_table_miRNAs_subjects_no_novel.txt_DESeq2_Normalized_adjusted.tsv"
pdata <- "0_Sample_data_subjects.txt"

countData <- read.delim(count, row.names=1, check.names = F)
phenoData <- read.delim(pdata, check.names=FALSE)

phenoData <- dplyr::filter(phenoData, id_mirna != "NO")
phenoData <- dplyr::filter(phenoData, to_remove %in% c("No", "WMS"))
phenoData <- column_to_rownames(phenoData, "id_mirna")
phenoData <- dplyr::filter(phenoData, cohort %in% c("IIGM-INT"))

countData <- countData[, row.names(phenoData)]
countData <- countData[rowSums(countData) > 0, ]

miRNA <- c("hsa-miR-5189-5p", "hsa-miR-6855-3p", "hsa-miR-1292-5p", "hsa-miR-1270", "hsa-miR-4634")

miRNA <- c("hsa-miR-1-3p", "hsa-miR-106b-5p", "hsa-miR-4635", "hsa-miR-4701-3p", "hsa-miR-1285-3p")

cdata <-  data.frame(phenoData, t(countData[miRNA, ]), check.names=F)

cdata <- dplyr::select(cdata, c(miRNA, "class2"))
cdata <- gather(cdata, miRNA, levels, -c("class2"))

# Lesion
colorc<-rev(colorRampPalette(brewer.pal(11,"RdBu"))(12))[c(2,12,10,8)]

# Class
colorc<-rev(colorRampPalette(brewer.pal(11,"RdBu"))(12))[c(2,10,12)]

cdata$class2 <- factor(cdata$class2, levels=c("Control_INT", "LS_CRAn_INT_T0", "LS_CRAp_INT_T0"))

cdata <- dplyr::filter(cdata, class2 %in% c("Control_INT", "LS_CRAn_INT_T0", "LS_CRAp_INT_T0"))

ggplot(cdata, aes(y=log10(levels), x=class2, col=class2)) +
  geom_violin(trim=FALSE, alpha=1) +
  scale_shape_manual(values=c(0,1,32)) +
  geom_jitter(width = 0.25, aes(col=class2), alpha=0.5, size=1.5)+
  geom_boxplot(width = 0.1, fill="white", outlier.alpha = 0) + 
  facet_wrap(~miRNA, nrow=1, scales = "free_y") + 
  theme_bw() +
  geom_pwc(hide.ns = T, label="p.signif", vjust = 0.5)+
  labs(y="log10(Normalized levels)", x="Risk group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="right")+
  scale_color_manual(values = colorc)
