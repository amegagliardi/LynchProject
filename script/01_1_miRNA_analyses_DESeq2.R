library("DESeq2")
library("tidyverse")
library("RColorBrewer")

path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

# Table containing the count data
cdata <- "0_Count_table_miRNAs_subjects_no_novel.txt"

# Table containing the covariates
pdata <- "0_Sample_data_subjects.txt"

#### Read the input files
countData <- read.delim(cdata, row.names=1, check.names=FALSE)
phenoData <- read.delim(pdata, check.names=FALSE)

phenoData <- dplyr::filter(phenoData, id_mirna != "NO")
phenoData <- dplyr::filter(phenoData, to_remove %in% c("No", "WMS"))

#phenoData$age <- scale(phenoData$age)
row.names(phenoData) <- phenoData$id_mirna
  
# Identifier of the reference covariate
refcov <- "comp0b"

# List of covariate to add to the LRT model
covtomodel <- c("sex", "age")

#### Count pre-processing
countData <- countData[, row.names(phenoData)]

#### Identification of the levels of the reference covariate
cov_class <- levels(as.factor(phenoData[,refcov]))

## Adjust for confounding factors
library("sva")

adjusted <- sva::ComBat_seq(as.matrix(countData), batch=phenoData$cohort, group=phenoData$diagnosis)

#### DESeq2 analysis on all the dataset
# Definition of the models
fullmodel <- formula(paste("~", paste(refcov, paste(covtomodel, collapse = "+"), sep="+")))
reducedmodel <- formula(paste("~", paste(covtomodel, collapse="+")))

# Creation of a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = adjusted,
                                       colData = phenoData, 
                                       design = fullmodel)

# DESeq2 analysis
dds <- DESeq(dds, betaPrior=FALSE, parallel = T)

# # Count normalization
norm_counts <- counts(dds, normalized=T)
 
# Computation of the average and the median number of reads in each sample class
medt <- matrix(NA, ncol=length(cov_class), nrow=nrow(norm_counts))
colnames(medt) <- cov_class
 
for(i in 1:length(cov_class)){
   
   subp <- subset(phenoData, phenoData[,refcov] == cov_class[i])
   subc <- norm_counts[,which(colnames(norm_counts) %in% rownames(subp))]
   
   medt[,i] <- rowMedians(subc)
   
   colnames(medt)[i] <- paste0("Median ", cov_class[i])
 }
 
out <- medt
row.names(out) <- rownames(norm_counts)
 
# Write the normalized count table and the table reporting the average and the median read counts
write.table(data.frame("ID"=rownames(norm_counts), round(norm_counts,3), check.names=F), paste0(cdata, "_DESeq2_Normalized_adjusted.tsv"), sep="\t", quote=FALSE, row.names=F)
 
write.table(data.frame("ID"=rownames(out), out), paste0(cdata, "_Median_Average_adjusted.tsv"), sep="\t", quote=FALSE, row.names=F)

# comp0 - T0 LS NEG vs control

#ALL
res_INT <- results(dds, contrast = c("class", "LS_INT_T0", "Control_INT"))
res_RUN <- results(dds, contrast = c("class", "LS_RUN_T0", "Control_RUN"))
write.table(res_INT, "2.2_Table_DESeq_LS_vs_control_T0_INT.tsv", quote = F, sep = "\t")
write.table(res_RUN, "2.2_Table_DESeq_LS_vs_control_T0_RUN.tsv", quote = F, sep = "\t")

#NEG
res_INT <- results(dds, contrast = c("comp0", "IT_NEG", "Control_INT"))
res_RUN <- results(dds, contrast = c("comp0", "RUN_NEG", "Control_RUN"))
write.table(res_INT, "2.2_Table_DESeq_LS_NEG_vs_control_T0_INT.tsv", quote = F, sep = "\t")
write.table(res_RUN, "2.2_Table_DESeq_LS_NEG_vs_control_T0_RUN.tsv", quote = F, sep = "\t")

#NEG_LP
res_INT <- results(dds, contrast = c("comp0", "IT_NEG-LP", "Control_INT"))
res_RUN <- results(dds, contrast = c("comp0", "RUN_NEG-LP", "Control_RUN"))
write.table(res_INT, "2.2_Table_DESeq_LS_NEG_LP_vs_control_T0_INT.tsv", quote = F, sep = "\t")
write.table(res_RUN, "2.2_Table_DESeq_LS_NEG_LP_vs_control_T0_RUN.tsv", quote = F, sep = "\t")

#POS
res_INT <- results(dds, contrast = c("comp0", "IT_POS", "Control_INT"))
res_RUN <- results(dds, contrast = c("comp0", "RUN_POS", "Control_RUN"))
write.table(res_INT, "2.2_Table_DESeq_LS_POS_vs_control_T0_INT.tsv", quote = F, sep = "\t")
write.table(res_RUN, "2.2_Table_DESeq_LS_POS_vs_control_T0_RUN.tsv", quote = F, sep = "\t")

# Comp1 - CRL+ vs CRL- (one sample)
res_comp1it <- results(dds, contrast = c("comp1", "IT_POS", "IT_NEG"))
write.table(res_comp1it, "2.2_Table_DESeq_Comp1_IT_POS_vs_NEG.tsv", quote = F, sep = "\t")

res_comp1run <- results(dds, contrast = c("comp1", "RUN_POS", "RUN_NEG"))
write.table(res_comp1run, "2.2_Table_DESeq_Comp1_RUN_POS_vs_NEG.tsv", quote = F, sep = "\t")

# Comp2 - CRL+ vs CRL- (all samples)
res_comp1it <- results(dds, contrast = c("comp2", "IT_POS", "IT_NEG"))
write.table(res_comp1it, "2.2_Table_DESeq_Comp2_IT_POS_vs_NEG.tsv", quote = F, sep = "\t")

res_comp1run <- results(dds, contrast = c("comp2", "RUN_POS", "RUN_NEG"))
write.table(res_comp1run, "2.2_Table_DESeq_Comp2_RUN_POS_vs_NEG.tsv", quote = F, sep = "\t")

# Comp T0 - T1
res_INT_POS <- results(dds, contrast = c("class0", "LS_INT_T0_POS", "Control_INT"))
res_INT_NEG <- results(dds, contrast = c("class0", "LS_INT_T0_NEG", "Control_INT"))
res_INT_T1_POS <- results(dds, contrast = c("class0", "LS_INT_T1_POS", "Control_INT"))
res_INT_T1_NEG <- results(dds, contrast = c("class0", "LS_INT_T1_NEG", "Control_INT"))

res_RUN_POS <- results(dds, contrast = c("class0", "LS_RUN_T0_POS", "Control_INT"))
res_RUN_NEG <- results(dds, contrast = c("class0", "LS_RUN_T0_NEG", "Control_INT"))
res_RUN_T1_POS <- results(dds, contrast = c("class0", "LS_RUN_T1_POS", "Control_INT"))
res_RUN_T1_NEG <- results(dds, contrast = c("class0", "LS_RUN_T1_NEG", "Control_INT"))

write.table(res_INT_POS, "2.2_Table_DESeq_LS_vs_control_T0_INT_POS.tsv", quote = F, sep = "\t")
write.table(res_INT_NEG, "2.2_Table_DESeq_LS_vs_control_T0_INT_NEG.tsv", quote = F, sep = "\t")
write.table(res_INT_T1_POS, "2.2_Table_DESeq_LS_vs_control_T1_INT_POS.tsv", quote = F, sep = "\t")
write.table(res_INT_T1_NEG, "2.2_Table_DESeq_LS_vs_control_T1_INT_NEG.tsv", quote = F, sep = "\t")

write.table(res_RUN_POS, "2.2_Table_DESeq_LS_vs_control_T0_RUN_POS.tsv", quote = F, sep = "\t")
write.table(res_RUN_NEG, "2.2_Table_DESeq_LS_vs_control_T0_RUN_NEG.tsv", quote = F, sep = "\t")
write.table(res_RUN_T1_POS, "2.2_Table_DESeq_LS_vs_control_T1_RUN_POS.tsv", quote = F, sep = "\t")
write.table(res_RUN_T1_NEG, "2.2_Table_DESeq_LS_vs_control_T1_RUN_NEG.tsv", quote = F, sep = "\t")

# Comp2
res_comp2it <- results(dds, contrast = c("comp2", "IT_POS", "IT_NEG"))
write.table(res_comp2it, "2.2_Table_DESeq_Comp2_IT_POS_vs_NEG.tsv", quote = F, sep = "\t")

res_comp2run <- results(dds, contrast = c("comp2", "RUN_POS", "RUN_NEG"))
write.table(res_comp2run, "2.2_Table_DESeq_Comp2_RUN_POS_vs_NEG.tsv", quote = F, sep = "\t")

# Comp3
res_comp3it <- results(dds, contrast = c("comp345", "IT_HGD", "IT_NEG"))
write.table(res_comp3it, "2.2_Table_DESeq_Comp3_IT_HGD_vs_NEG.tsv", quote = F, sep = "\t")

res_comp3run <- results(dds, contrast = c("comp345", "RUN_HGD", "RUN_NEG"))
write.table(res_comp3run, "2.2_Table_DESeq_Comp3_RUN_HGD_vs_NEG.tsv", quote = F, sep = "\t")

# Comp4
res_comp4it <- results(dds, contrast = c("comp345", "IT_LGD", "IT_NEG"))
write.table(res_comp4it, "2.2_Table_DESeq_Comp4_IT_LGD_vs_NEG.tsv", quote = F, sep = "\t")

res_comp4run <- results(dds, contrast = c("comp345", "RUN_LGD", "RUN_NEG"))
write.table(res_comp4run, "2.2_Table_DESeq_Comp4_RUN_LGD_vs_NEG.tsv", quote = F, sep = "\t")

# Comp5
res_comp5it <- results(dds, contrast = c("comp345", "IT_HGD", "IT_LGD"))
write.table(res_comp5it, "2.2_Table_DESeq_Comp5_IT_HGD_vs_LGD.tsv", quote = F, sep = "\t")

res_comp5run <- results(dds, contrast = c("comp345", "RUN_HGD", "RUN_LGD"))
write.table(res_comp5run, "2.2_Table_DESeq_Comp5_RUN_HGD_vs_LGD.tsv", quote = F, sep = "\t")

# Comp6
res_comp6it <- results(dds, contrast = c("comp6", "IT_POS", "IT_NEG"))
write.table(res_comp6it, "2.2_Table_DESeq_Comp6_IT_POS_vs_NEG.tsv", quote = F, sep = "\t")

res_comp6run <- results(dds, contrast = c("comp6", "RUN_POS", "RUN_NEG"))
write.table(res_comp6run, "2.2_Table_DESeq_Comp6_RUN_POS_vs_NEG.tsv", quote = F, sep = "\t")

# Comp7
res_comp7it <- results(dds, contrast = c("comp7", "IT_POS", "IT_NEG"))
write.table(res_comp7it, "2.2_Table_DESeq_Comp7_IT_POS_vs_NEG.tsv", quote = F, sep = "\t")

res_comp7run <- results(dds, contrast = c("comp7", "RUN_POS", "RUN_NEG"))
write.table(res_comp7run, "2.2_Table_DESeq_Comp7_RUN_POS_vs_NEG.tsv", quote = F, sep = "\t")

# Comp8
res_comp8it <- results(dds, contrast = c("comp8", "IT_NEG-LP", "IT_NEG"))
write.table(res_comp8it, "2.2_Table_DESeq_Comp8_IT_NEG-LP_vs_NEG.tsv", quote = F, sep = "\t")

res_comp8run <- results(dds, contrast = c("comp8", "RUN_NEG-LP", "RUN_NEG"))
write.table(res_comp8run, "2.2_Table_DESeq_Comp8_RUN_NEG-LP_vs_NEG.tsv", quote = F, sep = "\t")

# T0 LS vs control stratified by risk - merged RUN
res_INThr <- results(dds, contrast = c("class5", "LS_HR_INT_T0", "Control_INT"))
res_INTlr <- results(dds, contrast = c("class5", "LS_LR_INT_T0", "Control_INT"))
res_INThlr <- results(dds, contrast = c("class5", "LS_HR_INT_T0", "LS_LR_INT_T0"))

res_RUNhr <- results(dds, contrast = c("class5", "LS_HR_RUN_T0", "Control_RUN"))
res_RUNlr <- results(dds, contrast = c("class5", "LS_LR_RUN_T0", "Control_RUN"))
res_RUNhlr <- results(dds, contrast = c("class5", "LS_HR_RUN_T0", "LS_LR_RUN_T0"))

write.table(res_INThr, "2.2_Table_DESeq_HR_vs_control_T0_INT.tsv", quote = F, sep = "\t")
write.table(res_INTlr, "2.2_Table_DESeq_LR_vs_control_T0_INT.tsv", quote = F, sep = "\t")
write.table(res_RUNhr, "2.2_Table_DESeq_HR_vs_control_T0_RUN.tsv", quote = F, sep = "\t")
write.table(res_RUNlr, "2.2_Table_DESeq_LR_vs_control_T0_RUN.tsv", quote = F, sep = "\t")
write.table(res_RUNhlr, "2.2_Table_DESeq_HR_vs_LR_T0_RUN.tsv", quote = F, sep = "\t")
write.table(res_RUN2hr, "2.2_Table_DESeq_HR_vs_control_T0_RUN2.tsv", quote = F, sep = "\t")
write.table(res_RUN2lr, "2.2_Table_DESeq_LR_vs_control_T0_RUN2.tsv", quote = F, sep = "\t")

write.table(res_INThlr, "2.2_Table_DESeq_HR_vs_LR_T0_INT.tsv", quote = F, sep = "\t")
write.table(res_RUNhlr, "2.2_Table_DESeq_HR_vs_LR_T0_RUN.tsv", quote = F, sep = "\t")
write.table(res_RUN2hlr, "2.2_Table_DESeq_HR_vs_LR_T0_RUN2.tsv", quote = F, sep = "\t")

# T0 LS vs control stratified by lesion
res_INTn <- results(dds, contrast = c("class3", "LS_CRAn_INT_T0", "Control_INT"))
res_INTp <- results(dds, contrast = c("class2", "LS_CRAp_INT_T0", "Control_INT"))
res_INTnh <- results(dds, contrast = c("class3", "LS_CRAnh_INT_T0", "Control_INT"))

res_RUNn <- results(dds, contrast = c("class3", "LS_CRAn_RUN_T0", "Control_RUN"))
res_RUNp <- results(dds, contrast = c("class2", "LS_CRAp_RUN_T0", "Control_RUN"))
res_RUNnh <- results(dds, contrast = c("class3", "LS_CRAnh_RUN_T0", "Control_RUN"))

res_RUN2n <- results(dds, contrast = c("class3", "LS_CRAn_RUN2", "Control_RUN2"))
res_RUN2p <- results(dds, contrast = c("class2", "LS_CRAp_RUN2", "Control_RUN2"))
res_RUN2nh <- results(dds, contrast = c("class3", "LS_CRAnh_RUN2", "Control_RUN2"))

# T0 LS vs control stratified by lesion overall
res_INTn <- results(dds, contrast = c("class3b", "LS_CRAn_INT_T0", "Control_INT"))
res_INTp <- results(dds, contrast = c("class3b", "LS_CRAp_INT_T0", "Control_INT"))
res_RUNn <- results(dds, contrast = c("class3b", "LS_CRAn_RUN_T0", "Control_RUN"))
res_RUNp <- results(dds, contrast = c("class3b", "LS_CRAp_RUN_T0", "Control_RUN"))

# T0 LS CRL+ vs CRL- stratified by lesion
res_INTc <- results(dds, contrast = c("class3b", "LS_CRAp_INT_T0", "LS_CRAn_INT_T0"))
res_RUNc <- results(dds, contrast = c("class3b", "LS_CRAp_RUN_T0", "LS_CRAn_RUN_T0"))
res_RUN2c <- results(dds, contrast = c("class2", "LS_CRAp_RUN2", "LS_CRAn_RUN2"))

write.table(res_INTn, "2.2_Table_DESeq_CRAn_vs_control_T0_INT.tsv", quote = F, sep = "\t")
write.table(res_INTnh, "2.2_Table_DESeq_CRAnh_vs_control_T0_INT.tsv", quote = F, sep = "\t")
write.table(res_INTp, "2.2_Table_DESeq_CRAp_vs_control_T0_INT.tsv", quote = F, sep = "\t")
write.table(res_RUNn, "2.2_Table_DESeq_CRAn_vs_control_T0_RUN.tsv", quote = F, sep = "\t")
write.table(res_RUNnh, "2.2_Table_DESeq_CRAnh_vs_control_T0_RUN.tsv", quote = F, sep = "\t")
write.table(res_RUNp, "2.2_Table_DESeq_CRAp_vs_control_T0_RUN.tsv", quote = F, sep = "\t")

write.table(res_RUN2n, "2.2_Table_DESeq_CRAn_vs_control_T0_RUN2.tsv", quote = F, sep = "\t")
write.table(res_RUN2nh, "2.2_Table_DESeq_CRAnh_vs_control_T0_RUN2.tsv", quote = F, sep = "\t")
write.table(res_RUN2p, "2.2_Table_DESeq_CRAp_vs_control_T0_RUN2.tsv", quote = F, sep = "\t")

write.table(res_INTc, "2.2_Table_DESeq_CRAp_vs_CRAn_T0_INT.tsv", quote = F, sep = "\t")
write.table(res_RUNc, "2.2_Table_DESeq_CRAp_vs_CRAn_T0_RUN.tsv", quote = F, sep = "\t")
write.table(res_RUN2c, "2.2_Table_DESeq_CRAp_vs_CRAn_T0_RUN2.tsv", quote = F, sep = "\t")