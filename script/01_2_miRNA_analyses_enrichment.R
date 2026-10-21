library("RBiomirGS")
library("parallel")
library("tidyverse")
library("doParallel")

source("/Users/giulioferrero/Desktop/Software/0_R/mrnascan_custom.R")

path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

### Read gene set
path = "./../../../0_Data/Data_Annotations/msigdb_v2026.1.Hs_GMTs/"

libraries <- c("c5.go.bp.v2026.1")

### Download of validated miRNA targets
miRNA <- read.delim("~/Desktop/a")

rbiomirgs_mrnascan2(objTitle = "miRNA_target_val", 
                    mir = miRNA$ID,
                    sp = "hsa", 
                    queryType = "validated", 
                    parallelComputing = T, 
                    clusterType="PSOCK")

out_all <- data.frame(matrix(ncol=11))
names(out_all) <- c("GS", "converged", "loss", "gene.tested", "coef", "std.err", "t.value", "p.value", "adj.p.val", "Target_Gene", "Library")

for(i in 1:length(libraries)){
  
  ## Read gs files
  entrez_gmt <- paste0(path, libraries[i], ".Hs.entrez.gmt")
  symbol_gmt <- paste0(path, libraries[i], ".Hs.symbols.gmt")
  
  geneset <- rbiomirgs_gmt(file = entrez_gmt)
  geneset_symbol <- rbiomirgs_gmt(file = symbol_gmt)
  
  ### Logistic regression over gene sets - IWLS
  rbiomirgs_logistic(objTitle = "Table_Tested_iwls", 
                     mirna_DE = miRNA, 
                     var_mirnaName = "ID", 
                     var_mirnaFC = "log2FC", 
                     var_mirnaP = "pval", 
                     mrnalist = miRNA_target_val_mrna_entrez_list,
                     #mrnalist = gene,
                     gs_file = entrez_gmt, 
                     optim_method = "IWLS",
                     p.adj = "fdr", 
                     parallelComputing = T, 
                     clusterType="PSOCK")
  
  ### Retrieval of gene symbols
  gs_entrez <- geneset[as.character(Table_Tested_iwls_GS$GS)]
  gs_symbol <- geneset_symbol[as.character(Table_Tested_iwls_GS$GS)]
  
  Target_Gene <- rep(NA, nrow(Table_Tested_iwls_GS))
  
  for(j in 1:nrow(Table_Tested_iwls_GS)){
    Target_Gene[j] <- paste(gs_symbol[[j]][gs_entrez[[j]] %in% mrnascore$EntrezID], collapse=", ")
  }
  
  Table_Tested_iwls_GS <- cbind(Table_Tested_iwls_GS, Target_Gene)
  
  Library <- rep(libraries[i], nrow(Table_Tested_iwls_GS))
  
  Table_Tested_iwls_GS <- cbind(Table_Tested_iwls_GS, Library)
  
  out <- Table_Tested_iwls_GS
  
  out <- subset(out, out$adj.p.val < 0.1)
  
  if(ncol(out) == 11){
    out_all <- rbind(out_all, out)
  }  
}

out_all <- subset(out_all, out_all$adj.p.val < 0.05)
out_all <- subset(out_all, out_all$gene.tested > 2)

write.table(out_all, paste0("Table_Enrichment_all.tsv"), sep="\t", quote=F, row.names=F)

### Up
### Download of validated miRNA targets
miRNA <- read.delim("~/Desktop/a")
miRNA <- dplyr::filter(miRNA, log2FC > 0)

rbiomirgs_mrnascan2(objTitle = "miRNA_target_val", 
                    mir = miRNA$ID, 
                    sp = "hsa", 
                    queryType = "validated", 
                    parallelComputing = T, 
                    clusterType="PSOCK")

out_all <- data.frame(matrix(ncol=11))
names(out_all) <- c("GS", "converged", "loss", "gene.tested", "coef", "std.err", "t.value", "p.value", "adj.p.val", "Target_Gene", "Library")

for(i in 1:length(libraries)){
  
  ## Read gs files
  entrez_gmt <- paste0(path, libraries[i], ".Hs.entrez.gmt")
  symbol_gmt <- paste0(path, libraries[i], ".Hs.symbols.gmt")
  
  geneset <- rbiomirgs_gmt(file = entrez_gmt)
  geneset_symbol <- rbiomirgs_gmt(file = symbol_gmt)
  
  ### Logistic regression over gene sets - IWLS
  rbiomirgs_logistic(objTitle = "Table_Tested_iwls", 
                     mirna_DE = miRNA, 
                     var_mirnaName = "ID", 
                     var_mirnaFC = "log2FC", 
                     var_mirnaP = "pval", 
                     mrnalist = miRNA_target_val_mrna_entrez_list,
                     #mrnalist = gene,
                     gs_file = entrez_gmt, 
                     optim_method = "IWLS",
                     p.adj = "fdr", 
                     parallelComputing = T, 
                     clusterType="PSOCK")
  
  ### Retrieval of gene symbols
  gs_entrez <- geneset[as.character(Table_Tested_iwls_GS$GS)]
  gs_symbol <- geneset_symbol[as.character(Table_Tested_iwls_GS$GS)]
  
  Target_Gene <- rep(NA, nrow(Table_Tested_iwls_GS))
  
  for(j in 1:nrow(Table_Tested_iwls_GS)){
    Target_Gene[j] <- paste(gs_symbol[[j]][gs_entrez[[j]] %in% mrnascore$EntrezID], collapse=", ")
  }
  
  Table_Tested_iwls_GS <- cbind(Table_Tested_iwls_GS, Target_Gene)
  
  Library <- rep(libraries[i], nrow(Table_Tested_iwls_GS))
  
  Table_Tested_iwls_GS <- cbind(Table_Tested_iwls_GS, Library)
  
  out <- Table_Tested_iwls_GS
  
  out <- subset(out, out$adj.p.val < 0.1)
  
  if(ncol(out) == 11){
    out_all <- rbind(out_all, out)
  }  
}

#write.table(out_all, paste0("Table_Enrichment_", libraries[i], ".tsv"), sep="\t", quote=F, row.names=F)

out_all <- subset(out_all, out_all$adj.p.val < 0.05)
out_all <- subset(out_all, out_all$gene.tested > 2)

write.table(out_all, paste0("Table_Enrichment_up.tsv"), sep="\t", quote=F, row.names=F)

### Down
miRNA <- read.delim("~/Desktop/a")
miRNA <- filter(miRNA, log2FC < 0)

rbiomirgs_mrnascan2(objTitle = "miRNA_target_val", 
                    mir = miRNA$ID, 
                    sp = "hsa", 
                    queryType = "validated", 
                    parallelComputing = T, 
                    clusterType="PSOCK")

out_all <- data.frame(matrix(ncol=11))
names(out_all) <- c("GS", "converged", "loss", "gene.tested", "coef", "std.err", "t.value", "p.value", "adj.p.val", "Target_Gene", "Library")

for(i in 1:length(libraries)){
  
  ## Read gs files
  entrez_gmt <- paste0(path, libraries[i], ".Hs.entrez.gmt")
  symbol_gmt <- paste0(path, libraries[i], ".Hs.symbols.gmt")
  
  geneset <- rbiomirgs_gmt(file = entrez_gmt)
  geneset_symbol <- rbiomirgs_gmt(file = symbol_gmt)
  
  ### Logistic regression over gene sets - IWLS
  rbiomirgs_logistic(objTitle = "Table_Tested_iwls", 
                     mirna_DE = miRNA, 
                     var_mirnaName = "ID", 
                     var_mirnaFC = "log2FC", 
                     var_mirnaP = "pval", 
                     mrnalist = miRNA_target_val_mrna_entrez_list,
                     #mrnalist = gene,
                     gs_file = entrez_gmt, 
                     optim_method = "IWLS",
                     p.adj = "fdr", 
                     parallelComputing = F, 
                     clusterType="PSOCK")
  
  ### Retrieval of gene symbols
  gs_entrez <- geneset[as.character(Table_Tested_iwls_GS$GS)]
  gs_symbol <- geneset_symbol[as.character(Table_Tested_iwls_GS$GS)]
  
  Target_Gene <- rep(NA, nrow(Table_Tested_iwls_GS))
  
  for(j in 1:nrow(Table_Tested_iwls_GS)){
    Target_Gene[j] <- paste(gs_symbol[[j]][gs_entrez[[j]] %in% mrnascore$EntrezID], collapse=", ")
  }
  
  Table_Tested_iwls_GS <- cbind(Table_Tested_iwls_GS, Target_Gene)
  
  Library <- rep(libraries[i], nrow(Table_Tested_iwls_GS))
  
  Table_Tested_iwls_GS <- cbind(Table_Tested_iwls_GS, Library)
  
  out <- Table_Tested_iwls_GS
  
  out <- subset(out, out$adj.p.val < 0.1)
  
  if(ncol(out) == 11){
    out_all <- rbind(out_all, out)
  }  
}

out_all <- subset(out_all, out_all$adj.p.val < 0.05)
out_all <- subset(out_all, out_all$gene.tested > 2)

write.table(out_all, paste0("Table_Enrichment_down.tsv"), sep="\t", quote=F, row.names=F)