########################
### SIAMCAT MEASURES ###
########################

# 1. Setup ----------------------------------------------------------------

## 1.1 Libraries
library(tidyverse)
library(vegan)
library(ggpubr)
library(SIAMCAT)

## 1.2 Sources
WMS_dir <- "./WMS/Scripts/"
source(paste0(WMS_dir, "WMS_Custom_Functions.R"))
source(paste0(WMS_dir, "CRC_Lynch_Custom_Palettes.R"))

## seed for reproducibility
set.seed(1970)

# 2. Data -----------------------------------------------------------------

## 2.1 Loading
phenodata <- read_delim("WMS/data/0_metadata.tsv")
mpa_data <- read_delim("WMS/data/mpa_Dutch.tsv")

## 2.2 Filtering
i <- intersect(colnames(mpa_data)[-1], phenodata$SAMPLE_ID)

metadata <- phenodata %>% 
  dplyr::filter(cohort != "RUN2" & SAMPLE_ID %in% i) %>% 
  dplyr::mutate(diagnosis = ifelse(diagnosis == "LS", "Lynch", diagnosis),
                class = sub("IT_|_INT|RUN_|_RUN", "", class),
                comp0 = sub("IT_|_INT|RUN_|_RUN", "", comp0),
                comp1 = sub("IT_|_INT|RUN_|_RUN", "", comp1),
                comp2 = sub("IT_|_INT|RUN_|_RUN", "", comp2)) %>% 
  dplyr::select(SAMPLE_ID, diagnosis, time, class, comp0, comp1, comp2, age, sex) %>% 
  dplyr::rename("sample_id" = SAMPLE_ID) %>% 
  dplyr::filter(time %in% c("T0", "Control"))

mpa_species <- mpa_data %>%
  dplyr::select(clade_name, i[i %in% metadata$sample_id]) %>%
  select_taxon(taxa_id = "s") %>%
  ## Aligning mpa and metadata columns
  dplyr::select(ID, order(match(colnames(.), metadata$sample_id)))

mpa_species <- filter_nzv(mpa_species, scale_by_100 = 100)

if (any(colnames(mpa_species) != metadata$sample_id)) {
  stop("Samples in `mpa_species` are not aligned to metadata `sample_id`")
}

# 3. Data -----------------------------------------------------------------
comparison_columns <- c("diagnosis", "class", "comp0", "comp1", "comp2")

for (comp in comparison_columns) {
  if (comp == "diagnosis") {
    comparison_classes <- c("Control", "Lynch")
  } else if (comp == "class") {
    comparison_classes <- c("Control", "LS_T0")
  } else if (comp == "comp0") {
    comparison_classes <- c("Control", "NEG", "NEG-LP", "POS")
  } else {
    comparison_classes <- c("Control", "NEG", "POS")
  }
  columns_to_select  <- c("sample_id", comp, "age", "sex")

  siamcat_metadata <- metadata %>%
    dplyr::select(all_of(columns_to_select)) %>% 
    # dplyr::filter(str_detect(sample_id, sample_type)) %>% 
    dplyr::filter(.data[[comp]] %in% comparison_classes) %>% 
    tibble::column_to_rownames("sample_id") %>% 
    as.data.frame()

  siamcat_abundances <- mpa_species %>% 
    # column_to_rownames("ID") %>% 
    dplyr::select(all_of(rownames(siamcat_metadata))) %>% 
    as.matrix()

  dim(siamcat_metadata); dim(siamcat_abundances)

  ## Label creation
  pairs <- utils::combn(unique(as.character(comparison_classes)), 2)
  for (i in 1:ncol(pairs)) {
    pair <- pairs[, i]
    test <- paste0(comp ,"_", pair[2], "_vs_", pair[1])

    label.phenodata <- create.label(meta = siamcat_metadata,
                                    label = comp,
                                    case = pair[2],
                                    control = pair[1])

    ## Create Siamcat Object
    sc.obj <- siamcat(feat = siamcat_abundances,
                      meta = siamcat_metadata,
                      label = label.phenodata)

    ### Filtering
    sc.obj <- filter.features(sc.obj,
                              filter.method = 'abundance',
                              cutoff = 0.00001)

    ### Association test
    sc.obj <- check.associations(sc.obj,
                                formula = "feat ~ label + age + sex",
                                test = "lm",
                                alpha = 0.05,
                                mult.corr = "fdr",
                                pr.cutoff = 1e-6)

    association.table <- sc.obj@associations$assoc.results %>% 
      tibble::rownames_to_column("Feature_ID") %>% 
      tibble::as_tibble() %>% 
      dplyr::mutate(Test = test)

    association.plot(sc.obj, 
                    sort.by = "fc",
                    fn.plot = paste0("WMS/figures/SIAMCAT_AssociationPlot_s-level_", test, ".pdf"),
                    color.scheme = lesion_colors[pair[1:2]],
                    max.show = 25,
                    panels = c("fc", "prevalence"))

    write_delim(association.table, file = paste0("WMS/Results/SIAMCAT_AssociationTable_", test, ".txt"), delim = "\t", quote = "none")
  }  
}
