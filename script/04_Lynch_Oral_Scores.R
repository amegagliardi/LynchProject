########################
### SIAMCAT MEASURES ###
########################

# 1. Setup ----------------------------------------------------------------

## 1.1 Libraries
library(tidyverse)
library(vegan)
library(ggpubr)

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
oral_df <- read_delim(paste0(WMS_dir, "WMS_Oral_Species.txt"), delim = "\n")

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

if (any(colnames(mpa_species) != metadata$sample_id)) {
  stop("Samples in `mpa_species` are not aligned to metadata `sample_id`")
}

# 3. Compute oral scores --------------------------------------------------
oral_df <- oral_df %>%
  tidyr::separate_wider_delim(SGB_oral, delim = "|", names = c("ID", "SGB")) %>%
  dplyr::mutate(ID = str_replace_all(str_trim(ID), "_", " "),
                SGB = str_trim(SGB))

mpa_oral <- mpa_species %>% 
  dplyr::mutate(ID = str_remove_all(ID, "s__"),
                ID = str_replace_all(ID, "_", " "))

oral_scores <- mpa_oral %>%
  dplyr::filter(ID %in% oral_df$ID) %>%
  tidyr::pivot_longer(cols = -ID, names_to = "sample_id", values_to = "abundance") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(Oral_score_abundances = sum(abundance), 
                   Oral_score_species = sum(abundance > 0),.groups = "drop")

metadata <- metadata %>%
  dplyr::left_join(oral_scores, by = "sample_id")

##
comparison_columns <- c("diagnosis", "class", "comp0", "comp1", "comp2")
for (comp in comparison_columns) {
  if (comp == "diagnosis") {
    order_levels <- c("Control", "Lynch")
  } else if (comp == "class") {
    order_levels <- c("Control", "LS_T0")
  } else if (comp == "comp0") {
    order_levels <- c("Control", "NEG", "NEG-LP", "POS")
  } else {
    order_levels <- c("Control", "NEG", "POS")
  }

  oral_richness <- metadata %>% 
    dplyr::select(sample_id, !!sym(comp), time, Oral_score_species) %>% 
    dplyr::mutate(!!sym(comp) := ordered(.data[[comp]], levels = order_levels)) %>% 
    tidyr::drop_na(!!sym(comp)) %>% 
    ggplot(aes(x = .data[[comp]], y = Oral_score_species, fill = .data[[comp]])) +
    geom_boxplot(width = 0.4, outliers = FALSE) +
    geom_jitter(width = 0.3, alpha = 0.5) +
    geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = FALSE) +
    scale_fill_manual(values = lesion_colors, name = "Class") +
    labs(y = "Oral Score (Richness)", title = comp) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  oral_abundance <- metadata %>% 
    dplyr::select(sample_id, !!sym(comp), time, Oral_score_abundances) %>% 
    dplyr::mutate(!!sym(comp) := ordered(.data[[comp]], levels = order_levels)) %>% 
    tidyr::drop_na(!!sym(comp)) %>% 
    ggplot(aes(x = .data[[comp]], y = log10(Oral_score_abundances), fill = .data[[comp]])) +
    geom_boxplot(width = 0.4, outliers = FALSE) +
    geom_jitter(width = 0.3, alpha = 0.5) +
    geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = FALSE) +
    scale_fill_manual(values = lesion_colors, name = "Class") +
    labs(y = "Oral Score (Abundance (Log10))", title = comp) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  plots <- cowplot::plot_grid(oral_richness + guides(fill = "none") + theme(plot.title = element_blank()),
                              oral_abundance + guides(fill = "none") + theme(plot.title = element_blank()),
                              labels = "AUTO")

  title <- cowplot::ggdraw() +
    cowplot::draw_label(comp, fontface = "bold", size = 14)

  oral_plot <- cowplot::plot_grid(title, 
                                  plots,
                                  ncol = 1, 
                                  rel_heights = c(0.05, 1))

  ggsave(plot = oral_plot, filename = paste0("WMS/Figures/Oral_Scores_", comp, ".png"), width = 20, height = 20, units = "cm", dpi = 300)
  ggsave(plot = oral_richness, filename = paste0("WMS/Figures/Oral_Richness_", comp, ".png"), width = 20, height = 20, units = "cm", dpi = 300)
  ggsave(plot = oral_abundance, filename = paste0("WMS/Figures/Oral_Abundance_", comp, ".png"), width = 20, height = 20, units = "cm", dpi = 300)
}
