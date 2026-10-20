########################
### ECOLOGY MEASURES ###
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

## 2.2 Filtering
i <- intersect(colnames(mpa_data)[-1], phenodata$SAMPLE_ID)

metadata <- phenodata %>% 
  dplyr::filter(cohort != "RUN2" & SAMPLE_ID %in% i) %>% 
  dplyr::mutate(diagnosis = ifelse(diagnosis == "LS", "Lynch", diagnosis),
                class = sub("IT_|_INT|RUN_|_RUN", "", class),
                comp0 = sub("IT_|_INT|RUN_|_RUN", "", comp0),
                comp1 = sub("IT_|_INT|RUN_|_RUN", "", comp1),
                comp2 = sub("IT_|_INT|RUN_|_RUN", "", comp2)) %>% 
  dplyr::select(SAMPLE_ID, diagnosis, time, class, comp0, comp1, comp2) %>% 
  dplyr::rename("sample_id" = SAMPLE_ID) %>% 
  dplyr::filter(time %in% c("T0", "Control"))

mpa_species <- mpa_data %>% 
  dplyr::select(clade_name, i[i %in% metadata$sample_id]) %>%
  select_taxon(taxa_id = "s") %>%
  ## Aligning mpa and metadata columns
  dplyr::select(ID, order(match(colnames(.), metadata$sample_id)))

## 2.3 Manipulation
mpa_taxonomy <- mpa_data %>%
  dplyr::filter(stringr::str_detect(clade_name, "\\|t__")) %>%
  dplyr::transmute(ID = paste0(stringr::str_extract(clade_name, "s__[^|]+"), "_", stringr::str_extract(clade_name, "t__[^|]+")),
                   kingdom = stringr::str_match(clade_name, "(?:^|\\|)k__([^|]+)")[, 2],
                   phylum = stringr::str_match(clade_name, "(?:^|\\|)p__([^|]+)")[, 2],
                   class = stringr::str_match(clade_name, "(?:^|\\|)c__([^|]+)")[, 2],
                   order = stringr::str_match(clade_name, "(?:^|\\|)o__([^|]+)")[, 2],
                   family = stringr::str_match(clade_name, "(?:^|\\|)f__([^|]+)")[, 2],
                   genus = stringr::str_match(clade_name, "(?:^|\\|)g__([^|]+)")[, 2],
                   species = stringr::str_match(clade_name, "(?:^|\\|)s__([^|]+)")[, 2],
                   sgb = stringr::str_match(clade_name, "(?:^|\\|)t__([^|]+)")[, 2]) %>%
  dplyr::distinct(ID, .keep_all = TRUE) %>%
  select(-ID)

if (any(colnames(mpa_species)[-1] != metadata$sample_id)) {
  stop("Samples in `mpa_species` are not aligned to metadata `sample_id`")
}
comparison_columns <- c("diagnosis", "class", "comp0", "comp1", "comp2")

# 3. Ecology Measures -----------------------------------------------------
alpha_richness <- mpa_species %>% 
  tidyr::pivot_longer(cols = -ID, names_to = "sample_id", values_to = "abundance") %>% 
  dplyr::filter(abundance > 0) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise(alpha_richness = n())

alpha_shannon <- mpa_species %>% 
  tibble::column_to_rownames("ID") %>% 
  vegan::diversity(index = "shannon", MARGIN = 2) %>% 
  tibble::enframe(name = "sample_id", value = "alpha_shannon")

alpha_invSimpson <- mpa_species %>% 
  tibble::column_to_rownames("ID") %>% 
  vegan::diversity(index = "invsimpson", MARGIN = 2) %>% 
  tibble::enframe(name = "sample_id", value = "alpha_invSimpson")

alpha_diversity_df <- metadata %>% 
  dplyr::left_join(alpha_richness, by = "sample_id") %>% 
  dplyr::left_join(alpha_shannon, by = "sample_id") %>% 
  dplyr::left_join(alpha_invSimpson, by = "sample_id") %>% 
  dplyr::mutate(alpha_pielou = alpha_shannon / log(alpha_richness))

## 3.2 Beta diversity
beta_bray <- mpa_species %>% 
  tibble::column_to_rownames("ID") %>% 
  t() %>% 
  vegan::vegdist(method = "bray")

beta_pcoa <- cmdscale(beta_bray, eig = TRUE, k = 2)
pcoa_variance <- beta_pcoa$eig / sum(beta_pcoa$eig) * 100

beta_diversity_df <- beta_pcoa$points %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sample_id") %>% 
  dplyr::rename(PCoA1 = V1, PCoA2 = V2) %>% 
  dplyr::left_join(metadata, by = "sample_id")

beta_adonis <- adonis2(beta_bray ~ class, data = metadata)
beta_dispersion <- betadisper(beta_bray, group = metadata$class)
beta_anova <- anova(beta_dispersion)
beta_dispersion_perm <- permutest(beta_dispersion, permutations = 999)

beta_results <- purrr::map_dfr(comparison_columns, beta_test, beta_bray = beta_bray, metadata = metadata) %>% 
  dplyr::mutate( PERMANOVA_p_adj = p.adjust(PERMANOVA_p, method = "BH"))
write_delim(beta_results, "WMS/Results/Beta_Diversity_Test.tsv", quote = "none", delim = "\t")

pairwise_beta_results <- purrr::map_dfr(comparison_columns, beta_test_pairwise, beta_bray = beta_bray, metadata = metadata) %>% 
  dplyr::mutate( PERMANOVA_p_adj = p.adjust(PERMANOVA_p, method = "BH"))
write_delim(pairwise_beta_results, "WMS/Results/Beta_Diversity_pairwise-test.tsv", quote = "none", delim = "\t")

# 4. Statistics and plots -------------------------------------------------
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

  richness <- alpha_diversity_df %>% 
    dplyr::select(sample_id, !!sym(comp), time, alpha_richness) %>% 
    dplyr::filter(!!sym(comp) != "NotUsed") %>% 
    dplyr::mutate(!!sym(comp) := ordered(.data[[comp]], levels = order_levels)) %>% 
    tidyr::drop_na(!!sym(comp)) %>% 
    ggplot(aes(x = .data[[comp]], y = alpha_richness, fill = .data[[comp]])) +
    geom_boxplot(width = 0.4, outliers = FALSE) +
    geom_jitter(width = 0.3, alpha = 0.5) +
    geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = FALSE) +
    scale_fill_manual(values = lesion_colors, name = "Class") +
    labs(y = "Richness", title = comp) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  shannon <- alpha_diversity_df %>% 
    dplyr::select(sample_id, !!sym(comp), time, alpha_shannon) %>% 
    dplyr::filter(!!sym(comp) != "NotUsed") %>% 
    dplyr::mutate(!!sym(comp) := ordered(.data[[comp]], levels = order_levels)) %>% 
    tidyr::drop_na(!!sym(comp)) %>% 
    ggplot(aes(x = .data[[comp]], y = alpha_shannon, fill = .data[[comp]])) +
    geom_boxplot(width = 0.4, outliers = FALSE) +
    geom_jitter(width = 0.3, alpha = 0.5) +
    geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = FALSE) +
    scale_fill_manual(values = lesion_colors, name = "Class") +
    labs(y = "Shannon", title = comp) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  invSimp <- alpha_diversity_df %>% 
    dplyr::select(sample_id, !!sym(comp), time, alpha_invSimpson) %>% 
    dplyr::filter(!!sym(comp) != "NotUsed") %>% 
    dplyr::mutate(!!sym(comp) := ordered(.data[[comp]], levels = order_levels)) %>% 
    tidyr::drop_na(!!sym(comp)) %>% 
    ggplot(aes(x = .data[[comp]], y = alpha_invSimpson, fill = .data[[comp]])) +
    geom_boxplot(width = 0.4, outliers = FALSE) +
    geom_jitter(width = 0.3, alpha = 0.5) +
    geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = FALSE) +
    scale_fill_manual(values = lesion_colors, name = "Class") +
    labs(y = "InvSimp", title = comp) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  pielou <- alpha_diversity_df %>% 
    dplyr::select(sample_id, !!sym(comp), time, alpha_pielou) %>% 
    dplyr::filter(!!sym(comp) != "NotUsed") %>% 
    dplyr::mutate(!!sym(comp) := ordered(.data[[comp]], levels = order_levels)) %>% 
    tidyr::drop_na(!!sym(comp)) %>% 
    ggplot(aes(x = .data[[comp]], y = alpha_pielou, fill = .data[[comp]])) +
    geom_boxplot(width = 0.4, outliers = FALSE) +
    geom_jitter(width = 0.3, alpha = 0.5) +
    geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = FALSE) +
    scale_fill_manual(values = lesion_colors, name = "Class") +
    labs(y = "Evenness", title = comp) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  plots <- cowplot::plot_grid(richness + guides(fill = "none") + theme(plot.title = element_blank(), axis.text.x = element_blank()), 
                              shannon + guides(fill = "none") + theme(plot.title = element_blank(), axis.text.x = element_blank()),
                              invSimp + guides(fill = "none") + theme(plot.title = element_blank()), 
                              pielou + guides(fill = "none") + theme(plot.title = element_blank()), 
                              nrow = 2, ncol = 2, 
                              labels = "AUTO")
  title <- cowplot::ggdraw() +
    cowplot::draw_label(comp, fontface = "bold", size = 14)

  alpha_plot <- cowplot::plot_grid(title, 
                                  plots,
                                  ncol = 1, rel_heights = c(0.05, 1))

  ggsave(plot = alpha_plot, filename = paste0("WMS/Figures/Alpha_Diversity_", comp, ".png"), width = 20, height = 20, units = "cm", dpi = 300)
  ggsave(plot = richness, filename = paste0("WMS/Figures/Alpha_Diversity_Richness_", comp, ".png"), width = 20, height = 20, units = "cm", dpi = 300)
  ggsave(plot = shannon, filename = paste0("WMS/Figures/Alpha_Diversity_Shannon_", comp, ".png"), width = 20, height = 20, units = "cm", dpi = 300)
  ggsave(plot = invSimp, filename = paste0("WMS/Figures/Alpha_Diversity_InvSimp_", comp, ".png"), width = 20, height = 20, units = "cm", dpi = 300)
  ggsave(plot = pielou, filename = paste0("WMS/Figures/Alpha_Diversity_Evenness_", comp, ".png"), width = 20, height = 20, units = "cm", dpi = 300)

  ## 4.2 Beta diversity
  beta_plot <- beta_diversity_df %>% 
    # dplyr::filter(diagnosis != "NotUsed") %>% 
    dplyr::mutate(!!sym(comp) := ordered(.data[[comp]], levels = order_levels)) %>% 
    ggplot(aes(x = PCoA1, y = PCoA2, color = .data[[comp]])) +
    geom_point(size = 3.5) +
    scale_color_manual(values = lesion_colors, name = comp) +
    labs(x = paste0("PCoA-1 (", round(pcoa_variance[1], 1), "%)"), 
        y = paste0("PCoA-2 (", round(pcoa_variance[2], 1), "%)"),
        title = "Class") +
    theme_bw()

  ggsave(plot = beta_plot, filename = paste0("WMS/Figures/Beta_Diversity_", comp, ".png"), width = 20, height = 20, units = "cm", dpi = 300)
}
