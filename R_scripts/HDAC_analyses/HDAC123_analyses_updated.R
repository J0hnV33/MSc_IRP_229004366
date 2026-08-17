# Combined HDAC1/2/3 RNA-seq analysis
#
# Integrates the HDAC1/2 no-H3 dataset, HDAC3 dataset, and combined DESeq2
# result table for PCA, differential-expression, X-chromosome, and marker-gene
# analyses. Sections follow the existing execution order.
# Update the local input and output paths before running on another system.

# 1. Setup and result import

## 1.1 Packages 

library(DESeq2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(ggvenn)
library(ggVennDiagram)

## 1.2 Paths and output directories

base_dir <- "/home/jvk3/Desktop/HDAC_counts"

hdac12_counts_dir <- file.path(base_dir, "2_STAR_libinorm_counts")
hdac3_counts_dir  <- file.path(base_dir, "HDAC3_geneIDs")

all_res_dir <- file.path(base_dir, "All_Res_Tables")

combined_plot_dir <- file.path(base_dir, "Combined_HDAC1_2_3_outputs")
dir.create(combined_plot_dir, showWarnings = FALSE, recursive = TRUE)

## 1.3 Load the combined DESeq2 result table

all_res_hdac1_2_3_15fold <- read.csv(
  file.path(all_res_dir, "all_results_HDAC1_2_3_combined_15fold.csv"),
  stringsAsFactors = FALSE
)

all_res_hdac1_2_3_15fold <- all_res_hdac1_2_3_15fold %>%
  dplyr::mutate(
    contrast_short = factor(
      contrast_short,
      levels = c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
    ),
    experiment = factor(
      experiment,
      levels = c("HDAC12", "HDAC3")
    ),
    deg_status_15fold = factor(
      deg_status_15fold,
      levels = c("Down", "Not significant", "Up")
    ),
    CHR = as.character(CHR),
    START = suppressWarnings(as.numeric(START)),
    END = suppressWarnings(as.numeric(END))
  )

dim(all_res_hdac1_2_3_15fold)

all_res_hdac1_2_3_15fold %>%
  dplyr::count(experiment, contrast_short, deg_status_15fold)

# 2. Combined HDAC1/2/3 PCA

## 2.1 Define sample metadata

meta_pca <- tibble::tibble(
  sample_id = c(
    "A1", "A2", "A3",
    "H1", "H2",
    "B1", "B2", "B3",
    "A1", "A2", "B3",
    "B1", "B2", "A3"
  ),
  unique_sample = c(
    "HDAC12_A1", "HDAC12_A2", "HDAC12_A3",
    "HDAC12_H1", "HDAC12_H2",
    "HDAC12_B1", "HDAC12_B2", "HDAC12_B3",
    "HDAC3_A1", "HDAC3_A2", "HDAC3_B3",
    "HDAC3_B1", "HDAC3_B2", "HDAC3_A3"
  ),
  experiment = c(
    rep("HDAC12", 8),
    rep("HDAC3", 6)
  ),
  condition = c(
    "Scramble", "Scramble", "Scramble",
    "HDAC1_siRNA", "HDAC1_siRNA",
    "HDAC2_siRNA", "HDAC2_siRNA", "HDAC2_siRNA",
    "Scramble", "Scramble", "Scramble",
    "HDAC3_siRNA", "HDAC3_siRNA", "HDAC3_siRNA"
  ),
  rep = c(
    1, 2, 3,
    1, 2,
    1, 2, 3,
    1, 2, 3,
    1, 2, 3
  ),
  count_dir = c(
    rep(hdac12_counts_dir, 8),
    rep(hdac3_counts_dir, 6)
  )
) %>%
  dplyr::mutate(
    pca_group = dplyr::case_when(
      experiment == "HDAC12" & condition == "Scramble" ~ "HDAC12 Scramble",
      condition == "HDAC1_siRNA" ~ "HDAC1 KD",
      condition == "HDAC2_siRNA" ~ "HDAC2 KD",
      experiment == "HDAC3" & condition == "Scramble" ~ "HDAC3 Scramble",
      condition == "HDAC3_siRNA" ~ "HDAC3 KD",
      TRUE ~ condition
    ),
    sample_label = dplyr::case_when(
      pca_group == "HDAC12 Scramble" ~ paste0("HDAC12_Scr_rep", rep),
      pca_group == "HDAC1 KD" ~ paste0("HDAC1_KD_rep", rep),
      pca_group == "HDAC2 KD" ~ paste0("HDAC2_KD_rep", rep),
      pca_group == "HDAC3 Scramble" ~ paste0("HDAC3_Scr_rep", rep),
      pca_group == "HDAC3 KD" ~ paste0("HDAC3_KD_rep", rep),
      TRUE ~ unique_sample
    )
  ) %>%
  as.data.frame()

rownames(meta_pca) <- meta_pca$unique_sample

## 2.2 Import individual count files

read_count_for_pca <- function(directory, sample_id, unique_sample) {

  file <- list.files(
    directory,
    pattern = paste0("^", sample_id, "_.*_counts\\.txt$"),
    full.names = TRUE
  )

  file <- file[
    !grepl("bias|expression|distribution", basename(file), ignore.case = TRUE)
  ]

  if (length(file) != 1) {
    stop(paste("Problem finding count file for", unique_sample))
  }

  message(unique_sample, " -> ", basename(file))

  out <- read.delim(file, header = FALSE, stringsAsFactors = FALSE) %>%
    dplyr::select(gene_id = 1, count = 2) %>%
    dplyr::filter(!grepl("^__", gene_id)) %>%
    dplyr::mutate(
      gene_id = gsub("\\..*$", "", gene_id),
      count = as.numeric(count)
    ) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(
      count = sum(count, na.rm = TRUE),
      .groups = "drop"
    )

  colnames(out)[2] <- unique_sample

  return(out)
}

count_list_pca <- mapply(
  read_count_for_pca,
  directory = meta_pca$count_dir,
  sample_id = meta_pca$sample_id,
  unique_sample = meta_pca$unique_sample,
  SIMPLIFY = FALSE
)

## 2.3 Build the combined count matrix

counts_pca <- Reduce(
  function(x, y) dplyr::full_join(x, y, by = "gene_id"),
  count_list_pca
) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

counts_pca[is.na(counts_pca)] <- 0
storage.mode(counts_pca) <- "integer"

counts_pca <- counts_pca[, rownames(meta_pca)]

stopifnot(all(colnames(counts_pca) == rownames(meta_pca)))

## 2.4 Fit the PCA model

meta_pca$pca_group <- factor(
  meta_pca$pca_group,
  levels = c(
    "HDAC12 Scramble",
    "HDAC1 KD",
    "HDAC2 KD",
    "HDAC3 Scramble",
    "HDAC3 KD"
  )
)

dds_pca <- DESeqDataSetFromMatrix(
  countData = counts_pca,
  colData = meta_pca,
  design = ~ pca_group
)

keep <- rowSums(counts(dds_pca)) >= 10
dds_pca <- dds_pca[keep, ]

vsd_pca <- vst(dds_pca, blind = FALSE)

pca_combined_raw <- plotPCA(
  vsd_pca,
  intgroup = c("experiment", "pca_group"),
  ntop = 500,
  returnData = TRUE
)

percentVar_combined <- round(
  100 * attr(pca_combined_raw, "percentVar")
)

pca_combined_df <- as.data.frame(pca_combined_raw)
pca_combined_df$sample <- rownames(pca_combined_df)

pca_combined_df$sample_label <- meta_pca[
  pca_combined_df$sample,
  "sample_label"
]

## 2.5 Plot the combined PCA

pca_combined_plot <- ggplot(
  pca_combined_df,
  aes(
    x = PC1,
    y = PC2,
    colour = pca_group,
    shape = experiment,
    label = sample_label
  )
) +
  geom_point(size = 4, alpha = 0.9) +
  ggrepel::geom_text_repel(
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  theme_bw() +
  labs(
    title = "PCA of HDAC1/2 and HDAC3 samples",
    subtitle = "HDAC1/2 analysed without HDAC12_H3; PCA shown for sample-level overview only",
    x = paste0("PC1: ", percentVar_combined[1], "% variance"),
    y = paste0("PC2: ", percentVar_combined[2], "% variance"),
    colour = "Group",
    shape = "Experiment"
  )

pca_combined_plot

ggsave(
  filename = file.path(combined_plot_dir, "01_PCA_combined_HDAC1_2_3.png"),
  plot = pca_combined_plot,
  width = 9,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

# 3. Differential-expression summaries and target validation

## 3.1 Summarise DEG counts

deg_count_plot_df_combined <- dplyr::bind_rows(

  all_res_hdac1_2_3_15fold %>%
    dplyr::filter(deg_status_15fold != "Not significant") %>%
    dplyr::count(
      gene_set = "Genome-wide",
      contrast_short,
      deg_status_15fold,
      name = "n_genes"
    ),

  all_res_hdac1_2_3_15fold %>%
    dplyr::filter(
      CHR == "X",
      deg_status_15fold != "Not significant"
    ) %>%
    dplyr::count(
      gene_set = "Chromosome X",
      contrast_short,
      deg_status_15fold,
      name = "n_genes"
    )
) %>%
  dplyr::mutate(
    gene_set = factor(gene_set, levels = c("Chromosome X", "Genome-wide")),
    deg_status_15fold = factor(deg_status_15fold, levels = c("Down", "Up")),
    contrast_short = factor(contrast_short, levels = c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD"))
  )

deg_count_plot_df_combined

deg_count_plot_combined <- ggplot(
  deg_count_plot_df_combined,
  aes(
    x = contrast_short,
    y = n_genes,
    fill = deg_status_15fold
  )
) +
  geom_col(position = "dodge") +
  facet_wrap(~ gene_set, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "Down" = "red",
      "Up" = "forestgreen"
    ),
    labels = c(
      "Down" = "Downregulated",
      "Up" = "Upregulated"
    )
  ) +
  theme_bw() +
  labs(
    title = "Differentially expressed genes after HDAC knockdown",
    subtitle = "Separate DESeq2 analyses combined by rows; threshold: padj < 0.05 and |log2FC| ≥ log2(1.5)",
    x = "Knockdown condition",
    y = "Number of DEGs",
    fill = "Direction"
  )

deg_count_plot_combined

ggsave(
  filename = file.path(combined_plot_dir, "02_DEG_counts_combined_HDAC1_2_3.png"),
  plot = deg_count_plot_combined,
  width = 8,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

## 3.2 Validate HDAC knockdown at the expression level
# Uses PCA-matched colours and DESeq2 statistics.

### 3.2.1 Calculate normalised counts
dds_pca <- estimateSizeFactors(dds_pca)

norm_counts_combined <- counts(
  dds_pca,
  normalized = TRUE
)

### 3.2.2 Identify HDAC target-gene IDs
target_gene_ids_combined <- all_res_hdac1_2_3_15fold %>%
  dplyr::filter(SYMBOL %in% c("HDAC1", "HDAC2", "HDAC3")) %>%
  dplyr::select(gene_id, SYMBOL) %>%
  dplyr::distinct()

target_gene_ids_combined

target_expr_combined <- norm_counts_combined[
  target_gene_ids_combined$gene_id,
  ,
  drop = FALSE
] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  tidyr::pivot_longer(
    cols = -gene_id,
    names_to = "sample",
    values_to = "norm_count"
  ) %>%
  dplyr::left_join(
    target_gene_ids_combined,
    by = "gene_id"
  ) %>%
  dplyr::left_join(
    meta_pca %>%
      tibble::rownames_to_column("sample") %>%
      dplyr::select(sample, experiment, pca_group, sample_label),
    by = "sample"
  ) %>%
  dplyr::mutate(
    comparison_group = dplyr::case_when(
      SYMBOL == "HDAC1" & pca_group == "HDAC12 Scramble" ~ "Scramble",
      SYMBOL == "HDAC1" & pca_group == "HDAC1 KD" ~ "Knockdown",

      SYMBOL == "HDAC2" & pca_group == "HDAC12 Scramble" ~ "Scramble",
      SYMBOL == "HDAC2" & pca_group == "HDAC2 KD" ~ "Knockdown",

      SYMBOL == "HDAC3" & pca_group == "HDAC3 Scramble" ~ "Scramble",
      SYMBOL == "HDAC3" & pca_group == "HDAC3 KD" ~ "Knockdown",

      TRUE ~ NA_character_
    ),

    # Separate colour groups so each knockdown retains its PCA colour.
    colour_group = dplyr::case_when(
      SYMBOL == "HDAC1" & comparison_group == "Scramble" ~ "HDAC1 Scramble",
      SYMBOL == "HDAC1" & comparison_group == "Knockdown" ~ "HDAC1 KD",

      SYMBOL == "HDAC2" & comparison_group == "Scramble" ~ "HDAC2 Scramble",
      SYMBOL == "HDAC2" & comparison_group == "Knockdown" ~ "HDAC2 KD",

      SYMBOL == "HDAC3" & comparison_group == "Scramble" ~ "HDAC3 Scramble",
      SYMBOL == "HDAC3" & comparison_group == "Knockdown" ~ "HDAC3 KD",

      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(comparison_group)) %>%
  dplyr::mutate(
    comparison_group = factor(
      comparison_group,
      levels = c("Scramble", "Knockdown")
    ),
    SYMBOL = factor(
      SYMBOL,
      levels = c("HDAC1", "HDAC2", "HDAC3")
    )
  )

target_validation_cols <- c(
  "HDAC1 Scramble" = "#D95F5F",
  "HDAC1 KD"       = "#00BA38",

  "HDAC2 Scramble" = "#D95F5F",
  "HDAC2 KD"       = "#619CFF",

  "HDAC3 Scramble" = "#D95F5F",
  "HDAC3 KD"       = "#00BFC4"
)

### 3.2.3 Extract target-validation statistics

target_validation_stats <- all_res_hdac1_2_3_15fold %>%
  dplyr::filter(
    SYMBOL %in% c("HDAC1", "HDAC2", "HDAC3")
  ) %>%
  dplyr::mutate(
    contrast_short = as.character(contrast_short)
  ) %>%
  dplyr::filter(
    (SYMBOL == "HDAC1" & contrast_short == "HDAC1 KD") |
      (SYMBOL == "HDAC2" & contrast_short == "HDAC2 KD") |
      (SYMBOL == "HDAC3" & contrast_short == "HDAC3 KD")
  ) %>%
  dplyr::mutate(
    SYMBOL = factor(SYMBOL, levels = c("HDAC1", "HDAC2", "HDAC3")),
    sig_label = dplyr::case_when(
      is.na(padj) ~ "NA",
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    stat_label = paste0(
      "log2FC=", round(log2FoldChange, 2),
      "\n", sig_label
    )
  ) %>%
  dplyr::select(SYMBOL, log2FoldChange, padj, sig_label, stat_label)

### 3.2.4 Calculate label positions
target_label_positions <- target_expr_combined %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(
    y_max = max(norm_count, na.rm = TRUE),
    y_min = min(norm_count, na.rm = TRUE),
    y_range = y_max - y_min,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    y_pos = y_max + 0.12 * y_range
  )

target_validation_stats <- target_validation_stats %>%
  dplyr::left_join(
    target_label_positions,
    by = "SYMBOL"
  )

### 3.2.5 Plot target expression

target_validation_boxplot_combined <- ggplot(
  target_expr_combined,
  aes(
    x = comparison_group,
    y = norm_count,
    fill = colour_group
  )
) +
  geom_boxplot(
    colour = "black",
    outlier.shape = NA,
    width = 0.55,
    alpha = 0.78
  ) +
  geom_jitter(
    aes(colour = colour_group),
    width = 0.12,
    size = 2.5,
    alpha = 0.9
  ) +
  geom_text(
    data = target_validation_stats,
    aes(
      x = 1.5,
      y = y_pos,
      label = stat_label
    ),
    inherit.aes = FALSE,
    size = 3.5,
    lineheight = 0.9
  ) +
  facet_wrap(
    ~ SYMBOL,
    scales = "free_y"
  ) +
  scale_fill_manual(values = target_validation_cols) +
  scale_colour_manual(values = target_validation_cols) +
  theme_bw() +
  labs(
    title = "HDAC target expression after siRNA knockdown",
    subtitle = "Replicate-level DESeq2-normalised counts using matched scramble controls; DESeq2 statistics shown for target knockdown",
    x = NULL,
    y = "Normalised count",
    fill = NULL,
    colour = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

target_validation_boxplot_combined

ggsave(
  filename = file.path(
    combined_plot_dir,
    "03_HDAC_target_validation_boxplot_combined_coloured_stats.png"
  ),
  plot = target_validation_boxplot_combined,
  width = 9.5,
  height = 5,
  dpi = 300,
  bg = "white"
)

# 4. X-linked differential-expression analyses

## 4.1 Plot overlap of X-linked DEGs

x_deg_table_combined <- all_res_hdac1_2_3_15fold %>%
  dplyr::filter(
    CHR == "X",
    deg_status_15fold != "Not significant"
  )

x_deg_sets_combined <- list(
  "HDAC1 KD" = x_deg_table_combined %>%
    dplyr::filter(contrast_short == "HDAC1 KD") %>%
    dplyr::pull(gene_id) %>%
    unique(),

  "HDAC2 KD" = x_deg_table_combined %>%
    dplyr::filter(contrast_short == "HDAC2 KD") %>%
    dplyr::pull(gene_id) %>%
    unique(),

  "HDAC3 KD" = x_deg_table_combined %>%
    dplyr::filter(contrast_short == "HDAC3 KD") %>%
    dplyr::pull(gene_id) %>%
    unique()
)

# Use PCA-matched colours for the knockdown groups.
venn_kd_cols <- c(
  "HDAC1 KD" = "#00BA38",  # green
  "HDAC2 KD" = "#619CFF",  # blue
  "HDAC3 KD" = "#00BFC4"   # cyan
)

venn_x_combined <- ggvenn(
  x_deg_sets_combined,
  fill_color = unname(venn_kd_cols[names(x_deg_sets_combined)]),
  fill_alpha = 0.45,
  stroke_size = 0.8,
  stroke_color = "black",
  set_name_color = unname(venn_kd_cols[names(x_deg_sets_combined)]),
  set_name_size = 5,
  text_size = 5,
  text_color = "black",
  show_percentage = FALSE
) +
  labs(
    title = "Overlap of X-linked DEGs after HDAC knockdown",
    subtitle = "Threshold: padj < 0.05 and |log2FC| ≥ log2(1.5)"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = margin(20, 40, 20, 40)
  )

venn_x_combined

ggsave(
  filename = file.path(combined_plot_dir, "05_Xlinked_DEG_venn_combined.png"),
  plot = venn_x_combined,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

## 4.2 Build the X-linked DEG overlap table

x_overlap_gene_list_combined <- x_deg_table_combined %>%
  dplyr::mutate(
    contrast_id = dplyr::case_when(
      contrast_short == "HDAC1 KD" ~ "HDAC1",
      contrast_short == "HDAC2 KD" ~ "HDAC2",
      contrast_short == "HDAC3 KD" ~ "HDAC3",
      TRUE ~ as.character(contrast_short)
    ),
    gene_plot_label = dplyr::case_when(
      !is.na(SYMBOL) & SYMBOL != "" ~ SYMBOL,
      !is.na(gene_label) & gene_label != "" ~ gene_label,
      TRUE ~ gene_id
    )
  ) %>%
  dplyr::select(
    gene_id,
    gene_plot_label,
    SYMBOL,
    GENE_BIOTYPE,
    contrast_id,
    deg_status_15fold,
    log2FoldChange,
    padj
  ) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(
    names_from = contrast_id,
    values_from = c(deg_status_15fold, log2FoldChange, padj),
    values_fill = list(
      deg_status_15fold = "Not significant",
      log2FoldChange = NA_real_,
      padj = NA_real_
    )
  ) %>%
  dplyr::mutate(
    in_HDAC1 = deg_status_15fold_HDAC1 != "Not significant",
    in_HDAC2 = deg_status_15fold_HDAC2 != "Not significant",
    in_HDAC3 = deg_status_15fold_HDAC3 != "Not significant",
    n_contrasts = in_HDAC1 + in_HDAC2 + in_HDAC3,
    overlap_group = dplyr::case_when(
      in_HDAC1 & in_HDAC2 & in_HDAC3 ~ "Shared by all 3",
      in_HDAC1 & in_HDAC2 & !in_HDAC3 ~ "HDAC1 + HDAC2 only",
      in_HDAC1 & !in_HDAC2 & in_HDAC3 ~ "HDAC1 + HDAC3 only",
      !in_HDAC1 & in_HDAC2 & in_HDAC3 ~ "HDAC2 + HDAC3 only",
      in_HDAC1 & !in_HDAC2 & !in_HDAC3 ~ "HDAC1 only",
      !in_HDAC1 & in_HDAC2 & !in_HDAC3 ~ "HDAC2 only",
      !in_HDAC1 & !in_HDAC2 & in_HDAC3 ~ "HDAC3 only",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_contrasts),
    overlap_group,
    gene_plot_label
  )

x_overlap_gene_list_combined %>%
  dplyr::count(overlap_group)

shared_x_genes_combined <- x_overlap_gene_list_combined %>%
  dplyr::filter(n_contrasts >= 2) %>%
  dplyr::select(
    overlap_group,
    gene_id,
    gene_plot_label,
    SYMBOL,
    GENE_BIOTYPE,
    deg_status_15fold_HDAC1,
    deg_status_15fold_HDAC2,
    deg_status_15fold_HDAC3,
    log2FoldChange_HDAC1,
    log2FoldChange_HDAC2,
    log2FoldChange_HDAC3,
    padj_HDAC1,
    padj_HDAC2,
    padj_HDAC3
  ) %>%
  dplyr::arrange(overlap_group, gene_plot_label)

shared_x_genes_combined

shared_all3_x_gene <- x_overlap_gene_list_combined %>%
  dplyr::filter(overlap_group == "Shared by all 3") %>%
  dplyr::select(
    gene_id,
    gene_plot_label,
    SYMBOL,
    GENE_BIOTYPE,
    deg_status_15fold_HDAC1,
    deg_status_15fold_HDAC2,
    deg_status_15fold_HDAC3,
    log2FoldChange_HDAC1,
    log2FoldChange_HDAC2,
    log2FoldChange_HDAC3,
    padj_HDAC1,
    padj_HDAC2,
    padj_HDAC3
  )

shared_all3_x_gene

## 4.3 Plot shared X-linked genes as a log2FC heatmap

shared_x_plot_df <- x_overlap_gene_list_combined %>%
  dplyr::filter(n_contrasts >= 2) %>%
  dplyr::select(
    gene_plot_label,
    overlap_group,
    log2FoldChange_HDAC1,
    log2FoldChange_HDAC2,
    log2FoldChange_HDAC3
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("log2FoldChange_"),
    names_to = "contrast",
    values_to = "log2FC"
  ) %>%
  dplyr::mutate(
    contrast = dplyr::recode(
      contrast,
      "log2FoldChange_HDAC1" = "HDAC1 KD",
      "log2FoldChange_HDAC2" = "HDAC2 KD",
      "log2FoldChange_HDAC3" = "HDAC3 KD"
    )
  )

# Order genes by overlap group and then mean absolute effect size.
gene_order_shared <- shared_x_plot_df %>%
  dplyr::group_by(gene_plot_label, overlap_group) %>%
  dplyr::summarise(mean_abs_log2FC = mean(abs(log2FC), na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(overlap_group, dplyr::desc(mean_abs_log2FC)) %>%
  dplyr::pull(gene_plot_label)

shared_x_plot_df$gene_plot_label <- factor(
  shared_x_plot_df$gene_plot_label,
  levels = rev(unique(gene_order_shared))
)

shared_x_heatmap <- ggplot(
  shared_x_plot_df,
  aes(
    x = contrast,
    y = gene_plot_label,
    fill = log2FC
  )
) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(
    aes(label = ifelse(is.na(log2FC), "", round(log2FC, 2))),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "grey90",
    name = "log2FC"
  ) +
  theme_bw() +
  labs(
    title = "Shared X-linked DEGs across HDAC knockdowns",
    subtitle = "padj < 0.05 and |log2FC| ≥ log2(1.5)",
    x = "Knockdown condition",
    y = "Shared X-linked genes"
  )

shared_x_heatmap

ggsave(
  filename = file.path(combined_plot_dir, "shared_Xlinked_DEGs_heatmap.png"),
  plot = shared_x_heatmap,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

## 4.4 Plot X-linked DEGs by chromosome position

xist_start_bp <- 73817774
xist_end_bp   <- 73852754
xist_mid_mb   <- mean(c(xist_start_bp, xist_end_bp)) / 1e6

x_position_plot_df_combined <- x_deg_table_combined %>%
  dplyr::filter(
    !is.na(START),
    !is.na(END),
    is.finite(START),
    is.finite(END),
    !is.na(log2FoldChange),
    !is.na(padj)
  ) %>%
  dplyr::mutate(
    gene_mid_mb = ((START + END) / 2) / 1e6,
    distance_to_xist_mb = abs(gene_mid_mb - xist_mid_mb),
    gene_plot_label = dplyr::case_when(
      !is.na(SYMBOL) & SYMBOL != "" ~ SYMBOL,
      !is.na(gene_label) & gene_label != "" ~ gene_label,
      TRUE ~ gene_id
    )
  )

closest_xic_genes_combined <- x_position_plot_df_combined %>%
  dplyr::filter(gene_plot_label != "XIST") %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::slice_min(distance_to_xist_mb, n = 5, with_ties = FALSE) %>%
  dplyr::ungroup()

xic_label_df_combined <- x_position_plot_df_combined %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::summarise(
    y_pos = max(log2FoldChange, na.rm = TRUE) + 0.35,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    x_pos = xist_mid_mb,
    label = "XIST/XIC"
  )

x_position_plot_combined <- ggplot(
  x_position_plot_df_combined,
  aes(
    x = gene_mid_mb,
    y = log2FoldChange,
    colour = deg_status_15fold
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = xist_mid_mb, linetype = "dashed", colour = "black") +
  geom_point(size = 2.3, alpha = 0.85) +
  ggrepel::geom_text_repel(
    data = closest_xic_genes_combined,
    aes(label = gene_plot_label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.2,
    show.legend = FALSE
  ) +
  geom_text(
    data = xic_label_df_combined,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    angle = 90,
    vjust = -0.2,
    size = 3
  ) +
  facet_wrap(~ contrast_short, ncol = 1) +
  scale_colour_manual(
    values = c(
      "Down" = "red",
      "Up" = "forestgreen"
    ),
    labels = c(
      "Down" = "Downregulated",
      "Up" = "Upregulated"
    )
  ) +
  theme_bw() +
  labs(
    title = "Position of changing X-linked genes relative to XIST/XIC",
    subtitle = "Dashed vertical line marks the XIST locus; nearest significant X-linked DEGs are labelled",
    x = "Position on chromosome X (Mb)",
    y = "Shrunken log2 fold change",
    colour = "Direction"
  )

x_position_plot_combined

ggsave(
  filename = file.path(combined_plot_dir, "06_Xlinked_gene_position_XIST_XIC_combined.png"),
  plot = x_position_plot_combined,
  width = 9,
  height = 9,
  dpi = 300,
  bg = "white"
)

## 4.5 Plot the top 20 X-linked DEGs per contrast

top20_x_genes_combined <- x_deg_table_combined %>%
  dplyr::mutate(
    abs_log2FC = abs(log2FoldChange),
    neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)),
    ranking_score = abs_log2FC * neg_log10_padj,
    gene_plot_label = dplyr::case_when(
      !is.na(SYMBOL) & SYMBOL != "" ~ SYMBOL,
      !is.na(gene_label) & gene_label != "" ~ gene_label,
      TRUE ~ gene_id
    )
  ) %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::slice_max(
    order_by = ranking_score,
    n = 20,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

top20_x_plot_combined <- ggplot(
  top20_x_genes_combined,
  aes(
    x = log2FoldChange,
    y = reorder(gene_plot_label, ranking_score),
    colour = deg_status_15fold,
    size = neg_log10_padj
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.9) +
  facet_wrap(~ contrast_short, scales = "free_y") +
  scale_colour_manual(
    values = c(
      "Down" = "red",
      "Up" = "forestgreen"
    ),
    labels = c(
      "Down" = "Downregulated",
      "Up" = "Upregulated"
    )
  ) +
  theme_bw() +
  labs(
    title = "Top X-linked DEGs after HDAC knockdown",
    subtitle = "Ranked by -log10(padj) × absolute shrunken log2FC",
    x = "Shrunken log2 fold change",
    y = "Gene",
    colour = "Direction",
    size = "-log10(padj)"
  )

top20_x_plot_combined

ggsave(
  filename = file.path(combined_plot_dir, "07_Top20_Xlinked_DEGs_combined.png"),
  plot = top20_x_plot_combined,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

## 4.6 Plot the directional X-linked DEG heatmap

x_sig_genes_combined <- x_deg_table_combined %>%
  dplyr::pull(gene_id) %>%
  unique()

x_direction_heatmap_df_combined <- all_res_hdac1_2_3_15fold %>%
  dplyr::filter(
    gene_id %in% x_sig_genes_combined,
    contrast_short %in% c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
  ) %>%
  dplyr::mutate(
    gene_plot_label = dplyr::case_when(
      !is.na(SYMBOL) & SYMBOL != "" ~ SYMBOL,
      !is.na(gene_label) & gene_label != "" ~ gene_label,
      TRUE ~ gene_id
    ),
    is_significant = deg_status_15fold != "Not significant",
    neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)),
    ranking_score = abs(log2FoldChange) * neg_log10_padj
  )

# Retain the top 60 X-linked genes by maximum ranking score.
top_x_heatmap_genes <- x_direction_heatmap_df_combined %>%
  dplyr::group_by(gene_id, gene_plot_label) %>%
  dplyr::summarise(
    max_score = max(ranking_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::slice_max(max_score, n = 60, with_ties = FALSE) %>%
  dplyr::pull(gene_id)

x_direction_heatmap_top_df <- x_direction_heatmap_df_combined %>%
  dplyr::filter(gene_id %in% top_x_heatmap_genes)

gene_order_heatmap <- x_direction_heatmap_top_df %>%
  dplyr::group_by(gene_id, gene_plot_label) %>%
  dplyr::summarise(
    max_abs_lfc = max(abs(log2FoldChange), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(max_abs_lfc)) %>%
  dplyr::pull(gene_plot_label)

x_direction_heatmap_top_df <- x_direction_heatmap_top_df %>%
  dplyr::mutate(
    gene_plot_label = factor(gene_plot_label, levels = rev(gene_order_heatmap)),
    contrast_short = factor(contrast_short, levels = c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD"))
  )

x_direction_heatmap_plot_combined <- ggplot(
  x_direction_heatmap_top_df,
  aes(
    x = contrast_short,
    y = gene_plot_label,
    fill = log2FoldChange
  )
) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_point(
    data = x_direction_heatmap_top_df %>%
      dplyr::filter(is_significant),
    aes(x = contrast_short, y = gene_plot_label),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.6,
    fill = "black",
    colour = "black"
  ) +
  scale_fill_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0
  ) +
  theme_bw() +
  labs(
    title = "Directional changes in X-linked DEGs after HDAC knockdown",
    subtitle = "Top 60 X-linked DEGs; tile colour shows shrunken log2FC and black dots mark significant genes",
    x = NULL,
    y = "X-linked genes",
    fill = "Shrunken\nlog2FC"
  ) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 11)
  )

x_direction_heatmap_plot_combined

ggsave(
  filename = file.path(combined_plot_dir, "08_Xlinked_DEG_direction_heatmap_combined.png"),
  plot = x_direction_heatmap_plot_combined,
  width = 7,
  height = 11,
  dpi = 300,
  bg = "white"
)

# 5. Genome-wide volcano, MA, and overlap plots

dir.create(combined_plot_dir, showWarnings = FALSE, recursive = TRUE)

## 5.1 Prepare the genome-wide plotting table
plot_res_combined <- all_res_hdac1_2_3_15fold %>%
  dplyr::mutate(
    contrast_short = factor(
      contrast_short,
      levels = c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
    ),
    gene_label = dplyr::if_else(
      !is.na(SYMBOL) & SYMBOL != "",
      SYMBOL,
      gene_id
    ),
    neg_log10_padj = -log10(padj)
  )

### 5.1.1 Select the mean-expression measure
if ("baseMean" %in% colnames(plot_res_combined)) {
  plot_res_combined <- plot_res_combined %>%
    dplyr::mutate(mean_expr_for_ma = baseMean)
} else if (all(c("control_mean_norm_count", "treatment_mean_norm_count") %in% colnames(plot_res_combined))) {
  plot_res_combined <- plot_res_combined %>%
    dplyr::mutate(
      mean_expr_for_ma = rowMeans(
        cbind(control_mean_norm_count, treatment_mean_norm_count),
        na.rm = TRUE
      )
    )
} else {
  stop("Could not find baseMean or generic mean-count columns for the MA plot.")
}

### 5.1.2 Select HDAC and XIST labels
named_genes_combined <- plot_res_combined %>%
  dplyr::filter(
    (contrast_short == "HDAC1 KD" & SYMBOL %in% c("HDAC1", "XIST")) |
      (contrast_short == "HDAC2 KD" & SYMBOL %in% c("HDAC2", "XIST")) |
      (contrast_short == "HDAC3 KD" & SYMBOL %in% c("HDAC3", "XIST"))
  ) %>%
  dplyr::mutate(label_source = "HDAC/XIST") %>%
  dplyr::select(contrast_short, gene_id, gene_label, label_source)

### 5.1.3 Select the top eight DEGs per contrast
# Ranked by the smallest adjusted p-value among significant genes.
top8_genes_combined <- plot_res_combined %>%
  dplyr::filter(deg_status_15fold != "Not significant", !is.na(padj)) %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange)), .by_group = TRUE) %>%
  dplyr::slice_head(n = 8) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(label_source = "Top 8 DEG") %>%
  dplyr::select(contrast_short, gene_id, gene_label, label_source)

### 5.1.4 Combine the plot-label table
label_table_combined <- dplyr::bind_rows(
  named_genes_combined,
  top8_genes_combined
) %>%
  dplyr::distinct(contrast_short, gene_id, .keep_all = TRUE)

### 5.1.5 Add plot highlighting categories
plot_res_combined <- plot_res_combined %>%
  dplyr::left_join(
    label_table_combined,
    by = c("contrast_short", "gene_id", "gene_label")
  ) %>%
  dplyr::mutate(
    point_group = dplyr::case_when(
      label_source == "HDAC/XIST" ~ "HDAC/XIST",
      label_source == "Top 8 DEG" ~ "Top 8 DEG",
      deg_status_15fold != "Not significant" ~ "Other DEG",
      TRUE ~ "Not significant"
    ),
    point_group = factor(
      point_group,
      levels = c("Not significant", "Other DEG", "Top 8 DEG", "HDAC/XIST")
    )
  )

### 5.1.6 Export labelled genes
label_table_combined_export <- plot_res_combined %>%
  dplyr::filter(!is.na(label_source)) %>%
  dplyr::select(
    contrast_short,
    gene_id,
    SYMBOL,
    gene_label,
    log2FoldChange,
    padj,
    deg_status_15fold,
    label_source
  ) %>%
  dplyr::arrange(contrast_short, label_source, padj)

write.csv(
  label_table_combined_export,
  file.path(combined_plot_dir, "plot_labelled_genes_HDAC123.csv"),
  row.names = FALSE
)

## 5.2 Plot genome-wide differential expression

volcano_plot_combined <- ggplot(
  plot_res_combined %>% dplyr::filter(!is.na(padj)),
  aes(x = log2FoldChange, y = neg_log10_padj)
) +
  geom_point(
    data = plot_res_combined %>%
      dplyr::filter(!is.na(padj), point_group == "Not significant"),
    colour = "grey75",
    alpha = 0.6,
    size = 1.2
  ) +
  geom_point(
    data = plot_res_combined %>%
      dplyr::filter(!is.na(padj), point_group == "Other DEG"),
    colour = "#2C7FB8",
    alpha = 0.8,
    size = 1.4
  ) +
  geom_point(
    data = plot_res_combined %>%
      dplyr::filter(!is.na(padj), point_group == "Top 8 DEG"),
    colour = "#D95F02",
    alpha = 0.95,
    size = 2
  ) +
  geom_point(
    data = plot_res_combined %>%
      dplyr::filter(!is.na(padj), point_group == "HDAC/XIST"),
    colour = "black",
    alpha = 1,
    size = 2.4
  ) +
  ggrepel::geom_text_repel(
    data = plot_res_combined %>%
      dplyr::filter(!is.na(padj), !is.na(label_source)),
    aes(label = gene_label),
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  geom_vline(
    xintercept = c(-log2(1.5), log2(1.5)),
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    colour = "grey40"
  ) +
  facet_wrap(~ contrast_short, nrow = 1) +
  theme_bw() +
  labs(
    title = "Genome-wide volcano plots for HDAC knockdown",
    subtitle = "Thresholds: padj < 0.05 and |log2FC| ≥ log2(1.5); HDAC/XIST and top 8 DEGs labelled",
    x = "log2 fold change",
    y = expression(-log[10]("adjusted p-value"))
  )

volcano_plot_combined

ggsave(
  file.path(combined_plot_dir, "volcano_HDAC123_genomewide.png"),
  volcano_plot_combined,
  width = 12,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

## 5.3 Prepare the MA-plot data

ma_plot_df_combined <- all_res_hdac1_2_3_15fold %>%
  dplyr::mutate(
    contrast_short = factor(
      contrast_short,
      levels = c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
    ),
    gene_label = dplyr::if_else(
      !is.na(SYMBOL) & SYMBOL != "",
      SYMBOL,
      gene_id
    )
  )

### 5.3.1 Calculate mean expression
if ("baseMean" %in% colnames(ma_plot_df_combined)) {
  ma_plot_df_combined <- ma_plot_df_combined %>%
    dplyr::mutate(
      mean_expr_for_ma = baseMean,
      mean_expr_log10 = log10(mean_expr_for_ma + 1)
    )
} else {
  stop("baseMean column not found. Check colnames(all_res_hdac1_2_3_15fold).")
}

### 5.3.2 Select the top eight DEGs
top8_ma_genes <- ma_plot_df_combined %>%
  dplyr::filter(
    deg_status_15fold != "Not significant",
    !is.na(padj)
  ) %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange)), .by_group = TRUE) %>%
  dplyr::slice_head(n = 8) %>%
  dplyr::ungroup() %>%
  dplyr::select(contrast_short, gene_id) %>%
  dplyr::mutate(is_top8 = TRUE)

### 5.3.3 Define plotting categories
ma_plot_df_combined <- ma_plot_df_combined %>%
  dplyr::left_join(
    top8_ma_genes,
    by = c("contrast_short", "gene_id")
  ) %>%
  dplyr::mutate(
    is_top8 = dplyr::if_else(is.na(is_top8), FALSE, is_top8),
    point_group = dplyr::case_when(
      SYMBOL %in% c("HDAC1", "HDAC2", "HDAC3", "XIST") ~ "HDAC / XIST",
      is_top8 ~ "Top 8 DEG",
      deg_status_15fold == "Up" ~ "Upregulated",
      deg_status_15fold == "Down" ~ "Downregulated",
      TRUE ~ "Non-Significant"
    ),
    point_group = factor(
      point_group,
      levels = c(
        "Non-Significant",
        "Downregulated",
        "Upregulated",
        "Top 8 DEG",
        "HDAC / XIST"
      )
    )
  ) %>%
  dplyr::arrange(point_group)

### 5.3.4 Select MA-plot labels
ma_label_df_combined <- ma_plot_df_combined %>%
  dplyr::filter(
    point_group %in% c("HDAC / XIST", "Top 8 DEG")
  )

## 5.4 Plot genome-wide expression changes

ma_plot_combined_coloured <- ggplot(
  ma_plot_df_combined,
  aes(
    x = mean_expr_log10,
    y = log2FoldChange,
    colour = point_group
  )
) +
  geom_hline(
    yintercept = c(-log2(1.5), log2(1.5)),
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_point(
    size = 1.8,
    alpha = 0.8
  ) +
  ggrepel::geom_text_repel(
    data = ma_label_df_combined,
    aes(label = gene_label),
    colour = "black",
    size = 3.5,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ contrast_short,
    nrow = 1
  ) +
  scale_colour_manual(
    values = c(
      "Non-Significant" = "grey75",
      "Downregulated" = "red3",
      "Upregulated" = "forestgreen",
      "Top 8 DEG" = "dodgerblue3",
      "HDAC / XIST" = "black"
    )
  ) +
  theme_bw() +
  labs(
    title = "Genome-wide MA plots for HDAC knockdown",
    subtitle = "Upregulated genes in green, downregulated genes in red, top 8 DEGs in blue and HDAC/XIST genes in black",
    x = expression(log[10]("mean normalised expression + 1")),
    y = "log2 fold change",
    colour = "Gene group"
  )

ma_plot_combined_coloured

ggsave(
  file.path(
    combined_plot_dir,
    "MA_HDAC123_genomewide.png"
  ),
  ma_plot_combined_coloured,
  width = 12,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

## 5.5 Plot genome-wide DEG overlap

# Use the same HDAC colours as the combined PCA.
hdac_venn_cols <- c(
  "HDAC1 KD" = "#00BA38",  # green
  "HDAC2 KD" = "#619CFF",  # blue
  "HDAC3 KD" = "#00BFC4"   # cyan
)

genomewide_deg_sets_combined <- list(
  "HDAC1 KD" = plot_res_combined %>%
    dplyr::filter(
      contrast_short == "HDAC1 KD",
      deg_status_15fold != "Not significant"
    ) %>%
    dplyr::pull(gene_id) %>%
    unique(),

  "HDAC2 KD" = plot_res_combined %>%
    dplyr::filter(
      contrast_short == "HDAC2 KD",
      deg_status_15fold != "Not significant"
    ) %>%
    dplyr::pull(gene_id) %>%
    unique(),

  "HDAC3 KD" = plot_res_combined %>%
    dplyr::filter(
      contrast_short == "HDAC3 KD",
      deg_status_15fold != "Not significant"
    ) %>%
    dplyr::pull(gene_id) %>%
    unique()
)

venn_genomewide_combined <- ggvenn(
  genomewide_deg_sets_combined,
  fill_color = unname(hdac_venn_cols[names(genomewide_deg_sets_combined)]),
  fill_alpha = 0.55,
  stroke_color = "black",
  stroke_size = 0.6,
  set_name_color = "black",
  set_name_size = 5,
  text_color = "black",
  text_size = 5,
  show_percentage = FALSE
) +
  labs(
    title = "Overlap of genome-wide DEGs after HDAC knockdown",
    subtitle = "Threshold: padj < 0.05 and |log2FC| ≥ log2(1.5)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11),
    plot.margin = margin(20, 40, 20, 40)
  )

venn_genomewide_combined

ggsave(
  file.path(combined_plot_dir, "Venn_genomewide.png"),
  venn_genomewide_combined,
  width = 12,
  height = 5.5,
  dpi = 300,
  bg = "white"
)
colnames(meta_pca)
meta_pca

# 6. HDAC3-only PCA

## 6.1 Subset HDAC3 metadata
meta_hdac3_pca <- meta_pca %>%
  dplyr::filter(
    experiment == "HDAC3"
  ) %>%
  dplyr::mutate(
    unique_sample = as.character(unique_sample),
    condition = factor(
      condition,
      levels = c("Scramble", "HDAC3_siRNA")
    ),
    pca_group = factor(
      pca_group,
      levels = c("HDAC3 Scramble", "HDAC3 KD")
    )
  )

### 6.1.1 Check sample identifiers
meta_hdac3_pca$unique_sample
colnames(counts_pca)

## 6.2 Subset the HDAC3 count matrix
counts_hdac3_pca <- counts_pca[
  ,
  meta_hdac3_pca$unique_sample,
  drop = FALSE
]

## 6.3 Match metadata to the count matrix
meta_hdac3_pca <- meta_hdac3_pca %>%
  dplyr::slice(
    match(colnames(counts_hdac3_pca), unique_sample)
  )

meta_hdac3_pca <- as.data.frame(meta_hdac3_pca)
rownames(meta_hdac3_pca) <- meta_hdac3_pca$unique_sample

### 6.3.1 Check dimensions and sample order
dim(counts_hdac3_pca)
dim(meta_hdac3_pca)

stopifnot(ncol(counts_hdac3_pca) == nrow(meta_hdac3_pca))
stopifnot(all(colnames(counts_hdac3_pca) == rownames(meta_hdac3_pca)))

## 6.4 Build the DESeq2 object and run VST

dds_hdac3_pca <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_hdac3_pca,
  colData = meta_hdac3_pca,
  design = ~ condition
)

### 6.4.1 Filter low-count genes
keep_hdac3_pca <- rowSums(DESeq2::counts(dds_hdac3_pca)) >= 10
dds_hdac3_pca <- dds_hdac3_pca[keep_hdac3_pca, ]

### 6.4.2 Apply the variance-stabilising transformation
vsd_hdac3_pca <- DESeq2::vst(
  dds_hdac3_pca,
  blind = TRUE
)

### 6.4.3 Extract PCA data
pca_hdac3_df <- DESeq2::plotPCA(
  vsd_hdac3_pca,
  intgroup = "pca_group",
  returnData = TRUE
)

percentVar_hdac3 <- round(100 * attr(pca_hdac3_df, "percentVar"))

## 6.5 Plot the HDAC3-only PCA

hdac3_pca_plot <- ggplot(
  pca_hdac3_df,
  aes(
    x = PC1,
    y = PC2,
    colour = pca_group,
    label = name
  )
) +
  geom_point(
    size = 3.5
  ) +
  ggrepel::geom_text_repel(
    size = 3,
    max.overlaps = Inf
  ) +
  theme_bw() +
  labs(
    title = "PCA of HDAC3 RNA-seq samples",
    subtitle = "VST-transformed raw counts; DESeq2 normalisation",
    x = paste0("PC1: ", percentVar_hdac3[1], "% variance"),
    y = paste0("PC2: ", percentVar_hdac3[2], "% variance"),
    colour = "Condition"
  )

hdac3_pca_plot

ggsave(
  filename = file.path(
    combined_plot_dir,
    "PCA_HDAC3_only.png"
  ),
  plot = hdac3_pca_plot,
  width = 7,
  height = 6,
  dpi = 300,
  bg = "white"
)

# 7. Keshet-highlighted pluripotency and differentiation genes
## 7.1 Normalised-expression heatmap after HDAC knockdown

keshet_hdac_plot_dir <- file.path(
  combined_plot_dir,
  "Keshet_pluripotency_gene_heatmap_HDAC"
)

dir.create(
  keshet_hdac_plot_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


### 7.1.1 Define the focused Keshet gene panel

keshet_gene_panel <- tibble::tribble(
  ~SYMBOL,    ~role_group,                         ~role_note,
  "POU5F1",   "ESC / pluripotency regulators",      "Core pluripotency",
  "SOX2",     "ESC / pluripotency regulators",      "Core pluripotency",
  "NANOG",    "ESC / pluripotency regulators",      "Core pluripotency",
  "PRDM14",   "ESC / pluripotency regulators",      "Pluripotency regulator",
  "SALL4",    "ESC / pluripotency regulators",      "Pluripotency maintenance",
  "ETV4",     "ESC / pluripotency regulators",      "Ambiguous pluripotency regulatory effect",
  "ZNF281",   "ESC / pluripotency regulators",      "Ambiguous pluripotency regulatory effect",

  "SMARCC1",  "Pluripotency attenuators",           "Attenuator / formative-state related",
  "TGIF1",    "Pluripotency attenuators",           "Attenuator of primed pluripotency GRN",
  "MIS18BP1", "Pluripotency attenuators",           "Putative attenuator",

  "ZNF462",   "Lineage-associated regulators",      "Mesendoderm-associated regulator",
  "REST",     "Lineage-associated regulators",      "Neural-associated regulator"
)

### 7.1.2 Calculate DESeq2-normalised counts

dds_pca <- DESeq2::estimateSizeFactors(dds_pca)

norm_counts_combined <- DESeq2::counts(
  dds_pca,
  normalized = TRUE
)

gene_annot_for_keshet_hdac <- all_res_hdac1_2_3_15fold %>%
  dplyr::select(gene_id, SYMBOL) %>%
  dplyr::distinct() %>%
  dplyr::filter(SYMBOL %in% keshet_gene_panel$SYMBOL)

missing_keshet_hdac_genes <- setdiff(
  keshet_gene_panel$SYMBOL,
  gene_annot_for_keshet_hdac$SYMBOL
)

missing_keshet_hdac_genes

### 7.1.3 Prepare metadata and matched comparisons

sample_meta_hdac <- meta_pca %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::select(sample, pca_group)

contrast_sample_map <- tibble::tribble(
  ~contrast_short, ~scramble_group,   ~kd_group,
  "HDAC1 KD",      "HDAC12 Scramble", "HDAC1 KD",
  "HDAC2 KD",      "HDAC12 Scramble", "HDAC2 KD",
  "HDAC3 KD",      "HDAC3 Scramble",  "HDAC3 KD"
)

### 7.1.4 Build the expression table

keshet_hdac_expr_long <- norm_counts_combined[
  gene_annot_for_keshet_hdac$gene_id,
  ,
  drop = FALSE
] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  tidyr::pivot_longer(
    cols = -gene_id,
    names_to = "sample",
    values_to = "norm_count"
  ) %>%
  dplyr::left_join(
    gene_annot_for_keshet_hdac,
    by = "gene_id"
  ) %>%
  dplyr::left_join(
    sample_meta_hdac,
    by = "sample"
  ) %>%
  dplyr::left_join(
    keshet_gene_panel,
    by = "SYMBOL"
  )

keshet_hdac_expr_mean <- keshet_hdac_expr_long %>%
  dplyr::group_by(SYMBOL, role_group, role_note, pca_group) %>%
  dplyr::summarise(
    mean_norm_count = mean(norm_count, na.rm = TRUE),
    mean_log2_count = log2(mean_norm_count + 1),
    .groups = "drop"
  )

keshet_hdac_heatmap_df <- contrast_sample_map %>%
  tidyr::pivot_longer(
    cols = c(scramble_group, kd_group),
    names_to = "comparison_type",
    values_to = "pca_group"
  ) %>%
  dplyr::mutate(
    condition_plot = dplyr::case_when(
      comparison_type == "scramble_group" ~ "Scramble",
      comparison_type == "kd_group" ~ "Knockdown"
    ),
    condition_plot = factor(
      condition_plot,
      levels = c("Scramble", "Knockdown")
    ),
    contrast_short = factor(
      contrast_short,
      levels = c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
    )
  ) %>%
  dplyr::left_join(
    keshet_hdac_expr_mean,
    by = "pca_group"
  ) %>%
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::group_by(contrast_short, SYMBOL) %>%
  dplyr::mutate(
    row_scaled_expr = as.numeric(scale(mean_log2_count)),
    row_scaled_expr = dplyr::if_else(
      is.na(row_scaled_expr),
      0,
      row_scaled_expr
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    role_group = factor(
      role_group,
      levels = c(
        "ESC / pluripotency regulators",
        "Pluripotency attenuators",
        "Lineage-associated regulators"
      )
    ),
    SYMBOL = factor(
      SYMBOL,
      levels = rev(keshet_gene_panel$SYMBOL)
    )
  )

### 7.1.5 Plot the normalised-expression heatmap

keshet_hdac_heatmap <- ggplot(
  keshet_hdac_heatmap_df,
  aes(
    x = condition_plot,
    y = SYMBOL,
    fill = row_scaled_expr
  )
) +
  geom_tile(
    colour = "grey85",
    linewidth = 0.3
  ) +
  facet_grid(
    role_group ~ contrast_short,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish
  ) +
  theme_bw() +
  labs(
    title = "Keshet-highlighted pluripotency and differentiation regulators after HDAC knockdown",
    subtitle = "Scramble vs knockdown comparisons; row-scaled log2 DESeq2-normalised counts",
    x = NULL,
    y = NULL,
    fill = "Row-scaled\nlog2 count"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

keshet_hdac_heatmap

ggsave(
  filename = file.path(
    keshet_hdac_plot_dir,
    "Keshet_pluripotency_differentiation_regulators_HDAC_heatmap.png"
  ),
  plot = keshet_hdac_heatmap,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)

readr::write_csv(
  keshet_hdac_log2fc_df,
  file.path(
    keshet_hdac_plot_dir,
    "Keshet_pluripotency_differentiation_regulators_HDAC_log2FC_heatmap_table.csv"
  )
)

## 7.2 DESeq2 log2FC heatmap after HDAC knockdown

keshet_hdac_plot_dir <- file.path(
  combined_plot_dir,
  "Keshet_pluripotency_gene_heatmap_HDAC"
)

dir.create(
  keshet_hdac_plot_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

### 7.2.1 Define the focused Keshet gene panel

keshet_gene_panel <- tibble::tribble(
  ~SYMBOL,    ~role_group,
  "POU5F1",   "ESC / pluripotency regulators",
  "SOX2",     "ESC / pluripotency regulators",
  "NANOG",    "ESC / pluripotency regulators",
  "PRDM14",   "ESC / pluripotency regulators",
  "SALL4",    "ESC / pluripotency regulators",
  "ETV4",     "ESC / pluripotency regulators",
  "ZNF281",   "ESC / pluripotency regulators",

  "SMARCC1",  "Pluripotency attenuators",
  "TGIF1",    "Pluripotency attenuators",
  "MIS18BP1", "Pluripotency attenuators",

  "ZNF462",   "Lineage-associated regulators",
  "REST",     "Lineage-associated regulators"
)

### 7.2.2 Prepare the DESeq2 log2FC table

keshet_hdac_log2fc_df <- all_res_hdac1_2_3_15fold %>%
  dplyr::filter(
    SYMBOL %in% keshet_gene_panel$SYMBOL,
    contrast_short %in% c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
  ) %>%
  dplyr::left_join(
    keshet_gene_panel,
    by = "SYMBOL"
  ) %>%
  dplyr::mutate(
    contrast_short = factor(
      contrast_short,
      levels = c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
    ),
    role_group = factor(
      role_group,
      levels = c(
        "ESC / pluripotency regulators",
        "Pluripotency attenuators",
        "Lineage-associated regulators"
      )
    ),
    SYMBOL = factor(
      SYMBOL,
      levels = rev(keshet_gene_panel$SYMBOL)
    ),
    is_sig = padj < 0.05 & abs(log2FoldChange) >= log2(1.5),
    log2FC_plot = pmax(
      pmin(log2FoldChange, 0.75),
      -0.75
    )
  )

# Check whether any selected genes are absent from the DESeq2 results.
missing_keshet_hdac_genes <- setdiff(
  keshet_gene_panel$SYMBOL,
  unique(keshet_hdac_log2fc_df$SYMBOL)
)

missing_keshet_hdac_genes

# Review the DESeq2 values for the selected genes.
keshet_hdac_log2fc_df %>%
  dplyr::select(
    SYMBOL,
    role_group,
    contrast_short,
    log2FoldChange,
    padj,
    deg_status_15fold,
    is_sig
  ) %>%
  dplyr::arrange(role_group, SYMBOL, contrast_short)

### 7.2.3 Plot the log2FC heatmap

keshet_hdac_log2fc_heatmap <- ggplot(
  keshet_hdac_log2fc_df,
  aes(
    x = contrast_short,
    y = SYMBOL,
    fill = log2FC_plot
  )
) +
  geom_tile(
    colour = "grey85",
    linewidth = 0.3
  ) +
  geom_point(
    data = keshet_hdac_log2fc_df %>%
      dplyr::filter(is_sig),
    aes(
      x = contrast_short,
      y = SYMBOL
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.8,
    fill = "black",
    colour = "black"
  ) +
  facet_grid(
    role_group ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish
  ) +
  theme_bw() +
  labs(
    title = "Keshet-highlighted pluripotency and differentiation regulators after HDAC knockdown",
    subtitle = "DESeq2 log2 fold change relative to matched Scramble controls; black dots indicate significant DEGs",
    x = "HDAC perturbation",
    y = NULL,
    fill = "log2FC"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text.y = element_text(face = "bold", angle = 0),
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.6, "lines")
  )

keshet_hdac_log2fc_heatmap

ggsave(
  filename = file.path(
    keshet_hdac_plot_dir,
    "Keshet_pluripotency_differentiation_regulators_HDAC_log2FC_heatmap.png"
  ),
  plot = keshet_hdac_log2fc_heatmap,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

readr::write_csv(
  keshet_hdac_log2fc_df,
  file.path(
    keshet_hdac_plot_dir,
    "Keshet_pluripotency_differentiation_regulators_HDAC_log2FC_heatmap_table.csv"
  )
)

# 8. Naive and primed marker analyses after HDAC knockdown

# The required packages are loaded again so this block can run independently.
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(readr)

## 8.1 Create the output folder

hdac_marker_plot_dir <- file.path(
  combined_plot_dir,
  "Naive_primed_marker_HDAC_boxplots"
)

dir.create(
  hdac_marker_plot_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

## 8.2 Define naive and primed marker genes

marker_genes <- tibble::tribble(
  ~SYMBOL,   ~marker_state,
  "DPPA5",   "Naive marker",
  "DNMT3L",  "Naive marker",
  "TFCP2L1", "Naive marker",
  "KLF5",    "Naive marker",

  "DNMT3B",  "Primed marker",
  "SALL2",   "Primed marker",
  "ZIC2",    "Primed marker",
  "PODXL",   "Primed marker"
)

## 8.3 Extract DESeq2-normalised counts

if (!exists("norm_counts_combined")) {
  norm_counts_combined <- DESeq2::counts(
    dds_pca,
    normalized = TRUE
  )
}

sample_meta_hdac <- meta_pca %>%
  as.data.frame()

if (!"sample" %in% colnames(sample_meta_hdac)) {
  sample_meta_hdac <- sample_meta_hdac %>%
    tibble::rownames_to_column("sample")
}

sample_meta_hdac <- sample_meta_hdac %>%
  dplyr::select(sample, pca_group)

## 8.4 Define matched HDAC comparisons

hdac_comparison_map <- tibble::tribble(
  ~contrast_short, ~pca_group,        ~hdac_comparison, ~condition_plot,
  "HDAC1 KD",      "HDAC12 Scramble", "HDAC1",          "Scramble",
  "HDAC1 KD",      "HDAC1 KD",        "HDAC1",          "KD",
  "HDAC2 KD",      "HDAC12 Scramble", "HDAC2",          "Scramble",
  "HDAC2 KD",      "HDAC2 KD",        "HDAC2",          "KD",
  "HDAC3 KD",      "HDAC3 Scramble",  "HDAC3",          "Scramble",
  "HDAC3 KD",      "HDAC3 KD",        "HDAC3",          "KD"
)

## 8.5 Build the marker-expression table

marker_gene_ids <- all_res_hdac1_2_3_15fold %>%
  dplyr::filter(SYMBOL %in% marker_genes$SYMBOL) %>%
  dplyr::group_by(SYMBOL, gene_id) %>%
  dplyr::summarise(
    max_baseMean = max(baseMean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::slice_max(
    order_by = max_baseMean,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

missing_marker_genes <- setdiff(
  marker_genes$SYMBOL,
  marker_gene_ids$SYMBOL
)

missing_marker_genes

hdac_marker_expr_df <- norm_counts_combined[
  marker_gene_ids$gene_id,
  ,
  drop = FALSE
] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  tidyr::pivot_longer(
    cols = -gene_id,
    names_to = "sample",
    values_to = "norm_count"
  ) %>%
  dplyr::left_join(
    marker_gene_ids %>% dplyr::select(gene_id, SYMBOL),
    by = "gene_id"
  ) %>%
  dplyr::left_join(
    marker_genes,
    by = "SYMBOL"
  ) %>%
  dplyr::left_join(
    sample_meta_hdac,
    by = "sample"
  ) %>%
  dplyr::left_join(
    hdac_comparison_map,
    by = "pca_group"
  ) %>%
  dplyr::filter(!is.na(contrast_short)) %>%
  dplyr::mutate(
    log2_norm_count = log2(norm_count + 1),
    hdac_comparison = factor(
      hdac_comparison,
      levels = c("HDAC1", "HDAC2", "HDAC3")
    ),
    condition_plot = factor(
      condition_plot,
      levels = c("Scramble", "KD")
    ),
    marker_state = factor(
      marker_state,
      levels = c("Naive marker", "Primed marker")
    ),
    SYMBOL = factor(
      SYMBOL,
      levels = marker_genes$SYMBOL
    )
  )

## 8.6 Extract marker-gene DESeq2 statistics

hdac_marker_stats <- all_res_hdac1_2_3_15fold %>%
  dplyr::filter(
    SYMBOL %in% marker_genes$SYMBOL,
    contrast_short %in% c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
  ) %>%
  dplyr::left_join(
    marker_genes,
    by = "SYMBOL"
  ) %>%
  dplyr::mutate(
    hdac_comparison = dplyr::case_when(
      contrast_short == "HDAC1 KD" ~ "HDAC1",
      contrast_short == "HDAC2 KD" ~ "HDAC2",
      contrast_short == "HDAC3 KD" ~ "HDAC3"
    ),

    hdac_comparison = factor(
      hdac_comparison,
      levels = c("HDAC1", "HDAC2", "HDAC3")
    ),

    marker_state = factor(
      marker_state,
      levels = c("Naive marker", "Primed marker")
    ),

    SYMBOL = factor(
      SYMBOL,
      levels = marker_genes$SYMBOL
    ),

    full_deg = !is.na(padj) &
      padj < 0.05 &
      abs(log2FoldChange) >= log2(1.5),

    sig_label = dplyr::case_when(
      full_deg & padj < 0.001 ~ "***",
      full_deg & padj < 0.01  ~ "**",
      full_deg & padj < 0.05  ~ "*",
      !full_deg & !is.na(padj) & padj < 0.05 ~ "†",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(
    SYMBOL,
    marker_state,
    contrast_short,
    hdac_comparison,
    log2FoldChange,
    padj,
    deg_status_15fold,
    full_deg,
    sig_label
  )

hdac_marker_stats

## 8.7 Calculate label positions

hdac_marker_label_pos <- hdac_marker_expr_df %>%
  dplyr::group_by(SYMBOL, marker_state, hdac_comparison) %>%
  dplyr::summarise(
    y_max = max(log2_norm_count, na.rm = TRUE),
    y_min = min(log2_norm_count, na.rm = TRUE),
    y_range = y_max - y_min,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    y_pos = y_max + dplyr::if_else(
      y_range == 0,
      0.25,
      0.18 * y_range
    )
  )

hdac_marker_stats_plot <- hdac_marker_stats %>%
  dplyr::left_join(
    hdac_marker_label_pos,
    by = c("SYMBOL", "marker_state", "hdac_comparison")
  )

## 8.7.5 Create facet labels for HDAC marker plots

marker_facet_levels <- c(
  "DPPA5\nNaive marker",
  "DNMT3L\nNaive marker",
  "TFCP2L1\nNaive marker",
  "KLF5\nNaive marker",
  "DNMT3B\nPrimed marker",
  "SALL2\nPrimed marker",
  "ZIC2\nPrimed marker",
  "PODXL\nPrimed marker"
)

hdac_marker_expr_df <- hdac_marker_expr_df %>%
  dplyr::mutate(
    facet_label = paste0(SYMBOL, "\n", marker_state),
    facet_label = factor(facet_label, levels = marker_facet_levels)
  )

hdac_marker_stats_plot <- hdac_marker_stats_plot %>%
  dplyr::mutate(
    facet_label = paste0(SYMBOL, "\n", marker_state),
    facet_label = factor(facet_label, levels = marker_facet_levels)
  )

## 8.8 Plot marker-expression boxplots

hdac_marker_expr_df <- hdac_marker_expr_df %>%
  dplyr::mutate(
    plot_group = dplyr::case_when(
      condition_plot == "Scramble" ~ "Scramble",
      condition_plot == "KD" & hdac_comparison == "HDAC1" ~ "HDAC1 KD",
      condition_plot == "KD" & hdac_comparison == "HDAC2" ~ "HDAC2 KD",
      condition_plot == "KD" & hdac_comparison == "HDAC3" ~ "HDAC3 KD"
    ),
    plot_group = factor(
      plot_group,
      levels = c("Scramble", "HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
    )
  )

marker_box_cols <- c(
  "Scramble" = "#D95F5F",  # orangey/brick
  "HDAC1 KD" = "#00BA38",  # green
  "HDAC2 KD" = "#619CFF",  # blue
  "HDAC3 KD" = "#00BFC4"   # cyan
)

hdac_marker_boxplot <- ggplot(
  hdac_marker_expr_df,
  aes(
    x = hdac_comparison,
    y = log2_norm_count
  )
) +
  geom_boxplot(
    aes(
      fill = plot_group,
      group = interaction(hdac_comparison, condition_plot)
    ),
    position = position_dodge(width = 0.75),
    width = 0.6,
    outlier.shape = NA,
    colour = "grey25",
    linewidth = 0.5
  ) +
  geom_point(
    aes(
      fill = plot_group,
      group = condition_plot
    ),
    shape = 21,
    colour = "black",
    stroke = 0.25,
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.75
    ),
    size = 1.8,
    alpha = 0.9
  ) +
  geom_text(
    data = hdac_marker_stats_plot,
    aes(
      x = hdac_comparison,
      y = y_pos,
      label = sig_label
    ),
    inherit.aes = FALSE,
    size = 3.2
  ) +
  facet_wrap(
    ~ facet_label,
    scales = "free_y",
    ncol = 4
  ) +
  scale_fill_manual(
    values = marker_box_cols,
    name = "Condition"
  ) +
  theme_bw() +
  labs(
    title = "Pluripotency marker expression after HDAC knockdown",
    subtitle = "* padj < 0.05 and |log2FC| ≥ log2(1.5); † padj < 0.05 only",
    x = "HDAC knockdown comparison",
    y = "log2(DESeq2-normalised count + 1)"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

hdac_marker_boxplot

ggsave(
  filename = file.path(
    hdac_marker_plot_dir,
    "Naive_primed_marker_expression_after_HDAC_knockdown_boxplots_clean.png"
  ),
  plot = hdac_marker_boxplot,
  width = 10,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

## 8.9 Naive-like marker score after HDAC knockdown

### 8.9.1 Scale expression within each marker gene
hdac_marker_score_df <- hdac_marker_expr_df %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::mutate(
    scaled_expr = as.numeric(scale(log2_norm_count))
  ) %>%
  dplyr::ungroup()

### 8.9.2 Calculate sample-level marker scores
hdac_marker_score_summary <- hdac_marker_score_df %>%
  dplyr::group_by(
    hdac_comparison,
    condition_plot,
    sample,
    marker_state
  ) %>%
  dplyr::summarise(
    mean_scaled_expr = mean(scaled_expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = marker_state,
    values_from = mean_scaled_expr
  ) %>%
  dplyr::mutate(
    naive_like_score = `Naive marker` - `Primed marker`
  )

hdac_marker_score_summary

### 8.9.3 Test the naive-like score

hdac_marker_score_stats <- hdac_marker_score_summary %>%
  dplyr::group_by(hdac_comparison) %>%
  dplyr::summarise(
    n_scramble = sum(condition_plot == "Scramble"),
    n_kd = sum(condition_plot == "KD"),
    median_scramble = median(naive_like_score[condition_plot == "Scramble"], na.rm = TRUE),
    median_kd = median(naive_like_score[condition_plot == "KD"], na.rm = TRUE),
    delta_median = median_kd - median_scramble,
    p_value = wilcox.test(
      naive_like_score ~ condition_plot,
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    padj = p.adjust(p_value, method = "BH")
  )

hdac_marker_score_stats
