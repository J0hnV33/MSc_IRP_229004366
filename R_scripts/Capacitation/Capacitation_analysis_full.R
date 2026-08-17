# Capacitation bulk RNA-seq analysis
#
# End-to-end workflow for sample QC, DESeq2 differential expression,
# X-chromosome analyses, HDAC and pluripotency-marker expression, and plots.
# Sections follow the existing execution order.
# Update the local input and output paths before running on another system.

# 1. Setup and data import 

## 1.1 Packages 

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(readr)
library(ggh4x)
library(ggVennDiagram)

## 1.2 Paths and output directories

# Root directory for all capacitation inputs and outputs.
base_dir <- "/home/jvk3/Desktop/HDAC_counts"

cap_dir <- file.path(base_dir, "Capacitation")

cap_tpm_dir <- file.path(cap_dir, "libinorm_normalised")
cap_counts_dir <- file.path(cap_dir, "libinorm_raw_counts")

cap_meta_file <- file.path(cap_dir, "capacitation_meta_table.csv")

cap_plot_dir <- file.path(cap_dir, "Capacitation_outputs")
dir.create(cap_plot_dir, showWarnings = FALSE, recursive = TRUE)

## 1.3 Read metadata

cap_meta <- read.csv(
  cap_meta_file,
  stringsAsFactors = FALSE
)

cap_meta

colnames(cap_meta)

dim(cap_meta)

## 1.4 Clean metadata

cap_meta <- cap_meta %>%
  dplyr::mutate(
    condition_clean = dplyr::case_when(
      condition == "Naive" ~ "Naive",
      condition == "Cap. Day 3" ~ "Day3",
      condition == "Cap. Day 7" ~ "Day7",
      condition == "Primed" ~ "Primed",
      TRUE ~ condition
    ),
    condition_clean = factor(
      condition_clean,
      levels = c("Naive", "Day3", "Day7", "Primed")
    ),
    sample = as.character(sample),
    filename = as.character(filename),
    raw_count_path = file.path(cap_counts_dir, filename),
    tpm_filename = gsub("_raw_counts\\.txt$", "_expression.txt", filename),
    tpm_path = file.path(cap_tpm_dir, tpm_filename)
  )

cap_meta

all(file.exists(cap_meta$raw_count_path))
all(file.exists(cap_meta$tpm_path))

## 1.5 Read and combine raw count files

read_cap_count_file <- function(file, sample_name) {

  message("Reading: ", basename(file))

  x <- read.delim(
    file,
    header = FALSE,
    stringsAsFactors = FALSE
  )

  x <- x %>%
    dplyr::select(
      gene_id = 1,
      count = 2
    ) %>%
    dplyr::filter(
      !grepl("^__", gene_id)
    ) %>%
    dplyr::mutate(
      gene_id = gsub("\\..*$", "", gene_id),
      count = as.numeric(count)
    ) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(
      count = sum(count, na.rm = TRUE),
      .groups = "drop"
    )

  colnames(x)[2] <- sample_name

  x
}

cap_count_list <- mapply(
  read_cap_count_file,
  file = cap_meta$raw_count_path,
  sample_name = cap_meta$sample,
  SIMPLIFY = FALSE
)

cap_counts <- Reduce(
  function(x, y) dplyr::full_join(x, y, by = "gene_id"),
  cap_count_list
) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

cap_counts[is.na(cap_counts)] <- 0
storage.mode(cap_counts) <- "integer"

cap_counts <- cap_counts[, cap_meta$sample]

stopifnot(all(colnames(cap_counts) == cap_meta$sample))

dim(cap_counts)
head(cap_counts[, 1:4])

# 2. Quality control and normalisation

## 2.1 Sample-level QC summary

cap_sample_qc <- tibble::tibble(
  sample = colnames(cap_counts),
  total_counts = colSums(cap_counts),
  detected_genes = colSums(cap_counts > 0)
) %>%
  dplyr::left_join(
    cap_meta %>%
      dplyr::select(sample, condition_clean, rep),
    by = "sample"
  )

cap_sample_qc

write.csv(
  cap_sample_qc,
  file.path(cap_plot_dir, "capacitation_sample_QC.csv"),
  row.names = FALSE
)

## 2.2 Library-size and detected-gene QC plots
# Compare total library sizes and detected-gene counts across conditions.
cap_qc_library_plot <- ggplot(
  cap_sample_qc,
  aes(
    x = condition_clean,
    y = total_counts
  )
) +
  geom_boxplot(fill = "white", colour = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2.5, colour = "black") +
  theme_bw() +
  labs(
    title = "Capacitation sample library sizes",
    x = "Condition",
    y = "Total raw counts"
  )

cap_qc_library_plot

ggsave(
  filename = file.path(cap_plot_dir, "01_capacitation_library_size_QC.png"),
  plot = cap_qc_library_plot,
  width = 6,
  height = 4.5,
  dpi = 300,
  bg = "white"
)

cap_qc_detected_plot <- ggplot(
  cap_sample_qc,
  aes(
    x = condition_clean,
    y = detected_genes
  )
) +
  geom_boxplot(fill = "white", colour = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2.5, colour = "black") +
  theme_bw() +
  labs(
    title = "Detected genes across capacitation samples",
    x = "Condition",
    y = "Number of detected genes"
  )

cap_qc_detected_plot

ggsave(
  filename = file.path(cap_plot_dir, "02_capacitation_detected_genes_QC.png"),
  plot = cap_qc_detected_plot,
  width = 6,
  height = 4.5,
  dpi = 300,
  bg = "white"
)

## 2.3 DESeq2 model and PCA

cap_coldata <- cap_meta %>%
  dplyr::select(
    sample,
    condition_clean,
    rep
  ) %>%
  as.data.frame()

rownames(cap_coldata) <- cap_coldata$sample

dds_cap <- DESeqDataSetFromMatrix(
  countData = cap_counts,
  colData = cap_coldata,
  design = ~ condition_clean
)

keep <- rowSums(counts(dds_cap)) >= 10
dds_cap <- dds_cap[keep, ]

dds_cap <- DESeq(dds_cap)

vsd_cap <- vst(dds_cap, blind = FALSE)

pca_cap_raw <- plotPCA(
  vsd_cap,
  intgroup = "condition_clean",
  ntop = 500,
  returnData = TRUE
)

percentVar_cap <- round(100 * attr(pca_cap_raw, "percentVar"))

pca_cap_df <- as.data.frame(pca_cap_raw)
pca_cap_df$sample <- rownames(pca_cap_df)

pca_cap_plot <- ggplot(
  pca_cap_df,
  aes(
    x = PC1,
    y = PC2,
    colour = condition_clean,
    label = sample
  )
) +
  geom_point(size = 4, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) +
  theme_bw() +
  labs(
    title = "PCA of capacitation RNA-seq samples",
    subtitle = "VST-transformed raw counts; DESeq2 normalisation",
    x = paste0("PC1: ", percentVar_cap[1], "% variance"),
    y = paste0("PC2: ", percentVar_cap[2], "% variance"),
    colour = "Condition"
  )

pca_cap_plot

ggsave(
  filename = file.path(cap_plot_dir, "03_capacitation_PCA.png"),
  plot = pca_cap_plot,
  width = 7,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

## 2.4 Gene annotation and unfiltered normalised counts

# Reuse the combined HDAC annotation as the gene-annotation source.
all_res_annot_source <- read.csv(
  file.path(
    base_dir,
    "All_Res_Tables",
    "all_results_HDAC1_2_3_combined_15fold.csv"
  ),
  stringsAsFactors = FALSE
)

gene_annot_cap <- all_res_annot_source %>%
  dplyr::mutate(
    gene_id = gsub("\\..*$", "", gene_id),
    CHR = as.character(CHR)
  ) %>%
  dplyr::select(
    gene_id,
    SYMBOL,
    CHR,
    START,
    END,
    GENE_BIOTYPE,
    gene_label
  ) %>%
  dplyr::distinct(gene_id, .keep_all = TRUE)

## 2.5 Unfiltered normalised counts for marker analyses
dds_cap_all_genes <- DESeqDataSetFromMatrix(
  countData = cap_counts,
  colData = cap_coldata,
  design = ~ condition_clean
)

dds_cap_all_genes <- estimateSizeFactors(dds_cap_all_genes)

cap_norm_counts_all_genes <- counts(
  dds_cap_all_genes,
  normalized = TRUE
)

## 2.6 Match count identifiers to gene annotation

cap_count_ids <- tibble::tibble(
  count_gene_id = rownames(cap_norm_counts_all_genes)
)

annot_by_ensembl <- cap_count_ids %>%
  dplyr::left_join(
    gene_annot_cap %>%
      dplyr::rename(count_gene_id = gene_id),
    by = "count_gene_id"
  )

annot_by_symbol <- cap_count_ids %>%
  dplyr::left_join(
    gene_annot_cap %>%
      dplyr::rename(count_gene_id = SYMBOL),
    by = "count_gene_id"
  )

cap_gene_annot_final <- dplyr::if_else(
  sum(!is.na(annot_by_ensembl$SYMBOL)) >= sum(!is.na(annot_by_symbol$gene_id)),
  TRUE,
  FALSE
)

if (cap_gene_annot_final) {

  cap_gene_annot_final <- annot_by_ensembl %>%
    dplyr::mutate(mapping_used = "Ensembl ID")

} else {

  cap_gene_annot_final <- annot_by_symbol %>%
    dplyr::mutate(
      SYMBOL = count_gene_id,
      mapping_used = "Gene symbol"
    )
}

cap_gene_annot_final %>%
  dplyr::count(mapping_used)

## 2.7 Build a long-format annotated count table

cap_norm_long_all_genes <- cap_norm_counts_all_genes %>%
  as.data.frame() %>%
  tibble::rownames_to_column("count_gene_id") %>%
  tidyr::pivot_longer(
    cols = -count_gene_id,
    names_to = "sample",
    values_to = "norm_count"
  ) %>%
  dplyr::left_join(
    cap_meta %>%
      dplyr::select(sample, condition_clean, rep),
    by = "sample"
  ) %>%
  dplyr::left_join(
    cap_gene_annot_final,
    by = "count_gene_id"
  )

# 3. HDAC expression during capacitation

## 3.1 Prepare HDAC expression data

hdac_genes <- c("HDAC1", "HDAC2", "HDAC3")

cap_condition_cols <- c(
  "Naive"  = "#F8766D",
  "Day3"   = "#7CAE00",
  "Day7"   = "#00BFC4",
  "Primed" = "#C77CFF"
)

hdac_cap_expr <- cap_norm_long_all_genes %>%
  dplyr::filter(
    SYMBOL %in% hdac_genes
  ) %>%
  dplyr::mutate(
    SYMBOL = factor(SYMBOL, levels = hdac_genes),
    condition_clean = factor(
      condition_clean,
      levels = c("Naive", "Day3", "Day7", "Primed")
    )
  )

hdac_cap_expr %>%
  dplyr::count(SYMBOL, condition_clean)

## 3.2 Add DESeq2 statistics

all_res_capacitation_15fold <- readr::read_csv(
  "/home/jvk3/Desktop/HDAC_counts/Capacitation/Capacitation_outputs/all_results_capacitation_15fold.csv",
  show_col_types = FALSE
)

if (exists("all_res_capacitation_15fold")) {
  cap_res_for_stats <- all_res_capacitation_15fold
} else if (exists("cap_res")) {
  cap_res_for_stats <- cap_res
} else {
  stop("Could not find all_res_capacitation_15fold or cap_res")
}

hdac_deseq_stats <- cap_res_for_stats %>%
  dplyr::filter(
    SYMBOL %in% hdac_genes,
    contrast_short %in% c("Day3 vs Naive", "Day7 vs Naive", "Primed vs Naive")
  ) %>%
  dplyr::mutate(
    SYMBOL = factor(SYMBOL, levels = hdac_genes),
    condition_clean = dplyr::case_when(
      contrast_short == "Day3 vs Naive" ~ "Day3",
      contrast_short == "Day7 vs Naive" ~ "Day7",
      contrast_short == "Primed vs Naive" ~ "Primed",
      TRUE ~ NA_character_
    ),
    condition_clean = factor(
      condition_clean,
      levels = c("Naive", "Day3", "Day7", "Primed")
    ),
    stat_label = dplyr::case_when(
      is.na(padj) ~ "padj = NA",
      padj < 0.001 ~ paste0("padj < 0.001\nlog2FC = ", round(log2FoldChange, 2)),
      TRUE ~ paste0("padj = ", signif(padj, 2), "\nlog2FC = ", round(log2FoldChange, 2))
    )
  ) %>%
  dplyr::filter(!is.na(condition_clean))

hdac_deseq_stats

## 3.3 Calculate label positions

hdac_label_positions <- hdac_cap_expr %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(
    y_max = max(norm_count, na.rm = TRUE),
    y_min = min(norm_count, na.rm = TRUE),
    y_range = y_max - y_min,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    y_pos = y_max + dplyr::if_else(
      y_range == 0,
      0.08 * y_max,
      0.12 * y_range
    )
  )

hdac_deseq_stats <- hdac_deseq_stats %>%
  dplyr::left_join(
    hdac_label_positions,
    by = "SYMBOL"
  )

hdac_deseq_stats

## 3.4 Plot HDAC expression with DESeq2 statistics

hdac_cap_plot_deseq_stats <- ggplot(
  hdac_cap_expr,
  aes(
    x = condition_clean,
    y = norm_count,
    fill = condition_clean
  )
) +
  geom_boxplot(
    colour = "black",
    outlier.shape = NA,
    alpha = 0.75,
    width = 0.65
  ) +
  geom_jitter(
    aes(colour = condition_clean),
    width = 0.12,
    size = 2.4,
    alpha = 0.9
  ) +
  geom_text(
    data = hdac_deseq_stats,
    aes(
      x = condition_clean,
      y = y_pos,
      label = stat_label
    ),
    inherit.aes = FALSE,
    size = 3,
    lineheight = 0.9
  ) +
  facet_wrap(
    ~ SYMBOL,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(values = cap_condition_cols) +
  scale_colour_manual(values = cap_condition_cols) +
  theme_bw() +
  labs(
    title = "HDAC1/2/3 expression during capacitation",
    subtitle = "DESeq2-normalised counts across Naive, Day3, Day7 and Primed states; DESeq2 statistics shown relative to Naive",
    x = "Cell state / timepoint",
    y = "DESeq2-normalised count",
    fill = "Condition",
    colour = "Condition"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

hdac_cap_plot_deseq_stats

ggsave(
  file.path(
    cap_plot_dir,
    "HDAC1_2_3_expression_during_capacitation_coloured_DESeq2_stats.png"
  ),
  hdac_cap_plot_deseq_stats,
  width = 13,
  height = 5.5,
  dpi = 300
)

readr::write_csv(
  hdac_deseq_stats,
  file.path(
    cap_plot_dir,
    "HDAC1_2_3_capacitation_DESeq2_stats_vs_Naive.csv"
  )
)

## 3.5 Calculate pairwise Wilcoxon tests against Naive

comparisons_vs_naive <- tibble::tibble(
  group1 = "Naive",
  group2 = c("Day3", "Day7", "Primed"),
  x_start = 1,
  x_end = c(2, 3, 4)
)

hdac_stat_labels <- hdac_cap_expr %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::group_modify(~ {

    gene_df <- .x

    stat_df <- comparisons_vs_naive %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        p_value = {
          x <- gene_df %>%
            dplyr::filter(condition_clean == group1) %>%
            dplyr::pull(norm_count)

          y <- gene_df %>%
            dplyr::filter(condition_clean == group2) %>%
            dplyr::pull(norm_count)

          if (length(x) >= 2 && length(y) >= 2) {
            wilcox.test(x, y, exact = FALSE)$p.value
          } else {
            NA_real_
          }
        }
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        p_adj = p.adjust(p_value, method = "BH"),
        p_label = dplyr::case_when(
          is.na(p_adj) ~ "NA",
          p_adj < 0.001 ~ "***",
          p_adj < 0.01 ~ "**",
          p_adj < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      )

    y_max <- max(gene_df$norm_count, na.rm = TRUE)
    y_min <- min(gene_df$norm_count, na.rm = TRUE)
    y_range <- y_max - y_min

    stat_df %>%
      dplyr::mutate(
        y_pos = y_max + dplyr::row_number() * 0.12 * y_range
      )
  }) %>%
  dplyr::ungroup()

hdac_stat_labels

## 3.6 Plot HDAC expression with pairwise statistics

hdac_stat_labels <- hdac_stat_labels %>%
  dplyr::filter(p_label != "ns")

hdac_cap_plot <- ggplot(
  hdac_cap_expr,
  aes(
    x = condition_clean,
    y = norm_count,
    fill = condition_clean
  )
) +
  geom_boxplot(
    colour = "black",
    outlier.shape = NA,
    alpha = 0.75,
    width = 0.65
  ) +
  geom_jitter(
    aes(colour = condition_clean),
    width = 0.12,
    size = 2.4,
    alpha = 0.9
  ) +
  geom_segment(
    data = hdac_stat_labels,
    aes(
      x = x_start,
      xend = x_end,
      y = y_pos,
      yend = y_pos
    ),
    inherit.aes = FALSE,
    colour = "black"
  ) +
  geom_segment(
    data = hdac_stat_labels,
    aes(
      x = x_start,
      xend = x_start,
      y = y_pos,
      yend = y_pos - 0.03 * y_pos
    ),
    inherit.aes = FALSE,
    colour = "black"
  ) +
  geom_segment(
    data = hdac_stat_labels,
    aes(
      x = x_end,
      xend = x_end,
      y = y_pos,
      yend = y_pos - 0.03 * y_pos
    ),
    inherit.aes = FALSE,
    colour = "black"
  ) +
  geom_text(
    data = hdac_stat_labels,
    aes(
      x = (x_start + x_end) / 2,
      y = y_pos,
      label = p_label
    ),
    inherit.aes = FALSE,
    vjust = -0.4,
    size = 3.5
  ) +
  facet_wrap(
    ~ SYMBOL,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(values = cap_condition_cols) +
  scale_colour_manual(values = cap_condition_cols) +
  theme_bw() +
  labs(
    title = "HDAC1/2/3 expression during capacitation",
    subtitle = "DESeq2-normalised counts across Naive, Day3, Day7 and Primed states; Wilcoxon tests vs Naive",
    x = "Cell state / timepoint",
    y = "DESeq2-normalised count",
    fill = "Condition",
    colour = "Condition"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

hdac_cap_plot

ggsave(
  filename = file.path(cap_plot_dir, "05_HDAC_expression_capacitation.png"),
  plot = hdac_cap_plot,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# 4. X-chromosome expression and differential expression

## 4.1 Calculate the X-to-autosome expression ratio

cap_norm_annotated <- cap_norm_counts_all_genes %>%
  as.data.frame() %>%
  tibble::rownames_to_column("count_gene_id") %>%
  dplyr::left_join(
    cap_gene_annot_final,
    by = "count_gene_id"
  )

x_count_ids_cap <- cap_norm_annotated %>%
  dplyr::filter(CHR == "X") %>%
  dplyr::pull(count_gene_id)

autosome_count_ids_cap <- cap_norm_annotated %>%
  dplyr::filter(CHR %in% as.character(1:22)) %>%
  dplyr::pull(count_gene_id)

length(x_count_ids_cap)
length(autosome_count_ids_cap)

x_auto_ratio_cap <- tibble::tibble(
  sample = colnames(cap_norm_counts_all_genes),
  mean_X_expression = colMeans(
    cap_norm_counts_all_genes[x_count_ids_cap, , drop = FALSE]
  ),
  mean_autosome_expression = colMeans(
    cap_norm_counts_all_genes[autosome_count_ids_cap, , drop = FALSE]
  )
) %>%
  dplyr::mutate(
    X_autosome_ratio = mean_X_expression / mean_autosome_expression
  ) %>%
  dplyr::left_join(
    cap_meta %>%
      dplyr::select(sample, condition_clean, rep),
    by = "sample"
  )

x_auto_ratio_cap

x_auto_ratio_plot <- ggplot(
  x_auto_ratio_cap,
  aes(
    x = condition_clean,
    y = X_autosome_ratio
  )
) +
  geom_boxplot(
    fill = "white",
    colour = "black",
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.5,
    colour = "black"
  ) +
  theme_bw() +
  labs(
    title = "X:autosome expression ratio during capacitation",
    subtitle = "Mean X-linked expression divided by mean autosomal expression",
    x = "Cell state / timepoint",
    y = "X:autosome expression ratio"
  )

x_auto_ratio_plot

ggsave(
  filename = file.path(cap_plot_dir, "06_X_autosome_ratio_capacitation.png"),
  plot = x_auto_ratio_plot,
  width = 6,
  height = 4.5,
  dpi = 300,
  bg = "white"
)

write.csv(
  x_auto_ratio_cap,
  file.path(cap_plot_dir, "X_autosome_ratio_capacitation.csv"),
  row.names = FALSE
)

## 4.2 Fit DESeq2 contrasts across capacitation stages

cap_coldata$condition_clean <- factor(
  cap_coldata$condition_clean,
  levels = c("Naive", "Day3", "Day7", "Primed")
)

colData(dds_cap)$condition_clean <- factor(
  colData(dds_cap)$condition_clean,
  levels = c("Naive", "Day3", "Day7", "Primed")
)

dds_cap <- DESeq(dds_cap)

resultsNames(dds_cap)

## 4.3 Extract Naive-referenced DESeq2 results

lfc_cutoff_15fold <- log2(1.5)
padj_cutoff <- 0.05

cap_contrast_table <- tibble::tibble(
  contrast_short = c(
    "Day3 vs Naive",
    "Day7 vs Naive",
    "Primed vs Naive"
  ),
  numerator = c(
    "Day3",
    "Day7",
    "Primed"
  ),
  denominator = c(
    "Naive",
    "Naive",
    "Naive"
  )
)
get_cap_contrast <- function(numerator, denominator, contrast_short) {

  res <- DESeq2::results(
    dds_cap,
    contrast = c("condition_clean", numerator, denominator),
    alpha = padj_cutoff
  )

  res_tb <- as.data.frame(res) %>%
    tibble::rownames_to_column("count_gene_id") %>%
    dplyr::left_join(
      cap_gene_annot_final,
      by = "count_gene_id"
    ) %>%
    dplyr::mutate(
      contrast_short = contrast_short,
      numerator = numerator,
      denominator = denominator,
      deg_status_15fold = dplyr::case_when(
        !is.na(padj) &
          padj < padj_cutoff &
          log2FoldChange >= lfc_cutoff_15fold ~ "Up",
        !is.na(padj) &
          padj < padj_cutoff &
          log2FoldChange <= -lfc_cutoff_15fold ~ "Down",
        TRUE ~ "Not significant"
      ),
      deg_status_15fold = factor(
        deg_status_15fold,
        levels = c("Down", "Not significant", "Up")
      ),
      contrast_short = factor(
        contrast_short,
        levels = cap_contrast_table$contrast_short
      )
    )

  res_tb
}

cap_res_list <- mapply(
  get_cap_contrast,
  numerator = cap_contrast_table$numerator,
  denominator = cap_contrast_table$denominator,
  contrast_short = cap_contrast_table$contrast_short,
  SIMPLIFY = FALSE
)

all_res_capacitation_15fold <- dplyr::bind_rows(cap_res_list)

dim(all_res_capacitation_15fold)

## 4.4 Summarise genome-wide and X-linked DEGs

cap_deg_summary <- all_res_capacitation_15fold %>%
  dplyr::count(
    contrast_short,
    deg_status_15fold,
    name = "n_genes"
  )

cap_x_deg_summary <- all_res_capacitation_15fold %>%
  dplyr::filter(CHR == "X") %>%
  dplyr::count(
    contrast_short,
    deg_status_15fold,
    name = "n_X_genes"
  )

cap_deg_summary

cap_x_deg_summary

## 4.5 Save DESeq2 result tables

write.csv(
  all_res_capacitation_15fold,
  file.path(cap_plot_dir, "all_results_capacitation_15fold.csv"),
  row.names = FALSE
)

write.csv(
  cap_deg_summary,
  file.path(cap_plot_dir, "capacitation_DEG_summary_genomewide_15fold.csv"),
  row.names = FALSE
)

write.csv(
  cap_x_deg_summary,
  file.path(cap_plot_dir, "capacitation_DEG_summary_chrX_15fold.csv"),
  row.names = FALSE
)

## 4.6 Build the X-linked DEG table

x_deg_cap_table <- all_res_capacitation_15fold %>%
  dplyr::filter(
    CHR == "X",
    deg_status_15fold != "Not significant"
  ) %>%
  dplyr::mutate(
    gene_plot_label = dplyr::case_when(
      !is.na(SYMBOL) & SYMBOL != "" ~ SYMBOL,
      !is.na(gene_label) & gene_label != "" ~ gene_label,
      TRUE ~ count_gene_id
    )
  )

x_deg_cap_table %>%
  dplyr::count(contrast_short, deg_status_15fold)

## 4.7 Plot X-linked DEG counts

cap_x_deg_count_plot <- ggplot(
  cap_x_deg_summary %>%
    dplyr::filter(deg_status_15fold != "Not significant"),
  aes(
    x = contrast_short,
    y = n_X_genes,
    fill = deg_status_15fold
  )
) +
  geom_col(
    position = "dodge",
    colour = "black"
  ) +
  theme_bw() +
  labs(
    title = "X-linked differentially expressed genes during capacitation",
    subtitle = "Threshold: padj < 0.05 and |log2FC| ≥ log2(1.5)",
    x = "Comparison",
    y = "Number of X-linked DEGs",
    fill = "Direction"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

cap_x_deg_count_plot

ggsave(
  filename = file.path(cap_plot_dir, "07_UPDATED_Xlinked_DEG_counts_capacitation.png"),
  plot = cap_x_deg_count_plot,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

## 4.8 Plot overlap of all X-linked DEGs

x_deg_sets_vs_naive <- list(
  "Day3 vs Naive" = x_deg_cap_table %>%
    dplyr::filter(contrast_short == "Day3 vs Naive") %>%
    dplyr::pull(count_gene_id) %>%
    unique(),

  "Day7 vs Naive" = x_deg_cap_table %>%
    dplyr::filter(contrast_short == "Day7 vs Naive") %>%
    dplyr::pull(count_gene_id) %>%
    unique(),

  "Primed vs Naive" = x_deg_cap_table %>%
    dplyr::filter(contrast_short == "Primed vs Naive") %>%
    dplyr::pull(count_gene_id) %>%
    unique()
)

venn_x_cap_vs_naive <- ggVennDiagram(
  x_deg_sets_vs_naive,
  label = "count"
) +
  theme_void() +
  labs(
    title = "Overlap of X-linked DEGs during capacitation",
    subtitle = "Comparisons relative to naïve cells"
  )

venn_x_cap_vs_naive

ggsave(
  filename = file.path(cap_plot_dir, "08_Xlinked_DEG_venn_vs_Naive.png"),
  plot = venn_x_cap_vs_naive,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

## 4.9 Plot overlap of upregulated X-linked DEGs

x_up_sets_vs_naive <- list(
  "Day3 vs Naive" = x_deg_cap_table %>%
    dplyr::filter(contrast_short == "Day3 vs Naive", deg_status_15fold == "Up") %>%
    dplyr::pull(count_gene_id) %>%
    unique(),

  "Day7 vs Naive" = x_deg_cap_table %>%
    dplyr::filter(contrast_short == "Day7 vs Naive", deg_status_15fold == "Up") %>%
    dplyr::pull(count_gene_id) %>%
    unique(),

  "Primed vs Naive" = x_deg_cap_table %>%
    dplyr::filter(contrast_short == "Primed vs Naive", deg_status_15fold == "Up") %>%
    dplyr::pull(count_gene_id) %>%
    unique()
)

venn_x_up_cap_vs_naive <- ggVennDiagram(
  x_up_sets_vs_naive,
  label = "count"
) +
  theme_void() +
  labs(
    title = "Overlap of upregulated X-linked genes during capacitation",
    subtitle = "Comparisons relative to naïve cells"
  )

venn_x_up_cap_vs_naive

ggsave(
  filename = file.path(cap_plot_dir, "09_Xlinked_upregulated_DEG_venn_vs_Naive.png"),
  plot = venn_x_up_cap_vs_naive,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

## 4.10 Plot overlap of downregulated X-linked DEGs

x_down_sets_vs_naive <- list(
  "Day3 vs Naive" = x_deg_cap_table %>%
    dplyr::filter(contrast_short == "Day3 vs Naive", deg_status_15fold == "Down") %>%
    dplyr::pull(count_gene_id) %>%
    unique(),

  "Day7 vs Naive" = x_deg_cap_table %>%
    dplyr::filter(contrast_short == "Day7 vs Naive", deg_status_15fold == "Down") %>%
    dplyr::pull(count_gene_id) %>%
    unique(),

  "Primed vs Naive" = x_deg_cap_table %>%
    dplyr::filter(contrast_short == "Primed vs Naive", deg_status_15fold == "Down") %>%
    dplyr::pull(count_gene_id) %>%
    unique()
)

venn_x_down_cap_vs_naive <- ggVennDiagram(
  x_down_sets_vs_naive,
  label = "count"
) +
  theme_void() +
  labs(
    title = "Overlap of downregulated X-linked genes during capacitation",
    subtitle = "Comparisons relative to naïve cells"
  )

venn_x_down_cap_vs_naive

ggsave(
  filename = file.path(cap_plot_dir, "10_Xlinked_downregulated_DEG_venn_vs_Naive.png"),
  plot = venn_x_down_cap_vs_naive,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

## 4.11 Plot the X-linked DEG heatmap

# Retain genes significant in at least one X-linked contrast.
x_sig_genes_cap <- x_deg_cap_table %>%
  dplyr::pull(count_gene_id) %>%
  unique()

length(x_sig_genes_cap)

x_heatmap_cap_df <- all_res_capacitation_15fold %>%
  dplyr::filter(
    count_gene_id %in% x_sig_genes_cap
  ) %>%
  dplyr::mutate(
    gene_plot_label = dplyr::case_when(
      !is.na(SYMBOL) & SYMBOL != "" ~ SYMBOL,
      !is.na(gene_label) & gene_label != "" ~ gene_label,
      TRUE ~ count_gene_id
    ),
    is_significant = deg_status_15fold != "Not significant",
    neg_log10_padj = dplyr::case_when(
      !is.na(padj) ~ -log10(pmax(padj, .Machine$double.xmin)),
      TRUE ~ 0
    ),
    ranking_score = abs(log2FoldChange) * neg_log10_padj,
    contrast_short = factor(
      contrast_short,
      levels = c(
        "Day3 vs Naive",
        "Day7 vs Naive",
        "Primed vs Naive",
        "Day7 vs Day3",
        "Primed vs Day7"
      )
    )
  )

# Select the top 80 genes by combined significance and effect size.
top_x_heatmap_genes_cap <- x_heatmap_cap_df %>%
  dplyr::group_by(count_gene_id, gene_plot_label) %>%
  dplyr::summarise(
    max_score = max(ranking_score, na.rm = TRUE),
    max_abs_log2FC = max(abs(log2FoldChange), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(max_score),
    dplyr::desc(max_abs_log2FC)
  ) %>%
  dplyr::slice_head(n = 80) %>%
  dplyr::pull(count_gene_id)

x_heatmap_top_cap_df <- x_heatmap_cap_df %>%
  dplyr::filter(
    count_gene_id %in% top_x_heatmap_genes_cap
  )

# Order genes so the strongest combined scores appear at the top.
gene_order_x_heatmap_cap <- x_heatmap_top_cap_df %>%
  dplyr::group_by(count_gene_id, gene_plot_label) %>%
  dplyr::summarise(
    max_score = max(ranking_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(max_score)
  ) %>%
  dplyr::pull(gene_plot_label)

x_heatmap_top_cap_df <- x_heatmap_top_cap_df %>%
  dplyr::mutate(
    gene_plot_label = factor(
      gene_plot_label,
      levels = rev(unique(gene_order_x_heatmap_cap))
    )
  )

xlinked_cap_heatmap_plot <- ggplot(
  x_heatmap_top_cap_df,
  aes(
    x = contrast_short,
    y = gene_plot_label,
    fill = log2FoldChange
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.35
  ) +
  geom_point(
    data = x_heatmap_top_cap_df %>%
      dplyr::filter(is_significant),
    aes(
      x = contrast_short,
      y = gene_plot_label
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.8,
    fill = "black",
    colour = "black"
  ) +
  scale_fill_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0,
    name = "log2FC"
  ) +
  theme_bw() +
  labs(
    title = "Directional changes in X-linked genes during capacitation",
    subtitle = "Top 80 X-linked DEGs; cut-off: padj < 0.05 and |log2FC| ≥ log2(1.5); black dots mark significant contrasts",
    x = "Comparison",
    y = "X-linked genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill = "white", colour = "black")
  )

xlinked_cap_heatmap_plot

ggsave(
  filename = file.path(cap_plot_dir, "11_Xlinked_DEG_heatmap_capacitation.png"),
  plot = xlinked_cap_heatmap_plot,
  width = 9,
  height = 12,
  dpi = 300,
  bg = "white"
)

write.csv(
  x_heatmap_top_cap_df,
  file.path(cap_plot_dir, "Xlinked_DEG_heatmap_top80_capacitation.csv"),
  row.names = FALSE
)

write.csv(
  x_heatmap_cap_df,
  file.path(cap_plot_dir, "Xlinked_DEG_heatmap_all_capacitation.csv"),
  row.names = FALSE
)

## 4.12 Plot the top X-linked DEGs by contrast

top_x_cap_genes <- x_deg_cap_table %>%
  dplyr::filter(
    !is.na(log2FoldChange),
    !is.na(padj)
  ) %>%
  dplyr::mutate(
    neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)),
    abs_log2FC = abs(log2FoldChange),
    ranking_score = abs_log2FC * neg_log10_padj,
    gene_plot_label = dplyr::case_when(
      !is.na(SYMBOL) & SYMBOL != "" ~ SYMBOL,
      !is.na(gene_label) & gene_label != "" ~ gene_label,
      TRUE ~ count_gene_id
    ),
    contrast_short = factor(
      contrast_short,
      levels = c(
        "Day3 vs Naive",
        "Day7 vs Naive",
        "Primed vs Naive",
        "Day7 vs Day3",
        "Primed vs Day7"
      )
    )
  ) %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::slice_max(
    order_by = ranking_score,
    n = 10,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

# Create a clean within-contrast gene order.

top_x_cap_genes <- top_x_cap_genes %>%
  dplyr::arrange(
    contrast_short,
    ranking_score
  ) %>%
  dplyr::mutate(
    gene_contrast_label = paste(
      contrast_short,
      gene_plot_label,
      sep = "___"
    ),
    gene_contrast_label = factor(
      gene_contrast_label,
      levels = unique(gene_contrast_label)
    )
  )

# Plot the selected genes.
top_x_cap_plot <- ggplot(
  top_x_cap_genes,
  aes(
    x = log2FoldChange,
    y = gene_contrast_label,
    colour = deg_status_15fold,
    size = neg_log10_padj
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_point(alpha = 0.9) +
  facet_wrap(
    ~ contrast_short,
    scales = "free_y"
  ) +
  scale_y_discrete(
    labels = function(x) sub("^.*___", "", x)
  ) +
  scale_colour_manual(
    values = c(
      "Up" = "forestgreen",
      "Down" = "red"
    ),
    labels = c(
      "Up" = "Upregulated",
      "Down" = "Downregulated"
    )
  ) +
  theme_bw() +
  labs(
    title = "Top X-linked differentially expressed genes during capacitation",
    subtitle = "Top 10 per contrast; cut-off: padj < 0.05 and |log2FC| ≥ log2(1.5)",
    x = "Log2 fold change",
    y = "X-linked gene",
    colour = "Direction",
    size = "-log10(padj)"
  ) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 0),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )

top_x_cap_plot

ggsave(
  filename = file.path(cap_plot_dir, "12_UPDATED_Top10_Xlinked_DEGs_capacitation.png"),
  plot = top_x_cap_plot,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

write.csv(
  top_x_cap_genes,
  file.path(cap_plot_dir, "Top20_Xlinked_DEGs_capacitation.csv"),
  row.names = FALSE
)

# 5. Pluripotency-marker and HDAC trajectory analyses

## 5.1 Naive and primed marker expression

marker_genes <- tibble::tibble(
  SYMBOL = c(
    "DPPA5", "DNMT3L", "TFCP2L1", "KLF5",
    "DNMT3B", "SALL2", "ZIC2", "PODXL"
  ),
  marker_state = c(
    rep("Naive marker", 4),
    rep("Primed marker", 4)
  )
)

### 5.1.1 Extract DESeq2-normalised counts
cap_norm_counts <- DESeq2::counts(
  dds_cap,
  normalized = TRUE
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("count_gene_id")

### 5.1.2 Add gene annotation
cap_norm_counts_annot <- cap_norm_counts %>%
  dplyr::left_join(
    cap_gene_annot_final,
    by = "count_gene_id"
  )

### 5.1.3 Check marker-gene availability
marker_check <- marker_genes %>%
  dplyr::left_join(
    cap_norm_counts_annot %>%
      dplyr::distinct(SYMBOL, count_gene_id),
    by = "SYMBOL"
  )

marker_check

### 5.1.4 Identify sample columns
sample_cols <- colnames(DESeq2::counts(dds_cap))

### 5.1.5 Build sample metadata
cap_sample_meta <- cap_coldata %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_id")

### 5.1.6 Check sample identifiers
sample_cols
cap_sample_meta$sample_id

setdiff(sample_cols, cap_sample_meta$sample_id)
setdiff(cap_sample_meta$sample_id, sample_cols)

marker_expr_df <- cap_norm_counts_annot %>%
  dplyr::filter(SYMBOL %in% marker_genes$SYMBOL) %>%
  dplyr::left_join(
    marker_genes,
    by = "SYMBOL"
  ) %>%
  dplyr::select(
    count_gene_id,
    SYMBOL,
    marker_state,
    dplyr::all_of(sample_cols)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_cols),
    names_to = "sample_id",
    values_to = "norm_count"
  ) %>%
  dplyr::left_join(
    cap_sample_meta %>%
      dplyr::select(sample_id, condition_clean),
    by = "sample_id"
  ) %>%
  dplyr::mutate(
    condition_clean = factor(
      condition_clean,
      levels = c("Naive", "Day3", "Day7", "Primed")
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

marker_expr_df %>%
  dplyr::count(SYMBOL, marker_state, condition_clean)

pluripotency_marker_boxplot <- ggplot(
  marker_expr_df,
  aes(
    x = condition_clean,
    y = norm_count
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    fill = "white",
    colour = "black"
  ) +
  geom_jitter(
    width = 0.12,
    size = 2,
    alpha = 0.8
  ) +
  facet_wrap(
    marker_state ~ SYMBOL,
    scales = "free_y",
    ncol = 4
  ) +
  theme_bw() +
  labs(
    title = "Naive and primed pluripotency marker expression during capacitation",
    subtitle = "DESeq2-normalised counts across naïve, capacitating and primed states",
    x = "Cell state / timepoint",
    y = "DESeq2-normalised count"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.text = element_text(face = "bold")
  )

pluripotency_marker_boxplot

ggsave(
  file.path(
    cap_plot_dir,
    "pluripotency_marker_expression_during_capacitation.png"
  ),
  pluripotency_marker_boxplot,
  width = 12,
  height = 7,
  dpi = 300,
  bg = "white"
)

## 5.2 HDAC and pluripotency-marker heatmap

cap_marker_heatmap_genes <- tibble::tibble(
  SYMBOL = c(
    "HDAC1", "HDAC2", "HDAC3",
    "DPPA5", "DNMT3L", "TFCP2L1", "KLF5",
    "DNMT3B", "SALL2", "ZIC2", "PODXL"
  ),
  gene_group = c(
    rep("HDAC", 3),
    rep("Naive marker", 4),
    rep("Primed marker", 4)
  )
)

sample_cols <- colnames(DESeq2::counts(dds_cap))

cap_sample_meta <- cap_coldata %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_id")

hdac_marker_expr_df <- cap_norm_counts_annot %>%
  dplyr::filter(SYMBOL %in% cap_marker_heatmap_genes$SYMBOL) %>%
  dplyr::left_join(
    cap_marker_heatmap_genes,
    by = "SYMBOL"
  ) %>%
  dplyr::select(
    count_gene_id,
    SYMBOL,
    gene_group,
    dplyr::all_of(sample_cols)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_cols),
    names_to = "sample_id",
    values_to = "norm_count"
  ) %>%
  dplyr::left_join(
    cap_sample_meta %>%
      dplyr::select(sample_id, condition_clean),
    by = "sample_id"
  ) %>%
  dplyr::mutate(
    condition_clean = factor(
      condition_clean,
      levels = c("Naive", "Day3", "Day7", "Primed")
    ),
    SYMBOL = factor(
      SYMBOL,
      levels = cap_marker_heatmap_genes$SYMBOL
    ),
    gene_group = factor(
      gene_group,
      levels = c("HDAC", "Naive marker", "Primed marker")
    )
  )

### 5.2.1 Calculate mean expression by timepoint
hdac_marker_mean_df <- hdac_marker_expr_df %>%
  dplyr::group_by(SYMBOL, gene_group, condition_clean) %>%
  dplyr::summarise(
    mean_norm_count = mean(norm_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::mutate(
    scaled_expression = as.numeric(scale(log2(mean_norm_count + 1)))
  ) %>%
  dplyr::ungroup()

hdac_marker_heatmap <- ggplot(
  hdac_marker_mean_df,
  aes(
    x = condition_clean,
    y = SYMBOL,
    fill = scaled_expression
  )
) +
  geom_tile(
    colour = "white"
  ) +
  facet_grid(
    gene_group ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  theme_bw() +
  labs(
    title = "HDAC and pluripotency marker expression during capacitation",
    subtitle = "Mean DESeq2-normalised counts, row-scaled within each gene",
    x = "Cell state / timepoint",
    y = "Gene",
    fill = "Scaled\nexpression"
  ) +
  theme(
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

hdac_marker_heatmap

ggsave(
  file.path(
    cap_plot_dir,
    "hdac_marker_heatmap.png"
  ),
  hdac_marker_heatmap,
  width = 12,
  height = 7,
  dpi = 300,
  bg = "white"
)

## 5.3 HDAC expression versus marker-programme scores

sample_cols <- colnames(DESeq2::counts(dds_cap))

cap_sample_meta <- cap_coldata %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_id")

### 5.3.1 Build marker-programme scores per sample

marker_score_df <- cap_norm_counts_annot %>%
  dplyr::filter(SYMBOL %in% marker_genes$SYMBOL) %>%
  dplyr::left_join(
    marker_genes,
    by = "SYMBOL"
  ) %>%
  dplyr::select(
    count_gene_id,
    SYMBOL,
    marker_state,
    dplyr::all_of(sample_cols)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_cols),
    names_to = "sample_id",
    values_to = "norm_count"
  ) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::mutate(
    gene_scaled_expr = as.numeric(scale(log2(norm_count + 1)))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sample_id, marker_state) %>%
  dplyr::summarise(
    scaled_expr = mean(gene_scaled_expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    feature = dplyr::case_when(
      marker_state == "Naive marker" ~ "Naive marker score",
      marker_state == "Primed marker" ~ "Primed marker score"
    ),
    feature_group = "Marker programme"
  ) %>%
  dplyr::select(sample_id, feature, feature_group, scaled_expr)

### 5.3.2 Build scaled HDAC expression per sample

hdac_scaled_expr_df <- cap_norm_counts_annot %>%
  dplyr::filter(SYMBOL %in% c("HDAC1", "HDAC2", "HDAC3")) %>%
  dplyr::select(
    count_gene_id,
    SYMBOL,
    dplyr::all_of(sample_cols)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_cols),
    names_to = "sample_id",
    values_to = "norm_count"
  ) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::mutate(
    scaled_expr = as.numeric(scale(log2(norm_count + 1)))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    feature = SYMBOL,
    feature_group = "HDAC expression"
  ) %>%
  dplyr::select(sample_id, feature, feature_group, scaled_expr)

### 5.3.3 Combine HDAC expression and marker scores

hdac_marker_score_df <- dplyr::bind_rows(
  hdac_scaled_expr_df,
  marker_score_df
) %>%
  dplyr::left_join(
    cap_sample_meta %>%
      dplyr::select(sample_id, condition_clean),
    by = "sample_id"
  ) %>%
  dplyr::mutate(
    condition_clean = factor(
      condition_clean,
      levels = c("Naive", "Day3", "Day7", "Primed")
    ),
    feature = factor(
      feature,
      levels = c(
        "HDAC1", "HDAC2", "HDAC3",
        "Naive marker score", "Primed marker score"
      )
    ),
    feature_group = factor(
      feature_group,
      levels = c("HDAC expression", "Marker programme")
    )
  )

### 5.3.4 Summarise expression by timepoint

hdac_marker_score_summary <- hdac_marker_score_df %>%
  dplyr::group_by(feature_group, feature, condition_clean) %>%
  dplyr::summarise(
    mean_scaled_expr = mean(scaled_expr, na.rm = TRUE),
    sd_scaled_expr = sd(scaled_expr, na.rm = TRUE),
    .groups = "drop"
  )

### 5.3.5 Plot the combined trajectories

hdac_marker_score_plot <- ggplot(
  hdac_marker_score_summary,
  aes(
    x = condition_clean,
    y = mean_scaled_expr,
    group = feature,
    colour = feature
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ feature_group,
    ncol = 1
  ) +
  theme_bw() +
  labs(
    title = "HDAC expression relative to naive and primed marker programmes",
    subtitle = "Mean row-scaled expression across naive, capacitating and primed states",
    x = "Cell state / timepoint",
    y = "Scaled expression",
    colour = "Gene / marker score"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.text = element_text(face = "bold")
  )

hdac_marker_score_plot

## 5.4 Combined HDAC and marker trajectory plot

### 5.4.1 Define genes to plot
traj_gene_info <- tibble::tribble(
  ~SYMBOL,   ~gene_group,
  "HDAC1",   "HDAC",
  "HDAC2",   "HDAC",
  "HDAC3",   "HDAC",
  "DPPA5",   "Naive marker",
  "DNMT3L",  "Naive marker",
  "TFCP2L1", "Naive marker",
  "KLF5",    "Naive marker",
  "DNMT3B",  "Primed marker",
  "SALL2",   "Primed marker",
  "ZIC2",    "Primed marker",
  "PODXL",   "Primed marker"
)

### 5.4.2 Prepare sample metadata
cap_sample_meta <- as.data.frame(SummarizedExperiment::colData(dds_cap)) %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::select(sample_id, condition_clean) %>%
  dplyr::mutate(
    condition_clean = factor(
      condition_clean,
      levels = c("Naive", "Day3", "Day7", "Primed")
    )
  )

### 5.4.3 Identify sample columns
sample_cols <- colnames(DESeq2::counts(dds_cap))

### 5.4.4 Build the long-format expression table
traj_expr_df <- cap_norm_counts_annot %>%
  dplyr::filter(SYMBOL %in% traj_gene_info$SYMBOL) %>%
  dplyr::select(count_gene_id, SYMBOL, dplyr::all_of(sample_cols)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_cols),
    names_to = "sample_id",
    values_to = "norm_count"
  ) %>%
  dplyr::left_join(traj_gene_info, by = "SYMBOL") %>%
  dplyr::left_join(cap_sample_meta, by = "sample_id")

### 5.4.5 Calculate mean expression by condition
traj_summary <- traj_expr_df %>%
  dplyr::group_by(gene_group, SYMBOL, condition_clean) %>%
  dplyr::summarise(
    mean_norm_count = mean(norm_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::mutate(
    scaled_expr = as.numeric(scale(log2(mean_norm_count + 1)))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    gene_group = factor(
      gene_group,
      levels = c("HDAC", "Naive marker", "Primed marker")
    )
  )

### 5.4.6 Prepare endpoint labels
traj_labels <- traj_summary %>%
  dplyr::filter(condition_clean == "Primed")

### 5.4.7 Plot the trajectories
marker_hdac_trajectory_plot <- ggplot(
  traj_summary,
  aes(
    x = condition_clean,
    y = scaled_expr,
    group = SYMBOL
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_line(
    colour = "grey60",
    linewidth = 0.8
  ) +
  geom_point(
    size = 3,
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 0.8
  ) +
  geom_text(
    data = traj_labels,
    aes(label = SYMBOL),
    hjust = -0.05,
    size = 3.8,
    inherit.aes = FALSE,
    x = 4.05,
    y = traj_labels$scaled_expr
  ) +
  facet_grid(gene_group ~ .) +
  theme_bw() +
  labs(
    title = "Chronological expression trajectories of HDACs and pluripotency markers",
    subtitle = "Naive to primed progression; mean DESeq2-normalised counts, row-scaled within each gene",
    x = "Cell state",
    y = "Row-scaled log2 mean normalised count"
  ) +
  coord_cartesian(clip = "off") +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.18))) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.text.y = element_text(face = "bold"),
    plot.margin = margin(10, 80, 10, 10)
  )

marker_hdac_trajectory_plot

## 5.5 Marker trajectories shown against each HDAC

### 5.5.1 Define marker genes
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

### 5.5.2 Define HDAC reference genes
hdac_genes <- tibble::tribble(
  ~SYMBOL,
  "HDAC1",
  "HDAC2",
  "HDAC3"
)

### 5.5.3 Prepare sample metadata
cap_sample_meta <- as.data.frame(SummarizedExperiment::colData(dds_cap)) %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::select(sample_id, condition_clean) %>%
  dplyr::mutate(
    condition_clean = factor(
      condition_clean,
      levels = c("Naive", "Day3", "Day7", "Primed")
    )
  )

### 5.5.4 Identify sample columns
sample_cols <- colnames(DESeq2::counts(dds_cap))

### 5.5.5 Build the combined expression table
traj_genes_all <- dplyr::bind_rows(
  marker_genes %>% dplyr::mutate(gene_type = "Marker"),
  hdac_genes %>% dplyr::mutate(marker_state = "HDAC", gene_type = "HDAC")
)

traj_expr_df <- cap_norm_counts_annot %>%
  dplyr::filter(SYMBOL %in% traj_genes_all$SYMBOL) %>%
  dplyr::select(count_gene_id, SYMBOL, dplyr::all_of(sample_cols)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_cols),
    names_to = "sample_id",
    values_to = "norm_count"
  ) %>%
  dplyr::left_join(
    traj_genes_all %>% dplyr::select(SYMBOL, marker_state, gene_type),
    by = "SYMBOL"
  ) %>%
  dplyr::left_join(cap_sample_meta, by = "sample_id")

### 5.5.6 Calculate and scale mean expression
traj_summary <- traj_expr_df %>%
  dplyr::group_by(SYMBOL, marker_state, gene_type, condition_clean) %>%
  dplyr::summarise(
    mean_norm_count = mean(norm_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::mutate(
    scaled_expr = as.numeric(scale(log2(mean_norm_count + 1)))
  ) %>%
  dplyr::ungroup()

### 5.5.7 Separate marker and HDAC data
marker_traj_df <- traj_summary %>%
  dplyr::filter(gene_type == "Marker") %>%
  tidyr::crossing(
    focal_hdac = factor(c("HDAC1", "HDAC2", "HDAC3"),
                        levels = c("HDAC1", "HDAC2", "HDAC3"))
  ) %>%
  dplyr::mutate(
    marker_state = factor(
      marker_state,
      levels = c("Naive marker", "Primed marker")
    )
  )

hdac_ref_df <- traj_summary %>%
  dplyr::filter(gene_type == "HDAC") %>%
  dplyr::rename(focal_hdac = SYMBOL) %>%
  dplyr::mutate(
    focal_hdac = factor(
      focal_hdac,
      levels = c("HDAC1", "HDAC2", "HDAC3")
    )
  )

### 5.5.8 Prepare endpoint labels
marker_label_df <- marker_traj_df %>%
  dplyr::filter(condition_clean == "Primed")

### 5.5.9 Plot marker trajectories
marker_state_cols <- c(
  "Naive marker"  = "#F8766D",
  "Primed marker" = "#C77CFF"
)

marker_traj_by_hdac_plot <- ggplot() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey55"
  ) +
  geom_line(
    data = marker_traj_df,
    aes(
      x = condition_clean,
      y = scaled_expr,
      group = SYMBOL,
      colour = marker_state
    ),
    linewidth = 0.8,
    alpha = 0.9
  ) +
  geom_point(
    data = marker_traj_df,
    aes(
      x = condition_clean,
      y = scaled_expr,
      colour = marker_state
    ),
    size = 2.3
  ) +
  geom_line(
    data = hdac_ref_df,
    aes(
      x = condition_clean,
      y = scaled_expr,
      group = focal_hdac,
      linetype = "HDAC trajectory"
    ),
    colour = "black",
    linewidth = 1
  ) +
  geom_point(
    data = hdac_ref_df,
    aes(
      x = condition_clean,
      y = scaled_expr
    ),
    colour = "black",
    size = 2.2
  ) +
  geom_text(
    data = marker_label_df,
    aes(
      x = 4.08,
      y = scaled_expr,
      label = SYMBOL,
      colour = marker_state
    ),
    hjust = 0,
    size = 3.4,
    show.legend = FALSE
  ) +
  facet_wrap(~ focal_hdac, nrow = 1) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.22))) +
  scale_colour_manual(values = marker_state_cols) +
  scale_linetype_manual(
    values = c("HDAC trajectory" = "longdash")
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  labs(
    title = "Trajectory of pluripotency markers across capacitation",
    subtitle = "Each panel shows the 8 marker genes",
    x = "Cell state",
    y = "Row-scaled log2 mean normalised count",
    colour = NULL,
    linetype = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.margin = margin(10, 100, 10, 10)
  )

marker_traj_by_hdac_plot

ggsave(
  file.path(
    cap_plot_dir,
    "pluripotency_marker_trajectory_by_HDAC_colorchange.png"
  ),
  marker_traj_by_hdac_plot,
  width = 14,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

## 5.6 Marker-expression plots with DESeq2 statistics

### 5.6.1 Define condition colours
cap_condition_cols <- c(
  "Naive"  = "#F8766D",
  "Day3"   = "#7CAE00",
  "Day7"   = "#00BFC4",
  "Primed" = "#C77CFF"
)

### 5.6.2 Define marker-type colours
marker_type_cols <- c(
  "Naive marker"  = "#F8766D",
  "Primed marker" = "#C77CFF"
)

### 5.6.3 Define facet order and labels

marker_facet_levels <- paste0(
  marker_genes$SYMBOL,
  "\n",
  marker_genes$marker_state
)

marker_expr_df <- marker_expr_df %>%
  dplyr::mutate(
    condition_clean = factor(
      condition_clean,
      levels = c("Naive", "Day3", "Day7", "Primed")
    ),
    marker_state = factor(
      marker_state,
      levels = c("Naive marker", "Primed marker")
    ),
    SYMBOL = factor(
      SYMBOL,
      levels = marker_genes$SYMBOL
    ),
    facet_label = paste0(
      as.character(SYMBOL),
      "\n",
      as.character(marker_state)
    ),
    facet_label = factor(
      facet_label,
      levels = marker_facet_levels
    )
  )

### 5.6.4 Extract marker-gene DESeq2 statistics

if (exists("all_res_capacitation_15fold")) {
  cap_res_for_stats <- all_res_capacitation_15fold
} else if (exists("cap_res")) {
  cap_res_for_stats <- cap_res
} else {
  stop("Could not find all_res_capacitation_15fold or cap_res")
}

marker_deseq_stats <- cap_res_for_stats %>%
  dplyr::filter(
    SYMBOL %in% marker_genes$SYMBOL,
    contrast_short %in% c("Day3 vs Naive", "Day7 vs Naive", "Primed vs Naive")
  ) %>%
  dplyr::left_join(
    marker_genes,
    by = "SYMBOL"
  ) %>%
  dplyr::mutate(
    condition_clean = dplyr::case_when(
      contrast_short == "Day3 vs Naive" ~ "Day3",
      contrast_short == "Day7 vs Naive" ~ "Day7",
      contrast_short == "Primed vs Naive" ~ "Primed",
      TRUE ~ NA_character_
    ),
    condition_clean = factor(
      condition_clean,
      levels = c("Naive", "Day3", "Day7", "Primed")
    ),
    marker_state = factor(
      marker_state,
      levels = c("Naive marker", "Primed marker")
    ),
    SYMBOL = factor(
      SYMBOL,
      levels = marker_genes$SYMBOL
    ),
    facet_label = paste0(
      as.character(SYMBOL),
      "\n",
      as.character(marker_state)
    ),
    facet_label = factor(
      facet_label,
      levels = marker_facet_levels
    ),
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
  dplyr::filter(!is.na(condition_clean))

### 5.6.5 Calculate label positions

marker_label_positions <- marker_expr_df %>%
  dplyr::group_by(SYMBOL, facet_label) %>%
  dplyr::summarise(
    y_max = max(norm_count, na.rm = TRUE),
    y_min = min(norm_count, na.rm = TRUE),
    y_range = y_max - y_min,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    y_pos = y_max + 0.15 * y_range
  )

marker_deseq_stats <- marker_deseq_stats %>%
  dplyr::left_join(
    marker_label_positions,
    by = c("SYMBOL", "facet_label")
  )

### 5.6.6 Define facet-strip colours

strip_backgrounds <- marker_genes %>%
  dplyr::mutate(
    facet_label = paste0(SYMBOL, "\n", marker_state),
    strip_fill = marker_type_cols[as.character(marker_state)]
  ) %>%
  dplyr::arrange(
    match(facet_label, marker_facet_levels)
  ) %>%
  dplyr::pull(strip_fill)

strip_text_cols <- rep("black", length(strip_backgrounds))

### 5.6.7 Plot marker expression and statistics

pluripotency_marker_boxplot_coloured_stats <- ggplot(
  marker_expr_df,
  aes(
    x = condition_clean,
    y = norm_count,
    fill = condition_clean
  )
) +
  geom_boxplot(
    colour = "black",
    outlier.shape = NA,
    alpha = 0.75,
    width = 0.65
  ) +
  geom_jitter(
    aes(colour = condition_clean),
    width = 0.12,
    size = 2,
    alpha = 0.85
  ) +
  geom_text(
    data = marker_deseq_stats,
    aes(
      x = condition_clean,
      y = y_pos,
      label = stat_label
    ),
    inherit.aes = FALSE,
    size = 2.4,
    lineheight = 0.85
  ) +
  ggh4x::facet_wrap2(
    ~ facet_label,
    scales = "free_y",
    ncol = 4,
    strip = ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = strip_backgrounds,
        colour = "black"
      ),
      text_x = ggh4x::elem_list_text(
        colour = strip_text_cols,
        face = "bold"
      )
    )
  ) +
  scale_fill_manual(values = cap_condition_cols) +
  scale_colour_manual(values = cap_condition_cols) +
  theme_bw() +
  labs(
    title = "Naive and primed pluripotency marker expression during capacitation",
    subtitle = "DESeq2-normalised counts across Naive, Day3, Day7 and Primed states; DESeq2 statistics shown relative to Naive",
    x = "Cell state / timepoint",
    y = "DESeq2-normalised count",
    fill = "Condition",
    colour = "Condition"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

pluripotency_marker_boxplot_coloured_stats

ggsave(
  file.path(
    cap_plot_dir,
    "pluripotency_marker_expression_during_capacitation_coloured_stats.png"
  ),
  pluripotency_marker_boxplot_coloured_stats,
  width = 18,
  height = 8,
  dpi = 300,
  bg = "white"
)

readr::write_csv(
  marker_deseq_stats,
  file.path(
    cap_plot_dir,
    "pluripotency_marker_DESeq2_stats_vs_Naive.csv"
  )
)

# 6. Diagnostic summaries
## 6.1 HDAC expression summaries
hdac_cap_expr %>%
  dplyr::group_by(SYMBOL, condition_clean) %>%
  dplyr::summarise(
    mean_norm_count = mean(norm_count, na.rm = TRUE),
    sd_norm_count = sd(norm_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(n = Inf)

## 6.2 Marker-expression summaries
marker_expr_df %>%
  dplyr::group_by(marker_state, SYMBOL, condition_clean) %>%
  dplyr::summarise(
    mean_norm_count = mean(norm_count, na.rm = TRUE),
    sd_norm_count = sd(norm_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(n = Inf)

## 6.3 X-linked DEG totals
cap_x_deg_summary %>%
  dplyr::filter(
    contrast_short %in% c("Day3 vs Naive", "Day7 vs Naive", "Primed vs Naive"),
    deg_status_15fold != "Not significant"
  ) %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::summarise(
    total_X_DEGs = sum(n_X_genes),
    .groups = "drop"
  )

cap_x_deg_summary_naive <- cap_x_deg_summary %>%
  dplyr::filter(
    as.character(contrast_short) %in% c(
      "Day3 vs Naive",
      "Day7 vs Naive",
      "Primed vs Naive"
    )
  ) %>%
  dplyr::arrange(contrast_short, deg_status_15fold)

tibble::as_tibble(cap_x_deg_summary_naive) %>%
  print(n = Inf)
cap_x_deg_totals_naive <- cap_x_deg_summary_naive %>%
  dplyr::filter(deg_status_15fold != "Not significant") %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::summarise(
    total_X_DEGs = sum(n_X_genes),
    .groups = "drop"
  )

tibble::as_tibble(cap_x_deg_totals_naive) %>%
  print(n = Inf)

## 6.4 HDAC DESeq2 statistics

hdac_cap_deseq_stats <- all_res_capacitation_15fold %>%
  dplyr::filter(
    SYMBOL %in% c("HDAC1", "HDAC2", "HDAC3"),
    as.character(contrast_short) %in% c(
      "Day3 vs Naive",
      "Day7 vs Naive",
      "Primed vs Naive"
    )
  ) %>%
  dplyr::mutate(
    significance_label = dplyr::case_when(
      !is.na(padj) &
        padj < 0.05 &
        abs(log2FoldChange) >= log2(1.5) ~ "Significant by DEG threshold",
      !is.na(padj) &
        padj < 0.05 ~ "padj < 0.05 only",
      TRUE ~ "Not significant"
    )
  ) %>%
  dplyr::select(
    contrast_short,
    SYMBOL,
    baseMean,
    log2FoldChange,
    pvalue,
    padj,
    deg_status_15fold,
    significance_label
  ) %>%
  dplyr::arrange(SYMBOL, contrast_short)

tibble::as_tibble(hdac_cap_deseq_stats) %>%
  print(n = Inf, width = Inf)

## 6.5 Marker-gene DESeq2 statistics

marker_genes_cap <- c(
  "DPPA5", "DNMT3L", "TFCP2L1", "KLF5",
  "DNMT3B", "SALL2", "ZIC2", "PODXL"
)

marker_state_lookup <- tibble::tibble(
  SYMBOL = marker_genes_cap,
  marker_state = c(
    rep("Naive marker", 4),
    rep("Primed marker", 4)
  )
)

marker_cap_deseq_stats <- all_res_capacitation_15fold %>%
  dplyr::filter(
    SYMBOL %in% marker_genes_cap,
    as.character(contrast_short) %in% c(
      "Day3 vs Naive",
      "Day7 vs Naive",
      "Primed vs Naive"
    )
  ) %>%
  dplyr::left_join(marker_state_lookup, by = "SYMBOL") %>%
  dplyr::mutate(
    significance_label = dplyr::case_when(
      !is.na(padj) &
        padj < 0.05 &
        abs(log2FoldChange) >= log2(1.5) ~ "Significant by DEG threshold",
      !is.na(padj) &
        padj < 0.05 ~ "padj < 0.05 only",
      TRUE ~ "Not significant"
    )
  ) %>%
  dplyr::select(
    marker_state,
    SYMBOL,
    contrast_short,
    baseMean,
    log2FoldChange,
    pvalue,
    padj,
    deg_status_15fold,
    significance_label
  ) %>%
  dplyr::arrange(marker_state, SYMBOL, contrast_short)

tibble::as_tibble(marker_cap_deseq_stats) %>%
  print(n = Inf, width = Inf)

marker_cap_summary_counts <- marker_cap_deseq_stats %>%
  dplyr::count(
    marker_state,
    contrast_short,
    deg_status_15fold
  ) %>%
  tidyr::complete(
    marker_state,
    contrast_short,
    deg_status_15fold = c("Down", "Not significant", "Up"),
    fill = list(n = 0)
  ) %>%
  dplyr::arrange(marker_state, contrast_short, deg_status_15fold)

tibble::as_tibble(marker_cap_summary_counts) %>%
  print(n = Inf)

marker_cap_deseq_stats %>%
  dplyr::arrange(padj) %>%
  dplyr::select(
    marker_state,
    SYMBOL,
    contrast_short,
    log2FoldChange,
    padj,
    deg_status_15fold
  ) %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

## 6.6 Trajectory direction summaries

traj_summary %>%
  dplyr::filter(SYMBOL %in% c(
    "HDAC1", "HDAC2", "HDAC3",
    "DPPA5", "DNMT3L", "TFCP2L1", "KLF5",
    "DNMT3B", "SALL2", "ZIC2", "PODXL"
  )) %>%
  dplyr::select(SYMBOL, gene_type, condition_clean, mean_norm_count, scaled_expr) %>%
  dplyr::arrange(gene_type, SYMBOL, condition_clean) %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

traj_direction_summary <- traj_summary %>%
  dplyr::filter(SYMBOL %in% c(
    "HDAC1", "HDAC2", "HDAC3",
    "DPPA5", "DNMT3L", "TFCP2L1", "KLF5",
    "DNMT3B", "SALL2", "ZIC2", "PODXL"
  )) %>%
  dplyr::filter(condition_clean %in% c("Naive", "Primed")) %>%
  dplyr::select(SYMBOL, gene_type, condition_clean, scaled_expr) %>%
  tidyr::pivot_wider(
    names_from = condition_clean,
    values_from = scaled_expr
  ) %>%
  dplyr::mutate(
    primed_minus_naive_scaled = Primed - Naive,
    trajectory_direction = dplyr::case_when(
      primed_minus_naive_scaled > 0 ~ "Higher in primed",
      primed_minus_naive_scaled < 0 ~ "Lower in primed",
      TRUE ~ "No change"
    )
  ) %>%
  dplyr::arrange(gene_type, SYMBOL)

tibble::as_tibble(traj_direction_summary) %>%
  print(n = Inf, width = Inf)

# 7. Genome-wide MA plots
# Comparisons: Day3 vs Naive, Day7 vs Naive, and Primed vs Naive.

## 7.1 Load DESeq2 results

# Reuse the in-memory result object when available; otherwise load the saved table.

if (exists("all_res_capacitation_15fold")) {

  cap_res <- all_res_capacitation_15fold

} else {

  cap_res_path <- "/home/jvk3/Desktop/HDAC_counts/Capacitation/Capacitation_DESeq2_outputs/all_res_capacitation_15fold.csv"

  cap_res <- readr::read_csv(
    cap_res_path,
    show_col_types = FALSE
  )
}

## 7.2 Prepare the MA-plot data

log2fc_cutoff <- log2(1.5)

cap_ma_df <- cap_res %>%
  dplyr::mutate(
    contrast_short = as.character(contrast_short),

    ## Standardise contrast order
    contrast_short = dplyr::case_when(
      contrast_short %in% c("Day3 vs Naive", "Day3_vs_Naive") ~ "Day3 vs Naive",
      contrast_short %in% c("Day7 vs Naive", "Day7_vs_Naive") ~ "Day7 vs Naive",
      contrast_short %in% c("Primed vs Naive", "Primed_vs_Naive") ~ "Primed vs Naive",
      TRUE ~ contrast_short
    ),

    contrast_short = factor(
      contrast_short,
      levels = c("Day3 vs Naive", "Day7 vs Naive", "Primed vs Naive")
    ),

    ## Recalculate DEG class, just to make sure it is consistent
    deg_status_15fold = dplyr::case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange >= log2fc_cutoff ~ "Up",
      !is.na(padj) & padj < 0.05 & log2FoldChange <= -log2fc_cutoff ~ "Down",
      TRUE ~ "Not significant"
    ),

    deg_status_15fold = factor(
      deg_status_15fold,
      levels = c("Down", "Not significant", "Up")
    ),

    log10_baseMean = log10(baseMean + 1)
  ) %>%
  dplyr::filter(
    !is.na(contrast_short),
    !is.na(baseMean),
    !is.na(log2FoldChange)
  )

## 7.3 Check DEG counts

cap_ma_df %>%
  dplyr::count(contrast_short, deg_status_15fold)

## 7.4 Plot genome-wide expression changes

p_cap_ma <- ggplot(
  cap_ma_df,
  aes(x = log10_baseMean, y = log2FoldChange)
) +
  geom_point(
    data = cap_ma_df %>% dplyr::filter(deg_status_15fold == "Not significant"),
    colour = "grey75",
    alpha = 0.45,
    size = 0.7
  ) +
  geom_point(
    data = cap_ma_df %>% dplyr::filter(deg_status_15fold == "Down"),
    colour = "red",
    alpha = 0.7,
    size = 0.8
  ) +
  geom_point(
    data = cap_ma_df %>% dplyr::filter(deg_status_15fold == "Up"),
    colour = "forestgreen",
    alpha = 0.7,
    size = 0.8
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "solid",
    colour = "grey40",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = c(-log2fc_cutoff, log2fc_cutoff),
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.4
  ) +
  facet_wrap(~ contrast_short, ncol = 1) +
  theme_bw() +
  labs(
    title = "Capacitation bulk RNA-seq MA plots",
    subtitle = "Differential expression relative to Naive; green = upregulated, red = downregulated",
    x = "log10 mean normalised expression",
    y = "log2 fold change"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p_cap_ma

## 7.5 Add X-linked gene labels and save

cap_ma_label_df <- cap_ma_df %>%
  dplyr::filter(
    deg_status_15fold %in% c("Up", "Down"),
    CHR == "X"
  ) %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::slice_max(
    order_by = abs(log2FoldChange),
    n = 5,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

xist_xact_labels <- cap_ma_df %>%
  dplyr::filter(
    SYMBOL %in% c("XIST", "XACT")
  )

cap_ma_label_df <- dplyr::bind_rows(
  cap_ma_label_df,
  xist_xact_labels
) %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::distinct(contrast_short, SYMBOL, .keep_all = TRUE)

p_cap_ma_labelled <- p_cap_ma +
  ggrepel::geom_text_repel(
    data = cap_ma_label_df,
    aes(label = SYMBOL),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.3,
    point.padding = 0.2
  )

p_cap_ma_labelled

ggsave(
  file.path(cap_plot_dir, "Capacitation_bulk_MA_plots_vs_Naive_labelled.png"),
  p_cap_ma_labelled,
  width = 8,
  height = 10,
  dpi = 300
)

# 8. Keshet-highlighted pluripotency and differentiation genes
## 8.1 Mean-expression heatmap across capacitation

keshet_plot_dir <- file.path(cap_plot_dir, "Keshet_pluripotency_gene_heatmap")
dir.create(keshet_plot_dir, showWarnings = FALSE, recursive = TRUE)

### 8.1.1 Define the focused gene panel

keshet_gene_panel <- tibble::tribble(
  ~SYMBOL,    ~role_group,                         ~role_note,
  "POU5F1",   "ESC / pluripotency regulators",      "Core pluripotency; also linked to lineage exit",
  "SOX2",     "ESC / pluripotency regulators",      "Core pluripotency; mesendoderm bias in Keshet",
  "NANOG",    "ESC / pluripotency regulators",      "Core pluripotency; ectoderm bias in Keshet",
  "PRDM14",   "ESC / pluripotency regulators",      "Pluripotency regulator; ectoderm bias in Keshet",
  "SALL4",    "ESC / pluripotency regulators",      "Pluripotency maintenance",
  "ETV4",     "ESC / pluripotency regulators",      "Ambiguous pluripotency regulatory effect",
  "ZNF281",   "ESC / pluripotency regulators",      "Ambiguous pluripotency regulatory effect",

  "SMARCC1",  "Pluripotency attenuators",           "Attenuator / formative-state related",
  "TGIF1",    "Pluripotency attenuators",           "Attenuator of primed pluripotency GRN",
  "MIS18BP1", "Pluripotency attenuators",           "Putative attenuator",

  "ZNF462",   "Lineage-associated regulators",      "Mesendoderm-associated regulator",
  "REST",     "Lineage-associated regulators",      "Neural-associated regulator"
)

keshet_gene_panel

### 8.1.2 Prepare sample metadata

cap_sample_meta <- as.data.frame(SummarizedExperiment::colData(dds_cap)) %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::select(sample_id, condition_clean) %>%
  dplyr::mutate(
    condition_clean = factor(
      condition_clean,
      levels = c("Naive", "Day3", "Day7", "Primed")
    )
  )

sample_cols <- colnames(DESeq2::counts(dds_cap))

### 8.1.3 Extract gene expression

keshet_expr_long <- cap_norm_counts_annot %>%
  dplyr::filter(SYMBOL %in% keshet_gene_panel$SYMBOL) %>%
  dplyr::select(count_gene_id, SYMBOL, dplyr::all_of(sample_cols)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_cols),
    names_to = "sample_id",
    values_to = "norm_count"
  ) %>%
  dplyr::group_by(SYMBOL, sample_id) %>%
  dplyr::summarise(
    norm_count = sum(norm_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(keshet_gene_panel, by = "SYMBOL") %>%
  dplyr::left_join(cap_sample_meta, by = "sample_id")

# Check whether any selected genes are absent from the expression data.
missing_keshet_genes <- setdiff(
  keshet_gene_panel$SYMBOL,
  unique(keshet_expr_long$SYMBOL)
)

missing_keshet_genes

### 8.1.4 Calculate and scale mean expression

keshet_heatmap_df <- keshet_expr_long %>%
  dplyr::filter(!is.na(condition_clean)) %>%
  dplyr::group_by(SYMBOL, role_group, role_note, condition_clean) %>%
  dplyr::summarise(
    mean_norm_count = mean(norm_count, na.rm = TRUE),
    mean_log2_count = log2(mean_norm_count + 1),
    .groups = "drop"
  ) %>%
  dplyr::group_by(SYMBOL) %>%
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

### 8.1.5 Plot the mean-expression heatmap

keshet_cap_heatmap <- ggplot(
  keshet_heatmap_df,
  aes(
    x = condition_clean,
    y = SYMBOL,
    fill = row_scaled_expr
  )
) +
  geom_tile(
    colour = "grey85",
    linewidth = 0.3
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
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  theme_bw() +
  labs(
    title = "Keshet-highlighted pluripotency and differentiation regulators during capacitation",
    subtitle = "Selected genes from the Keshet et al. mode-of-action summary; row-scaled log2 normalised counts",
    x = "Capacitation stage",
    y = NULL,
    fill = "Row-scaled\nlog2 count"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text.y = element_text(face = "bold", angle = 0),
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.6, "lines")
  )

keshet_cap_heatmap

ggsave(
  filename = file.path(
    keshet_cap_log2fc_dir,
    "Keshet_pluripotency_differentiation_regulators_capacitation_heatmap.png"
  ),
  plot = keshet_cap_heatmap,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

## 8.2 DESeq2 log2FC heatmap
# Summarise DESeq2 fold changes for the focused Keshet gene panel.

keshet_cap_log2fc_dir <- file.path(
  cap_plot_dir,
  "Keshet_pluripotency_gene_log2FC_heatmap"
)

dir.create(
  keshet_cap_log2fc_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

### 8.2.1 Define the focused gene panel

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

### 8.2.2 Select the capacitation DESeq2 result table

if (exists("all_res_capacitation_15fold")) {
  cap_res_for_keshet <- all_res_capacitation_15fold
} else if (exists("cap_res")) {
  cap_res_for_keshet <- cap_res
} else {
  all_res_capacitation_15fold <- readr::read_csv(
    "/home/jvk3/Desktop/HDAC_counts/Capacitation/Capacitation_outputs/all_results_capacitation_15fold.csv",
    show_col_types = FALSE
  )
  cap_res_for_keshet <- all_res_capacitation_15fold
}

### 8.2.3 Prepare the log2FC heatmap table

keshet_cap_log2fc_df <- cap_res_for_keshet %>%
  dplyr::filter(
    SYMBOL %in% keshet_gene_panel$SYMBOL,
    contrast_short %in% c(
      "Day3 vs Naive",
      "Day7 vs Naive",
      "Primed vs Naive"
    )
  ) %>%
  dplyr::left_join(
    keshet_gene_panel,
    by = "SYMBOL"
  ) %>%
  dplyr::mutate(
    contrast_short = factor(
      contrast_short,
      levels = c(
        "Day3 vs Naive",
        "Day7 vs Naive",
        "Primed vs Naive"
      )
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

    ## Full DEG threshold
    is_sig = !is.na(padj) &
      padj < 0.05 &
      abs(log2FoldChange) >= log2(1.5),

    ## Cap colour scale to make moderate changes visible
    log2FC_plot = pmax(
      pmin(log2FoldChange, 2),
      -2
    )
  )

# Check whether any selected genes are absent from the DESeq2 results.
missing_keshet_cap_genes <- setdiff(
  keshet_gene_panel$SYMBOL,
  unique(keshet_cap_log2fc_df$SYMBOL)
)

missing_keshet_cap_genes

# Review which selected genes pass the full DEG threshold.
keshet_cap_log2fc_df %>%
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

### 8.2.4 Plot the log2FC heatmap

keshet_cap_log2fc_heatmap <- ggplot(
  keshet_cap_log2fc_df,
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
    data = keshet_cap_log2fc_df %>%
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
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  theme_bw() +
  labs(
    title = "Keshet-highlighted pluripotency and differentiation regulators during capacitation",
    subtitle = "DESeq2 log2 fold change relative to Naive; black dots indicate padj < 0.05 and |log2FC| ≥ log2(1.5)",
    x = "Capacitation comparison",
    y = NULL,
    fill = "log2FC"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text.y = element_text(face = "bold", angle = 0),
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.6, "lines")
  )

keshet_cap_log2fc_heatmap

ggsave(
  filename = file.path(
    keshet_cap_log2fc_dir,
    "Keshet_pluripotency_differentiation_regulators_capacitation_log2FC_heatmap.png"
  ),
  plot = keshet_cap_log2fc_heatmap,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

readr::write_csv(
  keshet_cap_log2fc_df,
  file.path(
    keshet_cap_log2fc_dir,
    "Keshet_pluripotency_differentiation_regulators_capacitation_log2FC_heatmap_table.csv"
  )
)

# 9. Chromosome X position and XCI-category analyses
# Includes all tested X-linked genes with Tukiainen XCI-category annotation.
## 9.1 Create the output folder

cap_xchr_plot_dir <- file.path(
  cap_plot_dir,
  "X_chromosome_position_plots"
)

dir.create(
  cap_xchr_plot_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

## 9.2 Define the XIST/XIC position

xist_start_bp <- 73817774
xist_end_bp   <- 73852754
xist_mid_mb   <- mean(c(xist_start_bp, xist_end_bp)) / 1e6

## 9.3 Join Tukiainen XCI annotations

tukiainen_path <- "/home/jvk3/Desktop/HDAC_counts/XCI_classifications/XCI_classification_outputs/Tukiainen_XCI_categories_clean.csv"

tukiainen_annot <- readr::read_csv(
  tukiainen_path,
  show_col_types = FALSE
) %>%
  dplyr::select(SYMBOL, tukiainen_status_grouped) %>%
  dplyr::distinct()

cap_xchr_plot_df <- all_res_capacitation_15fold %>%
  dplyr::left_join(
    tukiainen_annot,
    by = "SYMBOL"
  ) %>%
  dplyr::filter(
    CHR == "X",
    contrast_short %in% c(
      "Day3 vs Naive",
      "Day7 vs Naive",
      "Primed vs Naive"
    ),
    !is.na(START),
    !is.na(END),
    is.finite(START),
    is.finite(END),
    !is.na(log2FoldChange)
  ) %>%
  dplyr::mutate(
    gene_mid_mb = ((START + END) / 2) / 1e6,

    xci_plot_group = dplyr::case_when(
      tukiainen_status_grouped == "Subject to XCI" ~ "Subject to XCI",
      tukiainen_status_grouped == "Escapee" ~ "Escapee",
      TRUE ~ "Other"
    ),

    xci_plot_group = factor(
      xci_plot_group,
      levels = c("Subject to XCI", "Escapee", "Other")
    ),

    contrast_short = factor(
      contrast_short,
      levels = c(
        "Day3 vs Naive",
        "Day7 vs Naive",
        "Primed vs Naive"
      )
    ),

    is_sig = !is.na(padj) &
      padj < 0.05 &
      abs(log2FoldChange) >= log2(1.5),

    gene_plot_label = dplyr::case_when(
      !is.na(SYMBOL) & SYMBOL != "" ~ SYMBOL,
      !is.na(gene_label) & gene_label != "" ~ gene_label,
      TRUE ~ gene_id
    )
  )

## 9.4 Select genes to label
# Top five significant genes and top three significant escapees per comparison.

### 9.4.1 Select the top significant X-linked genes
cap_xchr_top_sig_label_df <- cap_xchr_plot_df %>%
  dplyr::filter(is_sig) %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::slice_max(
    order_by = abs(log2FoldChange),
    n = 5,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

### 9.4.2 Select the top significant escapees
cap_xchr_escapee_label_df <- cap_xchr_plot_df %>%
  dplyr::filter(
    is_sig,
    xci_plot_group == "Escapee"
  ) %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::slice_max(
    order_by = abs(log2FoldChange),
    n = 3,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

### 9.4.3 Combine labels and remove duplicates
cap_xchr_label_df <- dplyr::bind_rows(
  cap_xchr_top_sig_label_df,
  cap_xchr_escapee_label_df
) %>%
  dplyr::distinct(
    contrast_short,
    gene_plot_label,
    .keep_all = TRUE
  )

## 9.5 Plot X-linked genes across chromosome X

cap_xchr_bulk_plot <- ggplot(
  cap_xchr_plot_df,
  aes(
    x = gene_mid_mb,
    y = log2FoldChange,
    colour = xci_plot_group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_vline(
    xintercept = xist_mid_mb,
    linetype = "dashed",
    colour = "black"
  ) +
  geom_point(
    aes(
      alpha = is_sig,
      size = is_sig
    )
  ) +
  ggrepel::geom_text_repel(
    data = cap_xchr_label_df,
    aes(label = gene_plot_label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.2,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = xist_mid_mb,
    y = Inf,
    label = "XIST/XIC",
    angle = 90,
    vjust = 1.4,
    hjust = 1.1,
    size = 3
  ) +
  facet_wrap(
    ~ contrast_short,
    ncol = 1
  ) +
  scale_colour_manual(
    values = c(
      "Subject to XCI" = "#2166AC",
      "Escapee" = "#D55E00",
      "Other" = "grey65"
    )
  ) +
  scale_alpha_manual(
    values = c(
      "FALSE" = 0.35,
      "TRUE" = 0.9
    ),
    guide = "none"
  ) +
  scale_size_manual(
    values = c(
      "FALSE" = 1.3,
      "TRUE" = 2.3
    ),
    guide = "none"
  ) +
  theme_bw() +
  labs(
    title = "X-linked gene expression changes across chromosome X during capacitation",
    subtitle = "All tested X-linked genes are shown; colours indicate Tukiainen XCI category grouping",
    x = "Position on chromosome X (Mb)",
    y = "DESeq2 log2 fold change relative to Naive",
    colour = "XCI category"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

cap_xchr_bulk_plot

## 9.6 Save the chromosome plot and table

ggsave(
  filename = file.path(
    cap_xchr_plot_dir,
    "Capacitation_bulk_Xchr_all_X_genes_Tukiainen_categories.png"
  ),
  plot = cap_xchr_bulk_plot,
  width = 9,
  height = 8,
  dpi = 300,
  bg = "white"
)

readr::write_csv(
  cap_xchr_plot_df,
  file.path(
    cap_xchr_plot_dir,
    "Capacitation_bulk_Xchr_all_X_genes_Tukiainen_categories_table.csv"
  )
)

## 9.7 Statistical summaries
cap_xchr_plot_df <- cap_xchr_plot_df %>%
  dplyr::mutate(
    distance_to_xist_mb = abs(gene_mid_mb - xist_mid_mb),
    deg_direction = dplyr::case_when(
      is_sig & log2FoldChange > 0 ~ "Up",
      is_sig & log2FoldChange < 0 ~ "Down",
      TRUE ~ "Not significant"
    ),
    deg_direction = factor(
      deg_direction,
      levels = c("Down", "Not significant", "Up")
    )
  )

cap_xchr_category_counts <- cap_xchr_plot_df %>%
  dplyr::count(
    contrast_short,
    xci_plot_group,
    deg_direction,
    name = "n_genes"
  ) %>%
  tidyr::complete(
    contrast_short,
    xci_plot_group,
    deg_direction,
    fill = list(n_genes = 0)
  ) %>%
  dplyr::group_by(contrast_short, xci_plot_group) %>%
  dplyr::mutate(
    total_in_category = sum(n_genes),
    percent_in_category = round((n_genes / total_in_category) * 100, 2)
  ) %>%
  dplyr::ungroup()

cap_xchr_category_counts

cap_xchr_category_fisher <- cap_xchr_plot_df %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::group_modify(~{
    test_table <- table(.x$xci_plot_group, .x$is_sig)

    fisher_result <- fisher.test(test_table)

    tibble::tibble(
      test = "XCI category vs DEG status",
      p_value = fisher_result$p.value
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    padj = p.adjust(p_value, method = "BH")
  )

cap_xchr_category_fisher

cap_xchr_one_vs_rest_fisher <- cap_xchr_plot_df %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::group_modify(~{

    purrr::map_dfr(
      c("Subject to XCI", "Escapee", "Other"),
      function(category_name) {

        n_sig_group <- sum(.x$xci_plot_group == category_name & .x$is_sig)
        n_not_group <- sum(.x$xci_plot_group == category_name & !.x$is_sig)
        n_sig_rest <- sum(.x$xci_plot_group != category_name & .x$is_sig)
        n_not_rest <- sum(.x$xci_plot_group != category_name & !.x$is_sig)

        test_matrix <- matrix(
          c(
            n_sig_group, n_not_group,
            n_sig_rest, n_not_rest
          ),
          nrow = 2,
          byrow = TRUE,
          dimnames = list(
            c(category_name, "Rest"),
            c("Significant", "Not significant")
          )
        )

        fisher_result <- fisher.test(test_matrix)

        tibble::tibble(
          xci_plot_group = category_name,
          n_sig_group = n_sig_group,
          n_total_group = n_sig_group + n_not_group,
          odds_ratio = unname(fisher_result$estimate),
          p_value = fisher_result$p.value
        )
      }
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::mutate(
    padj = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::ungroup()

cap_xchr_one_vs_rest_fisher

cap_xchr_distance_stats <- cap_xchr_plot_df %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::summarise(
    n_sig = sum(is_sig, na.rm = TRUE),
    n_not_sig = sum(!is_sig, na.rm = TRUE),
    median_distance_sig_mb = median(distance_to_xist_mb[is_sig], na.rm = TRUE),
    median_distance_not_sig_mb = median(distance_to_xist_mb[!is_sig], na.rm = TRUE),
    p_value = wilcox.test(distance_to_xist_mb ~ is_sig)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    padj = p.adjust(p_value, method = "BH")
  )

cap_xchr_distance_stats
