## HDAC1/2 PCA WITHOUT HDAC12_H3 outlier
## Separate HDAC1/2 experiment analysis


## 1. Libraries

library(DESeq2)
library(apeglm)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(ggVennDiagram)


## 2. Paths

base_dir <- "/home/jvk3/Desktop/HDAC_counts"

hdac12_counts_dir <- file.path(base_dir,"2_STAR_libinorm_counts")

plot_dir <- file.path(base_dir, "plots")

plot_dir_hdac12_noH3 <- file.path(plot_dir,"HDAC12_no_H3_outputs")

dir.create(plot_dir_hdac12_noH3, showWarnings = FALSE, recursive = TRUE)


## 3. Create HDAC1/2 metadata

meta_hdac12 <- tibble::tibble(
  sample_id = c("A1", "A2", "A3", "H1", "H2", "H3", "B1", "B2", "B3"),
  condition = c("Scramble", "Scramble", "Scramble","HDAC1_siRNA", "HDAC1_siRNA", "HDAC1_siRNA",
    "HDAC2_siRNA", "HDAC2_siRNA", "HDAC2_siRNA"
  ),
  rep = c(1, 2, 3, 1, 2, 3, 1, 2, 3),
  experiment = "HDAC12"
) %>%
  dplyr::mutate(
    unique_sample = paste0("HDAC12_", sample_id)
  ) %>%
  as.data.frame()

rownames(meta_hdac12) <- meta_hdac12$unique_sample

meta_hdac12


## 4. Read HDAC1/2 count files

read_count <- function(sample_id, unique_sample) {
  
  file <- list.files(
    hdac12_counts_dir,
    pattern = paste0("^", sample_id, "_.*_counts\\.txt$"),
    full.names = TRUE
  )
  
  if (length(file) != 1) {
    stop(paste("Problem finding count file for", sample_id))
  }
  
  message(sample_id, " -> ", basename(file))
  
  read.delim(file, header = FALSE) %>%
    dplyr::select(gene_id = 1, count = 2) %>%
    dplyr::filter(!grepl("^__", gene_id)) %>%
    dplyr::mutate(
      gene_id = gsub("\\..*$", "", gene_id),
      count = as.numeric(count)
    ) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(
      count = sum(count),
      .groups = "drop"
    ) %>%
    dplyr::rename(!!unique_sample := count)
}

count_list_hdac12 <- mapply(
  read_count,
  meta_hdac12$sample_id,
  meta_hdac12$unique_sample,
  SIMPLIFY = FALSE
)

counts_hdac12 <- Reduce(
  function(x, y) dplyr::full_join(x, y, by = "gene_id"),
  count_list_hdac12
) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

counts_hdac12[is.na(counts_hdac12)] <- 0
storage.mode(counts_hdac12) <- "integer"

counts_hdac12 <- counts_hdac12[, rownames(meta_hdac12)]

stopifnot(all(colnames(counts_hdac12) == rownames(meta_hdac12)))

dim(counts_hdac12)


## 5. Remove HDAC12_H3 and create clear sample labels


meta_hdac12_noH3 <- meta_hdac12 %>%
  as.data.frame() %>%
  dplyr::filter(unique_sample != "HDAC12_H3") %>%
  dplyr::mutate(
    condition = factor(
      condition,
      levels = c("Scramble", "HDAC1_siRNA", "HDAC2_siRNA")
    ),
    sample_label = dplyr::case_when(
      condition == "Scramble"    ~ paste0("Scramble_rep", rep),
      condition == "HDAC1_siRNA" ~ paste0("HDAC1_KD_rep", rep),
      condition == "HDAC2_siRNA" ~ paste0("HDAC2_KD_rep", rep),
      TRUE ~ unique_sample
    )
  )

rownames(meta_hdac12_noH3) <- meta_hdac12_noH3$unique_sample

counts_hdac12_noH3 <- counts_hdac12[, rownames(meta_hdac12_noH3)]

stopifnot(all(colnames(counts_hdac12_noH3) == rownames(meta_hdac12_noH3)))

table(meta_hdac12_noH3$condition)

meta_hdac12_noH3 %>%
  dplyr::select(unique_sample, sample_label, condition, rep)

## save for future  
saveRDS(counts_hdac12, file.path(base_dir, "counts_hdac12.rds"))

## 6. Sample count QC

sample_count_qc_hdac12_noH3 <- tibble::tibble(
  sample = colnames(counts_hdac12_noH3),
  total_counts = colSums(counts_hdac12_noH3),
  detected_genes = colSums(counts_hdac12_noH3 > 0)
) %>%
  dplyr::left_join(
    meta_hdac12_noH3 %>%
      tibble::rownames_to_column("sample") %>%
      dplyr::select(sample, sample_label, condition, rep),
    by = "sample"
  )

sample_count_qc_hdac12_noH3

write.csv(
  sample_count_qc_hdac12_noH3,
  file.path(plot_dir_hdac12_noH3, "sample_count_QC_HDAC12_noH3.csv"),
  row.names = FALSE
)

## 7. DESeq2: HDAC1/2 only, H3 removed

dds_hdac12_noH3 <- DESeqDataSetFromMatrix(
  countData = counts_hdac12_noH3,
  colData = meta_hdac12_noH3,
  design = ~ condition
)

keep <- rowSums(counts(dds_hdac12_noH3)) >= 10
dds_hdac12_noH3 <- dds_hdac12_noH3[keep, ]

dds_hdac12_noH3 <- DESeq(
  dds_hdac12_noH3,
  fitType = "local"
)

resultsNames(dds_hdac12_noH3)


## 8. PCA without H3

vsd_hdac12_noH3 <- vst(
  dds_hdac12_noH3,
  blind = FALSE
)

pca_hdac12_noH3_raw <- plotPCA(
  vsd_hdac12_noH3,
  intgroup = "condition",
  ntop = 500,
  returnData = TRUE
)

percentVar_hdac12_noH3 <- round(
  100 * attr(pca_hdac12_noH3_raw, "percentVar")
)

pca_hdac12_noH3_df <- as.data.frame(pca_hdac12_noH3_raw)

pca_hdac12_noH3_df$sample <- rownames(pca_hdac12_noH3_df)

pca_hdac12_noH3_df$sample_label <- meta_hdac12_noH3[
  pca_hdac12_noH3_df$sample,
  "sample_label"
]

pca_hdac12_noH3_df %>%
  dplyr::select(sample, sample_label, condition)

pca_hdac12_noH3_plot <- ggplot(
  pca_hdac12_noH3_df,
  aes(
    x = PC1,
    y = PC2,
    colour = condition,
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
    title = "PCA of HDAC1/2 experiment without HDAC12_H3",
    subtitle = "DESeq2 VST; top 500 most variable genes",
    x = paste0("PC1: ", percentVar_hdac12_noH3[1], "% variance"),
    y = paste0("PC2: ", percentVar_hdac12_noH3[2], "% variance"),
    colour = "Condition"
  )

pca_hdac12_noH3_plot

ggsave(
  filename = file.path(plot_dir_hdac12_noH3, "01_PCA_HDAC12_noH3_labelled.png"),
  plot = pca_hdac12_noH3_plot,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)
