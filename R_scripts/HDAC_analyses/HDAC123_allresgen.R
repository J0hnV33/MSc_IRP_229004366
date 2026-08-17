## HDAC3-only all_res table
## Separate HDAC3 experiment analysis

## 1. Libraries

library(DESeq2)
library(apeglm)
library(dplyr)
library(tidyr)
library(tibble)
library(biomaRt)

## 2. Paths

base_dir <- "/home/jvk3/Desktop/HDAC_counts"

hdac3_counts_dir <- file.path(
  base_dir,
  "HDAC3_geneIDs"
)

all_res_dir <- file.path(
  base_dir,
  "All_Res_Tables"
)

dir.create(all_res_dir, showWarnings = FALSE, recursive = TRUE)

## 3. Create HDAC3 metadata

meta_hdac3 <- tibble::tibble(
  sample_id = c("A1", "A2", "B3", "B1", "B2", "A3"),
  condition = c(
    "Scramble", "Scramble", "Scramble",
    "HDAC3_siRNA", "HDAC3_siRNA", "HDAC3_siRNA"
  ),
  rep = c(1, 2, 3, 1, 2, 3),
  experiment = "HDAC3"
) %>%
  dplyr::mutate(
    unique_sample = paste0("HDAC3_", sample_id),
    condition = factor(
      condition,
      levels = c("Scramble", "HDAC3_siRNA")
    ),
    sample_label = dplyr::case_when(
      condition == "Scramble" ~ paste0("HDAC3_Scramble_rep", rep),
      condition == "HDAC3_siRNA" ~ paste0("HDAC3_KD_rep", rep),
      TRUE ~ unique_sample
    )
  ) %>%
  as.data.frame()

rownames(meta_hdac3) <- meta_hdac3$unique_sample

meta_hdac3

## 4. Read HDAC3 count files

read_count <- function(sample_id, unique_sample) {
  
  file <- list.files(
    hdac3_counts_dir,
    pattern = paste0("^", sample_id, "_.*_counts\\.txt$"),
    full.names = TRUE
  )
  
  file <- file[
    !grepl("bias|expression|distribution", basename(file), ignore.case = TRUE)
  ]
  
  if (length(file) != 1) {
    stop(paste("Problem finding count file for", sample_id))
  }
  
  message(sample_id, " -> ", basename(file))
  
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

count_list_hdac3 <- mapply(
  read_count,
  sample_id = meta_hdac3$sample_id,
  unique_sample = meta_hdac3$unique_sample,
  SIMPLIFY = FALSE
)

counts_hdac3 <- Reduce(
  function(x, y) dplyr::full_join(x, y, by = "gene_id"),
  count_list_hdac3
) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

counts_hdac3[is.na(counts_hdac3)] <- 0
storage.mode(counts_hdac3) <- "integer"

counts_hdac3 <- counts_hdac3[, rownames(meta_hdac3)]

stopifnot(all(colnames(counts_hdac3) == rownames(meta_hdac3)))

dim(counts_hdac3)

## 5. Sample count QC

sample_count_qc_hdac3 <- tibble::tibble(
  sample = colnames(counts_hdac3),
  total_counts = colSums(counts_hdac3),
  detected_genes = colSums(counts_hdac3 > 0)
) %>%
  dplyr::left_join(
    meta_hdac3 %>%
      tibble::rownames_to_column("sample") %>%
      dplyr::select(sample, sample_label, condition, rep),
    by = "sample"
  )

sample_count_qc_hdac3

write.csv(
  sample_count_qc_hdac3,
  file.path(all_res_dir, "sample_count_QC_HDAC3.csv"),
  row.names = FALSE
)

## 6. DESeq2: HDAC3 only

dds_hdac3 <- DESeqDataSetFromMatrix(
  countData = counts_hdac3,
  colData = meta_hdac3,
  design = ~ condition
)

keep <- rowSums(counts(dds_hdac3)) >= 10
dds_hdac3 <- dds_hdac3[keep, ]

dds_hdac3 <- DESeq(
  dds_hdac3,
  fitType = "local"
)

resultsNames(dds_hdac3)


## 7. Annotate genes using biomaRt

gene_ids_hdac3 <- rownames(dds_hdac3)

mart <- biomaRt::useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"
)

get_bm_chunked <- function(values, chunk_size = 2000) {
  
  chunks <- split(values, ceiling(seq_along(values) / chunk_size))
  
  out <- lapply(seq_along(chunks), function(i) {
    
    message("Running biomaRt chunk ", i, " of ", length(chunks))
    
    Sys.sleep(0.5)
    
    biomaRt::getBM(
      attributes = c(
        "ensembl_gene_id",
        "hgnc_symbol",
        "entrezgene_id",
        "chromosome_name",
        "start_position",
        "end_position",
        "gene_biotype",
        "description"
      ),
      filters = "ensembl_gene_id",
      values = chunks[[i]],
      mart = mart
    )
  })
  
  dplyr::bind_rows(out)
}

annot_raw_hdac3 <- get_bm_chunked(gene_ids_hdac3)

first_non_empty <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA)
  return(x[1])
}

gene_annot_hdac3 <- annot_raw_hdac3 %>%
  dplyr::mutate(
    chromosome_name = as.character(chromosome_name)
  ) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(
    SYMBOL = first_non_empty(hgnc_symbol),
    ENTREZID = first_non_empty(as.character(entrezgene_id)),
    CHR = first_non_empty(chromosome_name),
    START = suppressWarnings(min(start_position, na.rm = TRUE)),
    END = suppressWarnings(max(end_position, na.rm = TRUE)),
    GENE_BIOTYPE = first_non_empty(gene_biotype),
    DESCRIPTION = first_non_empty(description),
    .groups = "drop"
  ) %>%
  dplyr::rename(gene_id = ensembl_gene_id)

annotation_check_hdac3 <- tibble::tibble(
  total_genes_in_dds = length(gene_ids_hdac3),
  annotated_genes = nrow(gene_annot_hdac3),
  genes_with_symbol = sum(!is.na(gene_annot_hdac3$SYMBOL)),
  x_linked_genes = sum(gene_annot_hdac3$CHR == "X", na.rm = TRUE),
  x_linked_genes_with_position = sum(
    gene_annot_hdac3$CHR == "X" &
      !is.na(gene_annot_hdac3$START),
    na.rm = TRUE
  )
)

annotation_check_hdac3

## 8. Extract shrunken DESeq2 result

res_hdac3_lfc <- lfcShrink(
  dds_hdac3,
  coef = "condition_HDAC3_siRNA_vs_Scramble",
  type = "apeglm"
)

lfc_cutoff_15fold <- log2(1.5)
padj_cutoff <- 0.05

all_res_hdac3_15fold <- as.data.frame(res_hdac3_lfc) %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::left_join(gene_annot_hdac3, by = "gene_id") %>%
  dplyr::mutate(
    analysis_group = "HDAC3_only",
    experiment = "HDAC3",
    scramble_group = "HDAC3_Scramble",
    contrast = "HDAC3_siRNA_vs_HDAC3_Scramble",
    contrast_short = "HDAC3 KD",
    outlier_status = "No_outlier_removed",
    model_design = "~ condition",
    gene_label = dplyr::case_when(
      !is.na(SYMBOL) & SYMBOL != "" ~ SYMBOL,
      TRUE ~ gene_id
    ),
    deg_status_15fold = dplyr::case_when(
      !is.na(padj) &
        padj < padj_cutoff &
        log2FoldChange >= lfc_cutoff_15fold ~ "Up",
      !is.na(padj) &
        padj < padj_cutoff &
        log2FoldChange <= -lfc_cutoff_15fold ~ "Down",
      TRUE ~ "Not significant"
    )
  ) %>%
  dplyr::arrange(padj)

## 9. Add mean normalised counts

norm_counts_hdac3 <- counts(
  dds_hdac3,
  normalized = TRUE
)

scramble_samples_hdac3 <- rownames(meta_hdac3)[
  meta_hdac3$condition == "Scramble"
]

hdac3_kd_samples <- rownames(meta_hdac3)[
  meta_hdac3$condition == "HDAC3_siRNA"
]

mean_norm_counts_hdac3 <- tibble::tibble(
  gene_id = rownames(norm_counts_hdac3),
  Scramble_mean_norm_count = rowMeans(
    norm_counts_hdac3[, scramble_samples_hdac3, drop = FALSE],
    na.rm = TRUE
  ),
  HDAC3_KD_mean_norm_count = rowMeans(
    norm_counts_hdac3[, hdac3_kd_samples, drop = FALSE],
    na.rm = TRUE
  )
)

all_res_hdac3_15fold <- all_res_hdac3_15fold %>%
  dplyr::left_join(
    mean_norm_counts_hdac3,
    by = "gene_id"
  ) %>%
  dplyr::mutate(
    control_mean_norm_count = Scramble_mean_norm_count,
    treatment_mean_norm_count = HDAC3_KD_mean_norm_count
  )

all_res_hdac3_15fold %>%
  dplyr::select(
    gene_id,
    SYMBOL,
    contrast_short,
    log2FoldChange,
    padj,
    control_mean_norm_count,
    treatment_mean_norm_count,
    deg_status_15fold
  ) %>%
  head()


## 10. Save HDAC3 all_res table


write.csv(
  all_res_hdac3_15fold,
  file.path(all_res_dir, "all_results_HDAC3_15fold.csv"),
  row.names = FALSE
)

write.csv(
  annotation_check_hdac3,
  file.path(all_res_dir, "annotation_check_HDAC3.csv"),
  row.names = FALSE
)

list.files(all_res_dir)

all_res_hdac3_15fold %>%
  dplyr::count(deg_status_15fold)

all_res_hdac3_15fold %>%
  dplyr::filter(CHR == "X") %>%
  dplyr::count(deg_status_15fold)

all_res_hdac3_15fold %>%
  dplyr::filter(SYMBOL == "HDAC3") %>%
  dplyr::select(
    SYMBOL,
    log2FoldChange,
    padj,
    control_mean_norm_count,
    treatment_mean_norm_count,
    deg_status_15fold
  )


## Check HDAC1/2 and HDAC3 all_res table compatibility


all_res_hdac12_noH3_15fold <- read.csv(
  file.path(all_res_dir, "all_results_HDAC12_noH3_15fold.csv"),
  stringsAsFactors = FALSE
)

all_res_hdac3_15fold <- read.csv(
  file.path(all_res_dir, "all_results_HDAC3_15fold.csv"),
  stringsAsFactors = FALSE
)

setdiff(colnames(all_res_hdac12_noH3_15fold), colnames(all_res_hdac3_15fold))
setdiff(colnames(all_res_hdac3_15fold), colnames(all_res_hdac12_noH3_15fold))


## Standardise HDAC1/2 and HDAC3 all_res tables before combining

all_res_hdac12_noH3_15fold_clean <- all_res_hdac12_noH3_15fold %>%
  dplyr::select(
    -dplyr::any_of(c(
      "HDAC1_KD_mean_norm_count",
      "HDAC2_KD_mean_norm_count",
      "HDAC3_KD_mean_norm_count"
    ))
  )

all_res_hdac3_15fold_clean <- all_res_hdac3_15fold %>%
  dplyr::select(
    -dplyr::any_of(c(
      "HDAC1_KD_mean_norm_count",
      "HDAC2_KD_mean_norm_count",
      "HDAC3_KD_mean_norm_count"
    ))
  )

setdiff(colnames(all_res_hdac12_noH3_15fold_clean), colnames(all_res_hdac3_15fold_clean))
setdiff(colnames(all_res_hdac3_15fold_clean), colnames(all_res_hdac12_noH3_15fold_clean))



## Combine HDAC1/2 and HDAC3 all_res tables

all_res_hdac1_2_3_15fold <- dplyr::bind_rows(
  all_res_hdac12_noH3_15fold_clean,
  all_res_hdac3_15fold_clean
)

write.csv(
  all_res_hdac1_2_3_15fold,
  file.path(all_res_dir, "all_results_HDAC1_2_3_combined_15fold.csv"),
  row.names = FALSE
)

all_res_hdac1_2_3_15fold %>%
  dplyr::count(experiment, contrast_short, deg_status_15fold)

all_res_hdac1_2_3_15fold %>%
  dplyr::filter(CHR == "X") %>%
  dplyr::count(experiment, contrast_short, deg_status_15fold)

combined_deg_summary <- all_res_hdac1_2_3_15fold %>%
  dplyr::count(
    experiment,
    contrast_short,
    deg_status_15fold,
    name = "n_genes"
  )

combined_x_deg_summary <- all_res_hdac1_2_3_15fold %>%
  dplyr::filter(CHR == "X") %>%
  dplyr::count(
    experiment,
    contrast_short,
    deg_status_15fold,
    name = "n_genes"
  )

write.csv(
  combined_deg_summary,
  file.path(all_res_dir, "combined_DEG_summary_genomewide_15fold.csv"),
  row.names = FALSE
)

write.csv(
  combined_x_deg_summary,
  file.path(all_res_dir, "combined_DEG_summary_chrX_15fold.csv"),
  row.names = FALSE
)

combined_deg_summary
combined_x_deg_summary
