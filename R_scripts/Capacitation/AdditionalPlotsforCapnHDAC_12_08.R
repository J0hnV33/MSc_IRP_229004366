## Capacitation bulk TPM heatmap
## Tukiainen Subject-to-XCI genes only


library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(tibble)
library(purrr)
library(ggplot2)
library(scales)


## 1. Paths


cap_bulk_tpm_dir <- "/home/jvk3/Desktop/HDAC_counts/All three datsets bulk and ASE/Capacitation/bulk/2_STAR_libinorm_counts/libinorm_normalised"

cap_bulk_meta_path <- "/home/jvk3/Desktop/HDAC_counts/All three datsets bulk and ASE/Capacitation/meta_table_cap_bulk.csv"

tukiainen_path <- "/home/jvk3/Desktop/HDAC_counts/XCI_classifications/XCI_classification_outputs/Tukiainen_XCI_categories_clean.csv"

cap_bulk_heatmap_out_dir <- "/home/jvk3/Desktop/HDAC_counts/All three datsets bulk and ASE/Capacitation/Bulk_Tukiainen_heatmap_outputs"

dir.create(cap_bulk_heatmap_out_dir, recursive = TRUE, showWarnings = FALSE)


## 2. Read metadata and Tukiainen annotation


read_meta_auto <- function(path) {
  meta <- readr::read_csv(path, show_col_types = FALSE)
  
  if (ncol(meta) == 1) {
    meta <- readr::read_tsv(path, show_col_types = FALSE)
  }
  
  meta
}

cap_bulk_meta <- read_meta_auto(cap_bulk_meta_path)

cap_bulk_meta
colnames(cap_bulk_meta)

cap_bulk_meta_clean <- cap_bulk_meta %>%
  dplyr::mutate(
    filename = trimws(filename),
    sample_code = stringr::str_extract(filename, "^[A-Z][0-9]"),
    rep = as.integer(rep),
    condition_clean = dplyr::case_when(
      stringr::str_detect(condition, regex("naive", ignore_case = TRUE)) ~ "Naive",
      stringr::str_detect(condition, regex("day.?3|d3", ignore_case = TRUE)) ~ "Day3",
      stringr::str_detect(condition, regex("day.?7|d7", ignore_case = TRUE)) ~ "Day7",
      stringr::str_detect(condition, regex("primed", ignore_case = TRUE)) ~ "Primed",
      TRUE ~ as.character(condition)
    ),
    condition_clean = factor(
      condition_clean,
      levels = c("Naive", "Day3", "Day7", "Primed")
    )
  )

cap_bulk_meta_clean %>%
  dplyr::select(sample, filename, sample_code, rep, condition, condition_clean) %>%
  dplyr::arrange(condition_clean, rep) %>%
  print(n = Inf, width = Inf)

tukiainen_annot <- readr::read_csv(tukiainen_path, show_col_types = FALSE)

dim(tukiainen_annot)
colnames(tukiainen_annot)
head(tukiainen_annot)


## 3. Standardise Tukiainen columns


# Adjust these names only if your Tukiainen file uses different column names.
tuk_symbol_col <- "SYMBOL"
tuk_status_col <- "tukiainen_status_grouped"

if (!tuk_symbol_col %in% colnames(tukiainen_annot)) {
  stop("SYMBOL column not found in Tukiainen annotation. Check colnames(tukiainen_annot).")
}

if (!tuk_status_col %in% colnames(tukiainen_annot)) {
  stop("tukiainen_status_grouped column not found. Check colnames(tukiainen_annot).")
}

tukiainen_subject_genes <- tukiainen_annot %>%
  dplyr::rename(
    SYMBOL = dplyr::all_of(tuk_symbol_col),
    tukiainen_status_grouped = dplyr::all_of(tuk_status_col)
  ) %>%
  dplyr::filter(tukiainen_status_grouped == "Subject to XCI") %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

nrow(tukiainen_subject_genes)

head(tukiainen_subject_genes)


## 4. List and check bulk TPM files


cap_bulk_files <- tibble::tibble(
  file_path = list.files(
    cap_bulk_tpm_dir,
    pattern = "expression.txt$",
    full.names = TRUE
  )
) %>%
  dplyr::mutate(
    file_name = basename(file_path),
    sample_code = stringr::str_extract(file_name, "^[A-Z][0-9]"),
    size_mb = file.info(file_path)$size / 1024^2
  ) %>%
  dplyr::arrange(file_name)

cap_bulk_files %>%
  dplyr::select(file_name, sample_code, size_mb) %>%
  print(n = Inf, width = Inf)

# Check files missing from metadata
cap_bulk_files %>%
  dplyr::anti_join(
    cap_bulk_meta_clean,
    by = c("file_name" = "filename")
  )

# Check metadata rows missing files
cap_bulk_meta_clean %>%
  dplyr::anti_join(
    cap_bulk_files,
    by = c("filename" = "file_name")
  )


## 5. Inspect one bulk file


bulk_preview <- readr::read_tsv(
  cap_bulk_files$file_path[1],
  show_col_types = FALSE,
  n_max = 20
)

colnames(bulk_preview)
head(bulk_preview)


## 6. Read bulk TPM expression files
## Correct parser for libinorm expression.txt format


read_cap_bulk_expression <- function(path) {
  
  df_raw <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    col_names = FALSE
  )
  
  # The first row contains the real column names:
  # Gene | count | length | RPM | RPKM | RPK | TPM
  df_clean <- df_raw %>%
    tidyr::separate(
      col = X3,
      into = c("length", "RPM", "RPKM", "RPK", "TPM"),
      sep = "\t",
      remove = TRUE,
      fill = "right"
    ) %>%
    dplyr::rename(
      SYMBOL = X1,
      count = X2
    ) %>%
    dplyr::filter(SYMBOL != "Gene") %>%
    dplyr::mutate(
      file_name = basename(path),
      SYMBOL = as.character(SYMBOL),
      count = as.numeric(count),
      length = as.numeric(length),
      RPM = as.numeric(RPM),
      RPKM = as.numeric(RPKM),
      RPK = as.numeric(RPK),
      TPM = as.numeric(TPM)
    ) %>%
    dplyr::select(file_name, SYMBOL, count, length, RPM, RPKM, RPK, TPM)
  
  df_clean
}

cap_bulk_expr_long <- purrr::map_dfr(
  cap_bulk_files$file_path,
  read_cap_bulk_expression
)

dim(cap_bulk_expr_long)
head(cap_bulk_expr_long)
summary(cap_bulk_expr_long$TPM)


## 7. Add metadata using sample_code


cap_bulk_expr_meta <- cap_bulk_expr_long %>%
  dplyr::mutate(
    sample_code = stringr::str_extract(file_name, "^[A-Z][0-9]")
  ) %>%
  dplyr::left_join(
    cap_bulk_meta_clean %>%
      dplyr::select(
        sample,
        metadata_filename = filename,
        sample_code,
        rep,
        condition,
        condition_clean
      ),
    by = "sample_code"
  ) %>%
  dplyr::mutate(
    condition_clean = factor(
      condition_clean,
      levels = c("Naive", "Day3", "Day7", "Primed")
    )
  )

# Check failed joins
cap_bulk_expr_meta %>%
  dplyr::filter(is.na(condition_clean)) %>%
  dplyr::distinct(file_name, sample_code)

# Check sample assignment
cap_bulk_expr_meta %>%
  dplyr::distinct(file_name, sample_code, sample, condition_clean, rep) %>%
  dplyr::arrange(condition_clean, rep) %>%
  print(n = Inf, width = Inf)


## 8. Filter to Tukiainen Subject-to-XCI genes


cap_subject_expr <- cap_bulk_expr_meta %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::inner_join(
    tukiainen_subject_genes %>%
      dplyr::select(SYMBOL, tukiainen_status_grouped),
    by = "SYMBOL"
  ) %>%
  dplyr::filter(!is.na(TPM))

cap_subject_expr %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_genes = dplyr::n_distinct(SYMBOL),
    n_samples = dplyr::n_distinct(sample)
  )

cap_subject_expr %>%
  dplyr::count(condition_clean)


## 9. Mean TPM per timepoint


cap_subject_mean_expr <- cap_subject_expr %>%
  dplyr::group_by(SYMBOL, condition_clean) %>%
  dplyr::summarise(
    mean_tpm = mean(TPM, na.rm = TRUE),
    mean_log2_tpm = mean(log2(TPM + 1), na.rm = TRUE),
    .groups = "drop"
  )

cap_subject_wide <- cap_subject_mean_expr %>%
  dplyr::select(SYMBOL, condition_clean, mean_log2_tpm) %>%
  tidyr::pivot_wider(
    names_from = condition_clean,
    values_from = mean_log2_tpm
  ) %>%
  dplyr::filter(
    !is.na(Naive),
    !is.na(Day3),
    !is.na(Day7),
    !is.na(Primed)
  )

dim(cap_subject_wide)
head(cap_subject_wide)


## 10. Select Subject-to-XCI genes by highest Naive TPM


top_n_genes <- 25
min_naive_tpm <- 1

cap_subject_wide_silencing <- cap_subject_mean_expr %>%
  dplyr::select(SYMBOL, condition_clean, mean_log2_tpm, mean_tpm) %>%
  tidyr::pivot_wider(
    names_from = condition_clean,
    values_from = c(mean_log2_tpm, mean_tpm)
  ) %>%
  dplyr::filter(
    !is.na(mean_log2_tpm_Naive),
    !is.na(mean_log2_tpm_Day3),
    !is.na(mean_log2_tpm_Day7),
    !is.na(mean_log2_tpm_Primed),
    mean_tpm_Naive >= min_naive_tpm
  ) %>%
  dplyr::mutate(
    naive_to_primed_change = mean_log2_tpm_Primed - mean_log2_tpm_Naive
  ) %>%
  dplyr::filter(
    naive_to_primed_change < 0
  ) %>%
  dplyr::arrange(naive_to_primed_change) %>%   # strongest decrease first
  dplyr::slice_head(n = top_n_genes)

cap_subject_wide_naive_top <- cap_subject_wide_silencing %>%
  dplyr::select(
    SYMBOL,
    Naive = mean_log2_tpm_Naive,
    Day3 = mean_log2_tpm_Day3,
    Day7 = mean_log2_tpm_Day7,
    Primed = mean_log2_tpm_Primed,
    naive_mean_tpm = mean_tpm_Naive,
    naive_to_primed_change
  )

## 11. Row-scale log2 TPM and order by Naive-to-Primed change


expr_mat <- cap_subject_wide_naive_top %>%
  dplyr::select(SYMBOL, Naive, Day3, Day7, Primed) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  as.matrix()

expr_mat_z <- t(scale(t(expr_mat)))

expr_mat_z[expr_mat_z > 2] <- 2
expr_mat_z[expr_mat_z < -2] <- -2

gene_order <- cap_subject_wide_naive_top %>%
  dplyr::mutate(
    naive_to_primed_change = Primed - Naive
  ) %>%
  dplyr::arrange(naive_to_primed_change) %>%
  dplyr::pull(SYMBOL)

heatmap_df_naive_top <- as.data.frame(expr_mat_z) %>%
  tibble::rownames_to_column("SYMBOL") %>%
  tidyr::pivot_longer(
    cols = c("Naive", "Day3", "Day7", "Primed"),
    names_to = "condition_clean",
    values_to = "row_z"
  ) %>%
  dplyr::mutate(
    condition_clean = factor(
      condition_clean,
      levels = c("Naive", "Day3", "Day7", "Primed")
    ),
    SYMBOL = factor(SYMBOL, levels = rev(gene_order))
  )


## 12. Plot Naive-high Subject-to-XCI heatmap


p_cap_subject_heatmap_naive_top <- ggplot(
  heatmap_df_naive_top,
  aes(x = condition_clean, y = SYMBOL, fill = row_z)
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  theme_bw() +
  labs(
    x = "Capacitation stage",
    y = "Tukiainen Subject-to-XCI genes",
    fill = "Row-scaled\nlog2(TPM + 1)",
    title = "Naive-high Subject-to-XCI gene expression across capacitation",
    subtitle = paste0(
      "Top ", top_n_genes,
      " Tukiainen Subject-to-XCI genes ranked by mean Naive TPM"
    )
  ) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

p_cap_subject_heatmap_naive_top

ggsave(
  filename = file.path(
    cap_bulk_heatmap_out_dir,
    "HeatmapofSubjecttoXCI.png"
  ),
  plot = p_cap_subject_heatmap_naive_top,
  width = 9,
  height = 8,
  dpi = 300
)
