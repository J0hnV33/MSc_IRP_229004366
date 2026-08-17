## XCI gene classification analysis
## Tukiainen / Marks annotation of HDAC and capacitation results

library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

base_dir <- "/home/jvk3/Desktop/HDAC_counts"

xci_class_dir <- file.path(base_dir, "XCI_classifications")

tukiainen_dir <- file.path(
  xci_class_dir,
  "Tukiainen_2017",
  "Tukiainen_2017_Supplementary_Tables",
  "nature24265-s3"
)

marks_dir <- file.path(
  xci_class_dir,
  "Marks_2015"
)

all_res_dir <- file.path(base_dir, "All_Res_Tables")
cap_dir <- file.path(base_dir, "Capacitation")
cap_plot_dir <- file.path(cap_dir, "Capacitation_outputs")

xci_output_dir <- file.path(xci_class_dir, "XCI_classification_outputs")
dir.create(xci_output_dir, showWarnings = FALSE, recursive = TRUE)



## Load HDAC and capacitation results


all_res_hdac1_2_3_15fold <- read.csv(
  file.path(all_res_dir, "all_results_HDAC1_2_3_combined_15fold.csv"),
  stringsAsFactors = FALSE
)

all_res_capacitation_15fold <- read.csv(
  file.path(cap_plot_dir, "all_results_capacitation_15fold.csv"),
  stringsAsFactors = FALSE
)

all_res_hdac1_2_3_15fold <- all_res_hdac1_2_3_15fold %>%
  dplyr::mutate(
    gene_id = gsub("\\..*$", "", gene_id),
    SYMBOL = as.character(SYMBOL),
    CHR = as.character(CHR)
  )

all_res_capacitation_15fold <- all_res_capacitation_15fold %>%
  dplyr::mutate(
    SYMBOL = as.character(SYMBOL),
    CHR = as.character(CHR)
  )


## Inspect Tukiainen Table S1


tukiainen_table1_file <- file.path(
  tukiainen_dir,
  "Suppl.Table.1.xlsx"
)

readxl::excel_sheets(tukiainen_table1_file)

tukiainen_table1_raw <- readxl::read_excel(
  tukiainen_table1_file,
  sheet = 1
)

dim(tukiainen_table1_raw)

colnames(tukiainen_table1_raw)

head(tukiainen_table1_raw)


## Clean Tukiainen Table S1


tukiainen_table1 <- readxl::read_excel(
  tukiainen_table1_file,
  sheet = 1,
  skip = 1
)

dim(tukiainen_table1)

colnames(tukiainen_table1)

head(tukiainen_table1)


## Extract Tukiainen XCI categories


tukiainen_xci_categories <- tukiainen_table1 %>%
  dplyr::rename(
    SYMBOL = `Gene name`,
    tukiainen_gene_id = `Gene ID`,
    tukiainen_chr = Chr,
    tukiainen_start = `Start position`,
    tukiainen_end = `End position`,
    tukiainen_biotype = `Transcript type`,
    tukiainen_status = `Combined XCI status`
  ) %>%
  dplyr::mutate(
    SYMBOL = as.character(SYMBOL),
    tukiainen_gene_id = gsub("\\..*$", "", tukiainen_gene_id),
    tukiainen_status = tolower(as.character(tukiainen_status)),
    tukiainen_status_grouped = dplyr::case_when(
      tukiainen_status == "inactive" ~ "Subject to XCI",
      tukiainen_status == "escape" ~ "Escapee",
      tukiainen_status == "variable escape" ~ "Variable escapee",
      TRUE ~ "Unclassified"
    ),
    tukiainen_status_grouped = factor(
      tukiainen_status_grouped,
      levels = c(
        "Subject to XCI",
        "Escapee",
        "Variable escapee",
        "Unclassified"
      )
    )
  ) %>%
  dplyr::select(
    SYMBOL,
    tukiainen_gene_id,
    tukiainen_chr,
    tukiainen_start,
    tukiainen_end,
    tukiainen_biotype,
    tukiainen_status,
    tukiainen_status_grouped
  ) %>%
  dplyr::filter(
    !is.na(SYMBOL),
    SYMBOL != ""
  ) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

tukiainen_xci_categories %>%
  dplyr::count(tukiainen_status_grouped)


## Add Tukiainen categories to HDAC and capacitation results


hdac_tukiainen_annotated <- all_res_hdac1_2_3_15fold %>%
  dplyr::left_join(
    tukiainen_xci_categories,
    by = "SYMBOL"
  ) %>%
  dplyr::mutate(
    tukiainen_status_grouped = dplyr::if_else(
      is.na(as.character(tukiainen_status_grouped)),
      "Unclassified",
      as.character(tukiainen_status_grouped)
    ),
    tukiainen_status_grouped = factor(
      tukiainen_status_grouped,
      levels = c(
        "Subject to XCI",
        "Escapee",
        "Variable escapee",
        "Unclassified"
      )
    )
  )

cap_tukiainen_annotated <- all_res_capacitation_15fold %>%
  dplyr::left_join(
    tukiainen_xci_categories,
    by = "SYMBOL"
  ) %>%
  dplyr::mutate(
    tukiainen_status_grouped = dplyr::if_else(
      is.na(as.character(tukiainen_status_grouped)),
      "Unclassified",
      as.character(tukiainen_status_grouped)
    ),
    tukiainen_status_grouped = factor(
      tukiainen_status_grouped,
      levels = c(
        "Subject to XCI",
        "Escapee",
        "Variable escapee",
        "Unclassified"
      )
    )
  )


## Check Tukiainen annotation overlap


hdac_tukiainen_annotated %>%
  dplyr::filter(CHR == "X") %>%
  dplyr::count(tukiainen_status_grouped)

cap_tukiainen_annotated %>%
  dplyr::filter(CHR == "X") %>%
  dplyr::count(tukiainen_status_grouped)


## Tukiainen categories among X-linked DEGs


hdac_x_tukiainen_deg_summary <- hdac_tukiainen_annotated %>%
  dplyr::filter(
    CHR == "X",
    deg_status_15fold != "Not significant"
  ) %>%
  dplyr::count(
    contrast_short,
    deg_status_15fold,
    tukiainen_status_grouped,
    name = "n_genes"
  )

cap_x_tukiainen_deg_summary <- cap_tukiainen_annotated %>%
  dplyr::filter(
    CHR == "X",
    deg_status_15fold != "Not significant"
  ) %>%
  dplyr::count(
    contrast_short,
    deg_status_15fold,
    tukiainen_status_grouped,
    name = "n_genes"
  )

hdac_x_tukiainen_deg_summary

cap_x_tukiainen_deg_summary


## Save Tukiainen-annotated result tables

write.csv(
  tukiainen_xci_categories,
  file.path(xci_output_dir, "Tukiainen_XCI_categories_clean.csv"),
  row.names = FALSE
)

write.csv(
  hdac_tukiainen_annotated,
  file.path(xci_output_dir, "HDAC_results_Tukiainen_annotated.csv"),
  row.names = FALSE
)

write.csv(
  cap_tukiainen_annotated,
  file.path(xci_output_dir, "Capacitation_results_Tukiainen_annotated.csv"),
  row.names = FALSE
)

write.csv(
  hdac_x_tukiainen_deg_summary,
  file.path(xci_output_dir, "HDAC_Xlinked_DEGs_Tukiainen_summary.csv"),
  row.names = FALSE
)

write.csv(
  cap_x_tukiainen_deg_summary,
  file.path(xci_output_dir, "Capacitation_Xlinked_DEGs_Tukiainen_summary.csv"),
  row.names = FALSE
)


## Check exact Tukiainen status labels


tukiainen_table1 %>%
  dplyr::count(`Combined XCI status`)


## Re-clean Tukiainen XCI categories


tukiainen_xci_categories <- tukiainen_table1 %>%
  dplyr::rename(
    SYMBOL = `Gene name`,
    tukiainen_gene_id = `Gene ID`,
    tukiainen_chr = Chr,
    tukiainen_start = `Start position`,
    tukiainen_end = `End position`,
    tukiainen_biotype = `Transcript type`,
    tukiainen_status = `Combined XCI status`
  ) %>%
  dplyr::mutate(
    SYMBOL = as.character(SYMBOL),
    tukiainen_gene_id = gsub("\\..*$", "", tukiainen_gene_id),
    tukiainen_status = tolower(trimws(as.character(tukiainen_status))),
    tukiainen_status_grouped = dplyr::case_when(
      tukiainen_status %in% c("inactive", "subject to xci", "subject") ~ "Subject to XCI",
      tukiainen_status %in% c("escape", "escapee") ~ "Escapee",
      grepl("variable", tukiainen_status) ~ "Variable escapee",
      TRUE ~ "Unclassified"
    ),
    tukiainen_status_grouped = factor(
      tukiainen_status_grouped,
      levels = c(
        "Subject to XCI",
        "Escapee",
        "Variable escapee",
        "Unclassified"
      )
    )
  ) %>%
  dplyr::select(
    SYMBOL,
    tukiainen_gene_id,
    tukiainen_chr,
    tukiainen_start,
    tukiainen_end,
    tukiainen_biotype,
    tukiainen_status,
    tukiainen_status_grouped
  ) %>%
  dplyr::filter(
    !is.na(SYMBOL),
    SYMBOL != ""
  ) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

tukiainen_xci_categories %>%
  dplyr::count(tukiainen_status_grouped)


## Re-annotate HDAC and capacitation results


hdac_tukiainen_annotated <- all_res_hdac1_2_3_15fold %>%
  dplyr::left_join(
    tukiainen_xci_categories,
    by = "SYMBOL"
  ) %>%
  dplyr::mutate(
    tukiainen_status_grouped = dplyr::if_else(
      is.na(as.character(tukiainen_status_grouped)),
      "Unclassified",
      as.character(tukiainen_status_grouped)
    ),
    tukiainen_status_grouped = factor(
      tukiainen_status_grouped,
      levels = c(
        "Subject to XCI",
        "Escapee",
        "Variable escapee",
        "Unclassified"
      )
    )
  )

cap_tukiainen_annotated <- all_res_capacitation_15fold %>%
  dplyr::left_join(
    tukiainen_xci_categories,
    by = "SYMBOL"
  ) %>%
  dplyr::mutate(
    tukiainen_status_grouped = dplyr::if_else(
      is.na(as.character(tukiainen_status_grouped)),
      "Unclassified",
      as.character(tukiainen_status_grouped)
    ),
    tukiainen_status_grouped = factor(
      tukiainen_status_grouped,
      levels = c(
        "Subject to XCI",
        "Escapee",
        "Variable escapee",
        "Unclassified"
      )
    )
  )


## Unique X-linked gene annotation overlap


hdac_tukiainen_annotated %>%
  dplyr::filter(CHR == "X") %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::count(tukiainen_status_grouped)

cap_tukiainen_annotated %>%
  dplyr::filter(CHR == "X") %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::count(tukiainen_status_grouped)


## Tukiainen categories among X-linked DEGs


hdac_x_tukiainen_deg_summary <- hdac_tukiainen_annotated %>%
  dplyr::filter(
    CHR == "X",
    deg_status_15fold != "Not significant"
  ) %>%
  dplyr::count(
    contrast_short,
    deg_status_15fold,
    tukiainen_status_grouped,
    name = "n_genes"
  )

cap_x_tukiainen_deg_summary <- cap_tukiainen_annotated %>%
  dplyr::filter(
    CHR == "X",
    deg_status_15fold != "Not significant"
  ) %>%
  dplyr::count(
    contrast_short,
    deg_status_15fold,
    tukiainen_status_grouped,
    name = "n_genes"
  )

hdac_x_tukiainen_deg_summary
cap_x_tukiainen_deg_summary


## Plot HDAC X-linked DEGs by Tukiainen XCI category


hdac_x_tukiainen_deg_summary_plot <- hdac_x_tukiainen_deg_summary %>%
  dplyr::mutate(
    contrast_short = factor(
      contrast_short,
      levels = c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
    ),
    deg_status_15fold = factor(
      deg_status_15fold,
      levels = c("Down", "Up")
    ),
    tukiainen_status_grouped = factor(
      tukiainen_status_grouped,
      levels = c(
        "Subject to XCI",
        "Escapee",
        "Variable escapee",
        "Unclassified"
      )
    )
  ) %>%
  tidyr::complete(
    contrast_short,
    deg_status_15fold,
    tukiainen_status_grouped,
    fill = list(n_genes = 0)
  )

hdac_tukiainen_count_plot <- ggplot(
  hdac_x_tukiainen_deg_summary_plot,
  aes(
    x = contrast_short,
    y = n_genes,
    fill = tukiainen_status_grouped
  )
) +
  geom_col(
    colour = "black",
    linewidth = 0.25
  ) +
  facet_wrap(
    ~ deg_status_15fold,
    nrow = 1
  ) +
  scale_fill_grey(
    start = 0.25,
    end = 0.85
  ) +
  theme_bw() +
  labs(
    title = "Tukiainen XCI categories among HDAC-responsive X-linked genes",
    subtitle = "X-linked DEGs defined using padj < 0.05 and |log2FC| ≥ log2(1.5)",
    x = "HDAC perturbation",
    y = "Number of X-linked DEGs",
    fill = "Tukiainen category"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )

hdac_tukiainen_count_plot

ggsave(
  filename = file.path(xci_output_dir, "HDAC_Xlinked_DEGs_Tukiainen_category_counts.png"),
  plot = hdac_tukiainen_count_plot,
  width = 9,
  height = 5,
  dpi = 300,
  bg = "white"
)


## HDAC X-linked DEG category proportions


hdac_x_tukiainen_deg_prop_plot <- hdac_x_tukiainen_deg_summary_plot %>%
  dplyr::group_by(
    contrast_short,
    deg_status_15fold
  ) %>%
  dplyr::mutate(
    total_degs = sum(n_genes),
    proportion = n_genes / total_degs
  ) %>%
  dplyr::ungroup()

hdac_tukiainen_prop_plot <- ggplot(
  hdac_x_tukiainen_deg_prop_plot,
  aes(
    x = contrast_short,
    y = proportion,
    fill = tukiainen_status_grouped
  )
) +
  geom_col(
    colour = "black",
    linewidth = 0.25
  ) +
  facet_wrap(
    ~ deg_status_15fold,
    nrow = 1
  ) +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  scale_fill_grey(
    start = 0.25,
    end = 0.85
  ) +
  theme_bw() +
  labs(
    title = "Relative composition of HDAC-responsive X-linked genes",
    subtitle = "Proportion of X-linked DEGs in each Tukiainen XCI category",
    x = "HDAC perturbation",
    y = "Proportion of X-linked DEGs",
    fill = "Tukiainen category"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )

hdac_tukiainen_prop_plot

ggsave(
  filename = file.path(xci_output_dir, "HDAC_Xlinked_DEGs_Tukiainen_category_proportions.png"),
  plot = hdac_tukiainen_prop_plot,
  width = 9,
  height = 5,
  dpi = 300,
  bg = "white"
)


## Plot capacitation X-linked DEGs by Tukiainen XCI category


cap_x_tukiainen_deg_summary_plot <- cap_x_tukiainen_deg_summary %>%
  dplyr::mutate(
    contrast_short = factor(
      contrast_short,
      levels = c(
        "Day3 vs Naive",
        "Day7 vs Naive",
        "Primed vs Naive",
        "Day7 vs Day3",
        "Primed vs Day7"
      )
    ),
    deg_status_15fold = factor(
      deg_status_15fold,
      levels = c("Down", "Up")
    ),
    tukiainen_status_grouped = factor(
      tukiainen_status_grouped,
      levels = c(
        "Subject to XCI",
        "Escapee",
        "Variable escapee",
        "Unclassified"
      )
    )
  ) %>%
  tidyr::complete(
    contrast_short,
    deg_status_15fold,
    tukiainen_status_grouped,
    fill = list(n_genes = 0)
  )

cap_tukiainen_count_plot <- ggplot(
  cap_x_tukiainen_deg_summary_plot,
  aes(
    x = contrast_short,
    y = n_genes,
    fill = tukiainen_status_grouped
  )
) +
  geom_col(
    colour = "black",
    linewidth = 0.25
  ) +
  facet_wrap(
    ~ deg_status_15fold,
    nrow = 1
  ) +
  scale_fill_grey(
    start = 0.25,
    end = 0.85
  ) +
  theme_bw() +
  labs(
    title = "Tukiainen XCI categories among capacitation-responsive X-linked genes",
    subtitle = "X-linked DEGs defined using padj < 0.05 and |log2FC| ≥ log2(1.5)",
    x = "Capacitation comparison",
    y = "Number of X-linked DEGs",
    fill = "Tukiainen category"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )

cap_tukiainen_count_plot

ggsave(
  filename = file.path(xci_output_dir, "Capacitation_Xlinked_DEGs_Tukiainen_category_counts.png"),
  plot = cap_tukiainen_count_plot,
  width = 11,
  height = 5,
  dpi = 300,
  bg = "white"
)


## Capacitation X-linked DEG category proportions


cap_x_tukiainen_deg_prop_plot <- cap_x_tukiainen_deg_summary_plot %>%
  dplyr::group_by(
    contrast_short,
    deg_status_15fold
  ) %>%
  dplyr::mutate(
    total_degs = sum(n_genes),
    proportion = n_genes / total_degs
  ) %>%
  dplyr::ungroup()

cap_tukiainen_prop_plot <- ggplot(
  cap_x_tukiainen_deg_prop_plot,
  aes(
    x = contrast_short,
    y = proportion,
    fill = tukiainen_status_grouped
  )
) +
  geom_col(
    colour = "black",
    linewidth = 0.25
  ) +
  facet_wrap(
    ~ deg_status_15fold,
    nrow = 1
  ) +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  scale_fill_grey(
    start = 0.25,
    end = 0.85
  ) +
  theme_bw() +
  labs(
    title = "Relative composition of capacitation-responsive X-linked genes",
    subtitle = "Proportion of X-linked DEGs in each Tukiainen XCI category",
    x = "Capacitation comparison",
    y = "Proportion of X-linked DEGs",
    fill = "Tukiainen category"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )

cap_tukiainen_prop_plot

ggsave(
  filename = file.path(xci_output_dir, "Capacitation_Xlinked_DEGs_Tukiainen_category_proportions.png"),
  plot = cap_tukiainen_prop_plot,
  width = 11,
  height = 5,
  dpi = 300,
  bg = "white"
)


## Tukiainen-annotated X-linked log2FC heatmaps


tukiainen_heatmap_dir <- file.path(
  xci_output_dir,
  "Tukiainen_log2FC_heatmaps"
)

dir.create(
  tukiainen_heatmap_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


## HDAC X-linked log2FC heatmap by Tukiainen category


hdac_tukiainen_heatmap_long <- hdac_tukiainen_annotated %>%
  dplyr::filter(
    CHR == "X",
    tukiainen_status_grouped %in% c("Subject to XCI", "Escapee")
  ) %>%
  dplyr::mutate(
    contrast_short = factor(
      contrast_short,
      levels = c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
    ),
    tukiainen_status_grouped = factor(
      tukiainen_status_grouped,
      levels = c(
        "Subject to XCI",
        "Escapee"
      )
    ),
    gene_plot_label = dplyr::case_when(
      !is.na(SYMBOL) & SYMBOL != "" ~ SYMBOL,
      !is.na(gene_label) & gene_label != "" ~ gene_label,
      TRUE ~ gene_id
    ),
    is_sig = deg_status_15fold != "Not significant",
    neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)),
    ranking_score = dplyr::if_else(
      is_sig,
      abs(log2FoldChange) * neg_log10_padj,
      0
    )
  )


## Select top HDAC-responsive X-linked genes per category
############################################################
## Group heatmap genes by dominant HDAC response
############################################################

top_n_per_response_group <- 10

hdac_gene_primary_response <- hdac_tukiainen_heatmap_long %>%
  dplyr::filter(
    is_sig,
    tukiainen_status_grouped %in% c("Subject to XCI", "Escapee")
  ) %>%
  dplyr::group_by(
    gene_plot_label,
    tukiainen_status_grouped
  ) %>%
  dplyr::slice_max(
    order_by = ranking_score,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    primary_response = paste0(
      as.character(contrast_short),
      " ",
      dplyr::case_when(
        deg_status_15fold == "Up" ~ "upregulated",
        deg_status_15fold == "Down" ~ "downregulated",
        TRUE ~ "other"
      )
    ),
    primary_response = factor(
      primary_response,
      levels = c(
        "HDAC1 KD upregulated",
        "HDAC1 KD downregulated",
        "HDAC2 KD upregulated",
        "HDAC2 KD downregulated",
        "HDAC3 KD upregulated",
        "HDAC3 KD downregulated"
      )
    ),
    tukiainen_short = dplyr::case_when(
      tukiainen_status_grouped == "Subject to XCI" ~ "S",
      tukiainen_status_grouped == "Escapee" ~ "E",
      TRUE ~ "U"
    ),
    gene_display_label = paste0(gene_plot_label, " [", tukiainen_short, "]")
  ) %>%
  dplyr::filter(!is.na(primary_response))

hdac_tukiainen_heatmap_genes <- hdac_gene_primary_response %>%
  dplyr::group_by(primary_response) %>%
  dplyr::slice_max(
    order_by = ranking_score,
    n = top_n_per_response_group,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    primary_response,
    dplyr::desc(ranking_score)
  )

hdac_tukiainen_heatmap_genes %>%
  dplyr::count(primary_response, tukiainen_status_grouped)

############################################################
## Prepare grouped heatmap plotting table
############################################################

hdac_gene_order <- hdac_tukiainen_heatmap_genes %>%
  dplyr::arrange(
    primary_response,
    dplyr::desc(ranking_score)
  ) %>%
  dplyr::pull(gene_display_label)

hdac_tukiainen_heatmap_plot_df <- hdac_tukiainen_heatmap_long %>%
  dplyr::filter(
    tukiainen_status_grouped %in% c("Subject to XCI", "Escapee")
  ) %>%
  dplyr::inner_join(
    hdac_tukiainen_heatmap_genes %>%
      dplyr::select(
        gene_plot_label,
        tukiainen_status_grouped,
        primary_response,
        tukiainen_short,
        gene_display_label
      ),
    by = c("gene_plot_label", "tukiainen_status_grouped")
  ) %>%
  dplyr::mutate(
    gene_display_label = factor(
      gene_display_label,
      levels = rev(unique(hdac_gene_order))
    ),
    log2FC_plot = pmax(
      pmin(log2FoldChange, 2.5),
      -2.5
    )
  )

## Plot HDAC Tukiainen heatmap

############################################################
## Plot grouped HDAC Tukiainen heatmap
############################################################

hdac_tukiainen_log2FC_heatmap_grouped <- ggplot(
  hdac_tukiainen_heatmap_plot_df,
  aes(
    x = contrast_short,
    y = gene_display_label,
    fill = log2FC_plot
  )
) +
  geom_tile(
    colour = "grey85",
    linewidth = 0.25
  ) +
  geom_point(
    data = hdac_tukiainen_heatmap_plot_df %>%
      dplyr::filter(is_sig),
    aes(
      x = contrast_short,
      y = gene_display_label
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.4,
    fill = "black",
    colour = "black"
  ) +
  facet_grid(
    primary_response ~ .,
    scales = "free_y",
    space = "free_y",
    labeller = ggplot2::label_wrap_gen(width = 24)
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish
  ) +
  theme_bw() +
  labs(
    title = "HDAC-responsive X-linked genes grouped by dominant knockdown response",
    subtitle = "Top X-linked DEGs per group; [S] = Subject to XCI, [E] = Escapee; black dots = significant DEGs",
    x = "HDAC perturbation",
    y = "X-linked gene",
    fill = "log2FC"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text.y = element_text(face = "bold", angle = 0),
    panel.spacing.y = unit(0.5, "lines"),
    panel.grid = element_blank()
  )

hdac_tukiainen_log2FC_heatmap_grouped

hdac_heatmap_height <- max(
  7,
  0.2 * dplyr::n_distinct(hdac_tukiainen_heatmap_plot_df$gene_display_label) + 2
)

ggsave(
  filename = file.path(
    tukiainen_heatmap_dir,
    "HDAC_Tukiainen_Xlinked_log2FC_heatmap_grouped_by_response.png"
  ),
  plot = hdac_tukiainen_log2FC_heatmap_grouped,
  width = 8,
  height = hdac_heatmap_height,
  dpi = 300,
  bg = "white"
)
write.csv(
  hdac_tukiainen_heatmap_plot_df,
  file.path(
    tukiainen_heatmap_dir,
    "HDAC_Tukiainen_Xlinked_log2FC_heatmap_table.csv"
  ),
  row.names = FALSE
)


## Capacitation X-linked log2FC heatmap by Tukiainen category


cap_tukiainen_heatmap_long <- cap_tukiainen_annotated %>%
  dplyr::filter(
    CHR == "X",
    tukiainen_status_grouped %in% c("Subject to XCI", "Escapee")
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
    tukiainen_status_grouped = factor(
      tukiainen_status_grouped,
      levels = c(
        "Subject to XCI",
        "Escapee",
        "Variable escapee",
        "Unclassified"
      )
    ),
    gene_plot_label = dplyr::case_when(
      !is.na(SYMBOL) & SYMBOL != "" ~ SYMBOL,
      !is.na(gene_label) & gene_label != "" ~ gene_label,
      TRUE ~ count_gene_id
    ),
    is_sig = deg_status_15fold != "Not significant",
    neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)),
    ranking_score = dplyr::if_else(
      is_sig,
      abs(log2FoldChange) * neg_log10_padj,
      0
    )
  )

top_n_per_tukiainen_category <- 20

## Remove very large log2FC outliers 
cap_outlier_log2fc_cutoff <- 3

## Tighter heatmap colour scale
cap_heatmap_colour_limit <- 1.5

## Select top capacitation-responsive X-linked genes per category

cap_tukiainen_heatmap_genes <- cap_tukiainen_heatmap_long %>%
  dplyr::filter(is_sig) %>%
  dplyr::group_by(
    gene_plot_label,
    tukiainen_status_grouped
  ) %>%
  dplyr::summarise(
    max_ranking_score = max(ranking_score, na.rm = TRUE),
    max_abs_log2FC = max(abs(log2FoldChange), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    max_abs_log2FC <= cap_outlier_log2fc_cutoff
  ) %>%
  dplyr::group_by(tukiainen_status_grouped) %>%
  dplyr::slice_max(
    order_by = max_ranking_score,
    n = top_n_per_tukiainen_category,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    tukiainen_status_grouped,
    dplyr::desc(max_ranking_score)
  )

cap_tukiainen_heatmap_genes %>%
  dplyr::count(tukiainen_status_grouped)


## Prepare capacitation heatmap plotting table


cap_gene_order <- cap_tukiainen_heatmap_genes %>%
  dplyr::arrange(
    tukiainen_status_grouped,
    dplyr::desc(max_ranking_score)
  ) %>%
  dplyr::pull(gene_plot_label)

cap_tukiainen_heatmap_plot_df <- cap_tukiainen_heatmap_long %>%
  dplyr::semi_join(
    cap_tukiainen_heatmap_genes,
    by = c("gene_plot_label", "tukiainen_status_grouped")
  ) %>%
  dplyr::mutate(
    gene_plot_label = factor(
      gene_plot_label,
      levels = rev(unique(cap_gene_order))
    ),
    log2FC_plot = pmax(
      pmin(log2FoldChange, cap_heatmap_colour_limit),
      -cap_heatmap_colour_limit
    )
  )


## Plot capacitation Tukiainen heatmap


cap_tukiainen_log2FC_heatmap <- ggplot(
  cap_tukiainen_heatmap_plot_df,
  aes(
    x = contrast_short,
    y = gene_plot_label,
    fill = log2FC_plot
  )
) +
  geom_tile(
    colour = "grey85",
    linewidth = 0.25
  ) +
  geom_point(
    data = cap_tukiainen_heatmap_plot_df %>%
      dplyr::filter(is_sig),
    aes(
      x = contrast_short,
      y = gene_plot_label
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.4,
    fill = "black",
    colour = "black"
  ) +
  facet_grid(
    tukiainen_status_grouped ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-cap_heatmap_colour_limit, cap_heatmap_colour_limit),
    oob = scales::squish
  ) +
  theme_bw() +
  labs(
    title = "Capacitation-responsive X-linked genes by Tukiainen XCI category",
    subtitle = "Top X-linked DEGs per category after excluding extreme log2FC outliers; black dots indicate significant DEGs",
    x = "Capacitation comparison",
    y = "X-linked gene",
    fill = "log2FC"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text.y = element_text(face = "bold", angle = 0),
    panel.spacing.y = unit(0.6, "lines")
  )

cap_tukiainen_log2FC_heatmap

cap_heatmap_height <- max(
  7,
  0.18 * dplyr::n_distinct(cap_tukiainen_heatmap_plot_df$gene_plot_label) + 2
)

ggsave(
  filename = file.path(
    tukiainen_heatmap_dir,
    "Capacitation_Tukiainen_SubjectToXCI_Escapee_log2FC_heatmap_no_outliers_rescaled.png"
  ),
  plot = cap_tukiainen_log2FC_heatmap,
  width = 9,
  height = cap_heatmap_height,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    tukiainen_heatmap_dir,
    "Capacitation_Tukiainen_Xlinked_log2FC_heatmap.png"
  ),
  plot = cap_tukiainen_log2FC_heatmap,
  width = 9,
  height = cap_heatmap_height,
  dpi = 300,
  bg = "white"
)

write.csv(
  cap_tukiainen_heatmap_plot_df,
  file.path(
    tukiainen_heatmap_dir,
    "Capacitation_Tukiainen_Xlinked_log2FC_heatmap_table.csv"
  ),
  row.names = FALSE
)

hdac_tukiainen_heatmap_genes %>%
  dplyr::count(tukiainen_status_grouped)

cap_tukiainen_heatmap_genes %>%
  dplyr::count(tukiainen_status_grouped)



## Tukiainen category enrichment among X-linked DEGs


run_tukiainen_enrichment <- function(df, gene_col) {
  
  tukiainen_categories_to_test <- c(
    "Subject to XCI",
    "Escapee",
    "Variable escapee"
  )
  
  x_df <- df %>%
    dplyr::filter(CHR == "X") %>%
    dplyr::mutate(
      gene_key = as.character(.data[[gene_col]]),
      contrast_short = as.character(contrast_short),
      deg_status_15fold = as.character(deg_status_15fold),
      tukiainen_status_grouped = as.character(tukiainen_status_grouped)
    ) %>%
    dplyr::filter(
      tukiainen_status_grouped %in% tukiainen_categories_to_test,
      !is.na(gene_key),
      gene_key != ""
    ) %>%
    dplyr::distinct(
      contrast_short,
      gene_key,
      .keep_all = TRUE
    )
  
  test_grid <- expand.grid(
    contrast_short = unique(x_df$contrast_short),
    deg_direction = c("Down", "Up"),
    tukiainen_status_grouped = tukiainen_categories_to_test,
    stringsAsFactors = FALSE
  )
  
  enrichment_list <- lapply(seq_len(nrow(test_grid)), function(i) {
    
    this_contrast <- test_grid$contrast_short[i]
    this_direction <- test_grid$deg_direction[i]
    this_category <- test_grid$tukiainen_status_grouped[i]
    
    sub_df <- x_df %>%
      dplyr::filter(contrast_short == this_contrast)
    
    is_deg_direction <- sub_df$deg_status_15fold == this_direction
    is_category <- sub_df$tukiainen_status_grouped == this_category
    
    a <- sum(is_deg_direction & is_category)
    b <- sum(is_deg_direction & !is_category)
    c <- sum(!is_deg_direction & is_category)
    d <- sum(!is_deg_direction & !is_category)
    
    fisher_mat <- matrix(
      c(a, b, c, d),
      nrow = 2,
      byrow = TRUE
    )
    
    fisher_res <- stats::fisher.test(fisher_mat)
    
    data.frame(
      contrast_short = this_contrast,
      deg_direction = this_direction,
      tukiainen_status_grouped = this_category,
      DEG_in_category = a,
      DEG_not_category = b,
      nonDEG_in_category = c,
      nonDEG_not_category = d,
      odds_ratio = unname(fisher_res$estimate),
      p_value = fisher_res$p.value,
      stringsAsFactors = FALSE
    )
  })
  
  dplyr::bind_rows(enrichment_list) %>%
    dplyr::mutate(
      padj = p.adjust(p_value, method = "BH"),
      enrichment_status = dplyr::case_when(
        padj < 0.05 & odds_ratio > 1 ~ "Enriched",
        padj < 0.05 & odds_ratio < 1 ~ "Depleted",
        TRUE ~ "Not significant"
      ),
      log2_odds_ratio = log2(odds_ratio)
    )
}

hdac_tukiainen_enrichment <- run_tukiainen_enrichment(
  df = hdac_tukiainen_annotated,
  gene_col = "gene_id"
)

cap_tukiainen_enrichment <- run_tukiainen_enrichment(
  df = cap_tukiainen_annotated,
  gene_col = "count_gene_id"
)

hdac_tukiainen_enrichment %>%
  dplyr::arrange(padj)

cap_tukiainen_enrichment %>%
  dplyr::arrange(padj)

hdac_tukiainen_enrichment %>%
  dplyr::filter(padj < 0.05)

cap_tukiainen_enrichment %>%
  dplyr::filter(padj < 0.05)

write.csv(
  hdac_tukiainen_enrichment,
  file.path(xci_output_dir, "HDAC_Tukiainen_category_enrichment.csv"),
  row.names = FALSE
)

write.csv(
  cap_tukiainen_enrichment,
  file.path(xci_output_dir, "Capacitation_Tukiainen_category_enrichment.csv"),
  row.names = FALSE
)



## Plot Tukiainen enrichment results


plot_tukiainen_enrichment <- function(enrichment_df, plot_title, output_file, width = 10, height = 5) {
  
  enrichment_plot_df <- enrichment_df %>%
    dplyr::mutate(
      tukiainen_status_grouped = factor(
        tukiainen_status_grouped,
        levels = c(
          "Subject to XCI",
          "Escapee",
          "Variable escapee"
        )
      ),
      deg_direction = factor(
        deg_direction,
        levels = c("Down", "Up")
      ),
      neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)),
      log2_odds_ratio_plot = dplyr::case_when(
        is.infinite(log2_odds_ratio) & log2_odds_ratio > 0 ~ 4,
        is.infinite(log2_odds_ratio) & log2_odds_ratio < 0 ~ -4,
        TRUE ~ log2_odds_ratio
      ),
      log2_odds_ratio_plot = pmax(
        pmin(log2_odds_ratio_plot, 4),
        -4
      ),
      significant = padj < 0.05
    )
  
  enrichment_plot <- ggplot(
    enrichment_plot_df,
    aes(
      x = log2_odds_ratio_plot,
      y = tukiainen_status_grouped
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey50"
    ) +
    geom_point(
      aes(
        size = neg_log10_padj,
        shape = significant
      ),
      fill = "grey70",
      colour = "black"
    ) +
    facet_grid(
      deg_direction ~ contrast_short
    ) +
    scale_shape_manual(
      values = c(
        "FALSE" = 21,
        "TRUE" = 24
      )
    ) +
    theme_bw() +
    labs(
      title = plot_title,
      subtitle = "Positive log2 odds ratio indicates enrichment; BH-adjusted Fisher's exact test",
      x = "log2 odds ratio",
      y = "Tukiainen XCI category",
      size = "-log10 adjusted p-value",
      shape = "padj < 0.05"
    ) +
    theme(
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  
  print(enrichment_plot)
  
  ggsave(
    filename = output_file,
    plot = enrichment_plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
  
  return(enrichment_plot)
}

hdac_tukiainen_enrichment_plot <- plot_tukiainen_enrichment(
  enrichment_df = hdac_tukiainen_enrichment,
  plot_title = "Enrichment of Tukiainen XCI categories among HDAC-responsive X-linked genes",
  output_file = file.path(xci_output_dir, "HDAC_Tukiainen_category_enrichment_dotplot.png"),
  width = 10,
  height = 5
)

cap_tukiainen_enrichment_plot <- plot_tukiainen_enrichment(
  enrichment_df = cap_tukiainen_enrichment,
  plot_title = "Enrichment of Tukiainen XCI categories among capacitation-responsive X-linked genes",
  output_file = file.path(xci_output_dir, "Capacitation_Tukiainen_category_enrichment_dotplot.png"),
  width = 12,
  height = 5
)

hdac_tukiainen_enrichment %>%
  dplyr::arrange(padj)

cap_tukiainen_enrichment %>%
  dplyr::arrange(padj)
