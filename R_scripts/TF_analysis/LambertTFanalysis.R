## Lambert et al. 2018 TF analysis
## Import Human TF list and prepare for HDAC comparison

## 1. Libraries

library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)

## 2. Paths

base_dir <- "/home/jvk3/Desktop/HDAC_counts"

all_res_dir <- file.path(base_dir, "All_Res_Tables")

lambert_dir <- file.path(base_dir, "Lambert_TF")

lambert_plot_dir <- file.path(base_dir, "Lambert_TF_outputs")

dir.create(lambert_plot_dir, showWarnings = FALSE, recursive = TRUE)


## 3. Check sheet names


lambert_file <- list.files(
  lambert_dir,
  pattern = "\\.xlsx$",
  full.names = TRUE
)

lambert_file


lambert_sheets <- readxl::excel_sheets(lambert_file)

lambert_sheets


## 4. Read Lambert Table S1

lambert_tf_raw <- readxl::read_excel(
  lambert_file,
  sheet = "Table S1. Related to Figure 1B",
  skip = 1
)

dim(lambert_tf_raw)

colnames(lambert_tf_raw)

head(lambert_tf_raw)

## 5. Clean Lambert Table S1

lambert_tf_clean <- lambert_tf_raw %>%
  dplyr::rename(
    lambert_ensembl_id = ID,
    lambert_symbol = Name,
    TF_family = DBD,
    is_TF_raw = `...4`,
    TF_assessment = `TF assessment`,
    binding_mode = `Binding mode`,
    motif_status = `Motif status`
  ) %>%
  dplyr::mutate(
    gene_id = gsub("\\..*$", "", lambert_ensembl_id),
    is_Lambert_TF = dplyr::case_when(
      is_TF_raw == "Yes" ~ TRUE,
      is_TF_raw == "No" ~ FALSE,
      TRUE ~ NA
    )
  ) %>%
  dplyr::select(
    gene_id,
    lambert_ensembl_id,
    lambert_symbol,
    is_Lambert_TF,
    TF_family,
    TF_assessment,
    binding_mode,
    motif_status
  ) %>%
  dplyr::distinct(gene_id, .keep_all = TRUE)

dim(lambert_tf_clean)

lambert_tf_clean %>%
  dplyr::count(is_Lambert_TF)

lambert_tf_clean %>%
  dplyr::filter(is_Lambert_TF == TRUE) %>%
  dplyr::count(TF_assessment, sort = TRUE)


## 6. Load combined HDAC all_res table

all_res_hdac1_2_3_15fold <- read.csv(
  file.path(all_res_dir, "all_results_HDAC1_2_3_combined_15fold.csv"),
  stringsAsFactors = FALSE
)

all_res_hdac1_2_3_15fold <- all_res_hdac1_2_3_15fold %>%
  dplyr::mutate(
    gene_id = gsub("\\..*$", "", gene_id),
    CHR = as.character(CHR)
  )

## 7. Join Lambert TF annotation to HDAC results

all_res_hdac1_2_3_TF_annotated <- all_res_hdac1_2_3_15fold %>%
  dplyr::left_join(
    lambert_tf_clean,
    by = "gene_id"
  ) %>%
  dplyr::mutate(
    is_Lambert_TF = dplyr::if_else(
      is.na(is_Lambert_TF),
      FALSE,
      is_Lambert_TF
    ),
    TF_label = dplyr::case_when(
      is_Lambert_TF == TRUE ~ "Lambert TF",
      TRUE ~ "Not Lambert TF"
    )
  )

dim(all_res_hdac1_2_3_TF_annotated)

all_res_hdac1_2_3_TF_annotated %>%
  dplyr::count(TF_label)


## 8. Check HDAC-responsive transcription factors

tf_deg_summary <- all_res_hdac1_2_3_TF_annotated %>%
  dplyr::filter(
    is_Lambert_TF == TRUE,
    deg_status_15fold != "Not significant"
  ) %>%
  dplyr::count(
    contrast_short,
    deg_status_15fold,
    name = "n_TF_DEGs"
  )

tf_deg_summary

xlinked_tf_deg_summary <- all_res_hdac1_2_3_TF_annotated %>%
  dplyr::filter(
    is_Lambert_TF == TRUE,
    CHR == "X",
    deg_status_15fold != "Not significant"
  ) %>%
  dplyr::count(
    contrast_short,
    deg_status_15fold,
    name = "n_Xlinked_TF_DEGs"
  )

xlinked_tf_deg_summary

## 9. Save Lambert-annotated HDAC table

write.csv(
  all_res_hdac1_2_3_TF_annotated,
  file.path(lambert_plot_dir, "all_results_HDAC1_2_3_Lambert_TF_annotated_15fold.csv"),
  row.names = FALSE
)

write.csv(
  lambert_tf_clean,
  file.path(lambert_plot_dir, "Lambert_TF_cleaned_TableS1.csv"),
  row.names = FALSE
)

write.csv(
  tf_deg_summary,
  file.path(lambert_plot_dir, "Lambert_TF_DEG_summary_genomewide.csv"),
  row.names = FALSE
)

write.csv(
  xlinked_tf_deg_summary,
  file.path(lambert_plot_dir, "Lambert_TF_DEG_summary_chrX.csv"),
  row.names = FALSE
)

## List genome-wide HDAC-responsive Lambert TFs

tf_deg_list <- all_res_hdac1_2_3_TF_annotated %>%
  dplyr::filter(
    is_Lambert_TF == TRUE,
    deg_status_15fold != "Not significant"
  ) %>%
  dplyr::select(
    contrast_short,
    gene_id,
    SYMBOL,
    CHR,
    TF_family,
    TF_assessment,
    motif_status,
    log2FoldChange,
    padj,
    control_mean_norm_count,
    treatment_mean_norm_count,
    deg_status_15fold
  ) %>%
  dplyr::arrange(
    contrast_short,
    deg_status_15fold,
    padj
  )

tf_deg_list


## List X-linked HDAC-responsive Lambert TFs

xlinked_tf_deg_list <- all_res_hdac1_2_3_TF_annotated %>%
  dplyr::filter(
    is_Lambert_TF == TRUE,
    CHR == "X",
    deg_status_15fold != "Not significant"
  ) %>%
  dplyr::select(
    contrast_short,
    gene_id,
    SYMBOL,
    CHR,
    START,
    END,
    TF_family,
    TF_assessment,
    motif_status,
    log2FoldChange,
    padj,
    control_mean_norm_count,
    treatment_mean_norm_count,
    deg_status_15fold
  ) %>%
  dplyr::arrange(
    contrast_short,
    deg_status_15fold,
    padj
  )

xlinked_tf_deg_list

write.csv(
  tf_deg_list,
  file.path(lambert_plot_dir, "Lambert_TF_DEG_list_genomewide.csv"),
  row.names = FALSE
)

write.csv(
  xlinked_tf_deg_list,
  file.path(lambert_plot_dir, "Lambert_TF_DEG_list_chrX.csv"),
  row.names = FALSE
)


## Plot Lambert TF DEG counts


tf_deg_summary_plot <- ggplot(
  tf_deg_summary,
  aes(
    x = contrast_short,
    y = n_TF_DEGs,
    fill = deg_status_15fold
  )
) +
  geom_col(position = "dodge", colour = "black") +
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
    title = "Differentially expressed transcription factors after HDAC knockdown",
    subtitle = "TFs defined using Lambert et al. 2018 Table S1",
    x = "Knockdown condition",
    y = "Number of TF DEGs",
    fill = "Direction"
  )

tf_deg_summary_plot

ggsave(
  filename = file.path(lambert_plot_dir, "Lambert_TF_DEG_counts_genomewide.png"),
  plot = tf_deg_summary_plot,
  width = 7,
  height = 5,
  dpi = 300,
  bg = "white"
)

## Plot X-linked Lambert TF DEG counts

xlinked_tf_plot_df <- all_res_hdac1_2_3_TF_annotated %>%
  as.data.frame() %>%
  dplyr::filter(
    .data$is_Lambert_TF == TRUE,
    .data$CHR == "X",
    .data$deg_status_15fold != "Not significant"
  ) %>%
  dplyr::mutate(
    gene_plot_label = dplyr::case_when(
      !is.na(.data$SYMBOL) & .data$SYMBOL != "" ~ .data$SYMBOL,
      TRUE ~ .data$gene_id
    ),
    contrast_short = factor(
      .data$contrast_short,
      levels = c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
    )
  )

xlinked_tf_plot_df

xlinked_tf_lollipop_plot <- ggplot(
  xlinked_tf_plot_df,
  aes(
    x = log2FoldChange,
    y = reorder(gene_plot_label, log2FoldChange),
    colour = contrast_short
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = log2FoldChange,
      y = gene_plot_label,
      yend = gene_plot_label
    ),
    linewidth = 0.7
  ) +
  geom_point(size = 3) +
  facet_wrap(
    ~ contrast_short,
    scales = "free_y"
  ) +
  theme_bw() +
  labs(
    title = "X-linked transcription factors altered after HDAC knockdown",
    subtitle = "Lambert et al. 2018 transcription factor catalogue",
    x = "Shrunken log2 fold change",
    y = "X-linked transcription factors",
    colour = "Knockdown"
  )

xlinked_tf_lollipop_plot

ggsave(
  filename = file.path(
    lambert_plot_dir,
    "Xlinked_Lambert_TFs_lollipop.png"
  ),
  plot = xlinked_tf_lollipop_plot,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)


## Top changing Lambert TFs after HDAC knockdown

top_tf_plot_df <- all_res_hdac1_2_3_TF_annotated %>%
  dplyr::filter(
    is_Lambert_TF == TRUE,
    deg_status_15fold != "Not significant",
    !is.na(log2FoldChange),
    !is.na(padj)
  ) %>%
  dplyr::mutate(
    gene_plot_label = dplyr::case_when(
      !is.na(SYMBOL) & SYMBOL != "" ~ SYMBOL,
      TRUE ~ gene_id
    ),
    neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)),
    ranking_score = abs(log2FoldChange) * neg_log10_padj,
    contrast_short = factor(
      contrast_short,
      levels = c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
    )
  ) %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::slice_max(
    order_by = ranking_score,
    n = 20,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

top_tf_plot <- ggplot(
  top_tf_plot_df,
  aes(
    x = log2FoldChange,
    y = reorder(gene_plot_label, ranking_score),
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
  theme_bw() +
  labs(
    title = "Top differentially expressed transcription factors after HDAC knockdown",
    subtitle = "TFs defined using Lambert et al. 2018 Table S1; ranked by -log10(padj) × absolute log2FC",
    x = "Shrunken log2 fold change",
    y = "Transcription factor",
    colour = "Direction",
    size = "-log10(padj)"
  )

top_tf_plot

ggsave(
  filename = file.path(lambert_plot_dir, "Top20_Lambert_TF_DEGs_per_HDAC.png"),
  plot = top_tf_plot,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

write.csv(
  top_tf_plot_df,
  file.path(lambert_plot_dir, "Top20_Lambert_TF_DEGs_per_HDAC.csv"),
  row.names = FALSE
)


## TF-family breakdown among HDAC-responsive TFs


tf_family_summary <- all_res_hdac1_2_3_TF_annotated %>%
  dplyr::filter(
    is_Lambert_TF == TRUE,
    deg_status_15fold != "Not significant"
  ) %>%
  dplyr::mutate(
    TF_family = dplyr::case_when(
      is.na(TF_family) | TF_family == "" ~ "Unknown/not specified",
      TRUE ~ TF_family
    ),
    contrast_short = factor(
      contrast_short,
      levels = c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
    )
  ) %>%
  dplyr::count(
    contrast_short,
    TF_family,
    deg_status_15fold,
    name = "n_TF_DEGs"
  )

tf_family_summary

top_tf_families <- tf_family_summary %>%
  dplyr::group_by(TF_family) %>%
  dplyr::summarise(
    total_TF_DEGs = sum(n_TF_DEGs),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(total_TF_DEGs)) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(TF_family)

tf_family_plot_df <- tf_family_summary %>%
  dplyr::mutate(
    TF_family_plot = dplyr::case_when(
      TF_family %in% top_tf_families ~ TF_family,
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::group_by(
    contrast_short,
    TF_family_plot,
    deg_status_15fold
  ) %>%
  dplyr::summarise(
    n_TF_DEGs = sum(n_TF_DEGs),
    .groups = "drop"
  )

tf_family_breakdown_plot <- ggplot(
  tf_family_plot_df,
  aes(
    x = TF_family_plot,
    y = n_TF_DEGs,
    fill = deg_status_15fold
  )
) +
  geom_col(
    position = "dodge",
    colour = "black"
  ) +
  facet_wrap(
    ~ contrast_short,
    scales = "free_y"
  ) +
  theme_bw() +
  labs(
    title = "TF-family composition of HDAC-responsive transcription factors",
    subtitle = "Top TF families among Lambert-defined TF DEGs",
    x = "TF family / DNA-binding domain class",
    y = "Number of TF DEGs",
    fill = "Direction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

tf_family_breakdown_plot

ggsave(
  filename = file.path(lambert_plot_dir, "Lambert_TF_family_breakdown.png"),
  plot = tf_family_breakdown_plot,
  width = 11,
  height = 6,
  dpi = 300,
  bg = "white"
)

write.csv(
  tf_family_summary,
  file.path(lambert_plot_dir, "Lambert_TF_family_summary_all_families.csv"),
  row.names = FALSE
)

write.csv(
  tf_family_plot_df,
  file.path(lambert_plot_dir, "Lambert_TF_family_summary_top10_plus_other.csv"),
  row.names = FALSE
)



## 11. doRothEA TF-target analysis


library(dorothea)
library(igraph)
library(tidygraph)
library(ggraph)

data(dorothea_hs, package = "dorothea")

regulons <- dorothea_hs

regulons_high <- regulons %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))

dim(regulons_high)

head(regulons_high)

## 12. Significant Lambert TF DEGs for doRothEA

tf_degs_for_dorothea <- all_res_hdac1_2_3_TF_annotated %>%
  dplyr::filter(
    is_Lambert_TF == TRUE,
    deg_status_15fold != "Not significant",
    !is.na(SYMBOL),
    SYMBOL != ""
  ) %>%
  dplyr::mutate(
    tf = SYMBOL,
    tf_gene_id = gene_id,
    tf_CHR = CHR,
    tf_START = START,
    tf_END = END,
    tf_log2FoldChange = log2FoldChange,
    tf_padj = padj,
    tf_direction = deg_status_15fold
  ) %>%
  dplyr::select(
    contrast_short,
    tf,
    tf_gene_id,
    tf_CHR,
    tf_START,
    tf_END,
    TF_family,
    TF_assessment,
    motif_status,
    tf_log2FoldChange,
    tf_padj,
    tf_direction
  ) %>%
  dplyr::distinct()

tf_degs_for_dorothea %>%
  dplyr::count(contrast_short, tf_direction)


## 13. Prepare target gene DE information


target_gene_info <- all_res_hdac1_2_3_TF_annotated %>%
  dplyr::filter(
    !is.na(SYMBOL),
    SYMBOL != ""
  ) %>%
  dplyr::mutate(
    target = SYMBOL,
    target_gene_id = gene_id,
    target_CHR = CHR,
    target_START = START,
    target_END = END,
    target_log2FoldChange = log2FoldChange,
    target_padj = padj,
    target_direction = deg_status_15fold,
    target_is_DEG = deg_status_15fold != "Not significant"
  ) %>%
  dplyr::select(
    contrast_short,
    target,
    target_gene_id,
    target_CHR,
    target_START,
    target_END,
    target_log2FoldChange,
    target_padj,
    target_direction,
    target_is_DEG
  ) %>%
  dplyr::distinct()

## 14. Link HDAC-responsive TFs to doRothEA targets

tf_target_edges_all <- tf_degs_for_dorothea %>%
  dplyr::inner_join(
    regulons_high,
    by = "tf"
  ) %>%
  dplyr::left_join(
    target_gene_info,
    by = c("contrast_short", "target")
  ) %>%
  dplyr::filter(
    !is.na(target_gene_id)
  ) %>%
  dplyr::mutate(
    tf_is_Xlinked = tf_CHR == "X",
    target_is_Xlinked = target_CHR == "X",
    x_relationship = dplyr::case_when(
      tf_is_Xlinked & target_is_Xlinked ~ "X-linked TF -> X-linked target",
      tf_is_Xlinked & !target_is_Xlinked ~ "X-linked TF -> non-X target",
      !tf_is_Xlinked & target_is_Xlinked ~ "Non-X TF -> X-linked target",
      TRUE ~ "Non-X TF -> non-X target"
    ),
    target_direction_clean = dplyr::case_when(
      target_direction == "Up" ~ "Up",
      target_direction == "Down" ~ "Down",
      TRUE ~ "Not significant"
    ),
    tf_target_direction = paste0(tf_direction, " TF / ", target_direction_clean, " target")
  )

dim(tf_target_edges_all)

tf_target_edges_all %>%
  dplyr::count(contrast_short, x_relationship)


## 15. X-associated TF-target links


tf_target_edges_Xassociated <- tf_target_edges_all %>%
  dplyr::filter(
    tf_is_Xlinked | target_is_Xlinked
  )

tf_target_edges_Xassociated %>%
  dplyr::count(contrast_short, x_relationship)

tf_target_edges_Xassociated_sig_targets <- tf_target_edges_Xassociated %>%
  dplyr::filter(
    target_is_DEG == TRUE
  )

tf_target_edges_Xassociated_sig_targets %>%
  dplyr::count(
    contrast_short,
    x_relationship,
    tf_direction,
    target_direction_clean
  )

## save lists
write.csv(
  tf_target_edges_all,
  file.path(lambert_plot_dir, "doRothEA_all_TF_target_edges_HDAC.csv"),
  row.names = FALSE
)

write.csv(
  tf_target_edges_Xassociated,
  file.path(lambert_plot_dir, "doRothEA_Xassociated_TF_target_edges_HDAC.csv"),
  row.names = FALSE
)

write.csv(
  tf_target_edges_Xassociated_sig_targets,
  file.path(lambert_plot_dir, "doRothEA_Xassociated_TF_target_edges_significant_targets_HDAC.csv"),
  row.names = FALSE
)


## 16. Summary plot of X-associated TF-target links


## Inspect X-associated significant TF-target links


colnames(tf_target_edges_Xassociated_sig_targets)

xassoc_sig_edge_detail <- tf_target_edges_Xassociated_sig_targets %>%
  dplyr::select(
    contrast_short,
    tf,
    target,
    confidence,
    mor,
    x_relationship,
    tf_CHR,
    target_CHR,
    tf_direction,
    tf_log2FoldChange,
    tf_padj,
    target_direction_clean,
    target_log2FoldChange,
    target_padj
  ) %>%
  dplyr::mutate(
    regulation_type = dplyr::case_when(
      mor == 1 ~ "Activating",
      mor == -1 ~ "Repressing",
      TRUE ~ "Unknown"
    )
  ) %>%
  dplyr::arrange(
    contrast_short,
    tf,
    target
  )

xassoc_sig_edge_detail


## Plot X-associated significant TF-target links as listed heatmap


xassoc_sig_edge_long <- xassoc_sig_edge_detail %>%
  dplyr::mutate(
    pair_label = paste0(
      tf, " \u2192 ", target,
      "\n",
      x_relationship,
      " | ", regulation_type,
      " | confidence ", confidence
    )
  ) %>%
  dplyr::select(
    contrast_short,
    pair_label,
    tf,
    target,
    tf_log2FoldChange,
    target_log2FoldChange,
    tf_direction,
    target_direction_clean
  ) %>%
  tidyr::pivot_longer(
    cols = c(tf_log2FoldChange, target_log2FoldChange),
    names_to = "node_type",
    values_to = "log2FC"
  ) %>%
  dplyr::mutate(
    node_type = dplyr::recode(
      node_type,
      "tf_log2FoldChange" = "TF",
      "target_log2FoldChange" = "Target"
    ),
    node_name = dplyr::case_when(
      node_type == "TF" ~ tf,
      node_type == "Target" ~ target
    ),
    node_type = factor(node_type, levels = c("TF", "Target"))
  )

xassoc_sig_edge_plot <- ggplot(
  xassoc_sig_edge_long,
  aes(
    x = node_type,
    y = pair_label,
    fill = log2FC
  )
) +
  geom_tile(
    colour = "black",
    linewidth = 0.4,
    height = 0.75
  ) +
  geom_text(
    aes(
      label = paste0(node_name, "\n", round(log2FC, 2))
    ),
    size = 3.2
  ) +
  scale_fill_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0,
    name = "log2FC"
  ) +
  facet_wrap(
    ~ contrast_short,
    scales = "free_y"
  ) +
  theme_bw() +
  labs(
    title = "X-associated TF-target links after HDAC knockdown",
    subtitle = "HDAC-responsive Lambert TFs linked to significant doRothEA A/B/C targets",
    x = NULL,
    y = "TF-target relationship"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )

xassoc_sig_edge_plot

ggsave(
  filename = file.path(
    lambert_plot_dir,
    "doRothEA_Xassociated_significant_TF_target_links_listed.png"
  ),
  plot = xassoc_sig_edge_plot,
  width = 9,
  height = 4.5,
  dpi = 300,
  bg = "white"
)


## 17. Function to plot TF-target networks


plot_tf_target_network <- function(edge_df, contrast_to_plot, output_name) {
  
  plot_edges <- edge_df %>%
    dplyr::filter(
      contrast_short == contrast_to_plot
    ) %>%
    dplyr::mutate(
      tf_node = paste0("TF_", tf),
      target_node = paste0("Target_", target)
    )
  
  if (nrow(plot_edges) == 0) {
    message("No edges to plot for ", contrast_to_plot)
    return(NULL)
  }
  
  tf_nodes <- plot_edges %>%
    dplyr::select(
      node_id = tf_node,
      label = tf,
      log2FoldChange = tf_log2FoldChange,
      direction = tf_direction
    ) %>%
    dplyr::mutate(
      node_type = "TF"
    ) %>%
    dplyr::distinct()
  
  target_nodes <- plot_edges %>%
    dplyr::select(
      node_id = target_node,
      label = target,
      log2FoldChange = target_log2FoldChange,
      direction = target_direction_clean
    ) %>%
    dplyr::mutate(
      node_type = "Target"
    ) %>%
    dplyr::distinct()
  
  nodes <- dplyr::bind_rows(tf_nodes, target_nodes) %>%
    dplyr::mutate(
      node_colour_group = dplyr::case_when(
        node_type == "TF" & direction == "Up" ~ "TF Up",
        node_type == "TF" & direction == "Down" ~ "TF Down",
        node_type == "Target" & direction == "Up" ~ "Target Up",
        node_type == "Target" & direction == "Down" ~ "Target Down",
        TRUE ~ "Target not significant"
      )
    )
  
  edges <- plot_edges %>%
    dplyr::select(
      from = tf_node,
      to = target_node,
      dplyr::any_of(c(
        "confidence",
        "mor",
        "likelihood",
        "x_relationship"
      ))
    )
  
  g <- igraph::graph_from_data_frame(
    d = edges,
    vertices = nodes,
    directed = TRUE
  )
  
  tg <- tidygraph::as_tbl_graph(g) %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(
      degree = tidygraph::centrality_degree()
    )
  
  set.seed(123)
  
  p <- ggraph::ggraph(tg, layout = "fr") +
    ggraph::geom_edge_link(
      alpha = 0.25,
      colour = "grey60",
      arrow = grid::arrow(length = grid::unit(2, "mm")),
      end_cap = ggraph::circle(2, "mm")
    ) +
    ggraph::geom_node_point(
      aes(
        size = degree,
        colour = node_colour_group
      ),
      alpha = 0.9
    ) +
    ggraph::geom_node_text(
      aes(
        label = ifelse(node_type == "TF", label, "")
      ),
      repel = TRUE,
      size = 3
    ) +
    scale_colour_manual(
      values = c(
        "TF Up" = "#D55E00",
        "TF Down" = "#0072B2",
        "Target Up" = "#E69F00",
        "Target Down" = "#56B4E9",
        "Target not significant" = "grey70"
      )
    ) +
    scale_size(range = c(2, 8)) +
    theme_void() +
    labs(
      title = paste0("X-associated TF-target network: ", contrast_to_plot),
      subtitle = "Nodes coloured by TF/target direction of change; doRothEA A/B/C regulons",
      colour = "Node direction",
      size = "Degree"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )
  
  print(p)
  
  ggsave(
    filename = file.path(lambert_plot_dir, output_name),
    plot = p,
    width = 10,
    height = 9,
    dpi = 300,
    bg = "white"
  )
  
  return(p)
}


## 18. Plot X-associated significant target networks


network_HDAC1_X <- plot_tf_target_network(
  edge_df = tf_target_edges_Xassociated_sig_targets,
  contrast_to_plot = "HDAC1 KD",
  output_name = "doRothEA_Xassociated_network_HDAC1_KD.png"
)

network_HDAC2_X <- plot_tf_target_network(
  edge_df = tf_target_edges_Xassociated_sig_targets,
  contrast_to_plot = "HDAC2 KD",
  output_name = "doRothEA_Xassociated_network_HDAC2_KD.png"
)

network_HDAC3_X <- plot_tf_target_network(
  edge_df = tf_target_edges_Xassociated_sig_targets,
  contrast_to_plot = "HDAC3 KD",
  output_name = "doRothEA_Xassociated_network_HDAC3_KD.png"
)


tf_deg_summary
xlinked_tf_deg_summary
xlinked_tf_deg_list
xlinked_tf_plot_df
tf_target_edges_Xassociated_sig_targets
xassoc_sig_edge_detail

tf_target_edges_Xassociated_sig_targets %>%
  dplyr::count(contrast_short, x_relationship)



## Heatmap of TFs and predicted Subject-to-XCI target genes
## Lambert TFs + doRothEA TF-target relationships


library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(stringr)

## Output directory
tf_heatmap_dir <- file.path(lambert_plot_dir, "TF_target_heatmaps")
dir.create(tf_heatmap_dir, showWarnings = FALSE, recursive = TRUE)


## TF-target edges where target is X-linked and Subject to XCI


tukiainen_clean_path <- "/home/jvk3/Desktop/HDAC_counts/XCI_classifications/XCI_classification_outputs/Tukiainen_XCI_categories_clean.csv"

tukiainen_target_annot <- readr::read_csv(
  tukiainen_clean_path,
  show_col_types = FALSE
) %>%
  dplyr::select(
    target = SYMBOL,
    target_tukiainen_status_grouped = tukiainen_status_grouped
  ) %>%
  dplyr::distinct()

tf_target_subjectX_edges <- tf_target_edges_Xassociated %>%
  dplyr::left_join(
    tukiainen_target_annot,
    by = "target"
  ) %>%
  dplyr::filter(
    target_CHR == "X",
    target_tukiainen_status_grouped == "Subject to XCI"
  ) %>%
  dplyr::mutate(
    expected_relationship = dplyr::case_when(
      mor == 1 ~ "Activating",
      mor == -1 ~ "Inhibitory",
      TRUE ~ "Unknown"
    )
  )

tf_target_subjectX_edges %>%
  dplyr::count(contrast_short, expected_relationship)

tf_target_subjectX_edges %>%
  dplyr::select(
    contrast_short,
    tf,
    target,
    tf_CHR,
    target_CHR,
    target_tukiainen_status_grouped,
    mor,
    expected_relationship,
    confidence,
    tf_log2FoldChange,
    target_log2FoldChange,
    tf_padj,
    target_padj
  ) %>%
  dplyr::distinct()



## Heatmap of TFs and predicted Subject-to-XCI targets
## Lambert TF catalogue + doRothEA TF-target relationships


library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

tf_heatmap_dir <- file.path(lambert_plot_dir, "TF_target_heatmaps")
dir.create(tf_heatmap_dir, showWarnings = FALSE, recursive = TRUE)


## 1. Keep candidate TF-target links


tf_target_subjectX_edges_for_plot <- tf_target_subjectX_edges %>%
  dplyr::mutate(
    contrast_short = factor(
      contrast_short,
      levels = c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
    ),
    expected_relationship = factor(
      expected_relationship,
      levels = c("Activating", "Inhibitory", "Unknown")
    ),
    
    ## Simpler pair label
    pair_label = paste0(tf, " → ", target)
  )

tf_target_subjectX_edges_for_plot %>%
  dplyr::count(contrast_short, expected_relationship, confidence)


## 2. Create rows for TFs and their target genes


tf_rows <- tf_target_subjectX_edges_for_plot %>%
  dplyr::transmute(
    contrast_short,
    pair_label,
    expected_relationship,
    confidence,
    SYMBOL = tf,
    role = "TF",
    mor,
    tf,
    target
  )

target_rows <- tf_target_subjectX_edges_for_plot %>%
  dplyr::transmute(
    contrast_short,
    pair_label,
    expected_relationship,
    confidence,
    SYMBOL = target,
    role = "Target",
    mor,
    tf,
    target
  )

tf_target_rows_for_heatmap <- dplyr::bind_rows(
  tf_rows,
  target_rows
) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    role = factor(
      role,
      levels = c("TF", "Target")
    ),
    
    ## Cleaner heatmap row label
    row_label = paste0(pair_label, " | ", role, ": ", SYMBOL)
  )

genes_for_tf_target_heatmap <- unique(tf_target_rows_for_heatmap$SYMBOL)

## 3. Extract DESeq2-normalised counts
## Have to run at least up until dds_pca in HDAC123_analyses_updated.R

dds_pca <- DESeq2::estimateSizeFactors(dds_pca)

norm_counts_combined <- DESeq2::counts(
  dds_pca,
  normalized = TRUE
)

gene_annot_for_tf_target_heatmap <- all_res_hdac1_2_3_15fold %>%
  dplyr::select(gene_id, SYMBOL) %>%
  dplyr::distinct() %>%
  dplyr::filter(SYMBOL %in% genes_for_tf_target_heatmap)

sample_meta_for_tf_target_heatmap <- meta_pca %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::select(sample, pca_group)

tf_target_expr_long <- norm_counts_combined[
  gene_annot_for_tf_target_heatmap$gene_id,
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
    gene_annot_for_tf_target_heatmap,
    by = "gene_id"
  ) %>%
  dplyr::left_join(
    sample_meta_for_tf_target_heatmap,
    by = "sample"
  )


## 4. Define matched Scramble and KD groups for each contrast


contrast_sample_map <- tibble::tribble(
  ~contrast_short, ~scramble_group,      ~kd_group,
  "HDAC1 KD",      "HDAC12 Scramble",    "HDAC1 KD",
  "HDAC2 KD",      "HDAC12 Scramble",    "HDAC2 KD",
  "HDAC3 KD",      "HDAC3 Scramble",     "HDAC3 KD"
)

tf_target_expr_mean <- tf_target_expr_long %>%
  dplyr::group_by(SYMBOL, pca_group) %>%
  dplyr::summarise(
    mean_norm_count = mean(norm_count, na.rm = TRUE),
    log2_mean_count = log2(mean_norm_count + 1),
    .groups = "drop"
  )

tf_target_heatmap_df <- tf_target_rows_for_heatmap %>%
  dplyr::left_join(
    contrast_sample_map,
    by = "contrast_short"
  ) %>%
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
    )
  ) %>%
  dplyr::left_join(
    tf_target_expr_mean,
    by = c("SYMBOL", "pca_group")
  ) %>%
  dplyr::group_by(row_label) %>%
  dplyr::mutate(
    row_scaled_expr = as.numeric(scale(log2_mean_count))
  ) %>%
  dplyr::ungroup()

tf_target_heatmap_df %>%
  dplyr::select(
    contrast_short,
    pair_label,
    role,
    SYMBOL,
    condition_plot,
    mean_norm_count,
    log2_mean_count,
    row_scaled_expr
  )


## 5. Plot TF-target count heatmap


tf_target_expression_heatmap <- ggplot(
  tf_target_heatmap_df,
  aes(
    x = condition_plot,
    y = row_label,
    fill = row_scaled_expr
  )
) +
  geom_tile(
    colour = "grey85",
    linewidth = 0.25
  ) +
  facet_grid(
    contrast_short + expected_relationship ~ .,
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
    title = "Expression of candidate TFs and predicted Subject-to-XCI target genes",
    subtitle = "TFs identified using the Lambert catalogue; predicted TF-target relationships from doRothEA",
    x = NULL,
    y = NULL,
    fill = "Row-scaled\nlog2 count"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.5, "lines")
  )

tf_target_expression_heatmap

ggsave(
  filename = file.path(
    tf_heatmap_dir,
    "TF_predicted_SubjectToXCI_target_expression_heatmap_simplified_labels.png"
  ),
  plot = tf_target_expression_heatmap,
  width = 14,
  height = max(
    6,
    0.25 * dplyr::n_distinct(tf_target_heatmap_df$row_label) + 2
  ),
  dpi = 300,
  bg = "white"
)
readr::write_csv(
  tf_target_heatmap_df,
  file.path(
    tf_heatmap_dir,
    "TF_predicted_SubjectToXCI_target_expression_heatmap_table.csv"
  )
)