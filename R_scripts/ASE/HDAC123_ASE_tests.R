## Combine HDAC1/2 and HDAC3 ASE Xi-fraction outputs


library(tidyverse)
library(data.table)

hdac12_out <- "/home/jvk3/Desktop/HDAC_counts/All three datsets bulk and ASE/HDAC1 and 2/ASE/ase_outputs_Xi_fraction_HDAC12"
hdac3_out  <- "/home/jvk3/Desktop/HDAC_counts/All three datsets bulk and ASE/HDAC3/ASE/ase_outputs_Xi_fraction_HDAC3"

combined_ase_out <- "/home/jvk3/Desktop/HDAC_counts/All three datsets bulk and ASE/Combined_HDAC_ASE_Xi_fraction_outputs"

dir.create(combined_ase_out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(combined_ase_out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(combined_ase_out, "plots"), recursive = TRUE, showWarnings = FALSE)


## 1. Combine chromosome-wide %Xi per replicate


rep_xi_hdac12 <- readr::read_csv(
  file.path(hdac12_out, "tables", "pct_Xi_per_replicate.csv"),
  show_col_types = FALSE
) %>%
  mutate(dataset = "HDAC1/2")

rep_xi_hdac3 <- readr::read_csv(
  file.path(hdac3_out, "tables", "pct_Xi_per_replicate.csv"),
  show_col_types = FALSE
) %>%
  mutate(dataset = "HDAC3")

rep_xi_combined <- bind_rows(rep_xi_hdac12, rep_xi_hdac3) %>%
  mutate(
    contrast_short = case_when(
      condition == "siHDAC1" ~ "HDAC1 KD",
      condition == "siHDAC2" ~ "HDAC2 KD",
      condition == "HDAC3_siRNA" ~ "HDAC3 KD",
      condition == "Scramble" & dataset == "HDAC1/2" ~ "Scramble HDAC1/2",
      condition == "Scramble" & dataset == "HDAC3" ~ "Scramble HDAC3",
      TRUE ~ as.character(condition)
    ),
    contrast_short = factor(
      contrast_short,
      levels = c(
        "Scramble HDAC1/2",
        "HDAC1 KD",
        "HDAC2 KD",
        "Scramble HDAC3",
        "HDAC3 KD"
      )
    )
  )

readr::write_csv(
  rep_xi_combined,
  file.path(combined_ase_out, "tables", "combined_pct_Xi_per_replicate.csv")
)

rep_xi_combined


## 2. Plot chromosome-wide Xi-derived transcription


p_rep_xi_combined <- ggplot(
  rep_xi_combined,
  aes(x = contrast_short, y = pct_Xi)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(shape = dataset), width = 0.12, size = 2.8) +
  theme_bw() +
  labs(
    x = NULL,
    y = "% reads from Xi allele",
    shape = "Dataset",
    title = "Chromosome-wide Xi-derived transcription after HDAC knockdown",
    subtitle = "Xi contribution calculated from allele-specific reads assigned as Xi/(Xi + Xa)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

p_rep_xi_combined

ggsave(
  file.path(combined_ase_out, "plots", "Combined_HDAC_pct_Xi_per_replicate.png"),
  p_rep_xi_combined,
  width = 7,
  height = 5,
  dpi = 300
)


## 3. Combine gene-level Xi shift tables


gene_xi_hdac12 <- readr::read_csv(
  file.path(hdac12_out, "tables", "xi_shift_gene_vs_control.csv"),
  show_col_types = FALSE
) %>%
  mutate(dataset = "HDAC1/2")

gene_xi_hdac3 <- readr::read_csv(
  file.path(hdac3_out, "tables", "xi_shift_gene_vs_control.csv"),
  show_col_types = FALSE
) %>%
  mutate(dataset = "HDAC3")

gene_xi_combined <- bind_rows(gene_xi_hdac12, gene_xi_hdac3) %>%
  mutate(
    contrast_short = case_when(
      comparison == "siHDAC1" ~ "HDAC1 KD",
      comparison == "siHDAC2" ~ "HDAC2 KD",
      comparison == "HDAC3_siRNA" ~ "HDAC3 KD",
      TRUE ~ as.character(comparison)
    ),
    contrast_short = factor(
      contrast_short,
      levels = c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
    )
  )

readr::write_csv(
  gene_xi_combined,
  file.path(combined_ase_out, "tables", "combined_xi_shift_gene_vs_control.csv")
)

gene_xi_combined %>%
  dplyr::count(contrast_short)


## 4. Combine SNP-level Xi shift tables


snp_xi_hdac12 <- readr::read_csv(
  file.path(hdac12_out, "tables", "xi_shift_snp_vs_control.csv"),
  show_col_types = FALSE
) %>%
  mutate(dataset = "HDAC1/2")

snp_xi_hdac3 <- readr::read_csv(
  file.path(hdac3_out, "tables", "xi_shift_snp_vs_control.csv"),
  show_col_types = FALSE
) %>%
  mutate(dataset = "HDAC3")

snp_xi_combined <- bind_rows(snp_xi_hdac12, snp_xi_hdac3) %>%
  mutate(
    contrast_short = case_when(
      comparison == "siHDAC1" ~ "HDAC1 KD",
      comparison == "siHDAC2" ~ "HDAC2 KD",
      comparison == "HDAC3_siRNA" ~ "HDAC3 KD",
      TRUE ~ as.character(comparison)
    ),
    contrast_short = factor(
      contrast_short,
      levels = c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
    ),
    position_mb = position / 1e6,
    sig = ifelse(!is.na(fdr) & fdr < 0.05, "FDR < 0.05", "Not significant")
  )

readr::write_csv(
  snp_xi_combined,
  file.path(combined_ase_out, "tables", "combined_xi_shift_snp_vs_control.csv")
)

snp_xi_combined %>%
  dplyr::count(contrast_short, sig)

snp_xi_combined %>%
  mutate(
    xi_shift_direction = case_when(
      is.na(delta_xi) ~ "Missing",
      delta_xi > 0 ~ "Increased Xi fraction",
      delta_xi < 0 ~ "Decreased Xi fraction",
      TRUE ~ "No change"
    )
  ) %>%
  count(contrast_short, sig, xi_shift_direction)


## ASE X chromosome positional plot: delta Xi fraction


XIST_START <- 73820651
XIST_END   <- 73852753

xist_band_df <- tibble::tibble(
  xmin = XIST_START / 1e6,
  xmax = XIST_END / 1e6,
  ymin = -Inf,
  ymax = Inf
)

snp_xi_plot_df <- snp_xi_combined %>%
  filter(contig %in% c("X", "chrX")) %>%
  mutate(
    position_mb = position / 1e6,
    xi_shift_group = case_when(
      is.na(fdr) | fdr >= 0.05 ~ "Not significant",
      delta_xi > 0 ~ "Increased Xi fraction",
      delta_xi < 0 ~ "Decreased Xi fraction",
      TRUE ~ "No change"
    ),
    xi_shift_group = factor(
      xi_shift_group,
      levels = c(
        "Increased Xi fraction",
        "Decreased Xi fraction",
        "Not significant",
        "No change"
      )
    )
  )

p_ase_xchr_delta_xi <- ggplot(
  snp_xi_plot_df,
  aes(x = position_mb, y = delta_xi, colour = xi_shift_group)
) +
  geom_rect(
    data = xist_band_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "purple",
    alpha = 0.12
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 2.2, alpha = 0.85) +
  facet_wrap(~ contrast_short, ncol = 1) +
  scale_colour_manual(
    values = c(
      "Increased Xi fraction" = "red",
      "Decreased Xi fraction" = "blue",
      "Not significant" = "grey70",
      "No change" = "grey40"
    )
  ) +
  theme_bw() +
  labs(
    x = "Chromosome X position (Mb)",
    y = expression(Delta~"Xi fraction"),
    colour = "ASE shift",
    title = "Allele-specific Xi fraction shifts across chromosome X after HDAC knockdown",
    subtitle = "Positive values indicate increased Xi-derived transcription relative to Scramble"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_ase_xchr_delta_xi

ggsave(
  file.path(
    combined_ase_out,
    "plots",
    "Combined_HDAC_ASE_delta_Xi_fraction_chrX.png"
  ),
  p_ase_xchr_delta_xi,
  width = 9,
  height = 8,
  dpi = 300
)


## Bulk-style ASE chromosome X plot using gene-level delta Xi


library(tidyverse)
library(ggrepel)

bulk_raw <- readr::read_csv(hdac_tukiainen_path, show_col_types = FALSE)

# Check column names first
colnames(bulk_raw)

bulk_annot <- bulk_raw[bulk_raw$CHR == "X", ] %>%
  dplyr::select(
    SYMBOL,
    contrast_short,
    CHR,
    START,
    END,
    tukiainen_status_grouped,
    bulk_deg_status = deg_status_15fold,
    bulk_log2FC = log2FoldChange,
    bulk_padj = padj
  ) %>%
  dplyr::distinct(SYMBOL, contrast_short, .keep_all = TRUE)

ase_gene_bulkstyle_df <- gene_xi_combined %>%
  dplyr::filter(!is.na(gene), gene != "intergenic") %>%
  tidyr::separate_rows(gene, sep = ";") %>%
  dplyr::rename(SYMBOL = gene) %>%
  dplyr::left_join(
    bulk_annot,
    by = c("SYMBOL", "contrast_short")
  ) %>%
  dplyr::filter(CHR == "X", !is.na(START), !is.na(delta_xi)) %>%
  dplyr::mutate(
    position_mb = START / 1e6,
    ase_direction = dplyr::case_when(
      fdr < 0.05 & delta_xi > 0 ~ "Increased Xi fraction",
      fdr < 0.05 & delta_xi < 0 ~ "Decreased Xi fraction",
      TRUE ~ "Not significant"
    ),
    direction_y = dplyr::case_when(
      delta_xi > 0 ~ 1,
      delta_xi < 0 ~ -1,
      TRUE ~ 0
    ),
    abs_delta_xi = abs(delta_xi),
    contrast_short = factor(
      contrast_short,
      levels = c("HDAC1 KD", "HDAC2 KD", "HDAC3 KD")
    ),
    bulk_deg_status = dplyr::case_when(
      is.na(bulk_deg_status) ~ "Not in bulk table",
      TRUE ~ as.character(bulk_deg_status)
    ),
    bulk_deg_status = factor(
      bulk_deg_status,
      levels = c("Up", "Down", "Not significant", "Not in bulk table")
    )
  )


## Plot ASE direction across chromosome X


XIST_START <- 73820651
XIST_END   <- 73852753

xist_band_df <- tibble::tibble(
  xmin = XIST_START / 1e6,
  xmax = XIST_END / 1e6,
  ymin = -1.35,
  ymax = 1.35
)

label_df_ase <- ase_gene_bulkstyle_df %>%
  dplyr::filter(fdr < 0.05) %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::arrange(dplyr::desc(abs_delta_xi), .by_group = TRUE) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()

p_ase_bulkstyle <- ggplot(
  ase_gene_bulkstyle_df,
  aes(
    x = position_mb,
    y = direction_y,
    colour = ase_direction,
    size = abs_delta_xi
  )
) +
  geom_rect(
    data = xist_band_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "purple",
    alpha = 0.12
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(
    alpha = 0.85,
    position = position_jitter(width = 0, height = 0.05)
  ) +
  geom_text_repel(
    data = label_df_ase,
    aes(label = SYMBOL),
    size = 3,
    max.overlaps = 50,
    show.legend = FALSE
  ) +
  facet_wrap(~ contrast_short, ncol = 1) +
  scale_y_continuous(
    breaks = c(-1, 0, 1),
    labels = c("Decreased Xi", "0", "Increased Xi"),
    limits = c(-1.35, 1.35)
  ) +
  scale_colour_manual(
    values = c(
      "Increased Xi fraction" = "red",
      "Decreased Xi fraction" = "blue",
      "Not significant" = "grey70"
    )
  ) +
  scale_size_continuous(range = c(2, 6)) +
  theme_bw() +
  labs(
    x = "Position on chromosome X, GRCh38",
    y = expression("Direction of " * Delta * " Xi/(Xa+Xi)"),
    colour = "ASE shift",
    size = "|Δ Xi fraction|",
    title = "Allele-specific Xi-fraction direction across chromosome X after HDAC knockdown",
    subtitle = "Y-axis shows direction of Xi-fraction change; point shape shows matched bulk expression status"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

p_ase_bulkstyle

ggsave(
  file.path(
    combined_ase_out,
    "plots",
    "Combined_HDAC_ASE_bulk_style_Xi_direction_chrX.png"
  ),
  p_ase_bulkstyle,
  width = 10,
  height = 8,
  dpi = 300
)

ase_gene_bulkstyle_df %>%
  dplyr::count(contrast_short, ase_direction, bulk_deg_status)


## Directionality test: are significant ASE shifts biased towards increased Xi fraction?


ase_direction_stats <- ase_gene_bulkstyle_df %>%
  dplyr::filter(
    ase_direction %in% c("Increased Xi fraction", "Decreased Xi fraction")
  ) %>%
  dplyr::group_by(contrast_short) %>%
  dplyr::summarise(
    n_increased = sum(ase_direction == "Increased Xi fraction"),
    n_decreased = sum(ase_direction == "Decreased Xi fraction"),
    n_total = n_increased + n_decreased,
    prop_increased = n_increased / n_total,
    binom_p_one_sided = binom.test(
      n_increased,
      n_total,
      p = 0.5,
      alternative = "greater"
    )$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj_BH = p.adjust(binom_p_one_sided, method = "BH")
  )

ase_direction_stats

readr::write_csv(
  ase_direction_stats,
  file.path(
    combined_ase_out,
    "tables",
    "ASE_directionality_binomial_test_by_HDAC.csv"
  )
)


## Combined directionality test across all HDAC knockdowns

ase_direction_stats_combined <- ase_gene_bulkstyle_df %>%
  dplyr::filter(
    ase_direction %in% c("Increased Xi fraction", "Decreased Xi fraction")
  ) %>%
  dplyr::summarise(
    n_increased = sum(ase_direction == "Increased Xi fraction"),
    n_decreased = sum(ase_direction == "Decreased Xi fraction"),
    n_total = n_increased + n_decreased,
    prop_increased = n_increased / n_total,
    binom_p_one_sided = binom.test(
      n_increased,
      n_total,
      p = 0.5,
      alternative = "greater"
    )$p.value
  )

ase_direction_stats_combined



