# XCI analysis for two or more ordered conditions.
# conditions controls plotting order, while control_condition is used to call
# Xa/Xi alleles and as the reference for all statistical comparisons. Each
# non-control condition is compared with the control rather than performing
# every possible pairwise comparison.

library(data.table)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(readxl)
library(GenomicRanges)
library(biomaRt)
library(VGAM)

## Configuration

base_dir <- "/home/jvk3/Desktop/HDAC_counts/All three datsets bulk and ASE/Capacitation/ASE"

ase_dir <- file.path(base_dir, "8_ASE_counts")

meta_path <- file.path(base_dir, "meta_table_cap_ase.csv")

out_dir <- file.path(base_dir, "ase_outputs_Xi_fraction_capacitation")

gene_annot_path <- "/home/jvk3/Desktop/HDAC_counts/XCI_classifications/XCI_classification_outputs/HDAC_results_Tukiainen_annotated.csv"

min_depth <- 10
alpha <- 0.05

goi <- c("XIST", "XACT", "ATRX", "THOC2", "HUWE1", "G6PD")

## Tukiainen cleaned annotation file
xci_ref_path <- "/home/jvk3/Desktop/HDAC_counts/XCI_classifications/XCI_classification_outputs/Tukiainen_XCI_categories_clean.csv"

## Capacitation condition order
conditions <- c("Naive", "Day3", "Day7", "Primed")

## Use Primed to call Xa/Xi, because it should have clearer established XCI
xa_call_condition <- "Primed"

## Compare all timepoints back to Naive
control_condition <- "Naive"

stopifnot("XA_CALL_CONDITION must be one of CONDITIONS" = xa_call_condition %in% conditions)
stopifnot("CONTROL_CONDITION must be one of CONDITIONS" = control_condition %in% conditions)

comparison_conditions <- setdiff(conditions, control_condition)

## Use the Tukiainen category matching your cleaned file
candidate_xci_statuses <- c("Subject to XCI")

xa_tiers <- data.table(
  tier      = c("high", "medium", "low"),
  min_skew  = c(0.90,   0.80,     0.75),
  min_depth = c(15,     20,       30)
)

tiers_for_main_analysis <- c("high", "medium")

## X-inactivation center (XIC) region, GRCh38. A tight, literature-defined
## XIC (from a mapped deletion study, converted from hg19) sits at roughly
## chrX:72.86-74.15 Mb, and the core XIST/TSIX/JPX gene cluster sits at
## chrX:73.79-74.07 Mb. The padded window below comfortably contains both,
## with a few Mb of margin - a reasonable choice for a shaded landmark
## region, though NOT the tight/minimal XIC boundary. Narrow it to
## ~72.9-74.15 (or 73.79-74.07 for just the gene cluster) if you want the
## shading to reflect the literature boundary precisely rather than a
## padded window.
## X-inactivation center (XIC) region, GRCh38, TIGHT boundary.
## Source: a mapped XIC deletion study defined the region at hg19
## chrX:72,080,568-73,367,054 (https://doi.org/10.1098/rstb.2016.0359).
## Converted to GRCh38 using the ~780kb offset implied by XIST's own
## hg19->hg38 shift (hg19 73,040,486-73,072,588 -> hg38 73,817,774-
## 73,852,754). This offset-based conversion is approximate, not a formal
## liftOver - if you need publication-grade precision, run the hg19
## coordinates through UCSC liftOver (https://genome.ucsc.edu/cgi-bin/hgLiftOver)
## and replace the values below with the exact output.
##
## For an even tighter option - just the core XIST/TSIX/JPX gene cluster,
## no flanking regulatory region - use XIC_start_mb <- 73.79, XIC_end_mb <- 74.07.
xic_start_mb <- 72.86
xic_end_mb   <- 74.15
xic_start <- xic_start_mb * 1e6
xic_end   <- xic_end_mb   * 1e6

dir.create(file.path(out_dir, "plots"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

##  1. Load ASE counts + metadata, QC filter 

meta <- fread(meta_path)
setnames(meta, make.names(names(meta)))

## Check required columns
stopifnot("Metadata needs a 'rep' column" = "rep" %in% names(meta))
stopifnot("Metadata needs a 'condition' column" = "condition" %in% names(meta))
stopifnot("Metadata needs a 'filename' column" = "filename" %in% names(meta))
stopifnot("Metadata needs a 'sample' column" = "sample" %in% names(meta))

## Standardise condition names if needed
meta[, condition := dplyr::case_when(
  condition %in% c("Naive", "naive") ~ "Naive",
  condition %in% c("Day3", "Cap. Day 3", "Cap_Day_3", "D3") ~ "Day3",
  condition %in% c("Day7", "Cap. Day 7", "Cap_Day_7", "D7") ~ "Day7",
  condition %in% c("Primed", "primed") ~ "Primed",
  TRUE ~ as.character(condition)
)]

stopifnot("Metadata condition values must all be in CONDITIONS" =
            all(unique(meta$condition) %in% conditions))

meta[, filepath := file.path(ase_dir, filename)]

## Check all files exist
missing_files <- meta[!file.exists(filepath)]
if (nrow(missing_files) > 0) {
  print(missing_files[, .(sample, condition, rep, filename, filepath)])
  stop("Some ASE count files listed in the metadata do not exist.")
}

all_ase <- rbindlist(
  lapply(seq_len(nrow(meta)), function(i) {
    dt <- fread(meta$filepath[i])
    if (!"sample" %in% names(dt)) dt[, sample := meta$sample[i]]
    dt
  }),
  fill = TRUE
)
all_ase <- merge(all_ase, meta, by = "sample", all.x = TRUE)
all_ase[, contig := gsub("^chr", "", contig)]
all_ase[, condition := factor(condition, levels = conditions)]

max_other_frac   <- 0.05
max_lowmapq_frac <- 0.20
all_ase[, other_frac   := otherBases   / rawDepth]
all_ase[, lowmapq_frac := lowMAPQDepth / rawDepth]
n_before <- nrow(all_ase)
all_ase <- all_ase[(is.na(other_frac)   | other_frac   <= max_other_frac) &
                     (is.na(lowmapq_frac) | lowmapq_frac <= max_lowmapq_frac)]
cat("Site-quality filter removed", n_before - nrow(all_ase), "of", n_before, "rows\n")

rep_ratio <- all_ase[totalCount >= min_depth,
                     .(ref_ratio = refCount / totalCount),
                     by = .(contig, position, sample, condition)]
rep_wide <- dcast(rep_ratio, contig + position + condition ~ sample, value.var = "ref_ratio")
sample_cols <- setdiff(names(rep_wide), c("contig", "position", "condition"))
cat("\nReplicate ref_ratio correlation (pairwise complete obs), all conditions pooled:\n")
print(round(cor(rep_wide[, ..sample_cols], use = "pairwise.complete.obs"), 3))

##  2. Pool replicates per condition, per SNP 
combined <- all_ase[, .(
  refCount   = sum(refCount, na.rm = TRUE),
  altCount   = sum(altCount, na.rm = TRUE),
  totalCount = sum(totalCount, na.rm = TRUE),
  variantID  = variantID[1], refAllele = refAllele[1], altAllele = altAllele[1]
), by = .(contig, position, condition)]
combined[, ref_ratio := refCount / totalCount]



## 3. Gene annotation from local HDAC annotated results
##    This avoids biomaRt/Ensembl HTTP errors

genes_cache <- file.path(out_dir, "tables", "local_gene_coordinates_cache.rds")

if (file.exists(genes_cache)) {

  genes <- readRDS(genes_cache)

} else {

  if (!file.exists(gene_annot_path)) {
    stop("GENE_ANNOT_PATH does not exist. Check the path to HDAC_results_Tukiainen_annotated.csv")
  }

  genes <- readr::read_csv(
    gene_annot_path,
    show_col_types = FALSE
  ) %>%
    dplyr::select(
      ensembl_gene_id = gene_id,
      hgnc_symbol = SYMBOL,
      chr = CHR,
      start = START,
      end = END
    ) %>%
    dplyr::filter(
      !is.na(hgnc_symbol),
      hgnc_symbol != "",
      !is.na(chr),
      !is.na(start),
      !is.na(end)
    ) %>%
    dplyr::mutate(
      chr = gsub("^chr", "", as.character(chr)),
      start = as.integer(start),
      end = as.integer(end),
      gene_biotype = NA_character_,
      strand = "*"
    ) %>%
    dplyr::distinct(
      ensembl_gene_id,
      hgnc_symbol,
      chr,
      start,
      end,
      .keep_all = TRUE
    ) %>%
    dplyr::filter(
      chr %in% unique(combined$contig)
    )

  saveRDS(genes, genes_cache)
}

genes <- as.data.frame(genes)

genes_gr <- GenomicRanges::makeGRangesFromDataFrame(
  genes,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE
)

snp_gr <- GRanges(seqnames = combined$contig,
                  ranges = IRanges(start = combined$position, end = combined$position),
                  mcols = combined)
colnames(mcols(snp_gr)) <- gsub("^mcols\\.", "", colnames(mcols(snp_gr)))
snp_gr$snp_id <- seq_along(snp_gr)

hits <- findOverlaps(snp_gr, genes_gr)
pairs <- data.table(snp_id = snp_gr$snp_id[queryHits(hits)],
                    gene   = genes_gr$hgnc_symbol[subjectHits(hits)])
gene_annot <- pairs[, .(gene = paste(unique(gene), collapse = ";"), n_genes_overlapping = .N), by = snp_id]

snp_dt <- as.data.table(mcols(snp_gr))
snp_annotations <- merge(snp_dt, gene_annot, by = "snp_id", all.x = TRUE)
snp_annotations[is.na(gene), `:=`(gene = "intergenic", n_genes_overlapping = 0L)]

if (file.exists(xci_ref_path)) {

  xci_lookup <- fread(xci_ref_path)
  setnames(xci_lookup, make.names(names(xci_lookup)))

  ## Adjust depending on the exact column names in your cleaned file
  if ("SYMBOL" %in% names(xci_lookup)) {
    setnames(xci_lookup, "SYMBOL", "gene")
  }

  if ("tukiainen_status_grouped" %in% names(xci_lookup)) {
    setnames(xci_lookup, "tukiainen_status_grouped", "XCI_status")
  }

  xci_lookup <- unique(
    xci_lookup[
      !is.na(gene) & gene != "",
      .(gene, XCI_status)
    ],
    by = "gene"
  )

  snp_annotations <- merge(
    snp_annotations,
    xci_lookup,
    by = "gene",
    all.x = TRUE
  )

  snp_annotations[is.na(XCI_status), XCI_status := "Not in reference"]

} else {
  warning("XCI reference file not found - XCI_status will be NA throughout.")
  snp_annotations[, XCI_status := NA_character_]
}

snp_annotations[, condition := factor(condition, levels = conditions)]

## 4. Wide table: one row per SNP, one column set per condition 
# dcast automatically names columns <valuevar>_<condition label> for
# however many condition levels are present - this is what lets the rest
# of the script generalize to N conditions without hardcoding names.
wide <- dcast(snp_annotations,
              contig + position + gene + XCI_status + variantID + refAllele + altAllele ~ condition,
              value.var = c("refCount", "altCount", "totalCount"))
cc_snp <- as.data.table(wide)

for (cond in conditions) {
  rc <- paste0("refCount_", cond); ac <- paste0("altCount_", cond)
  if (!rc %in% names(cc_snp)) next
  cc_snp[, (paste0("ratio_", cond)) := get(rc) / (get(rc) + get(ac))]
  cc_snp[, (paste0("depth_", cond)) := get(rc) + get(ac)]
}

test_condition_effect <- function(ref_A, alt_A, ref_B, alt_B) {
  dat <- data.frame(condition = factor(c("A", "B"), levels = c("A", "B")),
                    ref = c(ref_A, ref_B), n = c(ref_A + alt_A, ref_B + alt_B))
  fit <- tryCatch(vglm(cbind(ref, n - ref) ~ condition, betabinomial, data = dat, trace = FALSE),
                  error = function(e) NULL)
  if (is.null(fit)) return(c(p_value = NA_real_, logit_shift = NA_real_))
  co <- summary(fit)@coef3
  c(p_value = co[2, "Pr(>|z|)"], logit_shift = co[2, "Estimate"])
}

##  4b. Call Xa/Xi per SNP from the xa_call_condition baseline 
ctrl_ratio_col <- paste0("ratio_", xa_call_condition)
ctrl_depth_col <- paste0("depth_", xa_call_condition)
candidates <- cc_snp[XCI_status %in% candidate_xci_statuses]

baseline_calls <- rbindlist(lapply(seq_len(nrow(xa_tiers)), function(i) {
  tier <- xa_tiers$tier[i]; min_skew <- xa_tiers$min_skew[i]; min_depth <- xa_tiers$min_depth[i]
  ratio_vals <- candidates[[ctrl_ratio_col]]
  depth_vals <- candidates[[ctrl_depth_col]]
  hit <- candidates[depth_vals >= min_depth & (ratio_vals >= min_skew | ratio_vals <= (1 - min_skew))]
  if (nrow(hit) == 0) return(NULL)
  hit_ratio <- hit[[ctrl_ratio_col]]
  data.table(contig = hit$contig, position = hit$position, gene = hit$gene,
             ratio_control = hit_ratio, tier = tier,
             xa_allele = fifelse(hit_ratio >= min_skew, "ref", "alt"))
}))
baseline_calls[, tier_rank := factor(tier, levels = xa_tiers$tier)]
setorder(baseline_calls, contig, position, tier_rank)
baseline_calls[, tier_rank := NULL]
baseline_calls <- unique(baseline_calls, by = c("contig", "position"))

cat(sprintf("\nXa/Xi candidate calls by tier (baseline = %s):\n", xa_call_condition))
print(table(baseline_calls$tier))

# Replicate-consistency check within the phase-calling condition only
rep_check <- rep_wide[condition == xa_call_condition]
rep_check <- merge(rep_check, baseline_calls[, .(contig, position, xa_allele)], by = c("contig", "position"))
id_vars <- c("contig", "position", "condition", "xa_allele")
rep_long <- melt(rep_check, id.vars = id_vars,
                 measure.vars = setdiff(names(rep_check), id_vars),
                 variable.name = "sample", value.name = "ref_ratio")
rep_long <- rep_long[!is.na(ref_ratio)]
rep_long[, agrees := fifelse(xa_allele == "ref", ref_ratio > 0.5, ref_ratio < 0.5)]
snp_consistency <- rep_long[, .(n_checked = .N, n_agree = sum(agrees)), by = .(contig, position)]
snp_consistency[, fully_consistent := n_agree == n_checked & n_checked > 0]

baseline_calls <- merge(baseline_calls, snp_consistency, by = c("contig", "position"), all.x = TRUE)
xa_xi_map_full <- baseline_calls[fully_consistent == TRUE]
xa_xi_map_full[, xi_allele := fifelse(xa_allele == "ref", "alt", "ref")]
xa_xi_map_full <- merge(xa_xi_map_full,
                        unique(cc_snp[, .(contig, position, variantID, refAllele, altAllele)]),
                        by = c("contig", "position"))
xa_xi_map_full[, Xa_base := fifelse(xa_allele == "ref", refAllele, altAllele)]
xa_xi_map_full[, Xi_base := fifelse(xi_allele == "ref", refAllele, altAllele)]
fwrite(xa_xi_map_full, file.path(out_dir, "tables", "xa_xi_allele_calls_all_tiers.csv"))

xa_ref_frac <- mean(xa_xi_map_full$xa_allele == "ref")
xa_binom <- binom.test(sum(xa_xi_map_full$xa_allele == "ref"), nrow(xa_xi_map_full), p = 0.5)
cat(sprintf("\nref called Xa in %.1f%% of %d calls (expect ~50%% if genuine biology). Binomial p = %.3g\n",
            100 * xa_ref_frac, nrow(xa_xi_map_full), xa_binom$p.value))
if (xa_ref_frac > 0.7 | xa_ref_frac < 0.3) {
  warning("Xa calls skewed toward one allele class - check reference-mapping bias.")
}

xa_xi_map <- xa_xi_map_full[tier %in% tiers_for_main_analysis]
fwrite(xa_xi_map, file.path(out_dir, "tables", "xa_xi_allele_calls_main.csv"))
cat(sprintf("\nUsing %d SNPs (tiers: %s) for the main analysis.\n",
            nrow(xa_xi_map), paste(tiers_for_main_analysis, collapse = ", ")))

##  5. Reorient every condition's counts into Xa/Xi space 
cc_xaxi <- merge(cc_snp, xa_xi_map[, .(contig, position, gene, xa_allele, xi_allele, tier)],
                 by = c("contig", "position", "gene"))

for (cond in conditions) {
  rc <- paste0("refCount_", cond); ac <- paste0("altCount_", cond)
  if (!rc %in% names(cc_xaxi)) next
  cc_xaxi[, (paste0("xa_count_", cond)) := fifelse(xa_allele == "ref", get(rc), get(ac))]
  cc_xaxi[, (paste0("xi_count_", cond)) := fifelse(xa_allele == "ref", get(ac), get(rc))]
  cc_xaxi[, (paste0("xi_fraction_", cond)) :=
            get(paste0("xi_count_", cond)) / (get(paste0("xa_count_", cond)) + get(paste0("xi_count_", cond)))]
}

# Pairwise test: each comparison condition vs control_condition, SNP level.
# Results are stacked long-format with a `comparison` column rather than
# one hardcoded pair, so this works for any number of conditions.
xi_ctrl_col <- paste0("xi_count_", control_condition)
xa_ctrl_col <- paste0("xa_count_", control_condition)

cc_xaxi_long <- rbindlist(lapply(comparison_conditions, function(cond) {
  xi_col <- paste0("xi_count_", cond); xa_col <- paste0("xa_count_", cond)
  d <- cc_xaxi[!is.na(get(xi_ctrl_col)) & !is.na(get(xi_col)) &
                 (get(xa_ctrl_col) + get(xi_ctrl_col)) >= min_depth &
                 (get(xa_col) + get(xi_col)) >= min_depth]
  if (nrow(d) == 0) return(NULL)
  res <- d[, as.list(test_condition_effect(get(xi_ctrl_col), get(xa_ctrl_col), get(xi_col), get(xa_col))),
           by = .(contig, position, gene)]
  res[, comparison := cond]
  res[, xi_fraction_control := d[[paste0("xi_fraction_", control_condition)]]]
  res[, xi_fraction_condition := d[[paste0("xi_fraction_", cond)]]]
  res[, delta_xi := xi_fraction_condition - xi_fraction_control]
  res
}))
cc_xaxi_long[, fdr := p.adjust(p_value, "BH"), by = comparison]   # FDR controlled within each comparison
fwrite(cc_xaxi_long, file.path(out_dir, "tables", "xi_shift_snp_vs_control.csv"))

# Gene level: sum counts per gene per condition first, then same pairwise loop
gene_counts <- cc_xaxi[, c(list(n_snps = .N),
                           lapply(conditions, function(cond) {
                             xa <- paste0("xa_count_", cond); xi <- paste0("xi_count_", cond)
                             if (!xa %in% names(cc_xaxi)) return(NULL)
                             setNames(list(sum(get(xa), na.rm = TRUE), sum(get(xi), na.rm = TRUE)),
                                      c(xa, xi))
                           }) |> unlist(recursive = FALSE)),
                       by = gene]

cc_xaxi_gene_long <- rbindlist(lapply(comparison_conditions, function(cond) {
  xi_col <- paste0("xi_count_", cond); xa_col <- paste0("xa_count_", cond)
  d <- gene_counts[(get(xa_ctrl_col) + get(xi_ctrl_col)) >= min_depth &
                     (get(xa_col) + get(xi_col)) >= min_depth]
  if (nrow(d) == 0) return(NULL)
  res <- d[, as.list(test_condition_effect(get(xi_ctrl_col), get(xa_ctrl_col), get(xi_col), get(xa_col))),
           by = .(gene, n_snps)]
  res[, comparison := cond]
  res[, xi_fraction_control := d[[xi_ctrl_col]] / (d[[xa_ctrl_col]] + d[[xi_ctrl_col]])]
  res[, xi_fraction_condition := d[[xi_col]] / (d[[xa_col]] + d[[xi_col]])]
  res[, delta_xi := xi_fraction_condition - xi_fraction_control]
  res
}))
cc_xaxi_gene_long[, fdr := p.adjust(p_value, "BH"), by = comparison]
fwrite(cc_xaxi_gene_long, file.path(out_dir, "tables", "xi_shift_gene_vs_control.csv"))

##  6. Chromosome-wide %Xi per replicate, all conditions 
rep_xaxi <- merge(all_ase[paste(contig, position) %in% paste(xa_xi_map$contig, xa_xi_map$position)],
                  xa_xi_map[, .(contig, position, xa_allele)], by = c("contig", "position"))
rep_xaxi[, xa_count := fifelse(xa_allele == "ref", refCount, altCount)]
rep_xaxi[, xi_count := fifelse(xa_allele == "ref", altCount, refCount)]

rep_xi_summary <- rep_xaxi[, .(pct_Xi = 100 * sum(xi_count) / (sum(xa_count) + sum(xi_count))),
                           by = .(sample, condition, rep)]
rep_xi_summary[, condition := factor(condition, levels = conditions)]
fwrite(rep_xi_summary, file.path(out_dir, "tables", "pct_Xi_per_replicate.csv"))

rep_xi_wide <- dcast(rep_xi_summary, rep ~ condition, value.var = "pct_Xi")
pairwise_xi_tests <- rbindlist(lapply(comparison_conditions, function(cond) {
  if (!all(c(cond, control_condition) %in% names(rep_xi_wide))) return(NULL)
  complete <- rep_xi_wide[!is.na(get(cond)) & !is.na(get(control_condition))]
  if (nrow(complete) < 2) return(NULL)
  tt <- t.test(complete[[cond]], complete[[control_condition]], paired = TRUE)
  data.table(comparison = cond, control = control_condition,
             mean_diff = tt$estimate, p_value = tt$p.value)
}))
pairwise_xi_tests[, p_adj := p.adjust(p_value, "BH")]
cat("\nPaired t-tests (%Xi), each condition vs", control_condition, ":\n")
print(pairwise_xi_tests)
fwrite(pairwise_xi_tests, file.path(out_dir, "tables", "pct_Xi_pairwise_tests.csv"))


## Main Figures


##  Fig A: chromosome-wide %Xi, bar + replicate points, all conditions 
rep_xi_mean <- rep_xi_summary[, .(mean_pct_Xi = mean(pct_Xi), sd_pct_Xi = sd(pct_Xi)), by = condition]

# Annotate each non-control bar with its p-value vs control
label_df <- merge(rep_xi_mean, pairwise_xi_tests, by.x = "condition", by.y = "comparison", all.x = TRUE)
label_df[, label := fifelse(is.na(p_value), "", sprintf("p = %.3g", p_adj))]
label_df[, y_label := mean_pct_Xi + sd_pct_Xi + 0.05 * diff(range(rep_xi_summary$pct_Xi, na.rm = TRUE))]

p_A <- ggplot() +
  geom_col(data = rep_xi_mean, aes(x = condition, y = mean_pct_Xi, fill = condition),
           width = 0.6, alpha = 0.7) +
  geom_errorbar(data = rep_xi_mean, aes(x = condition, ymin = mean_pct_Xi - sd_pct_Xi,
                                        ymax = mean_pct_Xi + sd_pct_Xi), width = 0.15) +
  geom_jitter(data = rep_xi_summary, aes(x = condition, y = pct_Xi),
              width = 0.08, size = 2.5, color = "black") +
  geom_text(data = label_df, aes(x = condition, y = y_label, label = label), size = 3.2) +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(title = "Chromosome-wide Xi-derived transcription",
       subtitle = sprintf("Bars = mean \u00b1 SD, points = replicates | each condition vs control (%s), paired t-test, BH-adjusted",
                          control_condition),
       x = NULL, y = "% of reads from Xi allele", fill = NULL)
ggsave(file.path(out_dir, "plots", "FigA_pct_Xi_by_condition.png"), p_A, width = 6, height = 5)

# Shared styling for every faceted figure below - visible panel borders and
# real spacing between facets, so it's unambiguous where one condition's
# panel ends and the next begins.
facet_style <- theme(
  panel.spacing    = unit(1.3, "lines"),
  panel.border     = element_rect(color = "grey30", fill = NA, linewidth = 0.7),
  strip.background = element_rect(fill = "grey85", color = "grey30"),
  strip.text       = element_text(face = "bold", size = 10)
)

##  Fig B: gene-level volcano, faceted by comparison 
# FDR values of exactly/effectively 0 give -log10(FDR) = Inf, which pushes
# those genes off the top of the plot (or collapses the axis). Cap them at
# a data-driven ceiling instead, and mark capped points with a distinct
# shape (triangle) so they're visibly flagged as "at least this significant"
# rather than looking like an ordinary point sitting near the cap.
nonzero_fdr <- cc_xaxi_gene_long$fdr[cc_xaxi_gene_long$fdr > 0 & is.finite(cc_xaxi_gene_long$fdr)]
y_cap_b <- if (length(nonzero_fdr) > 0) ceiling(-log10(min(nonzero_fdr))) + 2 else 10

cc_xaxi_gene_long[, neg_log10_fdr := -log10(fdr)]
cc_xaxi_gene_long[, is_capped := !is.finite(neg_log10_fdr) | neg_log10_fdr > y_cap_b]
cc_xaxi_gene_long[, neg_log10_fdr_plot := pmin(neg_log10_fdr, y_cap_b)]
cc_xaxi_gene_long[!is.finite(neg_log10_fdr_plot), neg_log10_fdr_plot := y_cap_b]


p_B <- ggplot(cc_xaxi_gene_long, aes(x = delta_xi, y = neg_log10_fdr_plot)) +
  geom_point(data = cc_xaxi_gene_long[fdr >= alpha], aes(shape = is_capped),
             color = "grey75", size = 2, alpha = 0.5) +
  geom_point(data = cc_xaxi_gene_long[fdr < alpha], aes(color = delta_xi > 0, shape = is_capped),
             size = 2.2) +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "blue"),
                     labels = c(`TRUE` = "increased", `FALSE` = "decreased")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c(`FALSE` = "observed", `TRUE` = "capped (FDR \u2248 0)")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
  geom_hline(yintercept = y_cap_b, linetype = "dotted", color = "grey40") +
  facet_wrap(~comparison) +
  theme_minimal() + facet_style +
  labs(title = sprintf("Gene-level Xi shift vs control (%s)", control_condition),
       subtitle = sprintf("Points above the dotted line have FDR \u2248 0 and are capped at %d for display", y_cap_b),
       x = "\u0394 Xi fraction (condition - control)", y = "-log10(FDR), capped",
       color = "FDR < 0.05", shape = NULL)
ggsave(file.path(out_dir, "plots", "FigB_gene_volcano_by_comparison.png"), p_B,
       width = 4 + 3 * length(comparison_conditions), height = 5)




##  Fig C: genome-wide positional plot, faceted by comparison
cc_xaxi_long_X <- cc_xaxi_long[contig %in% c("X", "chrX")]
p_C <- ggplot(cc_xaxi_long_X, aes(x = position, y = delta_xi)) +
  geom_point(data = cc_xaxi_long_X[fdr >= alpha], color = "grey75", size = 1, alpha = 0.5) +
  geom_point(data = cc_xaxi_long_X[fdr < alpha], aes(color = delta_xi > 0), size = 1.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  annotate("rect", xmin = xic_start, xmax = xic_end, ymin = -Inf, ymax = Inf,
           fill = "blue", alpha = 0.12) +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "blue"),
                     labels = c(`TRUE` = "increased", `FALSE` = "decreased")) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = " Mb")) +
  facet_wrap(~comparison, ncol = 1) +
  theme_minimal() + facet_style +
  labs(title = sprintf("Xi shift across the X chromosome, vs control (%s)", control_condition),
       subtitle = "Shaded band = XIC region, GRCh38 | grey = FDR \u2265 0.05 | red = increased Xi fraction, blue = decreased (FDR < 0.05)",
       x = "Position", y = "\u0394 Xi fraction", color = NULL)
ggsave(file.path(out_dir, "plots", "FigC_genomewide_position_by_comparison.png"), p_C,
       width = 10, height = 3 * length(comparison_conditions))


##  Fig D: multi-gene, multi-replicate, multi-condition heatmap 
gene_rep <- merge(all_ase[paste(contig, position) %in% paste(xa_xi_map$contig, xa_xi_map$position)],
                  xa_xi_map[, .(contig, position, gene, xa_allele)], by = c("contig", "position"))
gene_rep[, xa_count := fifelse(xa_allele == "ref", refCount, altCount)]
gene_rep[, xi_count := fifelse(xa_allele == "ref", altCount, refCount)]
gene_rep_summary <- gene_rep[, .(xi_fraction = sum(xi_count) / (sum(xa_count) + sum(xi_count))),
                             by = .(gene, sample, condition, rep)]
gene_rep_summary[, condition := factor(condition, levels = conditions)]
gene_rep_summary[, condition_rep := paste0(condition, "_rep", rep)]

heat_wide <- dcast(gene_rep_summary, gene ~ condition_rep, value.var = "xi_fraction")
desired_cols <- c("gene", unlist(lapply(conditions, function(cond)
  paste0(cond, "_rep", sort(unique(gene_rep_summary[condition == cond]$rep))))))
setcolorder(heat_wide, intersect(desired_cols, names(heat_wide)))
heat_mat <- as.matrix(heat_wide[, -1]); rownames(heat_mat) <- heat_wide$gene
heat_mat <- heat_mat[complete.cases(heat_mat), , drop = FALSE]

if (nrow(heat_mat) >= 2) {
  pheatmap(heat_mat, cluster_cols = FALSE,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           main = "Xi fraction per gene, per replicate, all conditions",
           filename = file.path(out_dir, "plots", "FigD_gene_replicate_heatmap.png"),
           width = 6 + length(conditions), height = max(4, nrow(heat_mat) * 0.25))
} else {
  warning("Fewer than 2 genes have complete coverage across all conditions/replicates - Fig D skipped.")
}

##  Fig E: candidate gene spotlight, trajectory across ALL conditions 
# Kept as connected lines deliberately (unlike Fig A) - for a timecourse
# like the capacitation set, the per-replicate trajectory across ordered
# conditions is exactly the thing worth seeing, not just a 2-group split.
spotlight_genes <- unique(c(goi, head(cc_xaxi_gene_long[order(fdr)]$gene, 6)))
spotlight_genes <- intersect(spotlight_genes, gene_rep_summary$gene)

if (length(spotlight_genes) > 0) {
  spotlight_data <- gene_rep_summary[gene %in% spotlight_genes]
  complete_pairs <- spotlight_data[, .N, by = .(gene, rep)][N == length(conditions), .(gene, rep)]
  n_dropped <- nrow(unique(spotlight_data[, .(gene, rep)])) - nrow(complete_pairs)
  if (n_dropped > 0) cat(sprintf("Fig E: dropping %d gene x replicate pair(s) missing >=1 condition.\n", n_dropped))
  spotlight_data <- merge(spotlight_data, complete_pairs, by = c("gene", "rep"))

  p_E <- ggplot(spotlight_data, aes(x = condition, y = xi_fraction, group = rep)) +
    geom_line(color = "grey60") +
    geom_point(aes(color = condition), size = 2.5) +
    scale_color_viridis_d() +
    facet_wrap(~gene, scales = "free_y") +
    theme_minimal() + facet_style +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Xi fraction trajectory across conditions (genes of interest / top hits)",
         x = NULL, y = "Xi fraction", color = NULL)
  ggsave(file.path(out_dir, "plots", "FigE_candidate_gene_spotlight.png"), p_E,
         width = 8, height = ceiling(length(spotlight_genes) / 3) * 2.5)
} else {
  warning("None of the genes of interest have confident Xa/Xi calls - Fig E skipped.")
}

##  Fig F: total Xi shift per gene, one bar plot per comparison 
for (cond in comparison_conditions) {
  sig_genes_bar <- cc_xaxi_gene_long[comparison == cond & fdr < alpha]
  if (nrow(sig_genes_bar) == 0) {
    warning(sprintf("No genes pass FDR < %.2f for %s vs %s - Fig F skipped for this comparison.",
                    alpha, cond, control_condition))
    next
  }
  p_F <- ggplot(sig_genes_bar, aes(x = reorder(gene, delta_xi), y = delta_xi, fill = delta_xi)) +
    geom_col() + coord_flip() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    labs(title = sprintf("Total Xi shift per gene: %s vs %s (FDR < %.2f)", cond, control_condition, alpha),
         x = "Gene", y = "\u0394 Xi fraction", fill = "\u0394 Xi\nfraction")
  ggsave(file.path(out_dir, "plots", sprintf("FigF_gene_shift_barplot_%s_vs_%s.png", cond, control_condition)),
         p_F, width = 7, height = max(4, nrow(sig_genes_bar) * 0.25))
}

cat("\nDone. Outputs in:", normalizePath(out_dir), "\n")



## Violin plot from existing capacitation ASE output
## Uses: xi_shift_gene_vs_control.csv


library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

## Define output directory if needed
out_dir <- "/home/jvk3/Desktop/HDAC_counts/All three datsets bulk and ASE/Capacitation/ASE/ase_outputs_Xi_fraction_capacitation"

## Load gene-level Xi shift table
xi_gene <- data.table::fread(
  file.path(out_dir, "tables", "xi_shift_gene_vs_control.csv")
)

## Check columns
colnames(xi_gene)

## Convert comparison table into long format:
## Naive comes from xi_fraction_control
## Day3/Day7/Primed come from xi_fraction_condition
cap_xi_violin_df <- dplyr::bind_rows(
  xi_gene %>%
    dplyr::select(
      gene,
      n_snps,
      xi_fraction = xi_fraction_control
    ) %>%
    dplyr::mutate(condition = "Naive"),

  xi_gene %>%
    dplyr::select(
      gene,
      n_snps,
      condition = comparison,
      xi_fraction = xi_fraction_condition
    )
) %>%
  dplyr::mutate(
    condition_clean = dplyr::case_when(
      condition %in% c("Naive", "naive") ~ "Naive",
      condition %in% c("Day3", "Cap. Day 3", "Cap_Day_3") ~ "Day3",
      condition %in% c("Day7", "Cap. Day 7", "Cap_Day_7") ~ "Day7",
      condition %in% c("Primed", "primed") ~ "Primed",
      TRUE ~ as.character(condition)
    ),
    condition_clean = factor(
      condition_clean,
      levels = c("Naive", "Day3", "Day7", "Primed")
    )
  ) %>%
  dplyr::filter(
    !is.na(xi_fraction),
    !is.na(condition_clean)
  ) %>%
  dplyr::distinct(
    gene,
    condition_clean,
    .keep_all = TRUE
  )

## Check number of genes per timepoint
cap_xi_violin_df %>%
  dplyr::count(condition_clean)


## Plot violin: Xi/(Xa+Xi) across capacitation


p_cap_xi_violin <- ggplot(
  cap_xi_violin_df,
  aes(x = condition_clean, y = xi_fraction)
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.5,
    fill = "grey85",
    colour = "black"
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    colour = "black"
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.5,
    alpha = 0.6
  ) +
  theme_bw() +
  labs(
    title = "Gene-level Xi contribution across capacitation",
    subtitle = "Each point represents one informative X-linked gene",
    x = NULL,
    y = "Xi/(Xa+Xi)"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank()
  )

p_cap_xi_violin   

ggsave(
  file.path(
    out_dir,
    "plots",
    "Capacitation_gene_level_Xi_fraction_violin.png"
  ),
  p_cap_xi_violin,
  width = 7,
  height = 5,
  dpi = 300
)   



## Violin plot from existing capacitation ASE output
## Uses: xi_shift_gene_vs_control.csv
## Subset to 6 genes of interest


library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

## Define output directory
out_dir <- "/home/jvk3/Desktop/HDAC_counts/All three datsets bulk and ASE/Capacitation/ASE/ase_outputs_Xi_fraction_capacitation"

## Define genes of interest
genes_of_interest <- c("APOO", "CCDC22", "DLG3", "GRIPAP1", "HUWE1", "PDK3")

## Load gene-level Xi shift table
xi_gene <- data.table::fread(
  file.path(out_dir, "tables", "xi_shift_gene_vs_control.csv")
)

## Check columns
colnames(xi_gene)

## Check whether gene symbols match exactly
intersect(genes_of_interest, unique(xi_gene$gene))
setdiff(genes_of_interest, unique(xi_gene$gene))

## Convert comparison table into long format:
## Naive comes from xi_fraction_control
## Day3/Day7/Primed come from xi_fraction_condition
cap_xi_violin_df <- dplyr::bind_rows(
  xi_gene %>%
    dplyr::select(
      gene,
      n_snps,
      xi_fraction = xi_fraction_control
    ) %>%
    dplyr::mutate(condition = "Naive"),

  xi_gene %>%
    dplyr::select(
      gene,
      n_snps,
      condition = comparison,
      xi_fraction = xi_fraction_condition
    )
) %>%
  dplyr::mutate(
    condition_clean = dplyr::case_when(
      condition %in% c("Naive", "naive") ~ "Naive",
      condition %in% c("Day3", "Cap. Day 3", "Cap_Day_3") ~ "Day3",
      condition %in% c("Day7", "Cap. Day 7", "Cap_Day_7") ~ "Day7",
      condition %in% c("Primed", "primed") ~ "Primed",
      TRUE ~ as.character(condition)
    ),
    condition_clean = factor(
      condition_clean,
      levels = c("Naive", "Day3", "Day7", "Primed")
    )
  ) %>%
  dplyr::filter(
    !is.na(xi_fraction),
    !is.na(condition_clean)
  ) %>%
  dplyr::distinct(
    gene,
    condition_clean,
    .keep_all = TRUE
  ) %>%
  dplyr::filter(
    gene %in% genes_of_interest
  )

## Check number of genes per timepoint
cap_xi_violin_df %>%
  dplyr::count(condition_clean)

## Optional: inspect the actual values
cap_xi_violin_df %>%
  dplyr::arrange(gene, condition_clean)


## Plot violin: Xi/(Xa+Xi) across capacitation
## Now restricted to the 6 genes of interest


p_cap_xi_violin <- ggplot(
  cap_xi_violin_df,
  aes(x = condition_clean, y = xi_fraction)
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.4,
    fill = "grey85",
    colour = "black"
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    colour = "black"
  ) +
  geom_line(
    aes(group = gene),
    colour = "grey60",
    alpha = 0.6,
    linewidth = 0.3
  ) +
  geom_point(
    aes(colour = gene),
    size = 3,
    alpha = 0.9
  ) +
  theme_bw() +
  labs(
    title = "Gene-level Xi contribution across capacitation",
    subtitle = "Six selected informative X-linked genes",
    x = NULL,
    y = "Xi/(Xa+Xi)",
    colour = "Gene"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank()
  )

p_cap_xi_violin

ggsave(
  file.path(
    out_dir,
    "plots",
    "Capacitation_gene_level_Xi_fraction_violin_6genes.png"
  ),
  p_cap_xi_violin,
  width = 7,
  height = 5,
  dpi = 300
)

## Stats for capacitation Xi/(Xa+Xi) violin plot
## Paired by gene, each condition compared against Naive


cap_xi_violin_stats <- cap_xi_violin_df %>%
  dplyr::filter(
    condition_clean %in% c("Naive", "Day3", "Day7", "Primed")
  ) %>%
  dplyr::select(gene, condition_clean, xi_fraction) %>%
  tidyr::pivot_wider(
    names_from = condition_clean,
    values_from = xi_fraction
  ) %>%
  tidyr::drop_na(Naive) %>%
  dplyr::summarise(
    Day3_p = wilcox.test(Day3, Naive, paired = TRUE, exact = FALSE)$p.value,
    Day7_p = wilcox.test(Day7, Naive, paired = TRUE, exact = FALSE)$p.value,
    Primed_p = wilcox.test(Primed, Naive, paired = TRUE, exact = FALSE)$p.value
  ) %>%
  tidyr::pivot_longer(
    everything(),
    names_to = "comparison",
    values_to = "p_value"
  ) %>%
  dplyr::mutate(
    condition_clean = dplyr::case_when(
      comparison == "Day3_p" ~ "Day3",
      comparison == "Day7_p" ~ "Day7",
      comparison == "Primed_p" ~ "Primed"
    ),
    padj = p.adjust(p_value, method = "BH"),
    sig_label = dplyr::case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

cap_xi_violin_stats

cap_xi_violin_summary <- cap_xi_violin_df %>%
  dplyr::group_by(condition_clean) %>%
  dplyr::summarise(
    n_genes = dplyr::n_distinct(gene),
    median_xi_fraction = median(xi_fraction, na.rm = TRUE),
    q1 = quantile(xi_fraction, 0.25, na.rm = TRUE),
    q3 = quantile(xi_fraction, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

cap_xi_violin_summary   

genes_of_interest <- c("APOO", "CCDC22", "DLG3", "GRIPAP1", "HUWE1", "PDK3")

xi_gene %>%
  dplyr::filter(gene %in% genes_of_interest) %>%
  dplyr::select(
    gene,
    comparison,
    n_snps,
    xi_fraction_control,
    xi_fraction_condition,
    delta_xi,
    p_value,
    fdr
  ) %>%
  dplyr::arrange(comparison, fdr)




## Capacitation ASE X chromosome plot
## Informative X-linked genes + Tukiainen categories


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(readr)
library(stringr)
library(tibble)
library(data.table)


## 1. Paths


out_dir <- "/home/jvk3/Desktop/HDAC_counts/All three datsets bulk and ASE/Capacitation/ASE/ase_outputs_Xi_fraction_capacitation"

cap_ase_xchr_dir <- file.path(
  out_dir,
  "plots",
  "X_chromosome_position_plots"
)

dir.create(
  cap_ase_xchr_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


## 2. Load ASE gene-level Xi shift table


xi_gene <- data.table::fread(
  file.path(out_dir, "tables", "xi_shift_gene_vs_control.csv")
)

colnames(xi_gene)


## 3. Load local gene-coordinate cache


gene_coord_annot <- readRDS(
  file.path(out_dir, "tables", "local_gene_coordinates_cache.rds")
)

gene_coord_annot <- as.data.frame(gene_coord_annot)

colnames(gene_coord_annot)
head(gene_coord_annot)   


## 4. Standardise gene-coordinate columns


pick_col <- function(df, candidates, label) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    stop(
      paste0(
        "Could not find column for ", label,
        ". Available columns are: ",
        paste(colnames(df), collapse = ", ")
      )
    )
  }
  hit[1]
}

symbol_col <- pick_col(
  gene_coord_annot,
  c("SYMBOL", "gene", "external_gene_name", "hgnc_symbol", "gene_name"),
  "gene symbol"
)

chr_col <- pick_col(
  gene_coord_annot,
  c("CHR", "chr", "chromosome_name", "seqnames", "chromosome"),
  "chromosome"
)

start_col <- pick_col(
  gene_coord_annot,
  c("START", "start", "start_position", "gene_start", "start_position_bp"),
  "gene start"
)

end_col <- pick_col(
  gene_coord_annot,
  c("END", "end", "end_position", "gene_end", "end_position_bp"),
  "gene end"
)

gene_coord_annot_clean <- gene_coord_annot %>%
  dplyr::transmute(
    SYMBOL = stringr::str_squish(as.character(.data[[symbol_col]])),
    CHR = as.character(.data[[chr_col]]),
    START = as.numeric(.data[[start_col]]),
    END = as.numeric(.data[[end_col]])
  ) %>%
  dplyr::mutate(
    CHR = dplyr::case_when(
      CHR %in% c("23", "chrX", "X") ~ "X",
      TRUE ~ CHR
    )
  ) %>%
  dplyr::filter(
    !is.na(SYMBOL),
    SYMBOL != "",
    !is.na(START),
    !is.na(END)
  ) %>%
  dplyr::distinct()   


## 5. Load Tukiainen XCI categories


tukiainen_path <- "/home/jvk3/Desktop/HDAC_counts/XCI_classifications/XCI_classification_outputs/Tukiainen_XCI_categories_clean.csv"

tukiainen_annot <- readr::read_csv(
  tukiainen_path,
  show_col_types = FALSE
) %>%
  dplyr::select(SYMBOL, tukiainen_status_grouped) %>%
  dplyr::mutate(
    SYMBOL = stringr::str_squish(SYMBOL),
    tukiainen_status_grouped = stringr::str_squish(tukiainen_status_grouped)
  ) %>%
  dplyr::distinct()   


## 6. Build ASE X chromosome plotting table


xist_start_bp <- 73817774
xist_end_bp   <- 73852754
xist_mid_mb   <- mean(c(xist_start_bp, xist_end_bp)) / 1e6

cap_ase_xchr_plot_df <- xi_gene %>%
  as.data.frame() %>%
  dplyr::mutate(
    SYMBOL = stringr::str_squish(as.character(gene)),
    comparison = as.character(comparison)
  ) %>%
  dplyr::left_join(
    gene_coord_annot_clean,
    by = "SYMBOL"
  ) %>%
  dplyr::left_join(
    tukiainen_annot,
    by = "SYMBOL"
  ) %>%
  dplyr::filter(
    CHR == "X",
    comparison %in% c(
      "Day3",
      "Day7",
      "Primed",
      "Day3_vs_Naive",
      "Day7_vs_Naive",
      "Primed_vs_Naive",
      "Day3 vs Naive",
      "Day7 vs Naive",
      "Primed vs Naive"
    ),
    !is.na(START),
    !is.na(END),
    !is.na(delta_xi)
  ) %>%
  dplyr::mutate(
    gene_mid_mb = ((START + END) / 2) / 1e6,
    distance_to_xist_mb = abs(gene_mid_mb - xist_mid_mb),

    comparison_clean = dplyr::case_when(
      comparison %in% c("Day3", "Day3_vs_Naive", "Day3 vs Naive") ~ "Day3 vs Naive",
      comparison %in% c("Day7", "Day7_vs_Naive", "Day7 vs Naive") ~ "Day7 vs Naive",
      comparison %in% c("Primed", "Primed_vs_Naive", "Primed vs Naive") ~ "Primed vs Naive",
      TRUE ~ comparison
    ),

    comparison_clean = factor(
      comparison_clean,
      levels = c(
        "Day3 vs Naive",
        "Day7 vs Naive",
        "Primed vs Naive"
      )
    ),

    xci_plot_group = dplyr::case_when(
      tukiainen_status_grouped == "Subject to XCI" ~ "Subject to XCI",
      tukiainen_status_grouped == "Escapee" ~ "Escapee",
      TRUE ~ "Other"
    ),

    xci_plot_group = factor(
      xci_plot_group,
      levels = c("Subject to XCI", "Escapee", "Other")
    ),

    is_sig = !is.na(fdr) & fdr < 0.05,

    ase_direction = dplyr::case_when(
      is_sig & delta_xi > 0 ~ "Increased Xi fraction",
      is_sig & delta_xi < 0 ~ "Decreased Xi fraction",
      TRUE ~ "Not significant"
    )
  )   

cap_ase_xchr_plot_df %>%
  dplyr::count(comparison_clean, xci_plot_group)

cap_ase_xchr_plot_df %>%
  dplyr::count(comparison_clean, ase_direction)   


## 7. Select genes to label


cap_ase_xchr_label_df <- cap_ase_xchr_plot_df %>%
  dplyr::filter(is_sig) %>%
  dplyr::group_by(comparison_clean) %>%
  dplyr::slice_max(
    order_by = abs(delta_xi),
    n = 6,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

cap_ase_escapee_label_df <- cap_ase_xchr_plot_df %>%
  dplyr::filter(
    is_sig,
    xci_plot_group == "Escapee"
  ) %>%
  dplyr::group_by(comparison_clean) %>%
  dplyr::slice_max(
    order_by = abs(delta_xi),
    n = 3,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

cap_ase_xchr_label_df <- dplyr::bind_rows(
  cap_ase_xchr_label_df,
  cap_ase_escapee_label_df
) %>%
  dplyr::distinct(
    comparison_clean,
    SYMBOL,
    .keep_all = TRUE
  )   


## 8. Plot ASE X chromosome position plot


cap_ase_xchr_plot <- ggplot(
  cap_ase_xchr_plot_df,
  aes(
    x = gene_mid_mb,
    y = delta_xi,
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
    data = cap_ase_xchr_label_df,
    aes(label = SYMBOL),
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
    ~ comparison_clean,
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
    title = "Allele-specific Xi-fraction changes across chromosome X during capacitation",
    subtitle = "Informative X-linked genes are shown; colours indicate Tukiainen XCI category grouping",
    x = "Position on chromosome X (Mb)",
    y = expression(Delta~"Xi fraction relative to Naive"),
    colour = "XCI category"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

cap_ase_xchr_plot   

ggsave(
  filename = file.path(
    cap_ase_xchr_dir,
    "Capacitation_ASE_Xchr_informative_X_genes_Tukiainen_categories.png"
  ),
  plot = cap_ase_xchr_plot,
  width = 9,
  height = 8,
  dpi = 300,
  bg = "white"
)

readr::write_csv(
  cap_ase_xchr_plot_df,
  file.path(
    cap_ase_xchr_dir,
    "Capacitation_ASE_Xchr_informative_X_genes_Tukiainen_categories_table.csv"
  )
)   


## Stats for capacitation ASE X chromosome plot


library(dplyr)
library(tidyr)
library(purrr)
library(readr)


## 1. Basic summary of Xi-fraction shifts


cap_ase_xchr_summary <- cap_ase_xchr_plot_df %>%
  dplyr::group_by(comparison_clean) %>%
  dplyr::summarise(
    n_genes = dplyr::n(),
    n_sig = sum(is_sig, na.rm = TRUE),
    n_increased = sum(delta_xi > 0, na.rm = TRUE),
    n_decreased = sum(delta_xi < 0, na.rm = TRUE),
    n_sig_increased = sum(is_sig & delta_xi > 0, na.rm = TRUE),
    n_sig_decreased = sum(is_sig & delta_xi < 0, na.rm = TRUE),
    median_delta_xi = median(delta_xi, na.rm = TRUE),
    mean_delta_xi = mean(delta_xi, na.rm = TRUE),
    min_delta_xi = min(delta_xi, na.rm = TRUE),
    max_delta_xi = max(delta_xi, na.rm = TRUE),
    .groups = "drop"
  )

cap_ase_xchr_summary


## 2. Wilcoxon signed-rank test: is median ΔXi different from 0?


cap_ase_delta_wilcox <- cap_ase_xchr_plot_df %>%
  dplyr::group_by(comparison_clean) %>%
  dplyr::group_modify(~{
    test_df <- .x %>%
      dplyr::filter(!is.na(delta_xi))

    test_result <- wilcox.test(
      test_df$delta_xi,
      mu = 0,
      alternative = "two.sided",
      exact = FALSE
    )

    tibble::tibble(
      n_genes = nrow(test_df),
      median_delta_xi = median(test_df$delta_xi, na.rm = TRUE),
      p_value = test_result$p.value
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    padj = p.adjust(p_value, method = "BH")
  )

cap_ase_delta_wilcox   


## 3. Binomial directionality test
## Tests whether decreased Xi-fraction shifts are more frequent than expected by chance


cap_ase_direction_binom <- cap_ase_xchr_plot_df %>%
  dplyr::filter(delta_xi != 0) %>%
  dplyr::group_by(comparison_clean) %>%
  dplyr::summarise(
    n_decreased = sum(delta_xi < 0, na.rm = TRUE),
    n_increased = sum(delta_xi > 0, na.rm = TRUE),
    n_total = n_decreased + n_increased,
    p_value = binom.test(
      x = n_decreased,
      n = n_total,
      p = 0.5,
      alternative = "greater"
    )$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    padj = p.adjust(p_value, method = "BH")
  )

cap_ase_direction_binom   
   
