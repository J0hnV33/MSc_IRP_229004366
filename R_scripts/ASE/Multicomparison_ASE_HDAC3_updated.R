
# XCI analysis for two or more ordered conditions.
# conditions controls plotting order, while control_condition is used both to
# call Xa/Xi alleles and as the reference for statistical comparisons. Each
# non-control condition is compared with the control; all pairwise comparisons
# are not performed.

library(data.table)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(readxl)
library(GenomicRanges)
library(biomaRt)
library(VGAM)

## CONFIG
base_dir <- "/home/jvk3/Desktop/HDAC_counts/All three datsets bulk and ASE/HDAC3/ASE"

ase_dir <- file.path(base_dir, "8_ASE_counts")

meta_path <- file.path(base_dir, "meta_table_flipped_HDAC3_ase.csv")

out_dir <- file.path(base_dir, "ase_outputs_Xi_fraction_HDAC3")

min_depth <- 10
alpha <- 0.05

goi <- c("ATRX", "THOC2", "HUWE1", "G6PD")

tukiainen_csv <- "/home/jvk3/Desktop/HDAC_counts/XCI_classifications/XCI_classification_outputs/Tukiainen_XCI_categories_clean.csv"
# edith these two lines per dataset
conditions <- c("Scramble", "HDAC3_siRNA")   # plotting order
control_condition <- "Scramble"                             # calling baseline + comparison reference
# For the capacitation set, use instead:
# CONDITIONS        <- c("naive", "cap1", "cap2", "primed")
# CONTROL_CONDITION <- "primed"


stopifnot("CONTROL_CONDITION must be one of CONDITIONS" = control_condition %in% conditions)
comparison_conditions <- setdiff(conditions, control_condition)   # every non-reference condition

candidate_xci_statuses <- c("Subject to XCI")

xa_tiers <- data.table(
  tier      = c("high", "medium", "low"),
  min_skew  = c(0.90,   0.80,     0.75),
  min_depth = c(15,     20,       30)
)
tiers_for_main_analysis <- c("high", "medium")

xist_start <- 73820651   # GRCh38 - verify against your build
xist_end   <- 73852753

dir.create(file.path(out_dir, "plots"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

##  1. Load ASE counts + metadata, QC filter
meta <- fread(meta_path)

setnames(meta, make.names(names(meta)))

stopifnot("metadata needs a 'rep' column" = "rep" %in% names(meta))
stopifnot("metadata needs a 'condition' column" = "condition" %in% names(meta))
stopifnot("metadata needs a 'filename' column" = "filename" %in% names(meta))
stopifnot("metadata needs a 'sample' column" = "sample" %in% names(meta))

stopifnot("metadata condition values must all be in CONDITIONS" =
            all(unique(meta$condition) %in% conditions))

meta[, filepath := file.path(ase_dir, filename)]

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

##  3. Gene + XCI status annotation
genes_cache <- file.path(out_dir, "tables", "ensembl_genes_cache.rds")
if (file.exists(genes_cache)) {
  genes <- readRDS(genes_cache)
} else {
  mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
  genes <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype",
                   "chromosome_name", "start_position", "end_position", "strand"),
    filters = "chromosome_name", values = unique(combined$contig), mart = mart
  )
  setDT(genes)
  setnames(genes, c("chromosome_name", "start_position", "end_position"), c("chr", "start", "end"))
  genes[, strand := ifelse(strand == 1, "+", "-")]
  saveRDS(genes, genes_cache)
}
genes_gr <- makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)

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

if (file.exists(tukiainen_csv)) {
  
  xci_lookup <- fread(tukiainen_csv)
  
  xci_lookup <- unique(
    xci_lookup[, .(
      gene = SYMBOL,
      XCI_status = tukiainen_status_grouped
    )],
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
  warning("Tukiainen CSV not found - XCI_status will be NA throughout.")
  snp_annotations[, XCI_status := NA_character_]
}
snp_annotations[, condition := factor(condition, levels = conditions)]

## = 4. Wide table: one row per SNP, one column set per condition
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

##  4b. Call Xa/Xi per SNP from the CONTROL_CONDITION baseline
ctrl_ratio_col <- paste0("ratio_", control_condition)
ctrl_depth_col <- paste0("depth_", control_condition)
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

cat(sprintf("\nXa/Xi candidate calls by tier (baseline = %s):\n", control_condition))
print(table(baseline_calls$tier))

# Replicate-consistency check within the control condition only
rep_check <- rep_wide[condition == control_condition]
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

##  5. Reorient EVERY condition's counts into Xa/Xi space
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

# Pairwise test: each comparison condition vs CONTROL_CONDITION, SNP level.
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

## 
##  Main figures 
## 

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

##  Fig B: gene-level volcano, faceted by comparison 
p_B <- ggplot(cc_xaxi_gene_long, aes(x = delta_xi, y = -log10(fdr), color = fdr < alpha)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("grey60", "firebrick")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
  facet_wrap(~comparison) +
  theme_minimal() +
  labs(title = sprintf("Gene-level Xi shift vs control (%s)", control_condition),
       x = "\u0394 Xi fraction (condition - control)", y = "-log10(FDR)", color = "FDR < 0.05")
ggsave(file.path(out_dir, "plots", "FigB_gene_volcano_by_comparison.png"), p_B,
       width = 4 + 3 * length(comparison_conditions), height = 5)

## Fig C: genome-wide positional plot, faceted by comparison 
cc_xaxi_long_X <- cc_xaxi_long[contig %in% c("X", "chrX")]
p_C <- ggplot(cc_xaxi_long_X, aes(x = position, y = delta_xi, color = fdr < alpha)) +
  geom_point(alpha = 0.6, size = 1.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = (xist_start + xist_end) / 2, linetype = "dotted", color = "blue") +
  scale_color_manual(values = c("grey60", "firebrick")) +
  facet_wrap(~comparison, ncol = 1) +
  theme_minimal() +
  labs(title = sprintf("Xi shift across the X chromosome, vs control (%s)", control_condition),
       x = "Position (bp)", y = "\u0394 Xi fraction", color = "FDR < 0.05")
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
# like the capacitation set, the per-replicate trajectory ACROSS ordered
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
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Xi fraction trajectory across conditions (genes of interest / top hits)",
         x = NULL, y = "Xi fraction", color = NULL)
  ggsave(file.path(out_dir, "plots", "FigE_candidate_gene_spotlight.png"), p_E,
         width = 8, height = ceiling(length(spotlight_genes) / 3) * 2.5)
} else {
  warning("None of the genes of interest have confident Xa/Xi calls - Fig E skipped.")
}

## Fig F: total Xi shift per gene, one bar plot per comparison
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
