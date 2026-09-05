#!/usr/bin/env Rscript
#
# Differential binding under an externally supplied scale factor, in both
# edgeR and DESeq2, from one shared count matrix.
#
# THE PARAMETERIZATION TRAP
# -------------------------
# The two engines take an external scaling in different units, and the
# difference is silent: pass the same numbers to both and you get results that
# disagree by tens of percent with no warning, no error, and plausible-looking
# output.
#
#   edgeR   effective library size = lib.size * norm.factors
#           -> norm.factors must be set to (1/s) / lib.size, then centred
#              so their product is 1. norm.factors is a CORRECTION to depth.
#
#   DESeq2  normalized counts = counts / sizeFactors
#           -> sizeFactors must be set to (1/s) directly. sizeFactors is the
#              COMPLETE per-sample divisor, depth already included.
#
# Multiplying a DESeq2 size factor by library size, as if it were an edgeR
# norm factor, double-counts sequencing depth. In our hands that manufactured
# a ~40% apparent discrepancy between two engines that in fact agreed. The
# check at the bottom of this file fails loudly if the parameterization ever
# drifts back.
#
# s = the multiplicative scale applied to a sample's counts, from
# scale_factors.tsv. Larger s means "this sample's signal must be scaled up".

suppressPackageStartupMessages({
  library(edgeR); library(DESeq2); library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts",      type = "character"),
  make_option("--samplesheet", type = "character"),
  make_option("--scale-factors", type = "character", dest = "scale_factors"),
  make_option("--method",      type = "character"),
  make_option("--numerator",   type = "character"),
  make_option("--denominator", type = "character"),
  make_option("--fdr",         type = "double", default = 0.05),
  make_option("--control-target", type = "character", default = "IgG",
              dest = "control_target"),
  make_option("--out",         type = "character")
)))

## ---- effective-size helpers: the whole point of this script -----------------

# edgeR wants a depth CORRECTION whose product across samples is 1.
edger_norm_factors <- function(s, lib_size) {
  nf <- (1 / s) / lib_size
  nf / exp(mean(log(nf)))            # centre; edgeR requires prod(nf) == 1
}

# DESeq2 wants the COMPLETE divisor. Depth is already in it. Do not multiply
# this by lib.size.
deseq_size_factors <- function(s) {
  sf <- 1 / s
  sf / exp(mean(log(sf)))
}

## ---- load -------------------------------------------------------------------

cts_raw <- read.delim(opt$counts, check.names = FALSE)
regions <- sprintf("%s:%d-%d", cts_raw[[1]], cts_raw[[2]], cts_raw[[3]])
cts     <- as.matrix(cts_raw[, -(1:3), drop = FALSE])
rownames(cts) <- regions

meta <- read.csv(opt$samplesheet)
sf   <- read.delim(opt$scale_factors, row.names = 1, check.names = FALSE)

# Control libraries share a condition label with the samples they control, so
# selecting on condition alone would silently pull IgG into the test group.
# It happens not to reach the matrix today only because the counts step already
# dropped it; that is luck, not a guarantee, so exclude it explicitly here too.
in_contrast <- meta$condition %in% c(opt$numerator, opt$denominator)
is_control  <- meta$target == opt$control_target
keep_samples <- intersect(meta$sample[in_contrast & !is_control], colnames(cts))

dropped <- setdiff(meta$sample[in_contrast & is_control], keep_samples)
if (length(dropped)) message(sprintf("[%s] excluded %s control librar%s: %s",
    opt$method, opt$control_target, if (length(dropped) == 1) "y" else "ies",
    paste(dropped, collapse = ", ")))

if (length(keep_samples) < 4) {
  stop(sprintf(paste0("contrast %s vs %s has %d usable samples (need >= 4). ",
                      "Check that the samplesheet conditions match the ",
                      "contrasts file and that control libraries are labelled ",
                      "target=%s."),
               opt$numerator, opt$denominator, length(keep_samples),
               opt$control_target))
}

cts   <- cts[, keep_samples, drop = FALSE]
group <- factor(meta$condition[match(keep_samples, meta$sample)],
                levels = c(opt$denominator, opt$numerator))
s     <- sf[keep_samples, opt$method]
stopifnot(all(is.finite(s)), all(s > 0))

# Filter on raw counts only. Filtering after normalization would let the
# method under test choose its own regions and invalidate the comparison.
keep <- rowSums(cts) >= 10 * ncol(cts)
cts  <- cts[keep, , drop = FALSE]
message(sprintf("[%s] %s vs %s: %d samples, %d regions after filtering",
                opt$method, opt$numerator, opt$denominator, ncol(cts), nrow(cts)))

## ---- edgeR ------------------------------------------------------------------

y <- DGEList(counts = cts, group = group)
y$samples$norm.factors <- edger_norm_factors(s, y$samples$lib.size)
design <- model.matrix(~group)
y  <- estimateDisp(y, design)
et <- glmQLFTest(glmQLFit(y, design), coef = 2)
res_edger <- topTags(et, n = Inf, sort.by = "none")$table

## ---- DESeq2 -----------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(cts, data.frame(group = group, row.names = colnames(cts)),
                              design = ~group)
sizeFactors(dds) <- deseq_size_factors(s)
dds <- DESeq(dds, quiet = TRUE)
res_deseq <- as.data.frame(results(dds, contrast = c("group", opt$numerator, opt$denominator)))

## ---- parameterization invariant --------------------------------------------
#
# Both engines were handed the same s. Their effective per-sample scaling must
# therefore agree up to a single global constant. Comparing RATIOS to the first
# sample removes that constant. If this fails, one of the two helpers above has
# been changed to the wrong unit and every downstream number is wrong.

eff_edger <- 1 / (y$samples$lib.size * y$samples$norm.factors)
eff_deseq <- 1 / sizeFactors(dds)
rel <- (eff_edger / eff_edger[1]) / (eff_deseq / eff_deseq[1])
if (max(abs(rel - 1)) > 0.02) {
  stop(sprintf(paste0("PARAMETERIZATION MISMATCH: edgeR and DESeq2 effective scaling ",
                      "differ by up to %.1f%% after removing the global constant. ",
                      "edgeR norm.factors correct lib.size; DESeq2 sizeFactors ARE ",
                      "the divisor. See the header of this file."),
               100 * max(abs(rel - 1))))
}
message(sprintf("[%s] parameterization check passed (max deviation %.3f%%)",
                opt$method, 100 * max(abs(rel - 1))))

## ---- emit -------------------------------------------------------------------

out <- data.frame(
  region      = rownames(cts),
  method      = opt$method,
  contrast    = sprintf("%s_vs_%s", opt$numerator, opt$denominator),
  edger_lfc   = res_edger$logFC,
  edger_fdr   = res_edger$FDR,
  deseq_lfc   = res_deseq$log2FoldChange,
  deseq_fdr   = res_deseq$padj,
  check.names = FALSE
)
write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("[%s] edgeR sig %d up / %d down | DESeq2 sig %d up / %d down",
                opt$method,
                sum(out$edger_fdr < opt$fdr & out$edger_lfc > 0, na.rm = TRUE),
                sum(out$edger_fdr < opt$fdr & out$edger_lfc < 0, na.rm = TRUE),
                sum(out$deseq_fdr < opt$fdr & out$deseq_lfc > 0, na.rm = TRUE),
                sum(out$deseq_fdr < opt$fdr & out$deseq_lfc < 0, na.rm = TRUE)))
