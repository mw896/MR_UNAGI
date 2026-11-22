if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) invisible(Sys.setlocale("LC_CTYPE", ".UTF-8"))
library(data.table)
library(TwoSampleMR)
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 3L) stop("Usage: Rscript MR.R PROJECT_DIR PLINK_1.9 EUR_bfile_prefix")
endpoint <- "finngen_R12_C3_MYELOMA_MONOCLONAL"
ex <- file.path(a[1], "clean", "exposure")
hits <- as.data.table(read.delim(file.path(ex, "JAK1_JAK2_eqtlgen_p5e8_candidate_rows.tsv.gz")))
af <- as.data.table(read.delim(file.path(ex, "JAK1_JAK2_eqtlgen_af_for_clumped_snps.tsv.gz")))
out <- read.delim(file.path(a[1], "clean", "outcome", paste0(endpoint, "_for_JAK1_JAK2.tsv.gz")))
out <- format_data(out, type = "outcome", phenotype_col = "phenotype")

clump <- function(x) {
  input <- tempfile()
  fwrite(x[, .(SNP, P = Pvalue)], input, sep = "\t")
  cmd <- c("--bfile", shQuote(a[3]), "--clump", shQuote(input), "--clump-p1", "5e-8",
           "--clump-p2", "5e-8", "--clump-r2", "0.1", "--clump-kb", "1000", "--out", shQuote(input))
  status <- system2(a[2], cmd, stdout = FALSE)
  stopifnot(status == 0L, file.exists(paste0(input, ".clumped")))
  x[SNP %in% fread(paste0(input, ".clumped"), fill = TRUE)$SNP]
}

run_gene <- function(g) {
  x <- hits[GeneSymbol == g & grepl("^rs[0-9]+", SNP) & Pvalue > 0 & Pvalue <= 5e-8]
  x <- merge(clump(unique(x[order(Pvalue)], by = "SNP")), af, by = "SNP", all.x = TRUE)
  x[, eaf := fifelse(AssessedAllele == AlleleB, AlleleB_all,
                     fifelse(AssessedAllele == AlleleA, 1 - AlleleB_all, NA_real_))]
  x[, se := 1 / sqrt(2 * eaf * (1 - eaf) * (NrSamples + Zscore^2))]
  iv <- x[is.finite(se) & se > 0, .(SNP, beta = Zscore * se, se, eaf, pval = Pvalue,
            effect_allele = AssessedAllele, other_allele = OtherAllele, id = g, phenotype = g)]
  dat <- harmonise_data(format_data(as.data.frame(iv), type = "exposure", phenotype_col = "phenotype"), out, action = 2)
  f <- (dat$beta.exposure / dat$se.exposure)^2
  dat <- dat[dat$mr_keep & is.finite(f) & f > 10 , ]
  stopifnot(nrow(dat) > 0)
  method <- if (nrow(dat) == 1L) "mr_wald_ratio" else "mr_ivw"
  generate_odds_ratios(mr(dat, method_list = method))
}
res <- rbindlist(lapply(c("JAK1", "JAK2"), run_gene))
res[, q_within_endpoint := p.adjust(pval, "BH")]
print(res[, .(exposure, method, nsnp, or, or_lci95, or_uci95, pval, q_within_endpoint)])
