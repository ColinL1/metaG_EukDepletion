#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]; sheet <- args[2]
cov_files <- strsplit(args[3], "\\s+")[[1]]
cat_files <- strsplit(args[4], "\\s+")[[1]]

meta <- fread(sheet, select = c("sample","species","treatment","buffer"))
cat_tab <- rbindlist(lapply(cat_files, fread), fill = TRUE)
setkey(cat_tab, contig)

read_cov <- function(f) {
  s <- sub("\\..*$", "", basename(f))
  d <- fread(f, sep = "\t", header = TRUE)
  setnames(d, 1, "contig")
  # CoverM prefixes every metric column with the bam name; match on the suffix
  pick <- function(pat) { i <- grep(pat, names(d), ignore.case = TRUE)
                          if (length(i)) d[[i[1]]] else NA_real_ }
  data.table(sample = s, contig = d$contig,
             length          = pick("Length$"),
             count           = pick("Read Count$"),
             mean            = pick("Mean$"),
             trimmed_mean    = pick("Trimmed Mean$"),
             covered_fraction= pick("Covered Fraction$"),
             tpm             = pick("TPM$"))
}

cov <- rbindlist(lapply(cov_files, read_cov))
cov <- cat_tab[cov, on = "contig"]              # attach taxonomy
cov[, sample := i.sample]; cov[, i.sample := NULL]
cov <- meta[cov, on = "sample"]
cov[, mode := mode]

fwrite(cov, sprintf("contig_table.%s.tsv.gz", mode), sep = "\t")

# per-sample mapping effort, for the rarefaction denominator
summ <- cov[, .(n_contigs_ref   = .N,
                n_contigs_hit   = sum(count > 0, na.rm = TRUE),
                bacterial_reads = sum(count, na.rm = TRUE)),
            by = .(mode, sample, species, treatment, buffer)]
fwrite(summ, sprintf("mapping_summary.%s.tsv", mode), sep = "\t")

fam <- cov[, .(reads        = sum(count, na.rm = TRUE),
               tpm          = sum(tpm, na.rm = TRUE),
               n_contigs    = .N,
               n_contigs_det= sum(count > 0, na.rm = TRUE),
               max_cov_frac = max(covered_fraction, na.rm = TRUE),
               bp           = sum(length, na.rm = TRUE)),
           by = .(mode, sample, species, treatment, buffer, family)]
fwrite(fam, sprintf("family_table.%s.tsv", mode), sep = "\t")
