#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args   <- commandArgs(trailingOnly = TRUE)
infile <- args[1]; sample <- args[2]; outfile <- args[3]

x <- fread(infile, sep = "\t", header = TRUE, quote = "",
           check.names = TRUE, na.strings = c("", "NA"))

# CAT writes the header as "# contig"
setnames(x, 1, "contig")
x[, contig := sub("^#\\s*", "", contig)]
x <- x[contig != "contig" & !is.na(contig)]

# CAT rank columns are plain names ("Vibrionaceae"), occasionally suffixed
# with "*" when the bit-score support at that rank is below CAT's threshold
# (see the "lineage scores" header). Some CAT builds instead embed the score
# directly ("Vibrionaceae: 0.87"); both are handled here.
get_name  <- function(v) sub(":.*$", "", trimws(v))
get_score <- function(v) suppressWarnings(as.numeric(sub("^.*:\\s*", "", v)))
# the "*" sits before the score when one is embedded ("Vibrionaceae*: 0.87")
is_unsupported <- function(v) grepl("\\*\\s*(:\\s*[0-9.]+)?$", trimws(v))

pick <- function(dt, rank) if (rank %in% names(dt)) dt[[rank]] else NA_character_

out <- data.table(
  contig        = paste0(sample, "__", x$contig),
  sample        = sample,
  classification= if ("classification" %in% names(x)) x$classification else NA_character_,
  family        = get_name(pick(x, "family")),
  family_score  = get_score(pick(x, "family")),
  order         = get_name(pick(x, "order")),
  class         = get_name(pick(x, "class")),
  phylum        = get_name(pick(x, "phylum")),
  superkingdom  = get_name(pick(x, "superkingdom"))
)

# normalise the non-assignments (and below-threshold, "*"-marked calls) so
# they survive as one explicit category rather than fragmenting real families
bad <- c("NA", "no support", "not classified", "no taxid assigned", "*", "")
unsupported <- is_unsupported(pick(x, "family"))
out[family %in% bad | is.na(family) | unsupported, `:=`(family = "Unassigned_family",
                                          family_score = NA_real_)]

fwrite(out, outfile, sep = "\t")
