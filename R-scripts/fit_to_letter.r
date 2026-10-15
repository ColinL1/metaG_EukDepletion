#! /usr/bin/env Rscript
# ==============================================================================
# Fit the main-figure PDFs onto a US Letter portrait page (copies only).
#
# Reads  ../plots/pdf/Figure_{2,3,4}.pdf (run the figure scripts first)
# Writes ../plots/pdf/letter/Figure_{2,3,4}.pdf
# Each figure is scaled uniformly (vector, no re-plotting) so it fits inside the
# page minus MARGIN_MM on every side, and centred.
# Needs the command-line tools `pdfinfo` (poppler) and `pdfjam` (TeX Live).
# ==============================================================================
source("00_setup.r")

# ---- PARAMS -------------------------------------------------------------------------
FIGURES   <- c("Figure_2", "Figure_3", "Figure_4")
PAGE_IN   <- c(width = 8.5, height = 11)    # US Letter portrait
MARGIN_MM <- 10

for (tool in c("pdfinfo", "pdfjam")) {
    if (Sys.which(tool) == "") stop("`", tool, "` not found on PATH", call. = FALSE)
}

in_dir  <- file.path(paths$plots, "pdf")
out_dir <- file.path(in_dir, "letter")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
check_inputs(file.path(in_dir, paste0(FIGURES, ".pdf")))

usable_pt <- PAGE_IN * 72 - 2 * MARGIN_MM / 25.4 * 72

for (fig in FIGURES) {
    src <- file.path(in_dir, paste0(fig, ".pdf"))
    dst <- file.path(out_dir, paste0(fig, ".pdf"))

    info <- system2("pdfinfo", shQuote(src), stdout = TRUE)
    size <- as.numeric(str_match(grep("^Page size:", info, value = TRUE), "([0-9.]+) x ([0-9.]+)")[, 2:3])
    scale <- min(usable_pt / size)

    status <- system2("pdfjam", c(shQuote(src), "--papersize", shQuote(sprintf("{%sin,%sin}", PAGE_IN[1], PAGE_IN[2])),
                                  "--noautoscale", "true", "--scale", format(scale, digits = 6),
                                  "--outfile", shQuote(dst)),
                      stdout = FALSE, stderr = FALSE)
    if (status != 0) stop("pdfjam failed for ", src, call. = FALSE)
    message(sprintf("wrote %s (%.0f x %.0f pt scaled by %.3f)", dst, size[1], size[2], scale))
}
