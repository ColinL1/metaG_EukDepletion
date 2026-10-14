#! /usr/bin/env Rscript
# ==============================================================================
# Pipeline to run all the scripts in order, producing the figures and tables for the paper.
# Each script runs in its own fresh R session (separate Rscript process), so
# nothing left over from one script can affect the next.
# Runnable from any location:  Rscript /path/to/R-scripts/run_all.r
# ==============================================================================

# ---- locate R-scripts/ (this file's folder) and work from there ----------------
script_dir <- local({
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg) > 0) {
        dirname(normalizePath(sub("^--file=", "", file_arg[1])))        # Rscript run_all.r
    } else if (!is.null(sys.frames()[[1]]$ofile)) {
        dirname(normalizePath(sys.frames()[[1]]$ofile))                 # source("run_all.r")
    } else {
        stop("Cannot locate run_all.r; run it with Rscript or source() it.", call. = FALSE)
    }
})
setwd(script_dir)

scripts <- c(
    "01_stats_anova.r",
    "02_figure_2_qpcr.r",
    "03_figure_3_reads.r",
    "04_figure_4_family_intersections.r",
    "05_figure_S2_aip_families.r",
    "tables/table_counts.r",
    "tables/table_mags.r"
    ## Optional: refit the main figures onto Letter pages for submission
    # , "fit_to_letter.r"
)

rscript <- file.path(R.home("bin"), "Rscript")
for (s in scripts) {
    message("\n==== ", s, " ====")
    status <- system2(rscript, s)
    if (status != 0) stop(s, " failed (exit status ", status, "); later scripts were not run.", call. = FALSE)
}
message("\nAll scripts finished.")
