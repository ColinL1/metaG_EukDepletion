#! /usr/bin/env Rscript
# ==============================================================================
# Saturation curves (not in the paper): expected number of bacterial families
# vs. family-assigned bacterial reads, one curve per library (vegan::rarefy,
# exact). A plateau = richness saturated at that depth.
# Run from R-scripts/:  Rscript exploratory/saturation_curves.r
# Output: ../plots/extra/extra_01_saturation_curves.{png,pdf}
#         ../plots/extra/extra_01_saturation_endpoints.csv
# ==============================================================================
source("00_setup.r")

# ---- PARAMS -------------------------------------------------------------------------
MODES    <- c("persample", "pooled")   # quantification modes to plot (rows)
N_POINTS <- 40                         # depths evaluated per curve

grp_levels <- c("H2 - Fresh", "H2 - DESS", "F003 - Fresh", "F003 - DESS",
                "Acropora - DESS", "Porites - DESS", "Pocillopora - DESS")
quant_dir <- if (INPUT_SOURCE == "binac2") paths$cat_quant else paths$legacy_cat_quant

curve_for_mode <- function(mode) {
    contig_file <- file.path(quant_dir, mode, sprintf("contig_table.%s.tsv.gz", mode))
    check_inputs(contig_file)
    d <- read_tsv(contig_file, show_col_types = FALSE) %>%
        filter(superkingdom == "Bacteria", family != "Unassigned_family", !is.na(family), family != "")

    meta <- d %>%
        distinct(sample, species, treatment, buffer) %>%
        mutate(Trt      = factor(canon_extraction(treatment), levels = extraction_levels),
               Strain   = canon_species(species),
               Buffer   = recode(canon_buffer(buffer), PBS = "Fresh"),
               grouping = paste(Strain, Buffer, sep = " - "))

    counts <- d %>%
        group_by(sample, family) %>%
        summarise(n = sum(count), .groups = "drop") %>%
        pivot_wider(names_from = family, values_from = n, values_fill = 0) %>%
        column_to_rownames("sample") %>%
        as.matrix()

    map_dfr(rownames(counts), function(s) {
        x      <- counts[s, ]
        depths <- unique(round(c(seq(0, sum(x), length.out = N_POINTS), sum(x))))
        depths <- depths[depths > 0]
        tibble(sample = s, mode = mode, depth = depths, families = as.numeric(vegan::rarefy(x, depths)))
    }) %>%
        left_join(meta, by = "sample")
}

curves <- map_dfr(MODES, curve_for_mode) %>%
    filter(grouping %in% grp_levels) %>%
    mutate(grouping = factor(grouping, levels = grp_levels),
           mode = factor(mode, levels = MODES))

# curve endpoints = observed richness at full depth
ends <- curves %>% group_by(mode, sample) %>% slice_max(depth, n = 1) %>% ungroup()
dir.create(paths$plots_extra, recursive = TRUE, showWarnings = FALSE)
write_csv(ends %>% select(mode, sample, grouping, Trt, bacterial_reads = depth, families_observed = families),
          file.path(paths$plots_extra, "extra_01_saturation_endpoints.csv"))

p <- ggplot(curves, aes(depth / 1e6, families, group = sample, color = Trt)) +
    geom_line(linewidth = 0.7) +
    geom_point(data = ends, size = 1.6) +
    facet_grid(mode ~ grouping, scales = "free_x") +
    scale_color_manual(name = "Protocol", values = col_protocol) +
    labs(x = "Family-assigned bacterial reads (millions)", y = "Expected bacterial families") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1))

save_figure(p, "extra_01_saturation_curves", "extra", width = 16, height = 6, dpi = 300,
            formats = c("png", "pdf"), bg = c(png = "white", pdf = "white"))
