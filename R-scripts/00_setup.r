#! /usr/bin/env Rscript
# ==============================================================================
# 00_setup.r - shared configuration for every script in R-scripts/
#
# Sourced at the top of each figure / stats / table script:
#     source("00_setup.r")
# Scripts in subfolders are also run from R-scripts/ (e.g. Rscript tables/table_counts.r)
# All scripts must be run with the working directory set to R-scripts/
# (e.g. `cd R-scripts && Rscript 02_figure_2_qpcr.r`).
# ==============================================================================

# ---- input folders ------------------------------------------------------------
# Paths are relative to R-scripts/. The defaults match a fresh run from the
# repository root:
#   nextflow run main.nf ... --outdir results       -> results/
#   cd mag_recovery && ./run_cobinning.sh            -> mag_recovery/results_cobinning/
#   cd mag_recovery && ./run_coassembly.sh           -> mag_recovery/results_coassembly/
# If you used another --outdir or moved the folders, change these two lines.
PIPELINE_RESULTS <- file.path("..", "results")        # main.nf --outdir
MAG_RESULTS      <- file.path("..", "mag_recovery")   # folder holding results_cobinning/ and results_coassembly/

# ---- working directory check ----------------------------------------------------
if (basename(getwd()) != "R-scripts") {
    stop("Run from the R-scripts/ folder (current wd: ", getwd(), ")", call. = FALSE)
}
ROOT <- ".."
if (!interactive()) pdf(NULL)   # grob building opens a device; avoid a stray Rplots.pdf

invisible(capture.output(library(pacman)))
p_load(tidyverse, ggpubr)

# ==============================================================================
# Paths
# ==============================================================================
paths <- list(
    # main.nf outputs
    reads_stats   = file.path(PIPELINE_RESULTS, "reads_based", "stats", "stats_reads.tsv"),
    cat_quant     = file.path(PIPELINE_RESULTS, "coassembly_based", "quantification"),
    # nf-core/mag outputs (mag_recovery/)
    mags_co       = file.path(MAG_RESULTS, "results_coassembly", "GenomeBinning", "QC", "gunc_checkm_summary.tsv"),
    mags_non_co   = file.path(MAG_RESULTS, "results_cobinning",  "GenomeBinning", "QC", "gunc_checkm_summary.tsv"),
    # inputs not produced by a pipeline
    qpcr          = file.path(ROOT, "qPCR"),
    # outputs
    plots         = file.path(ROOT, "plots"),
    plots_supp    = file.path(ROOT, "plots", "supplementary"),
    plots_extra   = file.path(ROOT, "plots", "extra"),
    results       = file.path(ROOT, "results"),
    tables        = file.path(ROOT, "results", "tables")
)
paths$anova_csv <- file.path(paths$results, "anova_posthoc_contrasts.csv")
# stop early, naming every missing input
check_inputs <- function(...) {
    files <- c(...)
    missing <- files[!file.exists(files)]
    if (length(missing) > 0) {
        stop("Missing input file(s):\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
    }
    invisible(TRUE)
}

# ==============================================================================
# unify names per species, extraction protocol and buffer
# ==============================================================================
canon_species <- function(x) {
    recode(x,
        Ac = "Acropora", Acro = "Acropora", Ahemp = "Acropora",
        Po = "Pocillopora", Poci = "Pocillopora", Pver = "Pocillopora",
        Por = "Porites", Pr = "Porites", Pori = "Porites", Plob = "Porites",
        F3 = "F003")
}
canon_extraction <- function(x) {
    recode(x, B = "BenzoBT", Benzo = "BenzoBT", MBX = "BBMB", MBS = "SMB")
}
canon_buffer <- function(x) recode(x, P = "PBS", D = "DESS")

extraction_levels <- c("BT", "BenzoBT", "MB", "BBMB", "SMB")
coral_species <- c(Porites = "P. lobata", Acropora = "A. hemprichii", Pocillopora = "P. verrucosa")
protocol_full <- c(BT = "Blood & Tissue Kit", BenzoBT = "Benzonase + BT", MB = "Microbiome Kit",
                   BBMB = "Bead-beating + Microbiome Kit", SMB = "Spinning + Microbiome Kit")
grouping_levels <- c("H2 - PBS", "H2 - DESS", "F003 - PBS", "F003 - DESS",
                     "Acropora - DESS", "Porites - DESS", "Pocillopora - DESS")

# ==============================================================================
# Palettes
# ==============================================================================
col_protocol <- c(BT = "#FFB94F", BenzoBT = "#D97474", MB = "#9AE0EC", BBMB = "#87ACE9", SMB = "#9C80BA")
col_protocol_full <- setNames(col_protocol, protocol_full[names(col_protocol)])
col_qpcr <- c("Blood & Tissue Kit" = "#FFB94F", "Microbiome Kit" = "#9AE0EC", "Benzonase + BT" = "#D97474",
              "Sonication + Benzonase + BT" = "#F2A7C5", "PMA + BT" = "#82A37A", "Sonication + PMA + BT" = "#AFDD84")

# ==============================================================================
# Plot helpers
# ==============================================================================
# Facet labeller: coral genera -> italic species names. Used through
#   as_labeller(function(x) custom_labeller(x, ...), label_parsed)
#   style      "inline"  -> "A. hemprichii - DESS" on one line
#              "stacked" -> species over buffer (atop)
#   aip_raw    TRUE: return H2/F003 labels unchanged (inline only)
#   pbs_fresh  TRUE: show the PBS buffer as "Fresh"
custom_labeller <- function(labels, style = c("inline", "stacked"), aip_raw = FALSE, pbs_fresh = FALSE) {
    style <- match.arg(style)
    coral_name <- coral_species
    format_label <- function(label) {
        parts <- unlist(strsplit(as.character(label), " - ", fixed = TRUE))
        has_suffix <- length(parts) > 1
        if (pbs_fresh && has_suffix) parts[2] <- ifelse(parts[2] == "PBS", "Fresh", parts[2])
        aip <- parts[1] %in% c("H2", "F003")
        if (aip && aip_raw) return(label)
        if (style == "stacked") {
            if (aip) return(bquote(atop(plain(.(parts[1])), .(parts[2]))))
            sp <- if (parts[1] %in% names(coral_name)) coral_name[[parts[1]]] else parts[1]
            if (has_suffix) bquote(atop(italic(.(sp)), .(parts[2]))) else bquote(italic(.(sp)))
        } else if (parts[1] %in% names(coral_name)) {
            sp <- coral_name[[parts[1]]]
            if (has_suffix) bquote(italic(.(sp)) ~ "- " ~ .(parts[2])) else bquote(italic(.(sp)))
        } else if (aip) {
            if (has_suffix) bquote(.(parts[1]) ~ "- " ~ .(parts[2])) else bquote(.(parts[1]))
        } else {
            if (has_suffix) bquote(italic(.(parts[1])) ~ "- " ~ .(parts[2])) else bquote(italic(.(parts[1])))
        }
    }
    if (is.atomic(labels)) {
        labels <- sapply(labels, format_label)
    } else if (is.data.frame(labels)) {
        labels$value <- sapply(labels$value, format_label)
    }
    labels
}

extract_legend <- function(plot) {
    g <- ggplotGrob(plot)
    g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
}

# Save one figure as png/svg/pdf.
#   kind "main"  -> ../plots/{png,svg,pdf}/<name>.<ext>
#        "supp"  -> ../plots/supplementary/<name>.<ext>
#        "extra" -> ../plots/extra/<name>.<ext>
#   bg   background per format; NA keeps the theme's (transparent) background
save_figure <- function(plot, name, kind = c("main", "supp", "extra"), width, height,
                        formats = c("png", "svg", "pdf"),
                        bg = c(png = "white", svg = "white", pdf = "white"), dpi = 600, ...) {
    kind <- match.arg(kind)
    for (ext in formats) {
        dir <- switch(kind, main = file.path(paths$plots, ext), supp = paths$plots_supp, extra = paths$plots_extra)
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        args <- list(filename = file.path(dir, paste0(name, ".", ext)), plot = plot,
                     width = width, height = height, dpi = dpi, ...)
        if (!is.na(bg[[ext]])) args$bg <- bg[[ext]]
        do.call(ggsave, args)
        message("wrote ", args$filename)
    }
}

# PARAMS != DEFAULTS -> variant: goes to plots/extra/ with the changed params in the name
output_target <- function(name, params, defaults, kind) {
    changed <- names(defaults)[!mapply(identical, params[names(defaults)], defaults)]
    if (length(changed) == 0) return(list(name = name, kind = kind))
    tag <- paste0(changed, "-", vapply(params[changed], function(v) paste(format(v), collapse = "+"), ""),
                  collapse = "_")
    list(name = paste0("extra_", name, "__", gsub("[^A-Za-z0-9.+_-]", "", tag)), kind = "extra")
}

# ==============================================================================
# Table helpers (scripts in tables/; need flextable, officer, rempsyc loaded)
# ==============================================================================
# Group columns are printed only on the first row of each group, with a rule between groups.
style_table <- function(df, group_cols, left_cols, right_cols) {
    group <- do.call(paste, c(lapply(group_cols, function(g) as.character(df[[g]])), sep = "_"))
    first <- c(TRUE, group[-1] != group[-length(group)])
    ft <- nice_table(df, separate.header = TRUE) %>% align(align = "justify", part = "all")
    for (row in which(first)[-1]) {
        ft <- hline(ft, i = row - 1, border = fp_border(color = "black", width = 1))
    }
    formatters <- lapply(group_cols, function(g) {
        vals <- as.character(df[[g]])
        function(x) ifelse(first, vals, "")
    })
    ft <- do.call(set_formatter, c(list(ft), setNames(formatters, group_cols)))
    ft %>%
        align(j = left_cols, align = "left", part = "all") %>%
        align(j = right_cols, align = "right", part = "all") %>%
        autofit() %>%
        fontsize(size = 9, part = "all") %>%
        padding(padding.top = 3, padding.bottom = 3, part = "all") %>%
        line_spacing(space = 1.15, part = "all") %>%
        bold(part = "header") %>%
        width(j = left_cols, width = 1.2) %>%
        width(j = right_cols, width = 1.0)
}

# Coral species names in italics.
#   header labels "<genus> - <buffer>" -> "<italic species> - <buffer>"
#   body column j: genus values -> italic species names (Aiptasia strains stay upright)
# Column names cannot hold the species names: nice_table() splits headers on ".".
italic_species <- function(ft, j = NULL) {
    for (g in names(coral_species)) {
        for (col in grep(paste0("^", g, " - "), ft$col_keys, value = TRUE)) {
            ft <- compose(ft, part = "header", j = col,
                          value = as_paragraph(as_i(coral_species[[g]]), sub(g, "", col, fixed = TRUE)))
        }
    }
    if (!is.null(j)) {
        vals  <- as.character(ft$body$dataset[[j]])
        first <- c(TRUE, vals[-1] != vals[-length(vals)])   # style_table() prints a group label on its first row only
        for (i in which(first & vals %in% names(coral_species))) {
            g <- vals[i]
            ft <- compose(ft, i = i, j = j, part = "body", value = as_paragraph(as_i(coral_species[[g]])))
        }
    }
    bold(ft, part = "header")
}

save_table <- function(ft, name) {
    dir.create(paths$tables, recursive = TRUE, showWarnings = FALSE)
    file <- file.path(paths$tables, paste0(name, ".docx"))
    save_as_docx(ft, path = file)
    message("wrote ", file)
}

# ==============================================================================
# Loader: read-fraction table (host / symbiont / bacteria / ...)
# Returns one row per sample x mapping with columns
#   sample, species, replicate, extraction, buffer, grouping, mapping, num_seqs
#   mapping: Host, Symbiodiniaceae, Bacteria, Archaea, Other (viruses + non-mapped), total
# ==============================================================================
load_read_stats <- function() {
    .read_stats_pipeline() %>%
        mutate(species    = canon_species(species),
               extraction = canon_extraction(extraction),
               buffer     = canon_buffer(buffer),
               grouping   = paste(species, buffer, sep = " - "))
}

.read_stats_pipeline <- function() {
    check_inputs(paths$reads_stats)
    s <- read_tsv(paths$reads_stats, show_col_types = FALSE) %>%
        filter(str_detect(file, "(_R1|_1)\\.fq\\.gz$")) %>%
        mutate(category = case_when(
            str_detect(file, "_TRIM_paried_R1\\.fq\\.gz$") ~ "total",
            TRUE ~ str_match(file, "\\.([A-Za-z-]+)\\.(un)?mapped_1\\.fq\\.gz$")[, 2]))
    if (anyNA(s$category)) stop("unparsed file names in ", paths$reads_stats, ":\n  ",
                                paste(s$file[is.na(s$category)], collapse = "\n  "), call. = FALSE)
    s %>%
        mutate(mapping = recode(category, host = "Host", symbiodiniaceae = "Symbiodiniaceae",
                                bacteria = "Bacteria", archaea = "Archaea",
                                viruses = "Other", `non-mapped` = "Other")) %>%
        separate(sample, c("species", "replicate", "extraction", "buffer"), sep = "_", remove = FALSE) %>%
        group_by(sample, species, replicate, extraction, buffer, mapping) %>%
        summarise(num_seqs = sum(num_seqs), .groups = "drop")
}

# ==============================================================================
# Loader: CAT bacterial family table (Figure 4, Figure S2, exploratory CAT plots)
#
# 1. family x sample read counts from the per-contig coverM table
# 2. presence call: per-species drarefy() on TOTAL sequencing depth
#    (rarefy = TRUE), or any read at all (rarefy = FALSE)
# 3. dot/box size: full, non-rarefied, length-weighted coverage share per sample
# 4. top_n families per strain kept by name, the rest collapsed to "Other"
# Returns list(d, meta, fam, plot_dat)
# ==============================================================================
prepare_cat_families <- function(mode = "persample", rarefy = TRUE, presence_frac = 0.9,
                                 min_breadth = NULL, cov_method = "trimmed_mean", top_n = 35) {
    contig_file <- file.path(paths$cat_quant, mode, sprintf("contig_table.%s.tsv.gz", mode))
    check_inputs(contig_file)

    d <- read_tsv(contig_file, show_col_types = FALSE) %>%
        filter(superkingdom == "Bacteria", family != "Unassigned_family", !is.na(family), family != "")

    meta <- d %>%
        distinct(sample, species, treatment, buffer) %>%
        mutate(
            Treatment = unname(protocol_full[canon_extraction(treatment)]),
            Strain    = factor(canon_species(species), levels = c("F003", "H2", "Acropora", "Porites", "Pocillopora")),
            Buffer    = factor(recode(canon_buffer(buffer), PBS = "Fresh"), levels = c("Fresh", "DESS")),
            grouping  = paste(Strain, Buffer, sep = " - ")
        )

    # step 1 ---------------------------------------------------------------------
    counts <- d %>%
        group_by(sample, family) %>%
        summarise(n = sum(count), .groups = "drop") %>%
        pivot_wider(names_from = family, values_from = n, values_fill = 0) %>%
        column_to_rownames("sample") %>%
        as.matrix()
    species_of <- meta %>% distinct(sample, species) %>% deframe()

    # step 2 ---------------------------------------------------------------------
    if (rarefy) {
        total_reads <- cat_total_depth(meta)
        # "non-bacterial" bucket so each row sums to the library's TOTAL depth
        non_bacterial <- total_reads[rownames(counts)] - rowSums(counts)
        stopifnot(all(non_bacterial > 0))
        counts_ext <- cbind(counts, `.NonBacterial` = non_bacterial)

        detected <- matrix(NA_real_, nrow(counts), ncol(counts), dimnames = dimnames(counts))
        for (sp in unique(species_of)) {
            samples_sp    <- names(species_of)[species_of == sp]
            counts_ext_sp <- counts_ext[samples_sp, , drop = FALSE]
            prob_sp       <- vegan::drarefy(counts_ext_sp, min(rowSums(counts_ext_sp)))
            detected[samples_sp, ] <- prob_sp[, colnames(counts), drop = FALSE]
        }
        present <- detected >= presence_frac
        message("rarefaction depth per species (total reads, min = target):")
        for (sp in sort(unique(species_of))) {
            message("  ", sp, ": ", paste(sort(total_reads[names(species_of)[species_of == sp]]), collapse = " "))
        }
    } else {
        present <- counts > 0
    }

    if (!is.null(min_breadth)) {
        breadth <- d %>%
            group_by(sample, family) %>%
            summarise(cf = weighted.mean(covered_fraction, w = length), .groups = "drop") %>%
            pivot_wider(names_from = family, values_from = cf, values_fill = 0) %>%
            column_to_rownames("sample") %>% as.matrix()
        present <- present & (breadth[rownames(present), colnames(present)] >= min_breadth)
    }

    # step 3 ---------------------------------------------------------------------
    size <- d %>%
        group_by(sample, family) %>%
        summarise(cov = weighted.mean(.data[[cov_method]], w = length), .groups = "drop") %>%
        group_by(sample) %>%
        mutate(rel = cov / sum(cov)) %>%
        ungroup()

    present_long <- as_tibble(present, rownames = "sample") %>%
        pivot_longer(-sample, names_to = "family", values_to = "present")

    chk <- size %>% left_join(present_long, by = c("sample", "family")) %>% filter(rel > 0)
    message(sprintf("presence filter: %d / %d family-sample cells with coverage dropped",
                    sum(!chk$present, na.rm = TRUE), nrow(chk)))

    # step 4 ---------------------------------------------------------------------
    fam <- size %>%
        left_join(present_long, by = c("sample", "family")) %>%
        filter(present) %>%
        left_join(meta, by = "sample") %>%
        rename(Family = family)

    top_families <- fam %>%
        group_by(Strain, Family) %>%
        summarise(m = mean(rel), .groups = "drop") %>%
        arrange(Strain, desc(m)) %>%
        group_by(Strain) %>% slice_head(n = top_n) %>%
        ungroup() %>% select(Strain, Family)

    plot_dat <- fam %>%
        left_join(mutate(top_families, is_top = TRUE), by = c("Strain", "Family")) %>%
        mutate(Family = ifelse(!is.na(is_top), Family, "Other")) %>%
        select(-is_top) %>%
        group_by(Strain, Treatment, Buffer, grouping, sample, Family) %>%
        summarise(rel = sum(rel), .groups = "drop")
    family_order <- sort(setdiff(unique(plot_dat$Family), "Other"))
    plot_dat$Family <- factor(plot_dat$Family, levels = rev(c(family_order, "Other")))

    list(d = d, meta = meta, fam = fam, plot_dat = plot_dat)
}

# Total sequencing depth (R1 + R2 reads) per CAT sample, named by CAT sample id
cat_total_depth <- function(meta) {
    key <- function(sp, trt, buf) paste(canon_species(sp), canon_extraction(trt), canon_buffer(buf), sep = "_")

    check_inputs(paths$reads_stats)
    per_group <- read_tsv(paths$reads_stats, show_col_types = FALSE) %>%
        filter(str_detect(file, "_TRIM_paried_R1\\.fq\\.gz$")) %>%
        separate(sample, c("species", "replicate", "extraction", "buffer"), sep = "_") %>%
        mutate(key = key(species, extraction, buffer)) %>%
        group_by(key) %>%
        summarise(total = sum(num_seqs) * 2, .groups = "drop") %>%
        deframe()

    samples <- meta %>% distinct(sample, species, treatment, buffer) %>%
        mutate(key = key(species, treatment, buffer))
    unmatched <- samples$sample[!samples$key %in% names(per_group)]
    if (length(unmatched) > 0) {
        stop("No total sequencing depth found for CAT sample(s): ", paste(unmatched, collapse = ", "),
             "\n  Something is inconsistent between the CAT quantification",
             " and the read stats - check before going further.", call. = FALSE)
    }
    setNames(per_group[samples$key], samples$sample)
}
