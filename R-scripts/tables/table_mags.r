#! /usr/bin/env Rscript
# ==============================================================================
# MAG tables: number of medium (>=50% complete, <=10% contamination) and high
# (>=90% complete, <=5% contamination; MIMAG) quality MAGs (CheckM) per
# species x buffer x extraction protocol.
# Run from R-scripts/:  Rscript tables/table_mags.r
# Output: ../results/tables/table_mags_cobinning.docx   (manuscript table: Aiptasia Fresh and corals)
#         ../results/tables/table_mags_full.docx         (co-binning and co-assembly side by side)
# Protocol codes, Aiptasia buffer "Fresh", coral species in italics,
# "-" = sequenced but no MAG, "NA" = not sequenced.
# ==============================================================================
source("00_setup.r")
p_load(flextable, officer, rempsyc)

q_cols <- c("nMAGs", "nHQ-MAGs")
species_levels <- c("H2", "F003", "Acropora", "Pocillopora", "Porites")

# ---- data ---------------------------------------------------------------------
check_inputs(paths$mags_co, paths$mags_non_co)

# genome ids: <assembler>-<binner>-<strain>_<replicate>_<treatment>_<buffer>.<bin>           (co-binning)
#             <assembler>-<binner>-<group>-<strain>_<treatment>_<buffer>.<bin>              (co-assembly)
mags <- bind_rows(
    read_tsv(paths$mags_non_co, show_col_types = FALSE) %>%
        separate(genome, into = c("assembler", "binner", "temp"), sep = "-", remove = FALSE) %>%
        separate(temp, into = c("strain", "replicate", "treatment", "buffer_bin"), sep = "_") %>%
        mutate(origin = "Co-binning"),
    read_tsv(paths$mags_co, show_col_types = FALSE) %>%
        separate(genome, into = c("assembler", "binner", "group", "temp"), sep = "-", remove = FALSE) %>%
        separate(temp, into = c("strain", "treatment", "buffer_bin"), sep = "_") %>%
        mutate(origin = "Co-assembly")
) %>%
    separate(buffer_bin, into = c("buffer", "bin_id"), sep = "\\.") %>%
    mutate(strain    = canon_species(strain),
           treatment = canon_extraction(treatment),
           buffer    = recode(canon_buffer(buffer), PBS = "Fresh")) %>%
    filter(checkM.completeness >= 50, checkM.contamination <= 10) %>%
    mutate(completeness_cat = ifelse(checkM.completeness >= 90 & checkM.contamination <= 5, q_cols[2], q_cols[1]))

# counts per group; nMAGs includes the HQ MAGs
counts <- mags %>%
    count(strain, treatment, buffer, origin, completeness_cat) %>%
    pivot_wider(names_from = completeness_cat, values_from = n, values_fill = 0) %>%
    mutate(!!q_cols[1] := .data[[q_cols[1]]] + .data[[q_cols[2]]]) %>%
    pivot_longer(all_of(q_cols), names_to = "quality", values_to = "n") %>%
    transmute(Species = strain, Buffer = buffer, Extraction = treatment,
              col = paste(origin, quality, sep = "."), n = ifelse(n == 0, "-", as.character(n))) %>%
    pivot_wider(names_from = col, values_from = n)

# every protocol of every sequenced species x buffer; never-sequenced combinations are "NA"
not_sequenced <- tribble(
    ~Species, ~Buffer, ~Extraction,
    "H2",     "DESS",  "SMB",
    "F003",   "DESS",  "BenzoBT",
    "F003",   "DESS",  "MB",
    "F003",   "DESS",  "SMB")
value_cols <- c("Co-binning.nMAGs", "Co-binning.nHQ-MAGs", "Co-assembly.nMAGs", "Co-assembly.nHQ-MAGs")

full <- expand_grid(Species = species_levels, Buffer = c("Fresh", "DESS"), Extraction = extraction_levels) %>%
    filter(!(Species %in% names(coral_species) & Buffer == "Fresh")) %>%
    left_join(counts, by = c("Species", "Buffer", "Extraction")) %>%
    mutate(across(all_of(value_cols), ~ ifelse(is.na(.), "-", .))) %>%
    rows_update(not_sequenced %>% mutate(!!!setNames(rep(list("NA"), length(value_cols)), value_cols)),
                by = c("Species", "Buffer", "Extraction"), unmatched = "ignore") %>%
    mutate(Species = factor(Species, levels = species_levels),
           Buffer = factor(Buffer, levels = c("Fresh", "DESS")),
           Extraction = factor(Extraction, levels = extraction_levels)) %>%
    arrange(Species, Buffer, Extraction) %>%
    select(Species, Buffer, Extraction, all_of(value_cols))

# ---- table 1: co-binning, Aiptasia Fresh and corals ----------------------------------
table1 <- full %>%
    filter(!(Species %in% c("H2", "F003") & Buffer == "DESS")) %>%
    select(Species, Buffer, Extraction, nMAGs = `Co-binning.nMAGs`, `nHQ-MAGs` = `Co-binning.nHQ-MAGs`)
print(table1, n = Inf)
ft <- style_table(table1, group_cols = c("Species", "Buffer"), left_cols = c("Species", "Buffer", "Extraction"),
                  right_cols = q_cols)
save_table(italic_species(ft, j = "Species"), "table_mags_cobinning")

# ---- table 2: co-binning and co-assembly, every species x buffer ----------------------
table2 <- full
print(table2, n = Inf)
ft <- style_table(table2, group_cols = c("Species", "Buffer"), left_cols = c("Species", "Buffer", "Extraction"),
                  right_cols = value_cols) %>%
    align(i = 1, j = value_cols, align = "center", part = "header")   # "Co-binning" / "Co-assembly" spans
save_table(italic_species(ft, j = "Species"), "table_mags_full")
