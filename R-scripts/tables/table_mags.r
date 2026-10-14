#! /usr/bin/env Rscript
# ==============================================================================
# MAG tables: number of medium (>=50% complete) and high (>=90%) quality MAGs
# (CheckM contamination <= 10%) per species x buffer x extraction protocol.
# Run from R-scripts/:  Rscript tables/table_mags.r
# Output: ../results/tables/table_mags_cobinning.docx
#         ../results/tables/table_mags_full.docx      (co-binning and co-assembly side by side)
# ==============================================================================
source("00_setup.r")
p_load(flextable, officer, rempsyc)

extraction_long <- c(BT = "Blood & Tissue", BenzoBT = "Benzonase", MB = "Microbiome",
                     BBMB = "Microbiome & Bead-Beating", SMB = "Microbiome & Spinning")
q_cols <- c("≥50", "≥90")

# ---- data ---------------------------------------------------------------------
check_inputs(paths$mags_co, paths$mags_non_co)

# genome ids: <assembler>-<binner>-<strain>_<replicate>_<treatment>_<buffer>.<bin>           (co-binning)
#             <assembler>-<binner>-<group>-<strain>_<treatment>_<buffer>.<bin>              (co-assembly)
mags <- bind_rows(
    read_tsv(paths$mags_non_co, show_col_types = FALSE) %>%
        separate(genome, into = c("assembler", "binner", "temp"), sep = "-", remove = FALSE) %>%
        separate(temp, into = c("strain", "replicate", "treatment", "buffer_bin"), sep = "_") %>%
        mutate(origin = "co-binning"),
    read_tsv(paths$mags_co, show_col_types = FALSE) %>%
        separate(genome, into = c("assembler", "binner", "group", "temp"), sep = "-", remove = FALSE) %>%
        separate(temp, into = c("strain", "treatment", "buffer_bin"), sep = "_") %>%
        mutate(origin = "co-assembly")
) %>%
    separate(buffer_bin, into = c("buffer", "bin_id"), sep = "\\.") %>%
    mutate(strain    = canon_species(strain),
           treatment = extraction_long[canon_extraction(treatment)],
           buffer    = canon_buffer(buffer)) %>%
    filter(checkM.completeness >= 50, checkM.contamination <= 10) %>%
    mutate(completeness_cat = ifelse(checkM.completeness >= 90, q_cols[2], q_cols[1]))

# counts per group; ">=50" includes the ">=90" MAGs; 0 shown as "-"
summary_mags <- mags %>%
    count(strain, treatment, buffer, origin, completeness_cat) %>%
    pivot_wider(names_from = completeness_cat, values_from = n, values_fill = 0) %>%
    mutate(!!q_cols[1] := .data[[q_cols[1]]] + .data[[q_cols[2]]],
           across(all_of(q_cols), ~ ifelse(. == 0, "-", as.character(.)))) %>%
    transmute(Species = factor(strain, levels = c("H2", "F003", "Acropora", "Pocillopora", "Porites")),
              Buffer = factor(buffer, levels = c("PBS", "DESS")),
              Extraction = factor(treatment, levels = extraction_long),
              !!q_cols[1] := .data[[q_cols[1]]], !!q_cols[2] := .data[[q_cols[2]]], Origin = origin) %>%
    arrange(Species, Buffer, Extraction)

co_assembly_data <- summary_mags %>% filter(Origin == "co-assembly") %>% select(-Origin)
co_binning_data  <- summary_mags %>% filter(Origin == "co-binning") %>% select(-Origin)

# ---- table 1: co-binning ------------------------------------------------------------
table1 <- co_binning_data %>% mutate(across(everything(), ~ ifelse(is.na(.), "-", as.character(.))))
print(table1, n = Inf)
save_table(style_table(table1, group_cols = c("Species", "Buffer"), left_cols = c("Species", "Buffer", "Extraction"),
                       right_cols = q_cols), "table_mags_cobinning")

# ---- table 2: co-binning and co-assembly, every protocol listed -------------------------
# combinations that were never sequenced are shown as "NA" (vs "-" = sequenced, no MAG)
not_sequenced <- tribble(
    ~Species, ~Buffer, ~Extraction,
    "H2",     "DESS",  "Microbiome & Spinning",
    "F003",   "DESS",  "Benzonase",
    "F003",   "DESS",  "Microbiome",
    "F003",   "DESS",  "Microbiome & Spinning")

full <- full_join(co_assembly_data, co_binning_data, by = c("Species", "Buffer", "Extraction"),
                  suffix = c(".Co-assembly", ".Co-binning"))
table2 <- expand_grid(Species = levels(summary_mags$Species), Buffer = levels(summary_mags$Buffer),
                      Extraction = unname(extraction_long)) %>%
    filter(!(Species %in% c("Acropora", "Porites", "Pocillopora") & Buffer == "PBS")) %>%
    left_join(full %>% mutate(across(c(Species, Buffer, Extraction), as.character)),
              by = c("Species", "Buffer", "Extraction")) %>%
    mutate(across(ends_with(c(".Co-assembly", ".Co-binning")), ~ ifelse(is.na(.), "-", .))) %>%
    rows_update(not_sequenced %>% mutate(`≥50.Co-assembly` = "NA", `≥90.Co-assembly` = "NA",
                                         `≥50.Co-binning` = "NA", `≥90.Co-binning` = "NA"),
                by = c("Species", "Buffer", "Extraction"), unmatched = "ignore") %>%
    mutate(Species = factor(Species, levels = levels(summary_mags$Species)),
           Buffer = factor(Buffer, levels = levels(summary_mags$Buffer)),
           Extraction = factor(Extraction, levels = extraction_long)) %>%
    arrange(Species, Buffer, Extraction)
names(table2) <- gsub("^(≥[0-9]+)\\.(Co-assembly|Co-binning)$", "\\2.\\1", names(table2))
table2 <- table2 %>% select(Species, Buffer, Extraction, `Co-assembly.≥50`, `Co-assembly.≥90`, `Co-binning.≥50`, `Co-binning.≥90`)
print(table2, n = Inf)
save_table(style_table(table2, group_cols = c("Species", "Buffer"), left_cols = c("Species", "Buffer", "Extraction"),
                       right_cols = c("Co-assembly.≥50", "Co-assembly.≥90", "Co-binning.≥50", "Co-binning.≥90")),
           "table_mags_full")
