#! /usr/bin/env Rscript
# ==============================================================================
# Read-fraction tables: % of reads per mapping category (host, Symbiodiniaceae,
# bacteria, archaea, other) for every species x buffer x extraction protocol.
# Run from R-scripts/:  Rscript tables/table_counts.r
# Output: ../results/tables/table_counts_by_mapping.docx   (rows mapping x protocol,   cols species - buffer)
#         ../results/tables/table_counts_by_species.docx   (rows mapping x species,    cols protocol; Aiptasia PBS only)
# ==============================================================================
source("00_setup.r")
p_load(flextable, officer, rempsyc)

mapping_cols <- c("Host (%)", "Symbiodiniaceae (%)", "Bacteria (%)", "Archaea (%)", "Other (%)")
extraction_long <- c(BT = "Blood & Tissue", BenzoBT = "Benzonase", MB = "Microbiome",
                     BBMB = "Microbiome & Bead-Beating", SMB = "Microbiome & Spinning")
group_cols <- c("H2_PBS", "H2_DESS", "F003_PBS", "F003_DESS", "Acropora_DESS", "Pocillopora_DESS", "Porites_DESS")

# ---- data: % of total reads per category ----------------------------------------
reads <- load_read_stats(kaiju = FALSE)
if (INPUT_SOURCE == "legacy") {
    # legacy input: original (pre-Kaiju-rerun) counts, as table_counts.r always read them;
    # legacy categories did not cover archaea / viruses / unplaced reads: the remainder is reported as "Archaea"
    reads <- reads %>%
        pivot_wider(names_from = mapping, values_from = num_seqs) %>%
        mutate(Archaea = total - (Host + Symbiodiniaceae + Bacteria + Other)) %>%
        pivot_longer(c(total, Host, Symbiodiniaceae, Bacteria, Archaea, Other), names_to = "mapping", values_to = "num_seqs")
}

pct <- reads %>%
    group_by(species, extraction, buffer, mapping) %>%
    summarise(n = sum(num_seqs), .groups = "drop") %>%
    pivot_wider(names_from = mapping, values_from = n) %>%
    transmute(Species = species, Buffer = buffer,
              Extraction = factor(extraction_long[extraction], levels = extraction_long),
              `Host (%)`            = Host / total * 100,
              `Symbiodiniaceae (%)` = Symbiodiniaceae / total * 100,
              `Bacteria (%)`        = Bacteria / total * 100,
              `Archaea (%)`         = Archaea / total * 100,
              `Other (%)`           = Other / total * 100) %>%
    arrange(Species, Buffer, Extraction)

# ---- table A: rows mapping x protocol, columns species - buffer -------------------------
by_mapping <- pct %>%
    pivot_longer(all_of(mapping_cols), names_to = "Mapping", values_to = "Percentage") %>%
    pivot_wider(names_from = c(Species, Buffer), values_from = Percentage) %>%
    mutate(across(all_of(group_cols), ~ ifelse(is.na(.), "-", round(., 2))),
           Mapping = factor(Mapping, levels = mapping_cols)) %>%
    arrange(Mapping, Extraction) %>%
    select(Mapping, Extraction, all_of(group_cols))
names(by_mapping) <- gsub("_", " - ", names(by_mapping))
print(by_mapping, n = Inf)
save_table(style_table(by_mapping, group_cols = "Mapping", left_cols = c("Extraction", "Mapping"),
                       right_cols = gsub("_", " - ", group_cols)), "table_counts_by_mapping")

# ---- table B: rows mapping x species, columns protocol (Aiptasia PBS only) ----------------
by_species <- pct %>%
    filter(!(Species %in% c("H2", "F003") & Buffer == "DESS")) %>%
    pivot_longer(all_of(mapping_cols), names_to = "Mapping", values_to = "Percentage") %>%
    mutate(Mapping = factor(Mapping, levels = mapping_cols),
           Species = factor(Species, levels = c("H2", "F003", "Acropora", "Pocillopora", "Porites"))) %>%
    arrange(Mapping, Species, Buffer) %>%
    pivot_wider(names_from = Extraction, values_from = Percentage) %>%
    mutate(across(all_of(unname(extraction_long)), ~ ifelse(is.na(.), "-", round(., 2)))) %>%
    select(Mapping, Species, all_of(unname(extraction_long)))
print(by_species, n = Inf)
save_table(style_table(by_species, group_cols = "Mapping", left_cols = "Mapping", right_cols = unname(extraction_long)),
           "table_counts_by_species")
