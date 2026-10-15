#! /usr/bin/env Rscript
# ==============================================================================
# Read-fraction tables: % of reads per mapping category (host, Symbiodiniaceae,
# bacteria, archaea, other) per extraction protocol, for Aiptasia (Fresh) and corals (DESS).
# Run from R-scripts/:  Rscript tables/table_counts.r
# Output: ../results/tables/table_counts_by_mapping.docx   (manuscript table: rows mapping x protocol, cols species - buffer)
#         ../results/tables/table_counts_by_species.docx   (rows mapping x species, cols protocol)
# ==============================================================================
source("00_setup.r")
p_load(flextable, officer, rempsyc)

# Protocol codes (BT, BenzoBT, MB, BBMB, SMB), Aiptasia buffer "Fresh", coral species
# in italics, two fixed decimals, "-" = not sequenced.
mapping_cols <- c("Host (%)", "Symbiodiniaceae (%)", "Bacteria (%)", "Archaea (%)", "Other (%)")
group_cols <- c("H2_Fresh", "F003_Fresh", "Acropora_DESS", "Pocillopora_DESS", "Porites_DESS")
fmt_pct <- function(x) ifelse(is.na(x), "-", formatC(x, format = "f", digits = 2))

# ---- data: % of total reads per category ----------------------------------------
reads <- load_read_stats()

pct <- reads %>%
    group_by(species, extraction, buffer, mapping) %>%
    summarise(n = sum(num_seqs), .groups = "drop") %>%
    pivot_wider(names_from = mapping, values_from = n) %>%
    transmute(Species = species, Buffer = recode(buffer, PBS = "Fresh"),
              Extraction = factor(extraction, levels = extraction_levels),
              `Host (%)`            = Host / total * 100,
              `Symbiodiniaceae (%)` = Symbiodiniaceae / total * 100,
              `Bacteria (%)`        = Bacteria / total * 100,
              `Archaea (%)`         = Archaea / total * 100,
              `Other (%)`           = Other / total * 100) %>%
    arrange(Species, Buffer, Extraction)

# ---- table A: rows mapping x protocol, columns species - buffer (manuscript table) -------
by_mapping <- pct %>%
    pivot_longer(all_of(mapping_cols), names_to = "Mapping", values_to = "Percentage") %>%
    pivot_wider(names_from = c(Species, Buffer), values_from = Percentage) %>%
    mutate(across(all_of(group_cols), fmt_pct),
           Mapping = factor(Mapping, levels = mapping_cols)) %>%
    arrange(Mapping, Extraction) %>%
    select(Mapping, Protocol = Extraction, all_of(group_cols))
names(by_mapping) <- gsub("_", " - ", names(by_mapping))
print(by_mapping, n = Inf)
ft <- style_table(by_mapping, group_cols = "Mapping", left_cols = c("Protocol", "Mapping"),
                  right_cols = gsub("_", " - ", group_cols))
save_table(italic_species(ft), "table_counts_by_mapping")

# ---- table B: rows mapping x species, columns protocol (Aiptasia Fresh only) --------------
by_species <- pct %>%
    filter(!(Species %in% c("H2", "F003") & Buffer == "DESS")) %>%
    pivot_longer(all_of(mapping_cols), names_to = "Mapping", values_to = "Percentage") %>%
    mutate(Mapping = factor(Mapping, levels = mapping_cols),
           Species = factor(Species, levels = c("H2", "F003", "Acropora", "Pocillopora", "Porites"))) %>%
    arrange(Mapping, Species, Buffer) %>%
    pivot_wider(names_from = Extraction, values_from = Percentage) %>%
    mutate(across(all_of(extraction_levels), fmt_pct)) %>%
    select(Mapping, Species, all_of(extraction_levels))
print(by_species, n = Inf)
ft <- style_table(by_species, group_cols = "Mapping", left_cols = c("Mapping", "Species"), right_cols = extraction_levels)
save_table(italic_species(ft, j = "Species"), "table_counts_by_species")
