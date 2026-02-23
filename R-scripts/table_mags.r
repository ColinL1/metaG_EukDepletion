#! /usr/bin/env Rscript
cat("Loading libraries...")
library(pacman)
p_load(rempsyc, flextable, officer, broom, report, effectsize, tidyverse)

# --- Clear workspace ---
# rm(list = ls())

# %% Read the TSV files
cat("Reading MAG summary files...")
folder_path_coass <- "/Users/luigicolin/PhD/work_local/VSC_local/metaG_EukDepletion/input/gunc_checkm_summary_co.tsv"
folder_path <- "/Users/luigicolin/PhD/work_local/VSC_local/metaG_EukDepletion/input/gunc_checkm_summary_non_co.tsv"

mags_coass <- read_tsv(folder_path_coass)
mags <- read_tsv(folder_path)

# %% Parse co-binning data
mags_parsed <- mags %>%
    separate(genome, into = c("assembler", "binner", "temp"), sep = "-", remove = FALSE) %>%
    separate(temp, into = c("strain", "replicate", "treatment", "buffer_bin"), sep = "_", remove = FALSE) %>%
    separate(buffer_bin, into = c("buffer", "bin_id"), sep = "\\.", remove = FALSE) %>%
    mutate(origin = "co-binning")

# %% Parse co-assembly data  
mags_parsed_coass <- mags_coass %>%
    separate(genome, into = c("assembler", "binner", "group", "temp"), sep = "-", remove = FALSE) %>%
    separate(temp, into = c("strain", "treatment", "buffer_bin"), sep = "_", remove = FALSE) %>%
    separate(buffer_bin, into = c("buffer", "bin_id"), sep = "\\.", remove = FALSE) %>%
    mutate(origin = "co-assembly", replicate = NA)

# %% Combine datasets
combined_mags <- bind_rows(mags_parsed, mags_parsed_coass) %>%
    # Clean up treatment names
    mutate(
        treatment = case_when(
            treatment == "B" ~ "Benzonase",
            treatment == "BT" ~ "Blood & Tissue",
            treatment == "MB" ~ "Microbiome",
            treatment == "MBX" ~ "Microbiome & Bead-Beating",
            treatment == "MBS" ~ "Microbiome & Spinning",
            TRUE ~ treatment
        ),
        # Clean up strain names
        strain = case_when(
            strain == "Acro" ~ "Acropora",
            strain == "Por" ~ "Porites", 
            strain == "Poci" ~ "Pocillopora",
            TRUE ~ strain
        ),
        # Add buffer cleanup if needed
        buffer = case_when(
            buffer == "D" ~ "DESS",
            buffer == "P" ~ "PBS",
            TRUE ~ buffer
        )
    ) %>%
    # Filter for quality MAGs
    filter(checkM.completeness >= 50, checkM.contamination <= 10) %>%
    # Add completeness categories
    mutate(
        completeness_cat = ifelse(checkM.completeness >= 90, "≥90", "≥50"),
        type = ifelse(strain %in% c("H2", "F003"), "Aiptasia", "Coral")
    )

# %% Create summary data
combined_mags_summary <- combined_mags %>%
    group_by(strain, treatment, buffer, origin, completeness_cat) %>%
    summarise(bins_count = n(), .groups = "drop") %>%
    pivot_wider(
        names_from = completeness_cat,
        values_from = bins_count,
        values_fill = 0
    ) %>%
    # Sum ≥90 into ≥50 category
    mutate(`≥50` = `≥50` + `≥90`) %>%
    arrange(strain, buffer, treatment)

# %% Prepare data for tables
# Ensure all factor levels are consistent with table_counts.r
combined_mags_summary <- combined_mags_summary %>%
    mutate(
        strain = factor(strain, levels = rev(c("Porites", "Pocillopora", "Acropora", "F003", "H2"))),
        buffer = factor(buffer, levels = rev(c("DESS", "PBS"))),
        treatment = factor(treatment, levels = c("Blood & Tissue", "Benzonase", "Microbiome", "Microbiome & Bead-Beating", "Microbiome & Spinning"))
    ) %>%
    arrange(strain, buffer, treatment)

# Filter and prepare final dataset
final_mags_summary <- combined_mags_summary %>%
    # filter(strain != "F003") %>%  # Remove F003
    # filter(!(strain %in% c("Acropora", "Porites", "Pocillopora") & buffer == "PBS")) %>%  # Remove coral + PBS combinations
    # Replace 0s with "-" for better presentation
    mutate(
        `≥50` = ifelse(`≥50` == 0, "-", as.character(`≥50`)),
        `≥90` = ifelse(`≥90` == 0, "-", as.character(`≥90`))
    ) %>%
    select(strain, buffer, treatment, `≥50`, `≥90`, origin) %>%
    arrange(strain, buffer, treatment)

# Rename columns to match table_counts.r style
names(final_mags_summary) <- c("Species", "Buffer", "Extraction", "≥50", "≥90", "Origin")

cat("Data processing complete. Creating tables...")


# # # Summarize number of bins per strain, treatment, and buffer
# summary_mags <- mags_parsed %>%
#     filter(checkM.completeness >= 50) %>%
#     filter(checkM.contamination <= 10) %>%
#     group_by(strain, treatment, buffer) %>%
#     summarize(
#         bins_count = n(),
#         lineage_values = paste(checkM.lineage, collapse = "; "),
#         completeness_values = paste(checkM.completeness, collapse = "; "),
#         .groups = "drop"
#     )

# # Parse the "genome" column . # for coassembled 
# mags_parsed_coass <- mags_coass %>%
#     separate(genome, into = c("assembler", "binner", "group", "temp"), sep = "-", remove = FALSE) %>%
#     separate(temp, into = c("strain", "treatment", "buffer_bin"), sep = "_", remove = FALSE) %>%
#     separate(buffer_bin, into = c("buffer", "bin_id"), sep = "\\.", remove = FALSE)

# # Summarize number of bins per strain, treatment, and buffer
# summary_mags_coass <- mags_parsed_coass %>%
#   filter(checkM.completeness >= 50) %>%
#   filter(checkM.contamination <= 10) %>%
#   group_by(strain, treatment, buffer) %>%
#   summarize(
#     bins_count = n(),
#     lineage_values = paste(checkM.lineage, collapse = "; "),
#     completeness_values = paste(checkM.completeness, collapse = "; "),
#     .groups = "drop"
#   )

# summary_mags %>%
#     arrange(buffer, strain, treatment)

# summary_mags_coass %>%
#     arrange(buffer, strain, treatment)
# # Write the summary to a TSV file
# # write_tsv(summary_mags, paste0(folder_path,"gunc_checkm_summary_summary.tsv"))
# # write_tsv(summary_mags_coass, paste0(folder_path_coass,"gunc_checkm_summary_summary.tsv"))

# #add origin info column  
# summary_mags$origin <- "co-binning"
# mags_parsed$origin <- "co-binning"
# summary_mags_coass$origin <- "co-assembly"
# mags_parsed_coass$origin <- "co-assembly"

# # Combine the two summaries
# combined_summary <- bind_rows(summary_mags, summary_mags_coass)
# combined_mags <- bind_rows(mags_parsed, mags_parsed_coass)

# # Write the combined summary to a TSV file
# # write_tsv(combined_summary, paste0(folder_path,"gunc_checkm_summary_combined_summary.tsv"))

# combined_mags_small <- combined_mags %>%
#   group_by(strain, replicate, treatment, buffer, bin_id, origin) %>%
#   summarize(
#     completeness = checkM.completeness,
#     contamination = checkM.contamination,
#     lineage = checkM.lineage,
#     .groups = "drop"
#   ) %>%
#   ungroup()

# #clean up the data
# # combined_mags_small$buffer <- ifelse(combined_mags_small$buffer == "D", "DESS", combined_mags_small$buffer)
# # combined_mags_small$buffer <- ifelse(combined_mags_small$buffer == "P", "PBS", combined_mags_small$buffer)
# # unique(combined_mags_small$buffer)
# combined_mags_small$treatment <- ifelse(combined_mags_small$treatment == "B", "Benzonase", combined_mags_small$treatment)
# # combined_mags_small$treatment <- ifelse(combined_mags_small$treatment == "Benzo", "Benzonase", combined_mags_small$treatment)
# combined_mags_small$treatment <- ifelse(combined_mags_small$treatment == "BT", "Blood & Tissue Kit", combined_mags_small$treatment)
# combined_mags_small$treatment <- ifelse(combined_mags_small$treatment == "MB", "Microbiome Kit", combined_mags_small$treatment)
# combined_mags_small$treatment <- ifelse(combined_mags_small$treatment == "MBX", "Microbiome Kit + Bead Beating", combined_mags_small$treatment)
# combined_mags_small$treatment <- ifelse(combined_mags_small$treatment == "MBS", "Microbiome Kit + Spinning", combined_mags_small$treatment)
# unique(combined_mags_small$treatment)
# combined_mags_small$strain <- ifelse(combined_mags_small$strain == "Acro", "Acropora", combined_mags_small$strain)
# combined_mags_small$strain <- ifelse(combined_mags_small$strain == "Por", "Porites", combined_mags_small$strain)
# combined_mags_small$strain <- ifelse(combined_mags_small$strain == "Poci", "Pocillopora", combined_mags_small$strain)
# unique(combined_mags_small$strain)

# #add grouping info
# combined_mags_small$group <- ifelse(combined_mags_small$strain == "Acropora", "Coral", combined_mags_small$strain)
# combined_mags_small$group <- ifelse(combined_mags_small$strain == "Porites", "Coral", combined_mags_small$group)
# combined_mags_small$group <- ifelse(combined_mags_small$strain == "Pocillopora", "Coral", combined_mags_small$group)
# combined_mags_small$group <- ifelse(combined_mags_small$strain == "H2", "Aiptasia", combined_mags_small$group)
# combined_mags_small$group <- ifelse(combined_mags_small$strain == "F003", "Aiptasia", combined_mags_small$group)

# # factor the data
# combined_mags_small$buffer <- factor(combined_mags_small$buffer, levels = rev(c("DESS", "PBS")))
# combined_mags_small$treatment <- factor(combined_mags_small$treatment, levels = c("Benzonase", "Blood & Tissue Kit", "Microbiome Kit", "Microbiome Kit + Bead Beating", "Microbiome Kit + Spinning"))
# combined_mags_small$strain <- factor(combined_mags_small$strain, levels = rev(c("Acropora", "Porites", "Pocillopora", "F003", "H2")))

# # make lineage factors alphabetically
# #make column lineage_clean
# combined_mags_small$lineage_clean <- gsub("[gpscfko]__", "", gsub(" \\(UID[0-9]{4}\\)", "",combined_mags_small$lineage))
# combined_mags_small$lineage_clean <- gsub(" \\(UID[0-9]{3}\\)", "", combined_mags_small$lineage_clean)
# combined_mags_small$lineage_clean <- gsub(" \\(UID[0-9]{2}\\)", "", combined_mags_small$lineage_clean)
# combined_mags_small$lineage_clean <- gsub(" \\(UID[0-9]{1}\\)", "", combined_mags_small$lineage_clean)

# # sort combined_mags_small by lineage_clean
# # combined_mags_small$lineage_clean <- factor(combined_mags_small$lineage_clean, levels = rev(sort(unique(combined_mags_small$lineage_clean))))
# combined_mags_small$lineage <- factor(combined_mags_small$lineage, levels = rev(sort(unique(combined_mags_small$lineage))))

# # gsub("c__", "", unique(combined_mags_small$lineage))
# #  unique(combined_mags_small$lineage))
# # taxonomy <- unique(gsub("[gpscfko]__", "", gsub(" \\(UID[0-9]{4}\\)", "", unique(combined_mags_small$lineage))))
# # taxonomy <- unique(gsub(" \\(UID[0-9]{3}\\)", "", taxonomy))
# # taxonomy <- unique(gsub(" \\(UID[0-9]{2}\\)", "", taxonomy))
# # taxonomy <- unique(gsub(" \\(UID[0-9]{1}\\)", "", taxonomy))

# # unique(taxonomy)
# # Define custom colors for each taxonomy
# taxonomy_colors <- c(
#     "Bacteria" = "#1f77b4",
#     "Actinomycetales" = "#ff7f0e",
#     "Rhizobiales" = "#2ca02c",
#     "Rhodospirillales" = "#d62728",
#     "Gammaproteobacteria" = "#9467bd",
#     "algicola" = "#8c564b",
#     "Mollicutes" = "#e377c2",
#     "Rhodobacteraceae" = "#7f7f7f",
#     "root" = "#bcbd22",
#     "Archaea" = "#17becf",
#     "Cyanobacteria" = "#aec7e8",
#     "Alphaproteobacteria" = "#ffbb78",
#     "Flavobacteriaceae" = "#98df8a",
#     "Proteobacteria" = "#ff9896",
#     "Deltaproteobacteria" = "#c5b0d5",
#     "Vibrio" = "#c49c94",
#     "Burkholderiaceae" = "#f7b6d2",
#     "Bradyrhizobium" = "#c7c7c7",
#     "Mycoplasma" = "#dbdb8d",
#     "Cytophagales" = "#9edae5",
#     "Bacteroidetes" = "#ffbb78"
# )

# # Apply custom colors to the plot
# # Create a function to match taxonomy with colors
# get_taxonomy_color <- function(lineage, taxonomy_colors) {
#     for (taxon in names(taxonomy_colors)) {
#         if (grepl(taxon, lineage, ignore.case = FALSE)) {
#             return(taxonomy_colors[[taxon]])
#         }
#     }
#     return(NA)  # Return NA if no match is found
# }

# # Apply the function to get colors for each lineage
# combined_mags_small$lineage_color <- sapply(combined_mags_small$lineage, get_taxonomy_color, taxonomy_colors)

# # Convert taxon_color to a named vector
# taxon_color_named <- setNames(combined_mags_small$lineage_color, combined_mags_small$lineage)
# # p_load("lemon")
# # p_load("ggforce")

# combined_mags_small$grouping <- paste(combined_mags_small$strain, combined_mags_small$buffer, sep = " - ")

# unique(combined_mags_small$grouping)


# head(combined_mags_small)

# # Get only combined_mags_small with completeness > 50 and contamination < 10
# combined_mags_keep <- combined_mags_small %>%
#     filter(origin == "co-assembly") %>%
#     # filter(origin == "co-binning" ) %>%
#     filter(completeness >= 50, contamination <= 10)
# head(combined_mags_keep)

# # add column to mark if completeness is above 90
# combined_mags_keep$completeness_above_90 <- ifelse(combined_mags_keep$completeness >= 90, "\u226590", "\u226550")

# # get number of bins per strain, treatment, buffer, for completeness_above_90 and above 50
# %% Create Table 1: Co-assembly vs Co-binning comparison
# Create separate datasets for each origin
co_assembly_data <- final_mags_summary %>% filter(Origin == "co-assembly") %>% select(-Origin)
co_binning_data <- final_mags_summary %>% filter(Origin == "co-binning") %>% select(-Origin)

# Combine with origin as column suffix
comparison_data <- co_binning_data %>%
    # full_join(co_binning_data, by = c("Species", "Buffer", "Extraction"), suffix = c(".Co-assembly", ".Co-binning")) %>%
    # Replace NAs with "-"
    mutate(across(everything(), ~ifelse(is.na(.), "-", as.character(.)))) %>%
    mutate(
        Species = factor(Species, levels = c("H2", "F003", "Acropora", "Pocillopora", "Porites")),
        Buffer = factor(Buffer, levels = c("PBS", "DESS"))
    ) %>%
    arrange(Species, Buffer, Extraction)

# Create the comparison table
table1 <- nice_table(comparison_data, separate.header = TRUE) %>%
    align(align = "justify", part = "all")

# Add borders for species/buffer group changes
species_col <- comparison_data$Species
buffer_col <- comparison_data$Buffer
group_id <- paste(species_col, buffer_col, sep = "_")
group_changes <- which(c(TRUE, group_id[-1] != group_id[-length(group_id)]))
group_changes_filtered <- group_changes[group_changes > 1]

for(row in group_changes_filtered) {
    table1 <- table1 %>%
        flextable::hline(i = row - 1, border = fp_border(color = "black", width = 1))
}

# Apply the same formatting as table_counts.r
table1_clean <- table1 %>%
    flextable::set_formatter(
        Species = function(x) {
            species_col <- comparison_data$Species
            buffer_col <- comparison_data$Buffer
            result <- character(length(x))
            result[1] <- as.character(species_col[1])
            for(i in 2:length(x)) {
                if(species_col[i] != species_col[i-1] || buffer_col[i] != buffer_col[i-1]) {
                    result[i] <- as.character(species_col[i])
                } else {
                    result[i] <- ""
                }
            }
            return(result)
        },
        Buffer = function(x) {
            species_col <- comparison_data$Species
            buffer_col <- comparison_data$Buffer
            result <- character(length(x))
            result[1] <- as.character(buffer_col[1])
            for(i in 2:length(x)) {
                if(species_col[i] != species_col[i-1] || buffer_col[i] != buffer_col[i-1]) {
                    result[i] <- as.character(buffer_col[i])
                } else {
                    result[i] <- ""
                }
            }
            return(result)
        }
    ) %>%
    align(j = c("Species", "Buffer"), align = "left", part = "body") %>%
    align(j = c("≥50", "≥90"), align = "right", part = "all") %>%
    # Apply same formatting improvements as table_counts.r
    autofit() %>%
    fontsize(size = 9, part = "all") %>%
    padding(padding.top = 3, padding.bottom = 3, part = "all") %>%
    line_spacing(space = 1.15, part = "all") %>%
    bold(part = "header") %>%
    width(j = c("Species", "Buffer", "Extraction"), width = 1.2) %>%
    width(j = c("≥50", "≥90"), width = 1.0)

# %% Create Table 2: Full treatment table (similar to the last table in the original)

# Create separate datasets for each origin
co_assembly_data <- final_mags_summary %>% filter(Origin == "co-assembly") %>% select(-Origin)
co_binning_data <- final_mags_summary %>% filter(Origin == "co-binning") %>% select(-Origin)

# Combine with origin as column suffix
full_comparison_data <- co_binning_data %>%
    full_join(co_assembly_data, by = c("Species", "Buffer", "Extraction"), suffix = c(".Co-assembly", ".Co-binning")) %>%
    # Replace NAs with "-"
    mutate(across(everything(), ~ifelse(is.na(.), "-", as.character(.)))) %>%
    mutate(
        Species = factor(Species, levels = c("H2", "F003", "Acropora", "Pocillopora", "Porites")),
        Buffer = factor(Buffer, levels = c("PBS", "DESS"))
    ) %>%
    arrange(Species, Buffer, Extraction)


# Ensure all combinations exist
all_species <- unique(full_comparison_data$Species)
all_buffers <- unique(full_comparison_data$Buffer)
all_extractions <- c("Blood & Tissue", "Benzonase", "Microbiome", "Microbiome & Bead-Beating", "Microbiome & Spinning")
# names(final_mags_summary) <- c("Species", "Buffer", "Extraction", "≥50", "≥90", "origin")

full_combinations <- expand.grid(
    Species = all_species,
    Buffer = all_buffers,
    Extraction = all_extractions,
    stringsAsFactors = FALSE
) %>%
    # Apply same filters as original
    # filter(Species != "F003") %>%
    filter(!(Species %in% c("Acropora", "Porites", "Pocillopora") & Buffer == "PBS"))

# Use co-assembly data for this table (as in original)
full_table_data <- full_combinations %>%
    left_join(full_comparison_data, # %>% filter(origin == "co-binning") %>% select(-origin) 
              by = c("Species", "Buffer", "Extraction")) %>%
    mutate(
        `≥50.Co-assembly` = ifelse(is.na(`≥50.Co-assembly`), "-", `≥50.Co-assembly`),
        `≥90.Co-assembly` = ifelse(is.na(`≥90.Co-assembly`), "-", `≥90.Co-assembly`),
        `≥50.Co-binning` = ifelse(is.na(`≥50.Co-binning`), "-", `≥50.Co-binning`),
        `≥90.Co-binning` = ifelse(is.na(`≥90.Co-binning`), "-", `≥90.Co-binning`),
    ) %>%
    arrange(Species, Buffer, Extraction)

full_table_data <- full_table_data %>%
    mutate(
        across(
            c(`≥50.Co-assembly`, `≥90.Co-assembly`, `≥50.Co-binning`, `≥90.Co-binning`),
            ~ifelse(is.na(.), "-", as.character(.))
        ),
        # Set "NA" for the specified combinations only
        `≥50.Co-assembly` = ifelse(
            (Extraction == "Microbiome & Spinning" & Species == "H2" & Buffer == "DESS") |
            (Extraction == "Benzonase" & Species == "F003" & Buffer == "DESS") |
            (Extraction == "Microbiome" & Species == "F003" & Buffer == "DESS") |
            (Extraction == "Microbiome & Spinning" & Species == "F003" & Buffer == "DESS"),
            "NA",
            `≥50.Co-assembly`
        ),
        `≥90.Co-assembly` = ifelse(
            (Extraction == "Microbiome & Spinning" & Species == "H2" & Buffer == "DESS") |
            (Extraction == "Benzonase" & Species == "F003" & Buffer == "DESS") |
            (Extraction == "Microbiome" & Species == "F003" & Buffer == "DESS") |
            (Extraction == "Microbiome & Spinning" & Species == "F003" & Buffer == "DESS"),
            "NA",
            `≥90.Co-assembly`
        ),
        `≥50.Co-binning` = ifelse(
            (Extraction == "Microbiome & Spinning" & Species == "H2" & Buffer == "DESS") |
            (Extraction == "Benzonase" & Species == "F003" & Buffer == "DESS") |
            (Extraction == "Microbiome" & Species == "F003" & Buffer == "DESS") |
            (Extraction == "Microbiome & Spinning" & Species == "F003" & Buffer == "DESS"),
            "NA",
            `≥50.Co-binning`
        ),
        `≥90.Co-binning` = ifelse(
            (Extraction == "Microbiome & Spinning" & Species == "H2" & Buffer == "DESS") |
            (Extraction == "Benzonase" & Species == "F003" & Buffer == "DESS") |
            (Extraction == "Microbiome" & Species == "F003" & Buffer == "DESS") |
            (Extraction == "Microbiome & Spinning" & Species == "F003" & Buffer == "DESS"),
            "NA",
            `≥90.Co-binning`
        )
    )

# Invert the names at the dot for relevant columns
names(full_table_data) <- gsub(
    "^(≥[0-9]+)\\.(Co-assembly|Co-binning)$",
    "\\2.\\1",
    names(full_table_data)
    )

# Create the full table
table2 <- nice_table(full_table_data, separate.header = TRUE) %>%
    align(align = "justify", part = "all")

# Add borders for species/buffer group changes
species_col2 <- full_table_data$Species
buffer_col2 <- full_table_data$Buffer
group_id2 <- paste(species_col2, buffer_col2, sep = "_")
group_changes2 <- which(c(TRUE, group_id2[-1] != group_id2[-length(group_id2)]))
group_changes_filtered2 <- group_changes2[group_changes2 > 1]

for(row in group_changes_filtered2) {
    table2 <- table2 %>%
        flextable::hline(i = row - 1, border = fp_border(color = "black", width = 1))
}

# Apply formatting
table2_clean <- table2 %>%
    flextable::set_formatter(
        Species = function(x) {
            species_col <- full_table_data$Species
            buffer_col <- full_table_data$Buffer
            result <- character(length(x))
            result[1] <- as.character(species_col[1])
            for(i in 2:length(x)) {
                if(species_col[i] != species_col[i-1] || buffer_col[i] != buffer_col[i-1]) {
                    result[i] <- as.character(species_col[i])
                } else {
                    result[i] <- ""
                }
            }
            return(result)
        },
        Buffer = function(x) {
            species_col <- full_table_data$Species
            buffer_col <- full_table_data$Buffer
            result <- character(length(x))
            result[1] <- as.character(buffer_col[1])
            for(i in 2:length(x)) {
                if(species_col[i] != species_col[i-1] || buffer_col[i] != buffer_col[i-1]) {
                    result[i] <- as.character(buffer_col[i])
                } else {
                    result[i] <- ""
                }
            }
            return(result)
        }
    ) %>%
    align(j = c("Species", "Buffer"), align = "left", part = "body") %>%
    # align(j = c("≥50", "≥90"), align = "right", part = "all") %>%
    align(j = c("Co-assembly.≥50", "Co-assembly.≥90", "Co-binning.≥50", "Co-binning.≥90"), align = "right", part = "all") %>%

    # Apply same formatting improvements as table_counts.r
    autofit() %>%
    fontsize(size = 9, part = "all") %>%
    padding(padding.top = 3, padding.bottom = 3, part = "all") %>%
    line_spacing(space = 1.15, part = "all") %>%
    bold(part = "header") %>%
    width(j = c("Species", "Buffer", "Extraction"), width = 1.2) %>%
    width(j = c("Co-assembly.≥50", "Co-assembly.≥90", "Co-binning.≥50", "Co-binning.≥90"), width = 1.0)

# %% Print tables
cat("Printing Table 1: Co-assembly vs Co-binning comparison...\n")
print(table1_clean, preview = "docx")

cat("Printing Table 2: Full treatment table (Co-assembly only)...\n") 
print(table2_clean, preview = "docx")

cat("Tables generated successfully!\n")



