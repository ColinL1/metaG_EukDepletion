#! /usr/bin/env Rscript

cat("Loading libraries...\n")
invisible(capture.output((library(pacman))))
# p_load("tidyr", "tidyfast", "dplyr", "ggplot2", "ggpubr", "tidyverse", "scales", "RColorBrewer", "optparse", "pacman","httpgd", "grid", "stringr") # nolint: line_length_linter, cSpell. install lemon? 
# p_load("tidyr", "tidyfast", "dplyr", "tidyverse", "scales", "httpgd", "grid") # nolint: line_length_linter, cSpell. install lemon? 
p_load(rempsyc, flextable, officer, broom, report, effectsize, tidyverse)

# --- Clear workspace ---
# rm(list = ls())

# %% Read the file
data_reads <- read.table("<PROJECT_ROOT>/input/stats_reads_sanitised.txt", header = TRUE, sep = "\t")

#check that all file *_1.fq.gz are identical to *_2.fq.gz redo later down too
data_reads %>%
    filter(grepl("_1.fq.gz", file)) %>%
    mutate(file2 = gsub("_1.fq.gz", "_2.fq.gz", file)) %>%
    left_join(data_reads %>% filter(grepl("_2.fq.gz", file)), by = c("file2" = "file")) %>%
    filter(num_seqs.x != num_seqs.y)
# if empty then all good

# drop all *_2.fq.gz or *R2.fq.gz lines
data_reads <- data_reads %>% filter(!grepl("_2.fq.gz", file))
data_reads <- data_reads %>% filter(!grepl("R2.fq.gz", file))

#cleanup data_reads file column removing any path
data_reads$file <- basename(data_reads$file)

#copy file column to temp column
data_reads$temp <- data_reads$file

#change _assembly to .assembly
data_reads$temp <- gsub("_TRIM_paried_1.fq.gz", ".total", data_reads$temp)
data_reads$temp <- gsub("_TRIM_paried_R1.fq.gz", ".total", data_reads$temp)
data_reads$temp <- gsub("_L[0-9]\\.", "\\.", data_reads$temp)

#remove _EKDN230034839-1A_HG3F2DSX7 
data_reads$temp <- gsub("_EKDN230034839-1A_HG3F2DSX7", "", data_reads$temp)
data_reads$temp <- gsub("_EKDN230034842-1A_HG3WCDSX7", "", data_reads$temp)
data_reads$temp <- gsub("_EKDN230034838-1A_HG3F2DSX7", "", data_reads$temp)
data_reads$temp <- gsub("_EKDN230034843-1A_HFVN2DSX7", "", data_reads$temp)

data_reads$temp <- gsub("_TRIM_paried", "", data_reads$temp)
data_reads$temp <- gsub(".mapped_1.fq.gz|.unmapped_1.fq.gz", "", data_reads$temp)

#split file coluimn at the . to get sample name and mapping
data_reads <- data_reads %>% separate(temp, c("sample", "mapping"), sep = "\\.", remove = FALSE)

#split sample to get metadata using "_"
data_reads <- data_reads %>% separate(sample, c("species", "replicate", "extraction", "buffer"), sep = "_", remove = FALSE)

#remove temp column
data_reads$temp <- NULL

#sort columns in order sample species replicate extraction buffer mapping num_seqs sum_len min_len avg_len max_len file format type
data_reads <- data_reads[, c("sample", "species", "replicate", "extraction", "buffer", "mapping", "num_seqs", "sum_len", "min_len", "avg_len", "max_len", "file", "format", "type")]

# clean up species names if necessary
#"F003" "F3"   "H2"   "Ac"   "Acro" "Po"   "Poci" "Por"  "Pr" to full names
data_reads$species <- gsub("\\bF3\\b", "F003", data_reads$species)
data_reads$species <- gsub("\\bAcro\\b", "Acropora", data_reads$species)
data_reads$species <- gsub("\\bAc\\b", "Acropora", data_reads$species)
data_reads$species <- gsub("\\bPo\\b", "Pocillopora", data_reads$species)
data_reads$species <- gsub("\\bPoci\\b", "Pocillopora", data_reads$species)
data_reads$species <- gsub("\\bPor\\b", "Porites", data_reads$species)
data_reads$species <- gsub("\\bPr\\b", "Porites", data_reads$species)

# extractions methods  B   BT  MB  MBX MBS to Benzonase "Blood & Tissue" "Microbiome" "Microbiome & beadbeating" "Microbiome & Spinning"
data_reads$extraction <- gsub("\\bB\\b", "Benzonase", data_reads$extraction)
data_reads$extraction <- gsub("\\bBT\\b", "Blood & Tissue", data_reads$extraction)
data_reads$extraction <- gsub("\\bMB\\b", "Microbiome", data_reads$extraction)
data_reads$extraction <- gsub("\\bMBX\\b", "Microbiome & Bead-Beating", data_reads$extraction)
data_reads$extraction <- gsub("\\bMBS\\b", "Microbiome & Spinning", data_reads$extraction)

#P to PBS and D to DESS
data_reads$buffer <- gsub("\\bP\\b", "PBS", data_reads$buffer)
data_reads$buffer <- gsub("\\bD\\b", "DESS", data_reads$buffer)

# fix mapping name
data_reads$mapping <- gsub("\\bscleractina\\b", "Host", data_reads$mapping)
data_reads$mapping <- gsub("\\baiptasia\\b", "Host", data_reads$mapping)
data_reads$mapping <- gsub("\\bsymbiodiniaceae\\b", "Symbiodiniaceae", data_reads$mapping)
data_reads$mapping <- gsub("\\bbacteria\\b", "Bacteria", data_reads$mapping)
data_reads$mapping <- gsub("\\bother\\b", "Other", data_reads$mapping)

data_reads$mapping <- factor(data_reads$mapping, levels = (rev(c("Host", "Symbiodiniaceae","Bacteria", "Other","total"))))

#check if for "total" sum_len is the same as the sum of all the other sum_len per each mapping 
data_reads %>%
    filter(mapping != "total") %>%
    group_by(sample) %>%
    summarise(total = sum(sum_len)) %>%
    left_join(data_reads %>% filter(mapping == "total"), by = "sample") %>%
    mutate(diff = total - sum_len) %>%
    filter(diff != 0)

#if empty then all good
#add origin information to two dfs # assembly df discarded
data_reads$origin <- "reads"

#merge two dfs
data <- data_reads

data["grouping"] <- paste(data$species, data$buffer, sep = " - ")
data$grouping <- factor(data$grouping, levels = c("H2 - PBS", "H2 - DESS", "F003 - PBS", "F003 - DESS", "Acropora - DESS", "Porites - DESS", "Pocillopora - DESS"))

#set data$species as ordered factor
data$species <- factor(data$species, levels = rev(c( "Porites", "Pocillopora", "Acropora" ,"F003", "H2")))
data$buffer <- factor(data$buffer, levels = rev(sort(unique(data$buffer))))

data$extraction <- factor(data$extraction, levels = (c("Blood & Tissue", "Benzonase", "Microbiome",  "Microbiome & Bead-Beating", "Microbiome & Spinning")))

## make extractions groups to sort sample
data$extraction_groups <- data$extraction
data$extraction_groups <- gsub("Microbiome & Bead-Beating", "2", data$extraction_groups)
data$extraction_groups <- gsub("Microbiome & Spinning", "2", data$extraction_groups)

data$extraction_groups <- gsub("Benzonase", "1", data$extraction_groups)
data$extraction_groups <- gsub("Blood & Tissue", "1", data$extraction_groups)
data$extraction_groups <- gsub("Microbiome", "1", data$extraction_groups)

data$extraction <- factor(data$extraction, levels = rev(c("Blood & Tissue", "Benzonase", "Microbiome",  "Microbiome & Bead-Beating", "Microbiome & Spinning")))
# data$extraction <- factor(data$extraction, levels = rev(c("Blood & Tissue", "Benzonase", "Microbiome", " ",  "Microbiome & Bead-Beating", "Microbiome & Spinning")))

## add number of replicate uniquer to each grouping and treatment to the df.
n_samples <- data %>%
    filter(origin == "reads") %>%
    group_by(grouping, extraction) %>%
    summarise(n_sample = n_distinct(replicate))

data <- left_join(data, n_samples, by = c("grouping", "extraction"))


head(data)

#extract archea from total

data <- data[,c("sample", "grouping", "type", "species", "replicate", "extraction", "buffer", "mapping", "num_seqs")]

# Step 1: Pivot wider so each mapping becomes a column
df_wide <- data %>%
  pivot_wider(names_from = mapping, values_from = num_seqs) # values_fill = 0

# Step 2: Compute "arch"
df_with_arch <- df_wide %>%
  mutate(Archaea = total - (Host + Symbiodiniaceae + Bacteria + Other))

# Step 3: Pivot back to long format (adding arch as new mapping)
df_final <- df_with_arch %>%
  pivot_longer(cols = c(total, Host, Symbiodiniaceae, Bacteria, Archaea, Other),
               names_to = "mapping",
               values_to = "num_seqs") %>%
  arrange(sample, mapping)  # Optional: sort by sample then mapping

# test_data <- df_final %>%
#     filter(mapping != "total")

data <- df_final

# flextable(data)

# Summarise data by species, extraction, buffer, and mapping for num_seqs
summary_data <- data %>%
    group_by(species, extraction, buffer, mapping) %>%
    summarise(total_num_seqs = sum(num_seqs), .groups = "drop")

print(summary_data)

summary_data_wide <- summary_data %>%
    pivot_wider(names_from = mapping, values_from = total_num_seqs)

print(summary_data_wide)
summary_data_wide <- summary_data_wide[,c("species", "buffer", "extraction", "Host", "Symbiodiniaceae", "Bacteria", "Archaea", "Other", "total")]
# flextable(summary_data_wide)
# Create a percentage-based version of summary_data_wide
summary_data_pct <- summary_data_wide %>%
    mutate(
        Host_pct = Host / total * 100,
        Symbiodiniaceae_pct = Symbiodiniaceae / total * 100,
        Bacteria_pct = Bacteria / total * 100,
        Archaea_pct = Archaea / total * 100,
        Other_pct = Other / total * 100
    )

# Optionally, select only the percentage columns and relevant metadata
summary_data_pct_only <- summary_data_pct %>%
    select(species, buffer, extraction, Host_pct, Symbiodiniaceae_pct, Bacteria_pct, Archaea_pct, Other_pct) %>% 
    arrange(species, buffer, extraction)

    # Ensure extraction is a factor with the desired order before arranging
    summary_data_pct_only$extraction <- factor(
        summary_data_pct_only$extraction,
        levels = c("Blood & Tissue", "Benzonase", "Microbiome", "Microbiome & Bead-Beating", "Microbiome & Spinning")
    )
    summary_data_pct_only <- summary_data_pct_only %>%
        arrange(species, buffer, extraction)

print(summary_data_pct_only)

names(summary_data_pct_only) <- c("Species", "Buffer", "Extraction", "Host (%)", "Symbiodiniaceae (%)", "Bacteria (%)", "Archaea (%)", "Other (%)")

names(summary_data_pct_only) <- c(
    "Species", 
    "Buffer", 
    "Extraction", 
    "Host (\u0025)", 
    "Symbiodiniaceae (\u0025)", 
    "Bacteria (\u0025)", 
    "Archaea (\u0025)", 
    "Other (\u0025)"
)

# Transpose count table 1 ( same as suplemetary table 1 for dess and pbs effects)
# Prepare data for transposition
transposed_data <- summary_data_pct_only %>%
    pivot_longer(
        cols = c("Host (\u0025)", "Symbiodiniaceae (\u0025)", "Bacteria (\u0025)", "Archaea (\u0025)", "Other (\u0025)"),
        names_to = "Mapping",
        values_to = "Percentage"
    ) %>%
    pivot_wider(
        names_from = c(Species, Buffer),
        values_from = Percentage
    ) %>%
    mutate(across(
        c("H2_PBS", "H2_DESS", "F003_PBS", "F003_DESS", "Acropora_DESS", "Pocillopora_DESS", "Porites_DESS"),
        ~ ifelse(is.na(.), "-", round(., 2))
    ))

print(transposed_data)

names(transposed_data) <- gsub("_", ".", names(transposed_data)) # Replace spaces with underscores for easier column referencing
# names(transposed_data) <- gsub("_", " - ", names(transposed_data)) # Replace spaces with underscores for easier column referencing
counts <- nice_table(transposed_data, separate.header = TRUE) %>%
  align(align = "justify", part = "all")

# Add borders for species/buffer group changes
extraction_col <- transposed_data$Extraction
# mapping_col <- transposed_data$Mapping

# Create a combined species-buffer identifier to find group changes
group_id <- extraction_col # paste(extraction_col, mapping_col, sep = "_")
group_changes <- which(c(TRUE, group_id[-1] != group_id[-length(group_id)]))

# Add top borders for group changes (exclude first row since it's already at top)
group_changes_filtered <- group_changes[group_changes > 1]
for(row in group_changes_filtered) {
  counts <- counts %>%
    flextable::hline(i = row - 1, border = fp_border(color = "black", width = 1))
}

counts_clean <- counts %>%
    flextable::set_formatter(
        Extraction = function(x) {
            # Get the original data to check for changes
            extraction_col <- transposed_data$Extraction
            # mapping_col <- transposed_data$Mapping
            result <- character(length(x))
            result[1] <- as.character(extraction_col[1])  # Always show first
            for(i in 2:length(x)) {
                # Show species name when either species OR buffer changes
                if(extraction_col[i] != extraction_col[i-1]) {
                    result[i] <- as.character(extraction_col[i])
                } else {
                    result[i] <- ""
                }
            }
            return(result)
        },
        Mapping = function(x) {
            # Get the original data to check for changes within species
            extraction_col <- transposed_data$Extraction
            # mapping_col <- transposed_data$Mapping
            result <- character(length(x))
            result[1] <- as.character(mapping_col[1])  # Always show first
            for(i in 2:length(x)) {
                if(extraction_col[i] != extraction_col[i-1] ) {
                    result[i] <- as.character(mapping_col[i])
                } else {
                    result[i] <- as.character(mapping_col[i])
                }
            }
            return(result)
        }
    ) %>%
    align(j = c("Extraction", "Mapping"), align = "left", part = "body") %>%
    align(j = c("H2 - PBS", "H2 - DESS", "F003 - PBS", "F003 - DESS" ,"Acropora - DESS" ,"Pocillopora - DESS", "Porites - DESS"), align = "right", part = "all") %>%
    # Improve table spacing and formatting
    autofit() %>%  # Automatically adjust column widths
    fontsize(size = 9, part = "all") %>%  # Set font size
    padding(padding.top = 3, padding.bottom = 3, part = "all") %>%  # Add vertical padding
    line_spacing(space = 1.15, part = "all") %>%  # Adjust line spacing
    # Make header bold
    bold(part = "header") %>%
    # Add some spacing between columns
    width(j = c("Extraction", "Mapping"), width = 1.2) %>%
    width(j = c("H2 - PBS", "H2 - DESS", "F003 - PBS", "F003 - DESS" ,"Acropora - DESS" ,"Pocillopora - DESS", "Porites - DESS"), width = 1.0) 

print(counts_clean, preview = "docx")



# Transpose count table 2 sorted by mapping and extraction after. alternative to transpose count table 1
# Prepare data for transposition
transposed_data <- summary_data_pct_only %>%
    pivot_longer(
        cols = c("Host (\u0025)", "Symbiodiniaceae (\u0025)", "Bacteria (\u0025)", "Archaea (\u0025)", "Other (\u0025)"),
        names_to = "Mapping",
        values_to = "Percentage"
    ) %>%
    pivot_wider(
        names_from = c(Species, Buffer),
        values_from = Percentage
    ) %>%
    mutate(across(
        c("H2_PBS", "H2_DESS", "F003_PBS", "F003_DESS", "Acropora_DESS", "Pocillopora_DESS", "Porites_DESS"),
        ~ ifelse(is.na(.), "-", round(., 2))
    )) %>% 
    mutate(
        Mapping = factor(
            Mapping,
            levels = c("Host (\u0025)", "Symbiodiniaceae (\u0025)", "Bacteria (\u0025)", "Archaea (\u0025)", "Other (\u0025)")
        ),
        Extraction = factor(
            Extraction,
            levels = c("Blood & Tissue", "Benzonase", "Microbiome", "Microbiome & Bead-Beating", "Microbiome & Spinning")
        )
    ) %>%
    arrange(Mapping, Extraction)

print(transposed_data)
transposed_data <- transposed_data[,c("Mapping", "Extraction", "H2_PBS", "H2_DESS", "F003_PBS", "F003_DESS", "Acropora_DESS", "Pocillopora_DESS", "Porites_DESS")] 


# names(transposed_data) <- gsub("_", ".", names(transposed_data)) # Replace spaces with underscores for easier column referencing
names(transposed_data) <- gsub("_", " - ", names(transposed_data)) # Replace spaces with underscores for easier column referencing
counts <- nice_table(transposed_data, separate.header = TRUE) %>%
  align(align = "justify", part = "all")

# Add borders for species/buffer group changes
# extraction_col <- transposed_data$Extraction
mapping_col <- transposed_data$Mapping

# Create a combined species-buffer identifier to find group changes
group_id <- mapping_col # paste(extraction_col, mapping_col, sep = "_")
group_changes <- which(c(TRUE, group_id[-1] != group_id[-length(group_id)]))

# Add top borders for group changes (exclude first row since it's already at top)
group_changes_filtered <- group_changes[group_changes > 1]
for(row in group_changes_filtered) {
  counts <- counts %>%
    flextable::hline(i = row - 1, border = fp_border(color = "black", width = 1))
}

counts_clean <- counts %>%
    flextable::set_formatter(
        Extraction = function(x) {
            # Get the original data to check for changes
            extraction_col <- transposed_data$Extraction
            mapping_col <- transposed_data$Mapping
            result <- character(length(x))
            result[1] <- as.character(extraction_col[1])  # Always show first
            for(i in 2:length(x)) {
                # Show species name when either species OR buffer changes
                if(extraction_col[i] != extraction_col[i-1]) {
                    result[i] <- as.character(extraction_col[i])
                } else {
                    result[i] <- as.character(extraction_col[i])
                }
            }
            return(result)
        },
        Mapping = function(x) {
            # Get the original data to check for changes within species
            # extraction_col <- transposed_data$Extraction
            mapping_col <- transposed_data$Mapping
            result <- character(length(x))
            result[1] <- as.character(mapping_col[1])  # Always show first
            for(i in 2:length(x)) {
                if(mapping_col[i] != mapping_col[i-1] ) {
                    result[i] <- as.character(mapping_col[i])
                } else {
                    result[i] <- ""
                }
            }
            return(result)
        }
    ) %>%
    align(j = c("Extraction", "Mapping"), align = "left", part = "body") %>%
    align(j = c("H2 - PBS", "H2 - DESS", "F003 - PBS", "F003 - DESS" ,"Acropora - DESS" ,"Pocillopora - DESS", "Porites - DESS"), align = "right", part = "all") %>%
    # Improve table spacing and formatting
    autofit() %>%  # Automatically adjust column widths
    fontsize(size = 9, part = "all") %>%  # Set font size
    padding(padding.top = 3, padding.bottom = 3, part = "all") %>%  # Add vertical padding
    line_spacing(space = 1.15, part = "all") %>%  # Adjust line spacing
    # Make header bold
    bold(part = "header") %>%
    # Add some spacing between columns
    width(j = c("Extraction", "Mapping"), width = 1.2) %>%
    width(j = c("H2 - PBS", "H2 - DESS", "F003 - PBS", "F003 - DESS" ,"Acropora - DESS" ,"Pocillopora - DESS", "Porites - DESS"), width = 1.0) 

print(counts_clean, preview = "docx")




# Transpose count table 3 ( same as suplemetary table 1 for dess and pbs effects)
# Prepare data for transposition: extraction methods as columns
transposed_data <- summary_data_pct_only %>%
    filter(!(Species %in% c("H2", "F003") & Buffer == "DESS")) %>%
    pivot_longer(
        cols = c("Host (\u0025)", "Symbiodiniaceae (\u0025)", "Bacteria (\u0025)", "Archaea (\u0025)", "Other (\u0025)"),
        names_to = "Mapping",
        values_to = "Percentage"
    ) %>%
    mutate(
        Mapping = factor(
            Mapping,
            levels = c("Host (\u0025)", "Symbiodiniaceae (\u0025)", "Bacteria (\u0025)", "Archaea (\u0025)", "Other (\u0025)")
        ),
        Species = factor(
            Species,
            levels = c("H2", "F003", "Acropora", "Pocillopora", "Porites")
        ),
        Buffer = factor(
            Buffer,
            levels = c("PBS", "DESS")
        )
    ) %>%
    arrange(Mapping, Species, Buffer) %>%
    pivot_wider(
        names_from = Extraction,
        values_from = Percentage
    ) %>%
    mutate(across(
        c("Blood & Tissue", "Benzonase", "Microbiome", "Microbiome & Bead-Beating", "Microbiome & Spinning"),
        ~ ifelse(is.na(.), "-", round(., 2))
    ))

print(transposed_data)
transposed_data <- transposed_data %>%
    select(Mapping, Species, `Blood & Tissue`, Benzonase, Microbiome, `Microbiome & Bead-Beating`, `Microbiome & Spinning`)
counts <- nice_table(transposed_data, separate.header = TRUE) %>%
  align(align = "justify", part = "all")

# Add borders for mapping group changes
mapping_col <- transposed_data$Mapping
group_id <- mapping_col
group_changes <- which(c(TRUE, group_id[-1] != group_id[-length(group_id)]))
group_changes_filtered <- group_changes[group_changes > 1]
for(row in group_changes_filtered) {
  counts <- counts %>%
    flextable::hline(i = row - 1, border = fp_border(color = "black", width = 1))
}

counts_clean <- counts %>%
    flextable::set_formatter(
        Mapping = function(x) {
            mapping_col <- transposed_data$Mapping
            result <- character(length(x))
            result[1] <- as.character(mapping_col[1])
            for(i in 2:length(x)) {
                if(mapping_col[i] != mapping_col[i-1]) {
                    result[i] <- as.character(mapping_col[i])
                } else {
                    result[i] <- ""
                }
            }
            return(result)
        }
    ) %>%
    align(j = "Mapping", align = "left", part = "body") %>%
    align(j = c("Blood & Tissue", "Benzonase", "Microbiome", "Microbiome & Bead-Beating", "Microbiome & Spinning"), align = "right", part = "all") %>%
    autofit() %>%
    fontsize(size = 9, part = "all") %>%
    padding(padding.top = 3, padding.bottom = 3, part = "all") %>%
    line_spacing(space = 1.15, part = "all") %>%
    bold(part = "header") %>%
    width(j = "Mapping", width = 1.2) %>%
    width(j = c("Blood & Tissue", "Benzonase", "Microbiome", "Microbiome & Bead-Beating", "Microbiome & Spinning"), width = 1.0)

print(counts_clean, preview = "docx")
