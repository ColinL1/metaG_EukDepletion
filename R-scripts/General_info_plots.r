#!/usr/bin/env Rscript

# echo status
cat("Loading libraries...\n")
lib_list <- c("tidyr", "tidyfast", "dplyr", "ggplot2", "ggpubr", "tidyverse", "scales", "RColorBrewer", "optparse", "pacman","httpgd", "grid") # nolint: line_length_linter, cSpell.
# invisible(
#     capture.output(
#         suppressMessages(
lapply(lib_list, require, quietly = TRUE, character.only = TRUE)


# #load tydifast with pacman (it installs it if not present)
# invisible(
#     capture.output(
#         suppressMessages(
#                 p_load("tidyfast"))))

# # Run "webplot" server (only necessary for live use)
hgd() #  http://127.0.0.1:38653/live?token=HraEkCgA

data <- read_csv("/home/colinl/metaG/Git/metaG_EukDepletion/tmp_json/report_with_legths_ONT_no_trim.csv")

data <- data[,c(2:5)]
# sequencing_stats <- read.table("/home/colinl/metaG/Git/metaG_methods/results/metaG_indonesia/mmseqs2_reports/ONT_reads_NR_lca.tsv", header = TRUE) # nolint: line_length_linter.
metadata <- read_csv("/home/colinl/metaG/Git/metaG_EukDepletion/sample_metadata_sheet.csv")

sample_data <- merge(data, metadata, by.y = "matching_index", by.x = "file_name")
names(sample_data)[names(sample_data) == 'Origin.y'] <- 'Location'
names(sample_data)[names(sample_data) == 'Origin.x'] <- 'Origin'

sample_data["sample"] <- paste(sample_data$Strain, sample_data$Replicate, sample_data$Location, sample_data$DATA_Type, sep = "-") 
sample_data["sample_short"] <- paste(sample_data$Replicate, sample_data$Location, sample_data$DATA_Type, sep = "-") 

sample_data["sample_short"] <- paste(sample_data$Replicate, sample_data$Location, sep = "-") 

sample_data["identifier"] <- paste(sample_data$Strain, sample_data$Replicate, sample_data$Treatment, sample_data$Buffer, sample_data$DATA_Type, sample_data$Location, sample_data$Key, sep = "_") 

sample_data_l <- sample_data[sample_data$Key == "length_reads",]
sample_data_c <- sample_data[sample_data$Key != "length_reads",]

overall_numbers_c <- sample_data_c %>%
        filter(Origin == "Overall")
sample_data_clean_c <- sample_data_c %>%
        filter(Origin != "Overall")

overall_numbers_c["sequenced_total"] <- overall_numbers_c$Value

overall_numbers_c <- overall_numbers_c[,c("identifier", "sequenced_total")]

sample_data_2_c <- merge(sample_data_clean_c, overall_numbers_c, by = "identifier")

sample_data_2_c["Total_percentage"] <- (sample_data_2_c$Value/sample_data_2_c$sequenced_total)*100

sample_df_c <- sample_data_2_c[,c("file_name", "Strain", "Replicate", "Treatment", "Buffer", "DATA_Type", "Location", "sample", "sample_short", "Value", "sequenced_total", "Total_percentage", "Origin", "Key")]

sample_df_c <- sample_df_c %>%
                filter(Treatment != "Blood&Tissue-AdaptiveSeq")

sample_df_c$Strain <- factor(sample_df_c$Strain, levels = rev(c("Porites", "Pocillopora", "Acropora", "F003", "H2")))
sample_df_c$Buffer <- factor(sample_df_c$Buffer, levels = rev(c("PBS", "DESS")))
sample_df_c$Origin <- factor(sample_df_c$Origin, levels = rev(c("Other", "aiptasiidae", "Corals", "Symbiodiniaceae", "Bacteria")))

# head(sample_df)
cols <- c("Symbiodiniaceae" = "#1b9e77", "Bacteria" = "#7570b3", "aiptasiidae" = "#FDBF6F", "Corals" = "#FDBF6F", "Other" = "#696969", "Overall" = "#e7298a")  # "Viruses" = "#e7298a", nolint: line_length_linter.


## Porites

Por_ont_count <- sample_df_c %>%
        filter(Strain == "Porites") %>%
        filter(DATA_Type == "FastQ-ONT") %>%
                filter(Key == "total_bases") %>%
        ggplot(aes(x = sample_short, y = Value, fill = Origin)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "Porites - reads ONT (bp count)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

Por_ont_count_reads <- sample_df_c %>%
    filter(Strain == "Porites") %>%
    filter(DATA_Type == "FastQ-ONT") %>%
    filter(Key != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Value, fill = Origin)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "Porites - reads ONT (reads count)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

Por_ont_percentage <- sample_df_c %>%
        filter(Strain == "Porites") %>%
        filter(DATA_Type == "FastQ-ONT") %>%
                filter(Key != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Total_percentage, fill = Origin)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "Porites - reads ONT (bp percentage)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# sample_data_l
sample_data_l$Strain <- factor(sample_data_l$Strain, levels = rev(c("Porites", "Pocillopora", "Acropora", "F003", "H2")))
sample_data_l$Buffer <- factor(sample_data_l$Buffer, levels = rev(c("PBS", "DESS")))
sample_data_l$Origin <- factor(sample_data_l$Origin, levels = rev(c("Overall", "Other", "aiptasiidae", "Corals", "Symbiodiniaceae", "Bacteria")))

Por_ont_length <- sample_data_l %>%
        filter(Strain == "Porites") %>%
        filter(DATA_Type == "FastQ-ONT") %>%
                filter(Key != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Value, fill = Origin)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "Porites - reads ONT (read length)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

ggarrange(Por_ont_percentage, Por_ont_count,Por_ont_count_reads, Por_ont_length, ncol = 1)


## Acropora

Acro_ont_count <- sample_df_c %>%
        filter(Strain == "Acropora") %>%
        filter(DATA_Type == "FastQ-ONT") %>%
                filter(Key == "total_bases") %>%
        ggplot(aes(x = sample_short, y = Value, fill = Origin)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "Acropora - reads ONT (bp count)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

Acro_ont_count_reads <- sample_df_c %>%
                filter(Strain == "Acropora") %>%
                filter(DATA_Type == "FastQ-ONT") %>%
                filter(Key != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Value, fill = Origin)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "Acropora - reads ONT (reads count)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

Acro_ont_percentage <- sample_df_c %>%
                filter(Strain == "Acropora") %>%
                filter(DATA_Type == "FastQ-ONT") %>%
                filter(Key != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Total_percentage, fill = Origin)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "Acropora - reads ONT (bp percentage)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# sample_data_l
sample_data_l$Strain <- factor(sample_data_l$Strain, levels = rev(c("Acropora", "Pocillopora", "Acropora", "F003", "H2")))
sample_data_l$Buffer <- factor(sample_data_l$Buffer, levels = rev(c("PBS", "DESS")))
sample_data_l$Origin <- factor(sample_data_l$Origin, levels = rev(c("Overall", "Other", "aiptasiidae", "Corals", "Symbiodiniaceae", "Bacteria")))

Acro_ont_length <- sample_data_l %>%
        filter(Strain == "Acropora") %>%
        filter(DATA_Type == "FastQ-ONT") %>%
                filter(Key != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Value, fill = Origin)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "Acropora - reads ONT (read length)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

ggarrange(Acro_ont_percentage, Acro_ont_count,Acro_ont_count_reads, Acro_ont_length, ncol = 1)


### Pocillopora

Pocil_ont_count <- sample_df_c %>%
        filter(Strain == "Pocillopora") %>%
        filter(DATA_Type == "FastQ-ONT") %>%
                filter(Key == "total_bases") %>%
        ggplot(aes(x = sample_short, y = Value, fill = Origin)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "Pocillopora - reads ONT (bp count)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

Pocil_ont_count_reads <- sample_df_c %>%
    filter(Strain == "Pocillopora") %>%
    filter(DATA_Type == "FastQ-ONT") %>%
    filter(Key != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Value, fill = Origin)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "Pocillopora - reads ONT (reads count)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

Pocil_ont_percentage <- sample_df_c %>%
        filter(Strain == "Pocillopora") %>%
        filter(DATA_Type == "FastQ-ONT") %>%
                filter(Key != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Total_percentage, fill = Origin)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "Pocillopora - reads ONT (bp percentage)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# sample_data_l
sample_data_l$Strain <- factor(sample_data_l$Strain, levels = rev(c("Pocillopora", "Pocillopora", "Pocillopora", "F003", "H2")))
sample_data_l$Buffer <- factor(sample_data_l$Buffer, levels = rev(c("PBS", "DESS")))
sample_data_l$Origin <- factor(sample_data_l$Origin, levels = rev(c("Overall", "Other", "aiptasiidae", "Corals", "Symbiodiniaceae", "Bacteria")))

Pocil_ont_length <- sample_data_l %>%
        filter(Strain == "Pocillopora") %>%
        filter(DATA_Type == "FastQ-ONT") %>%
                filter(Key != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Value, fill = Origin)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "Pocillopora - reads ONT (read length)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

ggarrange(Pocil_ont_percentage, Pocil_ont_count,Pocil_ont_count_reads, Pocil_ont_length, ncol = 1)

## H2

H2_ont_count <- sample_df_c %>%
        filter(Strain == "H2") %>%
        filter(DATA_Type == "FastQ-ONT") %>%
                filter(Key == "total_bases") %>%
        ggplot(aes(x = sample_short, y = Value, fill = Origin)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "H2 - reads ONT (bp count)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

H2_ont_count_reads <- sample_df_c %>%
    filter(Strain == "H2") %>%
    filter(DATA_Type == "FastQ-ONT") %>%
    filter(Key != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Value, fill = Origin)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "H2 - reads ONT (reads count)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

H2_ont_percentage <- sample_df_c %>%
        filter(Strain == "H2") %>%
        filter(DATA_Type == "FastQ-ONT") %>%
                filter(Key != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Total_percentage, fill = Origin)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "H2 - reads ONT (bp percentage)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# sample_data_l
sample_data_l$Strain <- factor(sample_data_l$Strain, levels = rev(c("H2", "H2", "H2", "F003", "H2")))
sample_data_l$Buffer <- factor(sample_data_l$Buffer, levels = rev(c("PBS", "DESS")))
sample_data_l$Origin <- factor(sample_data_l$Origin, levels = rev(c("Overall", "Other", "aiptasiidae", "Corals", "Symbiodiniaceae", "Bacteria")))

H2_ont_length <- sample_data_l %>%
        filter(Strain == "H2") %>%
        filter(DATA_Type == "FastQ-ONT") %>%
                filter(Key != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Value, fill = Origin)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "H2 - reads ONT (read length)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

ggarrange(H2_ont_percentage, H2_ont_count,H2_ont_count_reads, H2_ont_length, ncol = 1)


## F003

F003_ont_count <- sample_df_c %>%
        filter(Strain == "F003") %>%
        filter(DATA_Type == "FastQ-ONT") %>%
                filter(Key == "total_bases") %>%
        ggplot(aes(x = sample_short, y = Value, fill = Origin)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "F003 - reads ONT (bp count)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

F003_ont_count_reads <- sample_df_c %>%
    filter(Strain == "F003") %>%
    filter(DATA_Type == "FastQ-ONT") %>%
    filter(Key != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Value, fill = Origin)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "F003 - reads ONT (reads count)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

F003_ont_percentage <- sample_df_c %>%
        filter(Strain == "F003") %>%
        filter(DATA_Type == "FastQ-ONT") %>%
                filter(Key != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Total_percentage, fill = Origin)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "F003 - reads ONT (bp percentage)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# sample_data_l
sample_data_l$Strain <- factor(sample_data_l$Strain, levels = rev(c("F003", "F003", "F003", "F003", "F003")))
sample_data_l$Buffer <- factor(sample_data_l$Buffer, levels = rev(c("PBS", "DESS")))
sample_data_l$Origin <- factor(sample_data_l$Origin, levels = rev(c("Overall", "Other", "aiptasiidae", "Corals", "Symbiodiniaceae", "Bacteria")))

F003_ont_length <- sample_data_l %>%
        filter(Strain == "F003") %>%
        filter(DATA_Type == "FastQ-ONT") %>%
                filter(Key != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Value, fill = Origin)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "F003 - reads ONT (read length)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

ggarrange(F003_ont_percentage, F003_ont_count,F003_ont_count_reads, F003_ont_length, ncol = 1)



# Por_contigs <- sample_df %>%
#         filter(Strain == "Porites") %>%
#         filter(DATA_Type != "FastQ-ONT") %>%
#                 filter(Type != "total_bases") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         ggtitle( "Porites - contigs illumina") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# # Por_contigs

# Por_bp <- sample_df %>%
#         filter(Strain == "Porites") %>%
#         filter(DATA_Type == "FastQ-ONT") %>%
#                 filter(Type != "total_bases") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         ggtitle( "Porites - reads ONT (bp count)") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))


sample_data <- merge(data, metadata, by.y = "matching_index", by.x = "file_name")
data <- read.table("/home/colinl/metaG/Git/metaG_EukDepletion/results/illumina/bracken_reports/beta_div_matrix.txt", header = TRUE)
data <- data[,-1]

data
colnames(data) <- c("X1",  "X2",  "X3",  "X4",  "X5",  "X6",  "X7",  "X8",  "X9",
        "X10", "X11", "X12", "X13", "X14", "X15", "X16", "X17", "X18", "X19",
        "X20", "X21", "X22", "X23", "X24", "X25", "X26", "X27", "X28", "X29",
        "X30", "X31", "X32", "X33", "X34", "X35", "X36", "X37","X38", "X39",
        "X40", "X41", "X42", "X43", "X44", "X45", "X46", "X47", "X48", "X49",
        "X50", "X51", "X52", "X53", "X54")
# rownames(data)
# Fill in the upper triangle and mirror it to fill in the lower triangle
for (i in 1:length(data)) {
        for (j in i:n) {
        square_matrix[i, j] <- data[i, j]
        square_matrix[j, i] <-data[i, j]
        }
}
data_num <- matrix(as.numeric(square_matrix),    # Convert to numeric matrix
                        ncol = ncol(square_matrix))
data_num   

heatmap(data_num)
dist_matrix <- as.dist(data_num)
hist(dist_matrix)
# p_load("ggplot")
heatmap.2(as.dist, 
        trace = "none",   # Don't show dendrograms
        col = colorRampPalette(c("white", "blue"))(100),   # Choose a color palette
        key = TRUE,       # Show color key
        keysize = 1.0,    # Set color key size
        scale = "none",   # Don't scale the values
        margins = c(8, 8), # Set margins
        labRow = NA,      # Don't label rows
        labCol = NA)      # Don't label columns

heatmap(data.matrix(data))

heatmap(as.matrix(data_num), cluster_columns = FALSE, 
                column_names_gp = gpar(cex=0.8),
                row_hclust_width = unit(3, "cm"),
                show_row_names = TRUE)


ggcorr(data.matrix(data), method = c("everything", "pearson"))

data %>%
ggplot(aes(group1, group2)) +    # Create default ggplot2 heatmap
  geom_tile(aes(fill = values))

head(sample_df_c[sample_df_c$Origin == "Bacteria",])

# # ggarrange(Por_contigs, Por_bp, ncol = 1)

# Por_bp
heatmap(as.matrix(sample_df_c))

head(sample_df_c)
heatmap_df <- sample_df_c[, c("sample", "Total_percentage", "Origin")]

heatmap(data.matrix(sample_df_c))

### get bacteria and general stats as variable i.e per sample bacteria content, length, bp, count, buffer and other 

head(sample_data_l)
head(sample_data_2_c)
head(merge(sample_data_2_c, sample_data_l, by = "file_name.y", no.dups = T))

heatmap_df1 <- sample_data_2_c %>%
        filter(DATA_Type == "FastQ-ONT")

heatmap_df2 <- sample_data_l %>%
        filter(DATA_Type == "FastQ-ONT")

heatmap_df1$file_name.y == heatmap_df2$file_name.y

head(merge(heatmap_df1, heatmap_df2, by = "file_name.y", no.dups = T))

for (i in (1:length(heatmap_df1$sample))) {
        if (heatmap_df1[i,]$file_name.y == heatmap_df2[i,]$file_name.y) {
        print(heatmap_df2[i,]$Value)
        print(heatmap_df1[i,])
        }

}

i
i = 5
head(sample_data)

heatmap_df_b <- sample_data %>%
        filter(Origin == "Bacteria")
heatmap_df_c <- sample_data %>%
        filter(Origin != "Bacteria", Origin != "Symbiodiniaceae", Origin != "Other", Origin != "Overall")
heatmap_df_s <- sample_data %>%
        filter(Origin == "Symbiodiniaceae")

head(heatmap_df_b)
head(heatmap_df_c)
head(heatmap_df_s)


head(heatmap_df_b)

heatmap_df_b <- heatmap_df_b[,c("file_name","Origin","Key", "Value", "Strain", "Replicate", "Treatment","Buffer","DATA_Type", "Location")]

heatmap_df_b_l <- heatmap_df_b %>%
        filter(DATA_Type != "FastQ-ONT") %>%
        pivot_wider(names_from = Key, values_from = Value)

heatmap_df_b_l <- heatmap_df_b_l[, c("file_name", "Strain", "Treatment", "Buffer", "Location", "length_reads", "total_reads", "total_bases")]

heatmap(data.matrix(heatmap_df_b_l,rownames.force = FALSE), scale = "column", col = terrain.colors(256))


data.matrix(heatmap_df_b_l, rownames.force = TRUE)


# Sample data frame
data <- heatmap_df_b_l

# Columns to keep as non-numeric
non_numeric_cols <- c("file_name")

# Convert specified columns to factors
for (col in non_numeric_cols) {
  data[[col]] <- as.factor(data[[col]])
}

# Convert the entire data frame to matrix (including the factors)
matrix_with_factors <- data.matrix(data, rownames.force = FALSE)

# Convert factors back to their original values
for (col in non_numeric_cols) {
  matrix_with_factors[, col] <- as.character(data[[col]])
}

# Print the resulting matrix
head(matrix_with_factors)

heatmap(matrix_with_factors)
colnames(matrix_with_factors)


# Convert specified columns to factors
non_numeric_cols <- c("file_name")
for (col in non_numeric_cols) {
  data[[col]] <- as.factor(data[[col]])
}

# Convert the entire data frame to matrix (including the factors)
matrix_with_factors <- data.matrix(data)

# Set the row names of the matrix to the "file_name" column
rownames(matrix_with_factors) <- data$file_name

# Remove the "file_name" column from the matrix
matrix_with_factors <- matrix_with_factors[, -which(colnames(matrix_with_factors) == "file_name")]

# Print the resulting matrix
print(matrix_with_factors)


# Assuming matrix_with_factors contains non-numeric factors
# Convert the factors to numeric (if applicable)
numeric_matrix <- as.matrix(matrix_with_factors[, -1])  # Exclude the row names column

# Convert all matrix values to numeric
numeric_matrix <- as.numeric(numeric_matrix)

# Set the numeric values back to the matrix
matrix_with_numeric <- matrix(numeric_matrix, ncol = ncol(matrix_with_factors) - 1)

# Create a heatmap
heatmap(matrix_with_numeric)

p_load(GGally)
p_load(ggforce)


gall_test <- sample_data[,c("file_name", "Origin", "Key", "Value", "Strain", "Replicate", "Treatment", "Buffer", "DATA_Type", "Location")] %>%
        filter(DATA_Type == "FastQ-ONT") %>%
        pivot_wider(names_from = Key, values_from = Value)

gall_test$Origin <- gsub("Corals", "Host", gall_test$Origin)
gall_test$Origin <- gsub("aiptasiidae", "Host", gall_test$Origin)

gall_test2 <- gall_test[,c("file_name", "Origin", "Strain", "Treatment", "Buffer", "total_bases", "length_reads")] %>%
        pivot_wider(names_from = Origin, values_from = c(total_bases, length_reads))
names(gall_test)
gall_test$Strain <- as.factor(gall_test$Strain)
gall_test$Treatment <- as.factor(gall_test$Treatment)
gall_test$Buffer <- as.factor(gall_test$Buffer)
gall_test$Origin <- as.factor(gall_test$Origin)

# From the help page:
# data(heatmap_df_b_l, package = "reshape")
ggpairs(
        gall_test[,c("Strain","Origin", "Treatment", "total_bases", "total_reads", "length_reads")] %>%
                filter(Origin != "Overall") %>%
                filter(Treatment != "Blood&Tissue-AdaptiveSeq"),
        # ggally_T2[,c("Strain","Origin", "Treatment", "Total_percentage", "Key")],
        # g_df_3[,c("Strain","Origin", "Treatment", "Precentage", "Key")],
        # gall_test2[,c(2:4,5,6,9,7,10:11,14,12)],#("Strain", "Treatment", "length_reads", "total_reads","Origin")] #,
        upper = list(continuous = "cor", combo = "autopoint", discrete = "autopoint"),
        lower = list(continuous = "smooth", combo = "box", discrete = "box"),
        diag = list(continuous = "autopointDiag", combo = "autopointDiag", discrete = "autopointDiag"),
        ggplot2::aes(colour = Origin)
)

sample_data_2
ggcorr(gall_test2[,c(2:4,5,6,9,7,10:11,14,12)], method = c("everything", "pearson"))

ggally_T2 <- sample_df_c %>%
        filter(DATA_Type == "FastQ-ONT") %>%
        filter(Treatment != "Blood&Tissue-AdaptiveSeq")

ggally_T2$Origin <- gsub("Corals", "Host", ggally_T2$Origin)
ggally_T2$Origin <- gsub("aiptasiidae", "Host", ggally_T2$Origin)

g_df_3 <- ggally_T2 %>%
  group_by(Strain, Treatment, Buffer, Origin, Key, Total_percentage) %>%
  summarise(Precentage = mean(Total_percentage)) 

# now summarise to my variables
suppressMessages(bac_data <- bacteria_df %>%
  filter(class %in% top_spp$class) %>%
  group_by(queryID, strain, Treatment, X_names, class) %>%
  summarise(count = n()) %>%
  mutate(Frequency = count / sum(count)))


sample_data_bk <- sample_data

colnames(sample_data)
sample_data_test <- sample_data[,c("file_name", "Origin", "Key", "Value", "file_name.y", "Strain", "Replicate", "Treatment", "Buffer" , "DATA_Type", "Location", "sample" ,"sample_short")]
sample_data <- sample_data_test %>%
        pivot_wider(names_from = c(Origin, Key), values_from = Value)

sample_data$reads_coral <- sample_data$Corals_total_reads / sample_data$Overall_total_reads
sample_data$reads_bacteria <- sample_data$Bacteria_total_reads / sample_data$Overall_total_reads
sample_data$length_coral <- sample_data$Corals_length_reads / sample_data$Overall_length_reads
sample_data$length_bacteria <- sample_data$Bacteria_length_reads / sample_data$Overall_length_reads
#    %>%
#   pivot_longer(!religion, names_to = "income", values_to = "count")
test <- sample_data[, c("file_name", "file_name.y", "Strain", "Replicate", "Treatment", "Buffer" , "DATA_Type", "Location", "sample" ,"sample_short", "reads_coral", "reads_bacteria", "length_coral", "length_bacteria")]

test %>%
        pivot_longer(reads_coral:length_bacteria, #names_to = "Key",
        names_to = c("type", "origin"),
        names_pattern = "(.*)_(.*)",
        values_to = "ratio") %>%
ggplot(aes(type, origin)) +    # Create default ggplot2 heatmap
  geom_tile(aes(fill = ratio))


data_alphadiv_Sh <- read.table("/home/colinl/metaG/Git/metaG_EukDepletion/results/illumina/bracken_reports/alpha_diversity_Sh.txt", sep = "\t", col.names = c("Sample", "Shannon_div"))
data_alphadiv_Si <- read.table("/home/colinl/metaG/Git/metaG_EukDepletion/results/illumina/bracken_reports/alpha_diversity_Si.txt", sep = "\t", col.names = c("Sample", "Simpson_div"))
data_alphadiv_ISi <- read.table("/home/colinl/metaG/Git/metaG_EukDepletion/results/illumina/bracken_reports/alpha_diversity_ISi.txt", sep = "\t", col.names = c("Sample", "Simpson_reciprocal_div"))
data_alphadiv_F <- read.table("/home/colinl/metaG/Git/metaG_EukDepletion/results/illumina/bracken_reports/alpha_diversity_F.txt", sep = "\t", col.names = c("Sample", "Fisher_div"))
data_alphadiv_BP <- read.table("/home/colinl/metaG/Git/metaG_EukDepletion/results/illumina/bracken_reports/alpha_diversity_BP.txt", sep = "\t", col.names = c("Sample", "Berger_parker_div"))

plot_df_test <- cbind(data_alphadiv_Sh,data_alphadiv_Si, data_alphadiv_ISi, data_alphadiv_F, data_alphadiv_BP)
# dt_separate(data_alphadiv, Sample, into = )

# sample_data <- 
plot_df_test <- plot_df_test[,c("Sample","Shannon_div","Simpson_div","Simpson_reciprocal_div","Fisher_div","Berger_parker_div")]
plot_df_test <- merge(plot_df_test, metadata, by.y = "matching_index", by.x = "Sample")

plot_df_test$sample_lab <- paste(plot_df_test$Strain, plot_df_test$Replicate, sep = "_")
plot_df_test %>%
        # pivot_longer(reads_coral:length_bacteria, #names_to = "Key",
        # names_to = c("type", "origin"),
        # names_pattern = "(.*)_(.*)",
        # values_to = "ratio") %>%
ggplot(aes(sample_lab, Treatment)) +    # Create default ggplot2 heatmap
        geom_tile(aes(fill = Simpson_reciprocal_div)) +
        geom_text(aes(label = Simpson_reciprocal_div)) +
        scale_fill_gradient(low = "white", high = "#1b98e0")
