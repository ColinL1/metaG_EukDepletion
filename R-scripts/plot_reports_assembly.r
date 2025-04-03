#!/usr/bin/env Rscript

# echo status
cat("Loading libraries...\n")
lib_list <- c("tidyr", "tidyfast", "dplyr", "ggplot2", "ggpubr", "tidyverse", "scales", "RColorBrewer", "optparse", "pacman","httpgd") # nolint: line_length_linter, cSpell.
invisible(capture.output(
                suppressMessages(
                        lapply(lib_list, require, quietly = TRUE, character.only = TRUE))))

# #load tydifast with pacman (it installs it if not present)
# invisible(
#     capture.output(
#         suppressMessages(
#                 p_load("tidyfast"))))

# # Run "webplot" server (only necessary for live use)
hgd() # http://127.0.0.1:43617/live?token=GdpJ8LuW

# data <- read_csv("/home/colinl/metaG/Git/metaG_EukDepletion/results/ONT_non_trim_report_2.csv")  #### nanopore data not assembled
# data <- read_csv("/home/colinl/metaG/Git/metaG_EukDepletion/contigs_repots_kaiju.csv")
data <- read_csv("/home/colinl/metaG/Git/metaG_EukDepletion/results/illumina_V3/Kaiju_reads_corals/fastp/report.csv")

data <- data[,c(2:5)]
# data %>%ù
#         filter(file_name == "F003_00_MB_pass")

# sequencing_stats <- read.table("/home/colinl/metaG/Git/metaG_methods/results/metaG_indonesia/mmseqs2_reports/ONT_reads_NR_lca.tsv", header = TRUE) # nolint: line_length_linter.
metadata <- read_csv("/home/colinl/metaG/Git/metaG_EukDepletion/sample_metadata_sheet.csv")

sample_data <- merge(data, metadata, by.x = "file_name", by.y = "Index", all.x = TRUE)

sample_data["sample"] <- paste(sample_data$Strain, sample_data$Replicate, sample_data$From, sample_data$DATA_Type, sep = "-") 
sample_data["sample_short"] <- paste(sample_data$Replicate, sample_data$From, sample_data$DATA_Type, sep = "-") 

sample_data["sample_short"] <- paste(sample_data$Replicate, sample_data$From, sep = "-") 
sample_data["sample_short"] <- paste(sample_data$Strain, sample_data$Replicate, sample_data$Treatment, sample_data$Buffer, sep = "_") 
sample_data["sample_short"] <- paste(sample_data$Strain, sample_data$Replicate, sample_data$Buffer, sep = "-") 

sample_data["identifier"] <- paste(sample_data$Strain, sample_data$Replicate, sample_data$Treatment, sample_data$Buffer, sample_data$DATA_Type, sample_data$From, sample_data$Type, sep = "_") 

sample_data <- sample_data %>% 
                        filter(Type == "total_reads")

overall_numebrs <- sample_data %>%
        filter(Kingdom == "TOTAL")
sample_data_clean <- sample_data %>%
        filter(Kingdom != "TOTAL")

overall_numebrs["sequenced_total"] <- overall_numebrs$total

overall_numebrs <- overall_numebrs[,c("identifier", "sequenced_total")]

sample_data_2 <- merge(sample_data_clean, overall_numebrs, by = "identifier")

sample_data_2["Total_percentage"] <- (as.numeric(sample_data_2$total)/as.numeric(sample_data_2$sequenced_total))*100

sample_df <- sample_data_2[,c("file_name", "Strain", "Replicate", "Treatment", "Buffer", "DATA_Type", "From", "sample", "sample_short", "total", "sequenced_total", "Total_percentage", "Kingdom", "Type")]

sample_df <- sample_df %>%
                filter(Strain != "apo-H2") %>%
                filter(Treatment != "Blood&Tissue-filtered") %>%
                filter(Treatment != "Microbiome-filtereds") %>%
                filter(Treatment != "Blood&Tissue-AdaptiveSeq")

sample_df$Strain <- factor(sample_df$Strain, levels = rev(c("Porites", "Pocillopora", "Acropora", "F003", "H2")))
sample_df$Buffer <- factor(sample_df$Buffer, levels = rev(c("PBS", "DESS")))
sample_df$Kingdom <- factor(sample_df$Kingdom, levels = rev(c("Other", "aiptasiidae", "Corals", "Symbiodiniaceae", "Bacteria")))

# head(sample_df)
cols <- c("Symbiodiniaceae" = "#1b9e77", "Bacteria" = "#7570b3", "aiptasiidae" = "#FDBF6F", "Corals" = "#FDBF6F", "Other" = "#696969")  # "Viruses" = "#e7298a", nolint: line_length_linter.

for (i in unique(sample_df$Strain)) {
        print(paste("Plotting:", i, sep = " "))
        title <- paste(i, "Reads illumina", sep = " - ")
        print(sample_df %>%
                filter(Strain == i) %>%
                ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
                geom_bar(stat = "identity") + # , position="fill") +
                scale_fill_manual(values = cols) +
                theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
                facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
                theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
                strip.placement = "top", # move strip to outside
                strip.background = element_blank(), # remove strip background
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
                ggtitle(title) +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE)))

        ggsave(paste(i, "Reads illumina.svg", sep = " - "), device = "svg")
        ggsave(paste(i, "Reads illumina.png", sep = " - "), device = "png")

}


for (i in unique(sample_df$Strain)) {
        print(paste("Plotting:", i, sep = " "))

        print(ggarrange(sample_df %>%
                filter(Strain == i) %>%
                ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
                geom_bar(stat = "identity") + # , position="fill") +
                scale_fill_manual(values = cols) +
                theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
                facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
                theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
                strip.placement = "top", # move strip to outside
                strip.background = element_blank(), # remove strip background
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
                ggtitle(paste(i, "Reads illumina", sep = " - ")) +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE)), 
                
                sample_df %>%       
                filter(Strain == i) %>%
                filter(Kingdom != "Other") %>%
                ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
                geom_bar(stat = "identity", position = "fill") +
                scale_fill_manual(values = cols) +
                theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
                facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
                theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
                strip.placement = "top", # move strip to outside
                strip.background = element_blank(), # remove strip background
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
                ggtitle(paste(i, "Reads illumina (no other)", sep = " - ")) +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE)), ncol = 1))
        ggsave(paste(i, "Reads illumina.svg", sep = " - "), device = "svg")
        ggsave(paste(i, "Reads illumina.png", sep = " - "), device = "png")
}

## Porties # contigs

# por_contigs_df <- sample_df %>%
#         # filter(Treatment != "Blood&Tissue-filtered",
#         #         Treatment != "Microbiome-filtereds",
#         #         Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "Porites") %>%
#         filter(DATA_Type != "FastQ-ONT") %>%
#                 filter(Type != "total_bases") 

Por_contigs <- sample_df %>%
        filter(Strain == "Porites") %>%
        # filter(DATA_Type != "FastQ-ONT") %>%
                # filter(Type != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
        geom_bar(stat = "identity") + # , position="fill") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "Porites - Reads illumina") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# Por_contigs

Por_bp <- sample_df %>%
        filter(Strain == "Porites") %>%
        # filter(DATA_Type == "FastQ-ONT") %>%
                filter(Kingdom != "Other") %>%
        ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
        geom_bar(stat = "identity", position="fill") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "Porites - Reads illumina (no other)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# Por_bp
ggarrange(Por_contigs, Por_bp, ncol = 1)



### acropora 

Acro_contigs <- sample_df %>%
        filter(Strain == "Acropora") %>%
        # filter(DATA_Type != "FastQ-ONT") %>%
                # filter(Type != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
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
        ggtitle( "Acropora - Reads illumina") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# Acro_contigs

Acro_bp <- sample_df %>%
        filter(Strain == "Acropora") %>%
        # filter(DATA_Type == "FastQ-ONT") %>%
                # filter(Type != "total_bases") %>%
                filter(Kingdom != "Other") %>%
        ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
        geom_bar(stat = "identity", position="fill") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "Acropora - Reads illumina (no other)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# Acro_bp
ggarrange(Acro_contigs, Acro_bp, ncol = 1)

### pocillopora

Poci_contigs <- sample_df %>%
        filter(Strain == "Pocillopora") %>%
        # filter(DATA_Type != "FastQ-ONT") %>%
                # filter(Type != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
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
        ggtitle( "Pocillopora - Reads illumina") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# Poci_contigs

Poci_bp <- sample_df %>%
                filter(Strain == "Pocillopora") %>%
        # filter(DATA_Type == "FastQ-ONT") %>%
                # filter(Type != "total_bases") %>%
                filter(Kingdom != "Other") %>%
        ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
        geom_bar(stat = "identity", position="fill") +
        scale_fill_manual(values = cols) +
        theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        strip.placement = "top", # move strip to outside
        strip.background = element_blank(), # remove strip background
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        ggtitle( "Pocillopora - - Reads illumina (no other)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# Poci_bp
ggarrange(Poci_contigs, Poci_bp, ncol = 1)



### H2

H2_contigs <- sample_df %>%
        filter(Strain == "H2") %>%
        # filter(DATA_Type != "FastQ-ONT") %>%
                filter(Type != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
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
        ggtitle( "H2 - contigs illumina") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# H2_contigs

H2_bp <- sample_df %>%
        filter(Strain == "H2") %>%
        # filter(DATA_Type == "FastQ-ONT") %>%
                filter(Type != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
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
        ggtitle( "H2 - assembly ONT") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# H2_bp
# ggarrange(H2_contigs, H2_bp, ncol = 1)

### F003

F003_contigs <- sample_df %>%
        filter(Strain == "F003") %>%
        # filter(DATA_Type != "FastQ-ONT") %>%
                filter(Type != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
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
        ggtitle( "F003 - contigs illumina") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# F003_contigs

F003_bp <- sample_df %>%
        filter(Strain == "F003") %>%
        # filter(DATA_Type == "FastQ-ONT") %>%
                filter(Type != "total_bases") %>%
        ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
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
        ggtitle( "F003 - assembly ONT") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# F003_bp
# ggarrange(F003_contigs, F003_bp, ncol = 1)

pdf("plot_coral_kaiju_preliminary.pdf", paper = "a4r")
ggarrange(Por_contigs, Por_bp, ncol = 1)
ggarrange(Acro_contigs, Acro_bp, ncol = 1)
ggarrange(Poci_contigs, Poci_bp, ncol = 1)

# ggarrange(H2_contigs, H2_bp, ncol = 1)
# ggarrange(F003_contigs, F003_bp, ncol = 1)

dev.off()

# sample_df %>%
#         filter(file_name == "F003_00_MB_pass")
# # cols <- c("Eukaryota" = "#1b9e77", "Bacteria" = "#7570b3", "Archaea" = "#FDBF6F", "Viruses" = "#e7298a", "Unclassified" = "#696969", "Unknown" = "#353535") # nolint: line_length_linter.
# cols <- c("Symbiodiniaceae" = "#1b9e77", "Bacteria" = "#7570b3", "aiptasiidae" = "#FDBF6F", "Corals" = "#FDBF6F", "Other" = "#696969")  # "Viruses" = "#e7298a", nolint: line_length_linter.

# sample_df_plot <- sample_df %>%
#         filter(Type != "total_bases")
#         # filter(Type == "total_reads")

# sample_df_plot %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "Porites") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 3, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         ggtitle( "Porites - reads data bp count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))


# sample_df_plot %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "Pocillopora") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 2, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         ggtitle( "Pocillopora - reads data bp count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# sample_df_plot %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "Acropora") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 2, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         ggtitle( "Acropora - reads data bp count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# # 
# sample_df_plot %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter( Strain == "H2") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 4, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         ggtitle("H2 - reads data bp count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# sample_df_plot %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "F003") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 3, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         ggtitle( "F003 - reads data bp count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# sample_df_plot <- sample_df %>%
#         # filter(Type != "total_bases")
#         filter(Type == "total_reads")

# sample_df_plot %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "Porites") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 3, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         ggtitle( "Porites - reads data read count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))


# sample_df_plot %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "Pocillopora") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 2, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         ggtitle( "Pocillopora - reads data read count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# sample_df_plot %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "Acropora") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 2, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         ggtitle( "Acropora - reads data read count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# # aiptasiidae
# sample_df_plot %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter( Strain == "H2") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 4, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         ggtitle("H2 - reads data read count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# sample_df_plot %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "F003") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 3, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         ggtitle( "F003 - reads data read count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))


# # sample_df_plot %>%
# #         filter(Treatment != "Blood&Tissue-filtered",
# #                 Treatment != "Microbiome-filtereds",
# #                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
# #         filter(Strain == "F003") %>%
# #         filter(sample == "F003-0-LAB-FastQ-ONT")
# dev.off()





# # data <- read_csv("/home/colinl/metaG/Git/metaG_EukDepletion/results/ONT_report.csv")
# # data <- read_csv("/home/colinl/metaG/Git/metaG_EukDepletion/results/contigs_report.csv")

# # data <- read_csv("/home/colinl/metaG/Git/metaG_EukDepletion/results/ONT_non_trim_report.csv")
# data <- read_csv("/home/colinl/metaG/Git/metaG_EukDepletion/results/ONT_contigs_non_trim_report.csv")

# data <- data[,c(2:5)]

# data %>%
#         filter(file_name == "F003_00_MB_pass")
# # sequencing_stats <- read.table("/home/colinl/metaG/Git/metaG_methods/results/metaG_indonesia/mmseqs2_reports/ONT_reads_NR_lca.tsv", header = TRUE) # nolint: line_length_linter.
# metadata <- read_csv("/home/colinl/metaG/Git/metaG_EukDepletion/sample_metadata_sheet.csv")

# # data$file_name <- gsub("ID_BST_Por_06","Por_BST-6_MBX", data$file_name)
# # data$file_name <- gsub("ID_GT_Por_06","Por_GT-6_MBX", data$file_name)
# # data$file_name <- gsub("ID_GT_Por_03","Por_GT-3_MBX", data$file_name)
# # data$file_name <- gsub("ID_GT_Por_01","Por_GT-1_MBX", data$file_name)
# # data$file_name <- gsub("ID_BST_Por_04" ,"Por_BST-4_MBX", data$file_name)
# # data$file_name <- gsub("ID_BST_Por_03","Por_BST-3_MBX", data$file_name)
# # data$file_name <- gsub("ID_BST_Por_01","Por_BST-1_MBX", data$file_name)
# # data$file_name <- gsub("ID_BST_Por_07","Por_BST-7_MBX", data$file_name)
# # data$file_name <- gsub("ID_GT_Por_09" ,"Por_GT-9_MBX", data$file_name)
# # data$file_name <- gsub("ID_GT_Por_04","Por_GT-4_MBX", data$file_name)
# # data$file_name <- gsub("ID_BST_Por_05" ,"Por_BST-5_MBX", data$file_name)

# # data <- separate(data, file_name, into = c("strain", "replicate", "Treatment","Buffer"), sep = "_", remove = FALSE) # nolint: line_length_linter.
# # ### clean up treatment names ###
# # unique(data$strain)
# # unique(data$replicate)
# # unique(data$Treatment)
# # unique(data$Buffer)
# # ### Blood and tissue kit ###
# # data$Treatment <- gsub("BTa", "Blood & Tissue_Aposymbiotic", data$Treatment)
# # data$Treatment <- gsub("BTf", "Blood & Tissue_Filtered", data$Treatment)
# # data$Treatment <- gsub("BTAD", "Blood & Tissue_AdaptiveSeq", data$Treatment)
# # data$Treatment <- gsub("BT", "Blood & Tissue_Untreated", data$Treatment)
# # ### Microbiome kit ###
# # data$Treatment <- gsub("MBa", "Microbiome_Aposymbiotic", data$Treatment)
# # data$Treatment <- gsub("MBf", "Microbiome_Filtered", data$Treatment)
# # data$Treatment <- gsub("MBX", "Modified Microbiome_Untreated", data$Treatment)
# # data$Treatment <- gsub("MB", "Microbiome_Untreated", data$Treatment)

# # sample$queryID <- str_extract(sample$queryID, "(?:[^_]*_){4}")

# # str_extract(metadata$file_name, "(?:[^_]*_){4}")
# # data$file_name
# # length(data$file_name)
# # length(sample_data$file_name)

# sample_data <- merge(data, metadata, by.y = "matching_index", all.x = TRUE)
# sample_data %>%
#         filter(file_name == "F003_00_MB_pass")
# # table(is.na(sample_data))

# # head(sample_data)

# ### split by analysis type bp vs reads count 


# sample_data["sample"] <- paste(sample_data$Strain, sample_data$Replicate, sample_data$From, sample_data$DATA_Type, sep = "-") 
# sample_data["sample_short"] <- paste(sample_data$Replicate, sample_data$From, sample_data$DATA_Type, sep = "-") 

# sample_data["identifier"] <- paste(sample_data$Strain, sample_data$Replicate, sample_data$Treatment, sample_data$Buffer, sample_data$DATA_Type, sample_data$From, sample_data$Type, sep = "_") 

# head(sample_data)
# length(filter(sample_data, Kingdom == "TOTAL")$file_name)
# sample_data$Kingdom

# overall_numebrs <- sample_data %>%
#         filter(Kingdom == "TOTAL")
# sample_data_clean <- sample_data %>%
#         filter(Kingdom != "TOTAL")

# overall_numebrs %>%
#         filter(file_name == "F003_00_MB_pass")

# sample_data_clean %>%
#         filter(file_name == "F003_00_MB_pass")


# head(overall_numebrs)


# overall_numebrs["sequenced_total"] <- overall_numebrs$total

# overall_numebrs <- overall_numebrs[,c("identifier", "sequenced_total")]

# sample_data_2 <- merge(sample_data_clean, overall_numebrs, by = "identifier")

# sample_data_2 %>%
#         filter(file_name == "F003_00_MB_pass")

# # head(sample_data_2)
# # table(is.na(sample_data_2))
# # length(sample_data$file_name)
# # length(sample_data_clean$file_name)
# # length(sample_data_2$file_name)


# # overall_numebrs %>% 
# #         filter(identifier == "Acropora_14_Modified-microbiome_DESS_FastQ-ONT_Kaust_total_reads")

# # sample_data %>% 
# #         filter(identifier == "Acropora_14_Modified-microbiome_DESS_FastQ-ONT_Kaust_total_reads")

# # sample_data_2 %>% 
# #         filter(identifier == "Acropora_14_Modified-microbiome_DESS_FastQ-ONT_Kaust_total_reads")

# # sample_data_2 %>% 
# #         filter(identifier == "Acropora_14_Modified-microbiome_DESS_FastQ-ONT_Kaust_total_bases")

# # sample_data %>% 
# #         filter(identifier == "Acropora_14_Modified-microbiome_DESS_FastQ-ONT_Kaust_total_bases")

# # overall_numebrs %>% 
# #         filter(identifier == "Acropora_14_Modified-microbiome_DESS_FastQ-ONT_Kaust_total_bases")


# sample_data_2["Total_percentage"] <- (sample_data_2$total/sample_data_2$sequenced_total)*100
# # # 
# # sample_data["sample"] <- paste(sample_data$Strain, sample_data$Replicate, sample_data$From, sep = "_") 


# sample_df <- sample_data_2[,c("file_name", "Strain", "Replicate", "Treatment", "Buffer", "DATA_Type", "From", "sample", "sample_short", "total", "sequenced_total", "Total_percentage", "Kingdom", "Type")]

# sample_df <- sample_df %>%
#                 filter(Strain != "apo-H2") %>%
#                 filter(Treatment != "Blood&Tissue-filtered") %>%
#                 filter(Treatment != "Microbiome-filtereds") %>%
#                 filter(Treatment != "Blood&Tissue-AdaptiveSeq")
# #                 #"Blood&Tissue-filtered"    "Microbiome-filtereds"     "Blood&Tissue-AdaptiveSeq"
# # unique(sample_df$Strain)
# # unique(sample_df$Treatment)

# sample_df$Strain <- factor(sample_df$Strain, levels = rev(c("Porites", "Pocillopora", "Acropora", "F003", "H2")))

# sample_df %>%
#         filter(file_name == "F003_00_MB_pass")

# # head(sample_df)

# # sample_df %>%
# #     filter( Type != "total_bases") %>%
# #         ggplot(aes(x = file_name, y = Total_percentage, fill = Kingdom)) +
# #         geom_bar(stat = "identity") + 
# #         facet_wrap(Strain ~ Treatment + Buffer, scale = "free")

# # sample_df %>%
# #     filter( Type != "total_bases") %>%
# #     filter( Kingdom == "Bacteria") %>%
# #     filter( Strain == "Porites") %>%
# #     filter( DATA_Type != "Fasta-Contigs") %>%
# #         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
# #         geom_bar(stat = "identity") +
# #         facet_grid( Treatment ~ Strain + Buffer, scale = "free_x")

# # sample_df %>%
# #     filter( Type == "total_reads") %>%
# #     filter( Kingdom == "Bacteria") %>%
# #     filter( Strain == "Porites") %>%
# #     filter( DATA_Type != "Fasta-Contigs") %>%
# #         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
# #         geom_bar(stat = "identity") +
# #         facet_grid( Treatment ~ Strain + Buffer, scale = "free_x")

# # sample_df %>%
# #     # filter( Type != "total_bases") %>%
# #     # filter( Kingdom == "Bacteria") %>%
# #     filter(Treatment != "Blood&Tissue-filtered", Treatment != "Microbiome-filtereds", Treatment != "Blood&Tissue-AdaptiveSeq") %>%
# #     filter( Strain == "Porites") %>%
# #     filter( Buffer == "PBS") %>%
# #     filter( Strain != "apo-H2") %>%
# #     filter( DATA_Type != "Fasta-Contigs") %>%
# #         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
# #         geom_bar(stat = "identity") +
# #         facet_grid(Treatment ~ Type + Strain + Buffer, scale = "free_x")
# # #"Modified-microbiome", "Benzonase", "Blood&Tissue", "Microbiome", "Blood&Tissue-filtered", "Microbiome-filtereds", "Blood&Tissue-AdaptiveSeq"


# # sample_df %>%
# #     # filter( Type != "total_bases") %>%
# #     # filter( Kingdom == "Bacteria") %>%
# #     filter(Treatment != "Blood&Tissue-filtered",
# #         Treatment != "Microbiome-filtereds",
# #         Treatment != "Blood&Tissue-AdaptiveSeq") %>%
# #     # filter(Strain == "Porites") %>%
# #     # filter(Buffer != "PBS") %>%
# #     filter(Strain != "apo-H2") %>%
# #     filter(DATA_Type != "Fasta-Contigs") %>%
# #         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
# #         geom_bar(stat = "identity") +
# #         facet_grid(Treatment ~ Strain+ Type  + Buffer, scale = "free_x")

# # #"Modified-microbiome", "Benzonase", "Blood&Tissue", "Microbiome", "Blood&Tissue-filtered", "Microbiome-filtereds", "Blood&Tissue-AdaptiveSeq"

# # coral_df <- sample_df %>%
# #         filter(Strain != "apo-H2", Strain != "H2", Strain != "F003")

# # aip_df <- sample_df %>%
# #         filter(Strain != "Porites", Strain != "Acropora", Strain != "Pocillopora")

# # cols <- c("Eukaryota" = "#1b9e77", "Bacteria" = "#7570b3", "Archaea" = "#FDBF6F", "Viruses" = "#e7298a", "Unclassified" = "#696969", "Unknown" = "#353535") # nolint: line_length_linter.
# cols <- c("Symbiodiniaceae" = "#1b9e77", "Bacteria" = "#7570b3", "aiptasiidae" = "#FDBF6F", "Corals" = "#FDBF6F", "Other" = "#696969")  # "Viruses" = "#e7298a", nolint: line_length_linter.

# sample_df_plot <- sample_df %>%
#         filter(Type != "total_bases")
#         # filter(Type == "total_reads")
# sample_df_plot%>%
#         filter(file_name == "F003_00_MB_pass")

# # pdf("test001.pdf", paper = "a4r")
# # # pdf("Plot%03d.pdf")

# sample_df_plot %>%
#         # filter( Kingdom == "Bacteria") %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "Porites") %>%
#         # filter(Buffer != "PBS") %>%
#         # filter(DATA_Type == "Fasta-Contigs") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 3, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         # strip.text.x = element_blank(), # remove x-axis label for strip
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         # # axis.ticks.y = element_blank(),
#         # strip.text.y = element_text(angle = -90, hjust = 0.5, vjust = 1))
#         ggtitle( "Porites - contigs data bp count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))


# sample_df_plot %>%
#         # filter( Kingdom == "Bacteria") %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "Pocillopora") %>%
#         # filter(Buffer != "PBS") %>%
#         # filter(DATA_Type == "Fasta-Contigs") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 2, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         # strip.text.x = element_blank(), # remove x-axis label for strip
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         # # axis.ticks.y = element_blank(),
#         # strip.text.y = element_text(angle = -90, hjust = 0.5, vjust = 1))
#         ggtitle( "Pocillopora - contigs data bp count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# sample_df_plot %>%
#         # filter( Kingdom == "Bacteria") %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "Acropora") %>%
#         # filter(Buffer != "PBS") %>%
#         # filter(DATA_Type == "Fasta-Contigs") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 2, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         # strip.text.x = element_blank(), # remove x-axis label for strip
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         # # axis.ticks.y = element_blank(),
#         # strip.text.y = element_text(angle = -90, hjust = 0.5, vjust = 1))
#         ggtitle( "Acropora - contigs data bp count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# # aiptasiidae
# sample_df_plot %>%
#         # filter( Kingdom == "Bacteria") %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter( Strain == "H2") %>%
#         # filter(Buffer != "PBS") %>%
#         # filter(DATA_Type == "Fasta-Contigs") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 4, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         # strip.text.x = element_blank(), # remove x-axis label for strip
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         # # axis.ticks.y = element_blank(),
#         # strip.text.y = element_text(angle = -90, hjust = 0.5, vjust = 1))
#         ggtitle("H2 - contigs data bp count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# sample_df_plot %>%
#         # filter( Kingdom == "Bacteria") %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "F003") %>%
#         # filter(Buffer != "PBS") %>%
#         # filter(DATA_Type == "Fasta-Contigs") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 3, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         # strip.text.x = element_blank(), # remove x-axis label for strip
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         # # axis.ticks.y = element_blank(),
#         # strip.text.y = element_text(angle = -90, hjust = 0.5, vjust = 1))
#         ggtitle( "F003 - contigs data bp count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# sample_df_plot <- sample_df %>%
#         # filter(Type != "total_bases")
#         filter(Type == "total_reads")

# sample_df_plot %>%
#         # filter( Kingdom == "Bacteria") %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "Porites") %>%
#         # filter(Buffer != "PBS") %>%
#         # filter(DATA_Type == "Fasta-Contigs") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 3, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         # strip.text.x = element_blank(), # remove x-axis label for strip
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         # # axis.ticks.y = element_blank(),
#         # strip.text.y = element_text(angle = -90, hjust = 0.5, vjust = 1))
#         ggtitle( "Porites - contigs data read count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))


# sample_df_plot %>%
#         # filter( Kingdom == "Bacteria") %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "Pocillopora") %>%
#         # filter(Buffer != "PBS") %>%
#         # filter(DATA_Type == "Fasta-Contigs") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 2, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         # strip.text.x = element_blank(), # remove x-axis label for strip
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         # # axis.ticks.y = element_blank(),
#         # strip.text.y = element_text(angle = -90, hjust = 0.5, vjust = 1))
#         ggtitle( "Pocillopora - contigs data read count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# sample_df_plot %>%
#         # filter( Kingdom == "Bacteria") %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "Acropora") %>%
#         # filter(Buffer != "PBS") %>%
#         # filter(DATA_Type == "Fasta-Contigs") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 2, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         # strip.text.x = element_blank(), # remove x-axis label for strip
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         # # axis.ticks.y = element_blank(),
#         # strip.text.y = element_text(angle = -90, hjust = 0.5, vjust = 1))
#         ggtitle( "Acropora - contigs data read count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))



# # aiptasiidae
# sample_df_plot %>%
#         # filter( Kingdom == "Bacteria") %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter( Strain == "H2") %>%
#         # filter(Buffer != "PBS") %>%
#         # filter(DATA_Type == "Fasta-Contigs") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 4, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         # strip.text.x = element_blank(), # remove x-axis label for strip
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         # # axis.ticks.y = element_blank(),
#         # strip.text.y = element_text(angle = -90, hjust = 0.5, vjust = 1))
#         ggtitle("H2 - contigs data read count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

# sample_df_plot %>%
#         # filter( Kingdom == "Bacteria") %>%
#         filter(Treatment != "Blood&Tissue-filtered",
#                 Treatment != "Microbiome-filtereds",
#                 Treatment != "Blood&Tissue-AdaptiveSeq") %>%
#         filter(Strain == "F003") %>%
#         # filter(Buffer != "PBS") %>%
#         # filter(DATA_Type == "Fasta-Contigs") %>%
#         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values = cols) +
#         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
#         facet_wrap(Treatment ~  Buffer, ncol = 3, scale = "free_x") +
#         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
#         strip.placement = "top", # move strip to outside
#         strip.background = element_blank(), # remove strip background
#         # strip.text.x = element_blank(), # remove x-axis label for strip
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
#         # # axis.ticks.y = element_blank(),
#         # strip.text.y = element_text(angle = -90, hjust = 0.5, vjust = 1))
#         ggtitle( "F003 - contigs data read count") +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))


read.table("/home/colinl/NCCT/work/d7/3af6be92dec8ccbf4713e6333cf61e/metadata_ampliseq_2.txt", sep="\t", header=TRUE, row.names=1)