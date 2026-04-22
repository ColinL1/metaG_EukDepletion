#!/usr/bin/env Rscript

# echo status
cat("Loading libraries...\n")
lib_list <- c("tidyr", "tidyfast", "dplyr", "ggplot2", "ggpubr", "tidyverse", "scales", "RColorBrewer", "optparse", "pacman","httpgd", "grid") # nolint: line_length_linter, cSpell.
invisible(capture.output(
                suppressMessages(
                        lapply(lib_list, require, quietly = TRUE, character.only = TRUE))))

option_list = list(
        make_option(c("-f", "--file"), type = "character", default=NULL,
                help = "parsed fastp report as csv table", metavar = "character"),
        make_option(c("-m", "--metadata"), type = "character", default = NULL,
                help = "metadata csv file. Necessary to do plots.", metavar = "character"),
        make_option(c("-r", "--read_assembly"), type = "character", default = "reads",
                help = 'string indicating type of sequencing input. Either "reads" or "assembly', metavar = "character"),
        make_option(c("-s", "--seq_type"), type = "character", default = "illumina",
                help = 'string indicating type of sequencing platform. Either "illumina" or "ONT', metavar = "character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

####if missing file
if (is.null(opt$file)){
        print_help(opt_parser)
        stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

if (is.null(opt$metadata)){
        print_help(opt_parser)
        stop("missing metadata.n", call.=FALSE)
}

# hgd() #http://127.0.0.1:48805/live?token=HMuXZuDC
# opt$file <- "/home/colinl/metaG/Git/metaG_EukDepletion/work/e2/38d9a233feb95b408a332028b26824/report_2.csv"
# opt$metadata <- "/home/colinl/metaG/Git/metaG_EukDepletion/sample_metadata_sheet.csv"
invisible(capture.output(suppressMessages(data <- read_csv(opt$file))))
invisible(capture.output(suppressMessages(data <- read_csv("/home/colinl/metaG/Git/metaG_EukDepletion/work/21/733dac022ca9fc301fe429c5d6c400/report.csv"))))
data <- data[,c(2:5)]

invisible(capture.output(suppressMessages(metadata <- read_csv(opt$metadata))))
invisible(capture.output(suppressMessages(metadata <- read_csv("/home/colinl/metaG/Git/metaG_EukDepletion/sample_metadata_sheet.csv"))))

sample_data <- merge(data, metadata, by.x = "file_name", by.y = "Index", all.x = TRUE)

# sample_data<- sample_data %>%
#         pivot_wider(names_from = Type, values_from = total)

sample_data["sample"] <- paste(sample_data$Strain, sample_data$Replicate, sample_data$From, sample_data$DATA_Type, sep = "-") 
# sample_data["sample_short"] <- paste(sample_data$Replicate, sample_data$From, sample_data$DATA_Type, sep = "-") 

# sample_data["sample_short"] <- paste(sample_data$Replicate, sample_data$From, sep = "-") 
# sample_data["sample_short"] <- paste(sample_data$Strain, sample_data$Replicate, sample_data$Treatment, sample_data$Buffer, sep = "_") 
sample_data["sample_short"] <- paste(sample_data$Strain, sample_data$Replicate, sample_data$Buffer, sep = "-") 

sample_data["identifier"] <- paste(sample_data$Strain, sample_data$Replicate, sample_data$Treatment, sample_data$Buffer, sample_data$DATA_Type, sample_data$From, sample_data$Type, sep = "_") 

# sample_data <- sample_data %>% 
#                         filter(Type == "total_bases")

# overall_numebrs <- sample_data %>%
#         filter(Kingdom == "TOTAL")
# sample_data_clean <- sample_data %>%
#         filter(Kingdom != "TOTAL")

# overall_numebrs <- sample_data_clean %>% 
#         pivot_wider ( names_from = Kingdom, values_from = total)

# overall_numebrs["sequenced_total"] <- overall_numebrs$Symbiodiniaceae + overall_numebrs$Bacteria + overall_numebrs$aiptasiidae + overall_numebrs$Other

# overall_numebrs <- overall_numebrs[,c("identifier", "sequenced_total")]

# sample_data_2 <- merge(sample_data_clean, overall_numebrs, by = "identifier")

# sample_data_2["Total_percentage"] <- (as.numeric(sample_data_2$total)/as.numeric(sample_data_2$sequenced_total))*100

# sample_df <- sample_data_2[,c("file_name", "Strain", "Replicate", "Treatment", "Buffer", "DATA_Type", "From", "sample", "sample_short", "total", "sequenced_total", "Total_percentage", "Kingdom", "Type")]

sample_df <- sample_data[,c("file_name", "Strain", "Replicate", "Treatment", "Buffer", "DATA_Type", "From", "sample", "sample_short", "total", "Kingdom", "Type")]

# sample_df <- sample_df %>%
#                 filter(Strain != "apo-H2") %>%
#                 filter(Kingdom != TOTAL)
#                 filter(Treatment != "Blood&Tissue-filtered") %>%
#                 filter(Treatment != "Microbiome-filtereds") %>%
#                 filter(Treatment != "Blood&Tissue-AdaptiveSeq")

sample_df$Strain <- factor(sample_df$Strain, levels = rev(c("Porites", "Pocillopora", "Acropora", "F003", "H2")))
sample_df$Buffer <- factor(sample_df$Buffer, levels = rev(c("PBS", "DESS")))
sample_df$Kingdom <- factor(sample_df$Kingdom, levels = rev(c("Other", "Aiptasiidae", "Scleractina", "Symbiodiniaceae", "Bacteria", "TOTAL", "Trim")))

# head(sample_df)
cols <- c("Symbiodiniaceae" = "#1b9e77", "Bacteria" = "#7570b3", "aiptasiidae" = "#FDBF6F", "Scleractina" = "#FDBF6F", "Other" = "#696969")#,"TOTAL" = '#696969', "Trim" = '#e7298a')  # "Viruses" = "#e7298a", nolint: line_length_linter.

for (i in unique(sample_df$Strain)) {
        print(paste("Plotting:", i, sep = " "))
        mapping_df <- sample_df %>%
                filter(Strain == i) %>%
                filter(Type == "total_reads") %>%
                filter(Kingdom != "TOTAL", Kingdom != "Trim")
        
        stats_df <- sample_df %>%
                filter(Strain == i) %>%
                filter(Type == "total_reads") %>%
                filter(Kingdom != "Other", Kingdom != "Aiptasiidae", Kingdom != "Scleractina", Kingdom != "Symbiodiniaceae", Kingdom != "Bacteria")

        reads_plot <- mapping_df %>%
                ggplot(aes(x = sample_short, y = total, fill = Kingdom)) +
                geom_bar(stat = "identity") + #, position="fill") +
                scale_fill_manual(values = cols) +
                theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
                facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
                xlab("") +
                ylab("") +
                theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
                strip.placement = "top", # move strip to outside
                strip.background = element_blank(), # remove strip background
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
                ggtitle(paste(i, paste(opt$seq_type, opt$read_assembly, "(number of reads mapping)", sep = " "), sep = " - ")) +
                geom_point(data = stats_df, aes(fill = Kingdom, shape = Kingdom, col = Kingdom), size = 8, stroke = 1, show.legend = FALSE) + 
                scale_shape_manual(values = c(4, 3)) +
                geom_text(data = stats_df, aes( label = Kingdom), nudge_x = 0.30, nudge_y = 1, check_overlap = T, vjust = -1) + 
                scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))
        #  reads_plot

        bp_plot <- sample_df %>%       
                filter(Strain == i) %>%
                filter(Type != "total_reads") %>%
                filter(Kingdom != "TOTAL", Kingdom != "Trim") %>%
                ggplot(aes(x = sample_short, y = total, fill = Kingdom)) +
                geom_bar(stat = "identity", position = "fill") +
                scale_fill_manual(values = cols) +
                theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
                facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
                xlab("") +
                ylab("") +
                theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
                strip.placement = "top", # move strip to outside
                strip.background = element_blank(), # remove strip background
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
                ggtitle(paste(i, paste(opt$seq_type, opt$read_assembly, "(total base pairs mapping)", sep = " "), sep = " - ")) +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

        print(ggarrange(reads_plot, bp_plot, ncol = 1))
        invisible(capture.output(suppressMessages(ggsave(paste(i, opt$seq_type, paste(opt$read_assembly, ".svg",sep = ""), sep = "_"), device = "svg", width = 40, height = 30, units = "cm")))) 
        invisible(capture.output(suppressMessages(ggsave(paste(i, opt$seq_type, paste(opt$read_assembly, ".png", sep = ""), sep = "_"), device = "png", width = 40, height = 30, units = "cm"))))
}

for (i in unique(sample_df$Strain)) {
        print(paste("Plotting - TEST:", i, sep = " "))
        mapping_df <- sample_df %>%
                filter(Strain == i) %>%
                filter(Type == "total_reads") %>%
                filter(Kingdom != "TOTAL", Kingdom != "Trim")
        
        stats_df <- sample_df %>%
                filter(Strain == i) %>%
                filter(Type == "total_reads") %>%
                filter(Kingdom != "Other", Kingdom != "Aiptasiidae", Kingdom != "Scleractina", Kingdom != "Symbiodiniaceae", Kingdom != "Bacteria")

        reads_plot <- mapping_df %>%
                ggplot(aes(x = sample_short, y = total, fill = Kingdom)) +
                geom_bar(stat = "identity", position = "dodge2") +
                scale_fill_manual(values = cols) +
                theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
                facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
                xlab("") +
                ylab("") +
                theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
                strip.placement = "top", # move strip to outside
                strip.background = element_blank(), # remove strip background
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
                ggtitle(paste(i, paste(opt$seq_type, opt$read_assembly, "(number of reads mapping)", sep = " "), sep = " - ")) +
                geom_point(data = stats_df, aes(fill = Kingdom, shape = Kingdom, col = Kingdom), size = 8, stroke = 1, show.legend = FALSE) + 
                scale_shape_manual(values = c(4, 3)) +
                geom_text(data = stats_df, aes( label = Kingdom), nudge_x = 0.30, nudge_y = 1, check_overlap = T, vjust = -1) + 
                scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))
        #  reads_plot

        bp_plot <- sample_df %>%       
                filter(Strain == i) %>%
                filter(Type != "total_reads") %>%
                filter(Kingdom != "TOTAL", Kingdom != "Trim") %>%
                ggplot(aes(x = sample_short, y = total, fill = Kingdom)) +
                geom_bar(stat = "identity", position = "fill") +
                scale_fill_manual(values = cols) +
                theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
                facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
                xlab("") +
                ylab("") +
                theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
                strip.placement = "top", # move strip to outside
                strip.background = element_blank(), # remove strip background
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
                ggtitle(paste(i, paste(opt$seq_type, opt$read_assembly, "(total base pairs mapping)", sep = " "), sep = " - ")) +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE))

        print(ggarrange(reads_plot, bp_plot, ncol = 1))
        invisible(capture.output(suppressMessages(ggsave(paste("TEST", i, opt$seq_type, paste(opt$read_assembly, ".svg",sep = ""), sep = "_"), device = "svg", width = 40, height = 30, units = "cm")))) 
        invisible(capture.output(suppressMessages(ggsave(paste("TEST", i, opt$seq_type, paste(opt$read_assembly, ".png", sep = ""), sep = "_"), device = "png", width = 40, height = 30, units = "cm"))))
}
