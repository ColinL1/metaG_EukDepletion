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
                help = "metadata csv file. Necessary to do plots.", metavar = "character")
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

invisible(capture.output(suppressMessages(data <- read_csv(opt$file))))
data <- data[,c(2:5)]

invisible(capture.output(suppressMessages(metadata <- read_csv(opt$metadata))))

sample_data <- merge(data, metadata, by.x = "file_name", by.y = "Index", all.x = TRUE)

sample_data["sample"] <- paste(sample_data$Strain, sample_data$Replicate, sample_data$From, sample_data$DATA_Type, sep = "-") 
# sample_data["sample_short"] <- paste(sample_data$Replicate, sample_data$From, sample_data$DATA_Type, sep = "-") 

# sample_data["sample_short"] <- paste(sample_data$Replicate, sample_data$From, sep = "-") 
# sample_data["sample_short"] <- paste(sample_data$Strain, sample_data$Replicate, sample_data$Treatment, sample_data$Buffer, sep = "_") 
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
        print(ggarrange(sample_df %>%
                filter(Strain == i) %>%
                ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
                geom_bar(stat = "identity") + # , position="fill") +
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
                xlab("") +
                ylab("") +
                theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
                strip.placement = "top", # move strip to outside
                strip.background = element_blank(), # remove strip background
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
                ggtitle(paste(i, "Reads illumina (no other)", sep = " - ")) +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE)), ncol = 1))
        invisible(capture.output(suppressMessages(ggsave(paste(i, "Reads-illumina.svg", sep = "_"), device = "svg", width = 40, height = 30, units = "cm")))) 
        invisible(capture.output(suppressMessages(ggsave(paste(i, "Reads-illumina.png", sep = "_"), device = "png", width = 40, height = 30, units = "cm"))))
}
