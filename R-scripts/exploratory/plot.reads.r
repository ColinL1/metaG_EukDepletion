#! /usr/bin/env Rscript

cat("Loading libraries...\n")
lib_list <- c("tidyr", "tidyfast", "dplyr", "ggplot2", "ggpubr", "tidyverse", "scales", "RColorBrewer", "optparse", "pacman","httpgd", "grid") # nolint: line_length_linter, cSpell.
invisible(capture.output(
                suppressMessages(
                        lapply(lib_list, require, quietly = TRUE, character.only = TRUE))))


# Read the file
data <- read.table("/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/mapping/stats_reads_sanitised.txt", header = TRUE, sep = "\t")

# Print the first few rows of the data
head(data)

#print the data structure
str(data)

#check that all file *_1.fq.gz are identical to *_2.fq.gz redo later down too
data %>%
    filter(grepl("_1.fq.gz", file)) %>%
    mutate(file2 = gsub("_1.fq.gz", "_2.fq.gz", file)) %>%
    left_join(data %>% filter(grepl("_2.fq.gz", file)), by = c("file2" = "file")) %>%
    filter(num_seqs.x != num_seqs.y)

# if empty then all good
# drop all *_2.fq.gz or *R2.fq.gz lines
data <- data %>% filter(!grepl("_2.fq.gz", file))
data <- data %>% filter(!grepl("R2.fq.gz", file))

#cleanup data file column remoiving any path
data$file <- basename(data$file)

#copy file column to temp column
data$temp <- data$file

#change _assembly to .assembly
data$temp <- gsub("_TRIM_paried_1.fq.gz", ".total", data$temp)
data$temp <- gsub("_TRIM_paried_R1.fq.gz", ".total", data$temp)
data$temp <- gsub("_L[0-9]\\.", "\\.", data$temp)

#remove _EKDN230034839-1A_HG3F2DSX7 
data$temp <- gsub("_EKDN230034839-1A_HG3F2DSX7", "", data$temp)
data$temp <- gsub("_EKDN230034842-1A_HG3WCDSX7", "", data$temp)
data$temp <- gsub("_EKDN230034838-1A_HG3F2DSX7", "", data$temp)
data$temp <- gsub("_EKDN230034843-1A_HFVN2DSX7", "", data$temp)

data$temp <- gsub("_TRIM_paried", "", data$temp)
data$temp <- gsub(".mapped_1.fq.gz|.unmapped_1.fq.gz", "", data$temp)

#split file coluimn at the . to get sample name and mapping
data <- data %>% separate(temp, c("sample", "mapping"), sep = "\\.", remove = FALSE)

# #remove _trim from sample column
# data$sample <- gsub("_TRIM", "", data$sample)

#split sample to get metadata using "_"
data <- data %>% separate(sample, c("species", "replicate", "extraction", "buffer"), sep = "_", remove = FALSE)

#remove temp column
data$temp <- NULL

#sort columns in order sample species replicate extraction buffer mapping num_seqs sum_len min_len avg_len max_len file format type
data <- data[, c("sample", "species", "replicate", "extraction", "buffer", "mapping", "num_seqs", "sum_len", "min_len", "avg_len", "max_len", "file", "format", "type")]
head(data)

# clean up species names if necessary
#"F003" "F3"   "H2"   "Ac"   "Acro" "Po"   "Poci" "Por"  "Pr" to full names

data$species <- gsub("\\bF3\\b", "F003", data$species)
data$species <- gsub("\\bAcro\\b", "Acropora", data$species)
data$species <- gsub("\\bAc\\b", "Acropora", data$species)
data$species <- gsub("\\bPo\\b", "Pocillopora", data$species)
data$species <- gsub("\\bPoci\\b", "Pocillopora", data$species)
data$species <- gsub("\\bPor\\b", "Porites", data$species)
data$species <- gsub("\\bPr\\b", "Porites", data$species)

# extractions methods  B   BT  MB  MBX MBS to Benzonase "Blood & Tissue" "Microbiome" "Microbiome & beadbeating" "Microbiome & spinning"
data$extraction <- gsub("\\bB\\b", "Benzonase", data$extraction)
data$extraction <- gsub("\\bBT\\b", "Blood & Tissue", data$extraction)
data$extraction <- gsub("\\bMB\\b", "Microbiome", data$extraction)
data$extraction <- gsub("\\bMBX\\b", "Microbiome & bead-beating", data$extraction)
data$extraction <- gsub("\\bMBS\\b", "Microbiome & spinning", data$extraction)

#P to PBS and D to DESS
data$buffer <- gsub("\\bP\\b", "PBS", data$buffer)
data$buffer <- gsub("\\bD\\b", "DESS", data$buffer)

# "non-bacteria" "assembly" to  "other" "total" "sleractina" or "Aiptasia" to "host"
# data$mapping <- gsub("\\bnon-bacteria\\b", "other", data$mapping)
# data$mapping <- gsub("\\bassembly\\b", "total", data$mapping)   
data$mapping <- gsub("\\bscleractina\\b", "host", data$mapping)
data$mapping <- gsub("\\baiptasia\\b", "host", data$mapping)

#set buffer methods mapping as factors 
data$species <- factor(data$species, levels = rev(c("Porites", "Pocillopora", "Acropora", "F003", "H2")))
data$buffer <- factor(data$buffer, levels = rev(sort(unique(data$buffer))))
data$extraction <- factor(data$extraction, levels = (sort(unique(data$extraction))))
data$mapping <- factor(data$mapping, levels = (c("host", "symbiodiniaceae","bacteria", "other","total")))

#check if for "total" sum_len is the same as the sum of all the other sum_len per each mapping 
data %>%
    filter(mapping != "total") %>%
    group_by(sample) %>%
    summarise(total = sum(sum_len)) %>%
    left_join(data %>% filter(mapping == "total"), by = "sample") %>%
    mutate(diff = total - sum_len) %>%
    filter(diff != 0)

#if empty then all good

# call hgd() #http://127.0.0.1:34499/live?token=ifWngxKU

# Create a boxplot of the sum_len column with extractions and strain as facests 
data %>%
    filter(mapping != "total") %>%
    ggplot(aes(x = extraction, y = sum_len, fill = mapping)) +
        geom_boxplot() +
        facet_wrap(extraction~ species, ncol = 5, scale = "free_x") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Boxplot of sum_len column with extractions and strain as facests",
            x = "Extraction",
            y = "Sum Length")

#make stacked barplot of sum_len with mapping as fill and extractions and species as facest
data %>%
    filter(mapping != "total") %>%
    ggplot(aes(x = extraction, y = sum_len, fill = mapping)) +
        geom_bar(stat = "identity", position = "fill") +
        facet_wrap( ~ species, scale = "free_x" , ncol = 1) + #, ncol = 5, scale = "free_x") +
        theme_minimal() +
        coord_flip() +
        theme(axis.text.y = element_text(angle = 45, hjust = 1)) +
        labs(title = "Stacked Barplot of sum_len with mapping as fill and extractions and species as facest",
            x = "Extraction",
            y = "Sum Length")

#make stacked barplot of num_seqs with mapping as fill and extractions and species as facest
data %>%
    filter(mapping != "total") %>%
    ggplot(aes(x = extraction, y = num_seqs, fill = mapping)) +
        geom_bar(stat = "identity", position = "fill") +
        facet_wrap( ~ species, scale = "free_x" , ncol = 1) + #, ncol = 5, scale = "free_x") +
        theme_minimal() +
        coord_flip() +
        theme(axis.text.y = element_text(angle = 45, hjust = 1)) +
        labs(title = "Stacked Barplot of num_seqs with mapping as fill and extractions and species as facest",
            x = "Extraction",
            y = "Num seqs")
