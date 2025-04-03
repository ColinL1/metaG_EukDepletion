#!/usr/bin/env Rscript

# echo status
cat("Loading libraries...\n")
lib_list <- c("tidyr", "tidyfast", "dplyr", "ggplot2", "ggpubr", "ggh4x", "ggsignif", "tidyverse", "scales", "RColorBrewer", "randomcoloR", "optparse", "pacman","httpgd", "grid", "vegan")# nolint: line_length_linter, cSpell.
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
hgd() #http://127.0.0.1:36571/live?token=XPsbjBwv   http://127.0.0.1:36525/live?token=DGweDG5f


# data <- read_csv("/home/colinl/metaG/Git/metaG_EukDepletion/tmp_json/report_with_legths_ONT_no_trim.csv")

# data <- data[,c(2:5)]
# # sequencing_stats <- read.table("/home/colinl/metaG/Git/metaG_methods/results/metaG_indonesia/mmseqs2_reports/ONT_reads_NR_lca.tsv", header = TRUE) # nolint: line_length_linter.
# metadata <- read_csv("/home/colinl/metaG/Git/metaG_EukDepletion/sample_metadata_sheet.csv")

customized_read_tsv <- function(file){
  read_tsv(file, show_col_types = FALSE) %>% # , col_names = c("read_count", "taxid", "species", "genus", "family", "order", "class", "phylum", "superkingdom")
    mutate(filePath = file)
}

# sample <- list.files(path = "/home/colinl/databases/kaiju/input_contigs/tsv", full.names = TRUE) %>% # list all the files
#   lapply(customized_read_tsv) %>% # read them all in with our custom function
#   reduce(bind_rows) %>% # stack them all on top of each other
#   select(filePath, read_count, taxid, species, genus, family, order, class, phylum, clade, superkingdom )# %>% # select the correct columns

sample1 <- customized_read_tsv("/home/colinl/metaG/Git/metaG_EukDepletion/results/host_raw_reads_taxid.tsv")
sample2 <- customized_read_tsv("/home/colinl/metaG/Git/metaG_EukDepletion/results/zoox_raw_reads_taxid.tsv")

sample1$Treatment <- "Blood & Tissue supernatant"
sample2$Treatment <- "Blood & Tissue - zoox pellet"

sample <- rbind(sample1,sample2)
# head(sample)

sample$Strain <- "Porites"
sample$replicate <- "00"
sample$Buffer <- "NA"

sample["file"] <-  gsub("/home/colinl/metaG/Git/metaG_EukDepletion/results/", "", sample$filePath)

sample$file <- gsub(".tsv", "", sample$file)

# sample$file <- gsub("_D", "_DESS", sample$file)
# sample$file <- gsub("_DESSESS", "_DESS", sample$file)
# sample$file <- gsub("_P", "_PBS", sample$file)
# sample$file <- gsub("_PBSBS", "_PBS", sample$file)

# sample$file <- gsub("Ac_", "Acro_", sample$file)
# sample$file <- gsub("F3_", "F003_", sample$file)
# sample$file <- gsub("Pr_", "Por_", sample$file)
# sample$file <- gsub("Po_", "Poci_", sample$file)
# # sample$file <- gsub("_L1", "", sample$file)
# sample$file <- gsub("_L2", "", sample$file)
# sample$file <- gsub("_L3", "", sample$file)
# sample$file <- gsub("_L4", "", sample$file)

unique(sample$file)
# sample$file <- gsub("_P_", "_PBS_", sample$file)
# sample$file <- gsub("_P_", "_PBS_", sample$file)

# sample$file <- gsub("_DESS_EDM", "_rDESS_Illumina-reads_", sample$file)
# sample$file <- gsub("_PBS_EDM", "_rPBS_Illumina-reads_", sample$file)
# sample$file <- gsub("_DESS", "_cDESS_Illumina-contigs_", sample$file)
# sample$file <- gsub("_PBS", "_cPBS_Illumina-contigs_", sample$file)
# sample$file <- gsub("_pass", "_rDESS_ONT-reads_", sample$file)

# ### quick fix ID porties names
# sample$file <- gsub("ID_BST_Por_06","Por_BST-6_MBX_DESS_ONT-reads_", sample$file)
# sample$file <- gsub("ID_GT_Por_06","Por_GT-6_MBX_DESS_ONT-reads_", sample$file)
# sample$file <- gsub("ID_GT_Por_03","Por_GT-3_MBX_DESS_ONT-reads_", sample$file)
# sample$file <- gsub("ID_GT_Por_01","Por_GT-1_MBX_DESS_ONT-reads_", sample$file)
# sample$file <- gsub("ID_BST_Por_04" ,"Por_BST-4_MBX_DESS_ONT-reads_", sample$file)
# sample$file <- gsub("ID_BST_Por_03","Por_BST-3_MBX_DESS_ONT-reads_", sample$file)
# sample$file <- gsub("ID_BST_Por_01","Por_BST-1_MBX_DESS_ONT-reads_", sample$file)
# sample$file <- gsub("ID_BST_Por_07","Por_BST-7_MBX_DESS_ONT-reads_", sample$file)
# sample$file <- gsub("ID_GT_Por_09" ,"Por_GT-9_MBX_DESS_ONT-reads_", sample$file)
# sample$file <- gsub("ID_GT_Por_04","Por_GT-4_MBX_DESS_ONT-reads_", sample$file)
# sample$file <- gsub("ID_BST_Por_05" ,"Por_BST-5_MBX_DESS_ONT-reads_", sample$file)

# sample$file <- gsub("Por_fresh-MBX_01", "Por_01_MBX_PBS_ONT-reads_", sample$file)
# sample$file <- gsub("Por_fresh-MBX_02", "Por_02_MBX_PBS_ONT-reads_", sample$file)
# sample$file <- gsub("Por_fresh-MBX_03", "Por_03_MBX_PBS_ONT-reads_", sample$file)
# sample$file <- gsub("Por_01_fresh-BT", "Por_01_BT_PBS_ONT-reads_",sample$file)
# sample$file <- gsub("Por_02_fresh-BT", "Por_02_BT_PBS_ONT-reads_",sample$file)
# sample$file <- gsub("Por_03_fresh-BT", "Por_03_BT_PBS_ONT-reads_",sample$file)

# sample$file <- gsub("_ONT-reads__ONT_contigs", "_ONT-contigs_",sample$file)

# sample$file <- str_extract(sample$file, "(?:[^_]*_){5}")
# sample$file <- gsub("reads_", "reads", sample$file)
# sample$file <- gsub("contigs_", "contigs", sample$file)
# sample$file <- gsub("rDESS", "DESS", sample$file)
# sample$file <- gsub("cDESS", "DESS", sample$file)
# sample$file <- gsub("rPBS", "PBS", sample$file)
# sample$file <- gsub("cPBS", "PBS", sample$file)

# sample$file <- gsub("H2_00_BT_DESS_ONT-", "H2_00_BT_PBS_ONT-", sample$file)
# sample$file <- gsub("H2_00_MB_DESS_ONT-", "H2_00_MB_PBS_ONT-", sample$file)

# sample$file <- gsub("H2_01_MBX_DESS_ONT-", "H2_01_MBX_PBS_ONT-", sample$file)
# sample$file <- gsub("H2_02_MBX_DESS_ONT-", "H2_01_MBX_PBS_ONT-", sample$file)
# sample$file <- gsub("H2_03_MBX_DESS_ONT-", "H2_01_MBX_PBS_ONT-", sample$file)

# sample$file <- gsub("F003_00_BT_DESS_ONT-", "F003_00_BT_PBS_ONT-", sample$file)
# sample$file <- gsub("F003_00_MB_DESS_ONT-", "F003_00_MB_PBS_ONT-", sample$file)

# sample$file <- gsub("F003_01_MBX_DESS_ONT-", "F003_01_MBX_PBS_ONT-", sample$file)
# sample$file <- gsub("F003_02_MBX_DESS_ONT-", "F003_01_MBX_PBS_ONT-", sample$file)
# sample$file <- gsub("F003_03_MBX_DESS_ONT-", "F003_01_MBX_PBS_ONT-", sample$file)

# sample <- dt_separate(sample, file, into = c("Strain","replicate","Treatment","Buffer","Data_Type"), sep = "_", remove = FALSE) # nolint: line_length_linter.
# sample <- dt_separate(sample, file, into = c("Strain","replicate","Treatment","Buffer"), sep = "_", remove = FALSE) # nolint: line_length_linter.
# sample <- dt_separate(sample, taxon_name, into = c("superkingdom","phylum","class","order","family","genus","species"), sep = ";", remove = FALSE) # nolint: line_length_linter.

# sample_bk <- sample
sample <- sample_bk

sample <- sample[,c("file", "read_count", "taxid", "superkingdom", "clade", "phylum", "class", "order", "family", "genus", "species",  "Strain", "replicate", "Treatment", "Buffer")]

sample <- sample %>%
    # filter(taxid != "0") %>%
    filter(Treatment != "BTAD") %>%
    filter(taxid != "1170705")# is artificial vector sequence <- Cosmid vector pOJ436

# # write root as root
sample[sample$taxid == "1",]$species <- "Root"
sample[sample$taxid == "1",]$genus <- "Root"
sample[sample$taxid == "1",]$family <- "Root"
sample[sample$taxid == "1",]$order <- "Root"
sample[sample$taxid == "1",]$class <- "Root"
sample[sample$taxid == "1",]$phylum <- "Root"
sample[sample$taxid == "1",]$clade <- "Root"
sample[sample$taxid == "1",]$superkingdom <- "Root"

sample[sample$taxid == "0",]$species <- "Other - (non-bacteria)"
sample[sample$taxid == "0",]$genus <- "Other - (non-bacteria)"
sample[sample$taxid == "0",]$family <- "Other - (non-bacteria)"
sample[sample$taxid == "0",]$order <- "Other - (non-bacteria)"
sample[sample$taxid == "0",]$class <- "Other - (non-bacteria)"
sample[sample$taxid == "0",]$phylum <- "Other - (non-bacteria)"
sample[sample$taxid == "0",]$clade <- "Other - (non-bacteria)"
sample[sample$taxid == "0",]$superkingdom <- "Other - (non-bacteria)"

sample[sample$taxid == "131567",]$species <- "cellular organisms"
sample[sample$taxid == "131567",]$genus <- "cellular organisms"
sample[sample$taxid == "131567",]$family <- "cellular organisms"
sample[sample$taxid == "131567",]$order <- "cellular organisms"
sample[sample$taxid == "131567",]$class <- "cellular organisms"
sample[sample$taxid == "131567",]$phylum <- "cellular organisms"
sample[sample$taxid == "131567",]$clade <- "cellular organisms"
sample[sample$taxid == "131567",]$superkingdom <- "cellular organisms"

# sample[sample$taxid == "10292",]$family <- "Herpesviridae"
# sample[sample$taxid == "10292"]$order <- "Herpesvirales"
# sample[sample$taxid == "10292"]$class <- "Herviviricetes"
# sample[sample$taxid == "10292",]$phylum <- "Peploviricota"
# sample[sample$taxid == "10292"]$superkingdom <- "Viruses"

##TODO: pull taxonomy form higher rank when missing
# "superkingdom" "clade"        "phylum"       "class"        "order"        "family"       "genus"        "species"

sample[sample$clade == "unclassified",]$clade <- paste("Unclassified", sample[sample$clade == "unclassified",]$phylum, sep = " ")
sample[sample$phylum == "unclassified",]$phylum <- paste("Unclassified", sample[sample$phylum == "unclassified",]$class , sep = " ")
sample[sample$class == "unclassified",]$class <- paste("Unclassified", sample[sample$class == "unclassified",]$order, sep = " ")
sample[sample$order == "unclassified",]$order <- paste("Unclassified", sample[sample$order == "unclassified",]$family, sep = " ")
sample[sample$family == "unclassified",]$family <- paste("Unclassified", sample[sample$family == "unclassified",]$genus, sep = " ")
sample[sample$genus == "unclassified",]$genus <- paste("Unclassified", sample[sample$genus == "unclassified",]$species, sep = " ")

sample[sample$clade == "Unclassified unclassified",]$clade 
unique(sample$phylum)
unique(sample$class)
unique(sample$order)
unique(sample$family)
unique(sample$genus)
unique(sample$species)

sample$clade <- gsub("Unclassified unclassified", "unclassified", sample$clade)
sample$phylum <- gsub("Unclassified unclassified", "unclassified", sample$phylum)
sample$class <- gsub("Unclassified unclassified", "unclassified", sample$class)
sample$order <- gsub("Unclassified unclassified", "unclassified", sample$order)
sample$family <- gsub("Unclassified unclassified", "unclassified", sample$family)
sample$genus <- gsub("Unclassified unclassified", "unclassified", sample$genus)
sample$species <- gsub("Unclassified unclassified", "unclassified", sample$species)

sample$clade <- gsub("Unclassified Unclassified", "Unclassified", sample$clade)
sample$phylum <- gsub("Unclassified Unclassified", "Unclassified", sample$phylum)
sample$class <- gsub("Unclassified Unclassified", "Unclassified", sample$class)
sample$order <- gsub("Unclassified Unclassified", "Unclassified", sample$order)
sample$family <- gsub("Unclassified Unclassified", "Unclassified", sample$family)
sample$genus <- gsub("Unclassified Unclassified", "Unclassified", sample$genus)
sample$species <- gsub("Unclassified Unclassified", "Unclassified", sample$species)

# # write unclassified taxon as "Unclassified" + the higher taxonomy (for clade) when unclassified
sample[sample$clade == "unclassified",]$species <- paste("Unclassified",sample[sample$clade == "unclassified",]$superkingdom, sep = " ")
sample[sample$clade == "unclassified",]$genus <- paste("Unclassified",sample[sample$clade == "unclassified",]$superkingdom, sep = " ")
sample[sample$clade == "unclassified",]$family <- paste("Unclassified",sample[sample$clade == "unclassified",]$superkingdom, sep = " ")
sample[sample$clade == "unclassified",]$order <- paste("Unclassified",sample[sample$clade == "unclassified",]$superkingdom, sep = " ")
sample[sample$clade == "unclassified",]$class <- paste("Unclassified",sample[sample$clade == "unclassified",]$superkingdom, sep = " ")
sample[sample$clade == "unclassified",]$phylum <- paste("Unclassified",sample[sample$clade == "unclassified",]$superkingdom, sep = " ")
sample[sample$clade == "unclassified",]$clade <- paste("Unclassified",sample[sample$clade == "unclassified",]$superkingdom, sep = " ")
#sample[sample$clade == "unclassified",]$genus <- paste("Unclassified",sample[sample$clade == "unclassified",]$superkingdom, sep = " ")

# # write unclassified taxon as "unclassified" + the higher taxonomy (for phylum) when unclassified
sample[sample$phylum == "unclassified",]$species <- paste("Unclassified",sample[sample$phylum == "unclassified",]$clade, sep = " ")
sample[sample$phylum == "unclassified",]$genus <- paste("Unclassified",sample[sample$phylum == "unclassified",]$clade, sep = " ")
sample[sample$phylum == "unclassified",]$family <- paste("Unclassified",sample[sample$phylum == "unclassified",]$clade, sep = " ")
sample[sample$phylum == "unclassified",]$order <- paste("Unclassified",sample[sample$phylum == "unclassified",]$clade, sep = " ")
sample[sample$phylum == "unclassified",]$class <- paste("Unclassified",sample[sample$phylum == "unclassified",]$clade, sep = " ")
sample[sample$phylum == "unclassified",]$phylum <- paste("Unclassified",sample[sample$phylum == "unclassified",]$clade, sep = " ")

# # write unclassified taxon as "unclassified" + the higher taxonomy (for class) when unclassified
sample[sample$class == "unclassified",]$species <- paste("Unclassified",sample[sample$class == "unclassified",]$phylum, sep = " ")
sample[sample$class == "unclassified",]$genus <- paste("Unclassified",sample[sample$class == "unclassified",]$phylum, sep = " ")
sample[sample$class == "unclassified",]$family <- paste("Unclassified",sample[sample$class == "unclassified",]$phylum, sep = " ")
sample[sample$class == "unclassified",]$order <- paste("Unclassified",sample[sample$class == "unclassified",]$phylum, sep = " ")
sample[sample$class == "unclassified",]$class <- paste("Unclassified",sample[sample$class == "unclassified",]$phylum, sep = " ")

# # write unclassified taxon as "unclassified" + the higher taxonomy (for order) when unclassified
sample[sample$order == "unclassified",]$species <- paste("Unclassified",sample[sample$order == "unclassified",]$class, sep = " ")
sample[sample$order == "unclassified",]$genus <- paste("Unclassified",sample[sample$order == "unclassified",]$class, sep = " ")
sample[sample$order == "unclassified",]$family <- paste("Unclassified",sample[sample$order == "unclassified",]$class, sep = " ")
sample[sample$order == "unclassified",]$order <- paste("Unclassified",sample[sample$order == "unclassified",]$class, sep = " ")

# # write unclassified taxon as "unclassified" + the higher taxonomy (for family) when unclassified
sample[sample$family == "unclassified",]$species <- paste("Unclassified",sample[sample$family == "unclassified",]$order, sep = " ")
sample[sample$family == "unclassified",]$genus <- paste("Unclassified",sample[sample$family == "unclassified",]$order, sep = " ")
sample[sample$family == "unclassified",]$family <- paste("Unclassified",sample[sample$family == "unclassified",]$order, sep = " ")

# # write unclassified taxon as "unclassified" + the higher taxonomy (for genus) when unclassified
sample[sample$genus == "unclassified",]$species <- paste("Unclassified",sample[sample$genus == "unclassified",]$family, sep = " ")
sample[sample$genus == "unclassified",]$genus <- paste("Unclassified",sample[sample$genus == "unclassified",]$family, sep = " ")

# # write unclassified taxon as "unclassified" + the higher taxonomy (for species) when unclassified
sample[sample$species == "unclassified",]$species <- paste("Unclassified",sample[sample$species == "unclassified",]$genus, sep = " ")

sample[grep("Unclassified", sample$species),]$species
sample[grep("Unclassified", sample$species),]
sample[grep("unclassified", sample$genus),]$genus
sample[grep("unclassified", sample$genus),]
sample[grep("Unclassified", sample$family),]$family
sample[grep("unclassified", sample$family),]
sample[grep("unclassified", sample$species),]$species
sample[grep("unclassified", sample$genus),]$genus
sample[grep("unclassified", sample$family),]
sample[grep("unclassified", sample$order),]
sample[grep("unclassified", sample$class),]
sample[grep("unclassified", sample$clade),]
## add genus when you fix the pyhton code. 



sample[sample$phylum == "Unclassified unclassified",]

sample_bk2 <- sample
# sample<-sample_bk2

# sample$Treatment <- gsub("MBX", "Modified-microbiome kit", sample$Treatment)
# sample$Treatment <- gsub("MB", "Microbiome kit", sample$Treatment)
# sample$Treatment <- gsub("B", "Benzonase", sample$Treatment)
# sample$Treatment <- gsub("BenzonaseT", "Blood & Tissue", sample$Treatment)

# sample$Strain <- gsub("Por", "Porites", sample$Strain)
# sample$Strain <- gsub("Poci", "Pocillopora", sample$Strain)
# sample$Strain <- gsub("Acro", "Acropora", sample$Strain)

# sample$Sample <- paste(sample$Strain, sample$replicate, sample$Treatment, sep = "-")
sample$Sample <- paste(sample$replicate, sample$Treatment, sep = "-")
sample$Sample <- sample$file

# sample$Buffer <- gsub("DESSESS", "DESS", sample$Buffer)
# sample$Buffer <- gsub("PBSBS", "PBS", sample$Buffer)
sample

options(dplyr.summarise.inform = FALSE)
# p_load("randomcoloR")
## Plot by strain
# dev.off()
distinctColorPalette(16)
col <- distinctColorPalette(14)
# cellular organisms"     "Unclassified Eukaryota" "Unclassified Bacteria"  "Root"                   "Pseudomonadota"         "Planctomycetota"                               |
#  [7] "Other"                  "Mucoromycota"           "Cyanobacteriota"        "Basidiomycota"          "Bacteroidota"           "Bacillota"                                     |
# [13] "Ascomycota"             "Actinomycetota"                            
# col <- c("cellular organisms" = "#E88A61", "Unclassified Bacteria" =  "#7EE4CE", "Thermodesulfobacteriota" = "#DDA6C7", "Root" =  "#82DE8D", "Pseudomonadota" =  "#D27BD7", "Planctomycetota" = "#DA538E", "Other" =  "#DAD9DD", "Myxococcota" = "#A7A161", "Euryarchaeota" =  "#DAE2B2", "Cyanobacteriota" =  "#7656E1" , "Campylobacterota" = "#88E857", "Bacteroidota" = "#D63FDC", "Bacillota" =  "#7AC3E1", "Actinomycetota" =  "#8C8FDA", "Verrucomicrobiota" =  "#82848B")
# pdf(file = "/home/colinl/metaG/Git/metaG_EukDepletion/Top_phylum.pdf", onefile = TRUE, paper = "a4r")
for (i in c("Porites")) { #, "Pocillopora", "Acropora", "F003", "H2"
    print(paste("Plotting:", i, sep = " "))
      top_spp <- sample %>%
        filter(taxid != 0) %>%
        filter(Strain == i) %>%
          group_by(Sample, Strain, Treatment, Buffer, species) %>%
          summarise(count = sum(read_count)) %>%
          top_n(20, wt = count)

      bac_data <- sample %>%
        filter(Strain == i) %>%
          filter(species %in% top_spp$species) %>%
          group_by(Sample, Strain, Treatment, Buffer, species) %>%
          summarise(count = sum(read_count)) %>%
          mutate(Frequency = count / sum(count))

      other_data <- sample %>%
        filter(Strain == i) %>%
          filter(!species %in% top_spp$species) %>%
          group_by(Sample, Strain, Treatment, Buffer, species) %>%
          summarise(count = sum(read_count)) %>%
          mutate(Frequency = count / sum(count))
      other_data$species <- "Other"

      plot_bacteria_phy_df <- rbind(bac_data, other_data)

      plot_bacteria_phy_df <- plot_bacteria_phy_df %>%
        arrange(desc(species))
      plot_bacteria_phy_df <- plot_bacteria_phy_df %>%
        arrange(desc(Strain))
      plot_bacteria_phy_df <- plot_bacteria_phy_df %>%
        arrange(desc(Treatment))

      plot_bacteria_phy_df$Strain <- factor(plot_bacteria_phy_df$Strain,
        levels = rev(unique(plot_bacteria_phy_df$Strain)))

      plot_bacteria_phy_df$Treatment <- factor(plot_bacteria_phy_df$Treatment,
        levels = rev(unique(plot_bacteria_phy_df$Treatment)))
      # print(length(unique((plot_bacteria_phy_df$species))))
      # print((unique((plot_bacteria_phy_df$species))))
      n <- length(unique((plot_bacteria_phy_df$species)))
      col <- distinctColorPalette(n)

      print(plot_bacteria_phy_df %>%
        group_by(Sample, Strain, Treatment, Buffer, species) %>%
        summarise(Frequency = mean(Frequency)) %>%
          ggplot(aes(x = Sample,
            y = Frequency,
            fill = species)) +
          geom_bar(stat = "identity", position = "fill") +
          facet_nested(Buffer ~ Treatment,
            scale = "free") + #, labeller = label_both
          ylab(element_blank()) +
          scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
          labels = scales::percent_format(accuracy = 1)) +
          scale_fill_manual(values = col) +
          theme_bw(base_size = 16) +
          theme(strip.background = element_rect(color = "black", fill = NA),
          # strip.background = element_rect(color = "black"), # remove strip background
            # strip.text.y = element_blank(), # remove x-axis label for strip
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black")) +
          scale_x_discrete(labels = function(x)
            str_wrap(x, width = 1, whitespace_only = FALSE)) +
          ggtitle(paste("Top 10 bacteria species:", i, "(ONT kaiju refeq - anna)", sep = " "))
      )
}
# dev.off()
for (i in c("superkingdom", "clade", "phylum", "class", "order", "family", "genus", "species")) { #, "Pocillopora", "Acropora", "F003", "H2" "",
    print(paste("Plotting:", i, sep = " "))
      top_spp <- sample %>%
        filter(taxid != 0) %>%
        filter(Strain == "Porites") %>%
          group_by(Sample, Strain, Treatment, Buffer, !!sym(i)) %>%
          summarise(count = sum(read_count)) %>%
          top_n(15, wt = count)

      bac_data <- sample %>%
        filter(Strain == "Porites") %>%
          filter(!!sym(i) %in% top_spp[[i]]) %>%
          group_by(Sample, Strain, Treatment, Buffer, !!sym(i)) %>%
          summarise(count = sum(read_count)) %>%
          mutate(Frequency = count / sum(count))

      other_data <- sample %>%
        filter(Strain == "Porites") %>%
        filter(!(!!sym(i)) %in% top_spp[[i]]) %>%
          group_by(Sample, Strain, Treatment, Buffer, !!sym(i)) %>%
          summarise(count = sum(read_count)) %>%
          mutate(Frequency = count / sum(count))
      other_data[[i]] <- "Other"

      plot_bacteria_phy_df <- rbind(bac_data, other_data)

      plot_bacteria_phy_df <- plot_bacteria_phy_df %>%
        arrange(desc(i))
      plot_bacteria_phy_df <- plot_bacteria_phy_df %>%
        arrange(desc(Strain))
      plot_bacteria_phy_df <- plot_bacteria_phy_df %>%
        arrange(desc(Treatment))

      n <- length(unique((plot_bacteria_phy_df[[i]])))
      print(n)}

      col <- distinctColorPalette(16)


for (i in c("superkingdom", "clade", "phylum", "class", "order", "family", "genus", "species")) { #, "Pocillopora", "Acropora", "F003", "H2" "",
    print(paste("Plotting:", i, sep = " "))
      top_spp <- sample %>%
        filter(taxid != 0) %>%
        filter(Strain == "Porites") %>%
          group_by(Sample, Strain, Treatment, Buffer, !!sym(i)) %>%
          summarise(count = sum(read_count)) %>%
          top_n(10, wt = count)

      bac_data <- sample %>%
        filter(Strain == "Porites") %>%
          filter(!!sym(i) %in% top_spp[[i]]) %>%
          group_by(Sample, Strain, Treatment, Buffer, !!sym(i)) %>%
          summarise(count = sum(read_count)) %>%
          mutate(Frequency = count / sum(count))

      other_data <- sample %>%
        filter(Strain == "Porites") %>%
        filter(!(!!sym(i)) %in% top_spp[[i]]) %>%
          group_by(Sample, Strain, Treatment, Buffer, !!sym(i)) %>%
          summarise(count = sum(read_count)) %>%
          mutate(Frequency = count / sum(count))
      other_data[[i]] <- "Other"

      plot_bacteria_phy_df <- rbind(bac_data, other_data)

      plot_bacteria_phy_df <- plot_bacteria_phy_df %>%
        arrange(desc(i))
      plot_bacteria_phy_df <- plot_bacteria_phy_df %>%
        arrange(desc(Strain))
      plot_bacteria_phy_df <- plot_bacteria_phy_df %>%
        arrange(desc(Treatment))

      plot_bacteria_phy_df$Strain <- factor(plot_bacteria_phy_df$Strain,
        levels = rev(unique(plot_bacteria_phy_df$Strain)))

      plot_bacteria_phy_df$Treatment <- factor(plot_bacteria_phy_df$Treatment,
        levels = rev(unique(plot_bacteria_phy_df$Treatment)))
      # print(length(unique((plot_bacteria_phy_df$i))))
      # print((unique((plot_bacteria_phy_df$i))))
      n <- length(unique((plot_bacteria_phy_df[[i]])))
      # col <- distinctColorPalette(n)

      print(plot_bacteria_phy_df %>%
        group_by(Sample, Strain, Treatment, Buffer, !!sym(i)) %>%
        summarise(Frequency = mean(Frequency)) %>%
          ggplot(aes(x = Sample,
            y = Frequency,
            fill = !!sym(i))) +
          geom_bar(stat = "identity", position = "fill") +
          facet_nested(Buffer ~ Treatment,
            scale = "free") + #, labeller = label_both
          ylab(element_blank()) +
          scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
          labels = scales::percent_format(accuracy = 1)) +
          scale_fill_manual(values = col) +
          theme_bw(base_size = 16) +
          theme(strip.background = element_rect(color = "black", fill = NA),
          # strip.background = element_rect(color = "black"), # remove strip background
            # strip.text.y = element_blank(), # remove x-axis label for strip
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black")) +
          scale_x_discrete(labels = function(x)
            str_wrap(x, width = 1, whitespace_only = FALSE)) +
          ggtitle(paste("Top 10 bacteria", i, "Porites Anna (kaiju refesq)", sep = " "))
      )
}

# pdf(file = "/home/colinl/metaG/Git/metaG_EukDepletion/rarefaction_curve.pdf", onefile = TRUE, paper = "a4r")
for (i in (1)) {
    print(paste("Plotting:", i, sep = " "))
      vegan_df <- sample #%>%
        # filter(Data_Type == i)
      vegan_df <- vegan_df[,c("file", "species")]
      vegan_df <- vegan_df[!grep(pattern = "Unclassified",
        x = vegan_df$species),] %>% 
          filter(species != "Root", species != "Other - (non-bacteria)")
    print(
      rarecurve(as.matrix(table(vegan_df)), step = 1000, cex = 0.75, las = 1)
    )
}
# dev.off()



for (i in c("species","superkingdom", "clade", "phylum", "class", "order", "family", "genus", "species")) {
    print(paste("Plotting:", i, sep = " "))
      vegan_df <- sample# %>%
        # filter(Strain == i)

      vegan_df <- vegan_df[,c("file", i)]
      vegan_df <- vegan_df[!grepl(pattern = "Unclassified",
                    x = vegan_df[[i]]), ] %>%
                      filter(!!sym(i) != "Root",
                        !!sym(i) != "Other - (non-bacteria)")

      shannon <- diversity(as.matrix(table(vegan_df)), index = "shannon")

      shannon_df <- data.frame(file = rownames(as.data.frame(shannon)),
        shannon = as.data.frame(shannon)[,1])

      shannon_df$Treatment <- "X"
      shannon_df[shannon_df$file == "host_raw_reads_taxid",]$Treatment <- "Blood & Tissue supernatant"
      shannon_df[shannon_df$file == "zoox_raw_reads_taxid",]$Treatment <- "Blood & Tissue - zoox pellet"
      
      shannon_df$Strain <- "Porites"
      shannon_df$Buffer <- "NA"
      shannon_df$Strain <- factor(shannon_df$Strain,
        levels = rev(unique(shannon_df$Strain)))
    print(shannon_df)
    print(
      shannon_df %>%
        # filter(Strain == i) %>%
      ggplot(aes(x = Treatment, y = shannon)) +
        geom_boxplot(aes(x = Treatment, y = shannon, fill = Treatment)) +
        geom_jitter(color="black", size = 0.4, alpha = 0.9) +
        theme_bw(base_size = 16) +
          theme(legend.position = "bottom",
            strip.background = element_rect(color = NA, fill = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black")) +
          scale_x_discrete(labels = function(x)
            str_wrap(x, width = 1, whitespace_only = FALSE)) +
        ggtitle(paste("Shannon diversity :", i, sep = " ")) +
        xlab("") +
        facet_grid( Strain ~ Buffer, scale = "free")
    )
}

# pdf(file = "/home/colinl/metaG/Git/metaG_EukDepletion/shannon_index.pdf", onefile = TRUE, paper = "a4r")
for (i in unique(sample$Strain)) {
    print(paste("Plotting:", i, sep = " "))
      vegan_df <- sample %>%
        filter(Strain == i)

      vegan_df <- vegan_df[,c("file", "species")]
      vegan_df <- vegan_df[!grep(pattern = "Unclassified",
                    x = vegan_df$species), ] %>%
                      filter(species != "Root",
                        species != "Other - (non-bacteria)")

      shannon <- diversity(as.matrix(table(vegan_df)), index = "shannon")

      shannon_df <- data.frame(file = rownames(as.data.frame(shannon)),
        shannon = as.data.frame(shannon)[,1])

      # shannon_df <- dt_separate(shannon_df, file,
      #   into = c("Strain", "replicate", "Treatment", "Buffer"),
      #   sep = "_", remove = FALSE)

      # shannon_df$Treatment <- gsub("MBX",
      #   "Modified-microbiome kit",shannon_df$Treatment)
      # shannon_df$Treatment <- gsub("MB",
      #   "Microbiome kit", shannon_df$Treatment)
      # shannon_df$Treatment <- gsub("B",
      #   "Benzonase", shannon_df$Treatment)
      # shannon_df$Treatment <- gsub("BenzonaseT",
      #   "Blood & Tissue", shannon_df$Treatment)

      # shannon_df$Strain <- gsub("Por", "Porites", shannon_df$Strain)
      # shannon_df$Strain <- gsub("Poci", "Pocillopora", shannon_df$Strain)
      # shannon_df$Strain <- gsub("Acro", "Acropora", shannon_df$Strain)
      shannon_df$Treatment <- "X"
      shannon_df[shannon_df$file == "host_raw_reads_taxid",]$Treatment <- "Blood & Tissue supernatant"
      shannon_df[shannon_df$file == "zoox_raw_reads_taxid",]$Treatment <- "Blood & Tissue - zoox pellet"
      
      shannon_df$Strain <- "Porites"
      shannon_df$Buffer <- "NA"
      shannon_df$Strain <- factor(shannon_df$Strain,
        levels = rev(unique(shannon_df$Strain)))

    print(
      shannon_df %>%
        filter(Strain == i) %>%
      ggplot(aes(x = Treatment, y = shannon)) +
        geom_boxplot(aes(x = Treatment, y = shannon, fill = Treatment)) +
        geom_jitter(color="black", size = 0.4, alpha = 0.9) +
        theme_bw(base_size = 16) +
          theme(legend.position = "bottom",
            strip.background = element_rect(color = NA, fill = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black")) +
          scale_x_discrete(labels = function(x)
            str_wrap(x, width = 1, whitespace_only = FALSE)) +
        ggtitle(paste("Shannon diversity :", i, sep = " ")) +
        xlab("") +
        # geom_signif(
        #   comparisons = list(
        #     c("Modified-microbiome kit", "Blood & Tissue"),
        #     c("Modified-microbiome kit","Microbiome kit"),
        #     c("Blood & Tissue","Microbiome kit"),
        #     c("Blood & Tissue","Benzonase")
        #                 ),
        #   map_signif_level = TRUE
        # ) +
        facet_grid( Strain ~ Buffer, scale = "free")
    )
}
# dev.off()

####  coral plot shannon ######
      
      vegan_df <- sample %>%
        filter(Strain != "H2", Strain != "F003")

      vegan_df <- vegan_df[,c("file", "species")]
      vegan_df <- vegan_df[!grep(pattern = "Unclassified",
                    x = vegan_df$species), ] %>%
                      filter(species != "Root",
                        species != "Other - (non-bacteria)")

      shannon <- diversity(as.matrix(table(vegan_df)), index = "shannon")

      shannon_df <- data.frame(file = rownames(as.data.frame(shannon)),
        shannon = as.data.frame(shannon)[,1])

      shannon_df <- dt_separate(shannon_df, file,
        into = c("Strain", "replicate", "Treatment", "Buffer"),
        sep = "_", remove = FALSE)

      shannon_df$Treatment <- gsub("MBX",
        "Modified-microbiome kit",shannon_df$Treatment)
      shannon_df$Treatment <- gsub("MB",
        "Microbiome kit", shannon_df$Treatment)
      shannon_df$Treatment <- gsub("B",
        "Benzonase", shannon_df$Treatment)
      shannon_df$Treatment <- gsub("BenzonaseT",
        "Blood & Tissue", shannon_df$Treatment)

      shannon_df$Strain <- gsub("Por", "Porites", shannon_df$Strain)
      shannon_df$Strain <- gsub("Poci", "Pocillopora", shannon_df$Strain)
      shannon_df$Strain <- gsub("Acro", "Acropora", shannon_df$Strain)

      shannon_df$Strain <- factor(shannon_df$Strain,
        levels = rev(unique(shannon_df$Strain)))

    print(
      shannon_df %>%
        filter(Strain != "H2", Strain != "F003") %>%
      ggplot(aes(x = Treatment, y = shannon)) +
        geom_boxplot(aes(x = Treatment, y = shannon, fill = Treatment)) +
        geom_jitter(color="black", size = 0.4, alpha = 0.9) +
        theme_bw(base_size = 16) +
          theme(legend.position = "bottom",
            strip.background = element_rect(color = NA, fill = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black")) +
          scale_x_discrete(labels = function(x)
            str_wrap(x, width = 1, whitespace_only = FALSE)) +
        ggtitle(paste("Shannon diversity :", "corals", sep = " ")) +
        xlab("") +
        # geom_signif(
        #   comparisons = list(
        #     c("Modified-microbiome kit", "Blood & Tissue"),
        #     c("Modified-microbiome kit","Microbiome kit"),
        #     c("Blood & Tissue","Microbiome kit"),
        #     c("Blood & Tissue","Benzonase")
        #                 ),
        #   map_signif_level = TRUE
        # ) +
        facet_grid( Strain ~ Buffer, scale = "free")
    )


#### aip plot shannon #######
      vegan_df <- sample %>%
        filter(Strain != "Acropora", Strain != "Pocillopora", Strain != "Porites")

      vegan_df <- vegan_df[,c("file", "species")]
      vegan_df <- vegan_df[!grep(pattern = "Unclassified",
                    x = vegan_df$species), ] %>%
                      filter(species != "Root",
                        species != "Other - (non-bacteria)")

      shannon <- diversity(as.matrix(table(vegan_df)), index = "shannon")

      shannon_df <- data.frame(file = rownames(as.data.frame(shannon)),
        shannon = as.data.frame(shannon)[,1])

      shannon_df <- dt_separate(shannon_df, file,
        into = c("Strain", "replicate", "Treatment", "Buffer"),
        sep = "_", remove = FALSE)

      shannon_df$Treatment <- gsub("MBX",
        "Modified-microbiome kit",shannon_df$Treatment)
      shannon_df$Treatment <- gsub("MB",
        "Microbiome kit", shannon_df$Treatment)
      shannon_df$Treatment <- gsub("B",
        "Benzonase", shannon_df$Treatment)
      shannon_df$Treatment <- gsub("BenzonaseT",
        "Blood & Tissue", shannon_df$Treatment)

      shannon_df$Strain <- gsub("Por", "Porites", shannon_df$Strain)
      shannon_df$Strain <- gsub("Poci", "Pocillopora", shannon_df$Strain)
      shannon_df$Strain <- gsub("Acro", "Acropora", shannon_df$Strain)

      shannon_df$Strain <- factor(shannon_df$Strain,
        levels = rev(unique(shannon_df$Strain)))

        shannon_df$Strain <- factor(shannon_df$Strain,
          levels = rev(unique(shannon_df$Strain)))

      data_grp <- shannon_df[,
        c("file", "Strain", "replicate", "Treatment", "Buffer")]

      data_richness <-
        estimateR(as.matrix(table(vegan_df))) # calculate richness and Chao1 using vegan package

      data_evenness <- 
        diversity(as.matrix(table(vegan_df))) /
          log(specnumber(as.matrix(table(vegan_df)))) # calculate evenness index using vegan package

      data_shannon <- diversity(as.matrix(table(vegan_df)),
        index = "shannon") # calculate Shannon index using vegan package

      data_alphadiv <- cbind(data_grp,
        t(data_richness),
        data_shannon,
        data_evenness) # combine all indices in one data table

      data_alphadiv_tidy <-  data_alphadiv %>%
        # mutate(sample_id = data_alphadiv$file) %>%
        gather(key   = alphadiv_index,
              value = obs_values,
              -file, -Strain, -replicate, -Treatment, -Buffer)

    print(
      shannon_df %>%
        filter(Strain != "Acropora", Strain != "Pocillopora", Strain != "Porites") %>%
      ggplot(aes(x = Treatment, y = shannon)) +
        geom_boxplot(aes(x = Treatment, y = shannon, fill = Treatment)) +
        geom_jitter(color="black", size = 0.4, alpha = 0.9) +
        theme_bw(base_size = 16) +
          theme(legend.position = "bottom",
            strip.background = element_rect(color = NA, fill = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black")) +
          scale_x_discrete(labels = function(x)
            str_wrap(x, width = 1, whitespace_only = FALSE)) +
        ggtitle(paste("Shannon diversity :", "Aiptasia", sep = " ")) +
        xlab("") +
        # geom_signif(
        #   comparisons = list(
        #     c("Modified-microbiome kit", "Blood & Tissue"),
        #     c("Modified-microbiome kit","Microbiome kit"),
        #     c("Blood & Tissue","Microbiome kit"),
        #     c("Blood & Tissue","Benzonase")
        #                 ),
        #   map_signif_level = TRUE
        # ) +
        facet_grid( Strain ~ Buffer, scale = "free")
    )
    print(
      data_alphadiv_tidy %>%
      filter(alphadiv_index !=  "se.chao1",
          alphadiv_index != "S.chao1",
          alphadiv_index != "S.obs",
          alphadiv_index != "S.ACE",
          alphadiv_index != "se.ACE") %>%
        ggplot(aes(x = obs_values, y = Treatment)) +
        geom_boxplot(aes(fill = Treatment)) +
        facet_nested(alphadiv_index ~ Buffer,
          scale = "free") + #, labeller = label_both
            # (Data_Type ~ Buffer + alphadiv_index, scale = "free") +
        labs(title = i, x = "", y = "", tag = "B") +
        coord_flip() +
        geom_jitter(color = "black", size = 0.4, alpha = 0.9) 
    )



for (i in unique(sample$Strain)) {
    print(paste("Plotting:", i, sep = " "))
      vegan_df <- sample %>%
        filter(Strain == i)

      vegan_df <- vegan_df[,c("file", "species")]
      vegan_df <- vegan_df[!grep(pattern = "Unclassified",
                    x = vegan_df$species), ] %>%
                      filter(species != "Root",
                        species != "Other - (non-bacteria)")
      shannon <- diversity(as.matrix(table(vegan_df)), index = "shannon")

        shannon_df <- data.frame(file = rownames(as.data.frame(shannon)),
          shannon = as.data.frame(shannon)[,1])

        shannon_df <- dt_separate(shannon_df, file,
          into = c("Strain", "replicate", "Treatment", "Buffer"),
          sep = "_", remove = FALSE)


        shannon_df$Treatment <- gsub("MBX",
          "Modified-microbiome kit",shannon_df$Treatment)
        shannon_df$Treatment <- gsub("MB",
          "Microbiome kit", shannon_df$Treatment)
        shannon_df$Treatment <- gsub("B",
          "Benzonase", shannon_df$Treatment)
        shannon_df$Treatment <- gsub("BenzonaseT",
          "Blood & Tissue", shannon_df$Treatment)

        shannon_df$Strain <- gsub("Por", "Porites", shannon_df$Strain)
        shannon_df$Strain <- gsub("Poci", "Pocillopora", shannon_df$Strain)
        shannon_df$Strain <- gsub("Acro", "Acropora", shannon_df$Strain)

        shannon_df$Strain <- factor(shannon_df$Strain,
          levels = rev(unique(shannon_df$Strain)))

      data_grp <- shannon_df[,
        c("file", "Strain", "replicate", "Treatment", "Buffer")]

      data_richness <-
        estimateR(as.matrix(table(vegan_df))) # calculate richness and Chao1 using vegan package

      data_evenness <- 
        diversity(as.matrix(table(vegan_df))) /
          log(specnumber(as.matrix(table(vegan_df)))) # calculate evenness index using vegan package

      data_shannon <- diversity(as.matrix(table(vegan_df)),
        index = "shannon") # calculate Shannon index using vegan package

      data_alphadiv <- cbind(data_grp,
        t(data_richness),
        data_shannon,
        data_evenness) # combine all indices in one data table

      data_alphadiv_tidy <-  data_alphadiv %>%
        # mutate(sample_id = data_alphadiv$file) %>%
        gather(key   = alphadiv_index,
              value = obs_values,
              -file, -Strain, -replicate, -Treatment, -Buffer)

    print(
      data_alphadiv_tidy %>%
      filter(alphadiv_index !=  "se.chao1",
          alphadiv_index != "S.chao1",
          alphadiv_index != "S.obs",
          alphadiv_index != "S.ACE",
          alphadiv_index != "se.ACE") %>%
        ggplot(aes(x = obs_values, y = Treatment)) +
        geom_boxplot(aes(fill = Treatment)) +
        facet_nested(alphadiv_index ~ Buffer,
          scale = "free") + #, labeller = label_both
            # (Data_Type ~ Buffer + alphadiv_index, scale = "free") +
        labs(title = i, x = "", y = "", tag = "B") +
        coord_flip() +
        geom_jitter(color = "black", size = 0.4, alpha = 0.9) 
    )
}

unique((data_alphadiv_tidy$alphadiv_index))
pairs(data_alphadiv[,c("S.obs","S.chao1", "data_shannon", "data_evenness")])
i<- "Porites"

summary(aov(data_shannon ~ Treatment * Strain, data = data_alphadiv))

aov_test <- aov(data_shannon ~ Treatment * Strain, data = data_alphadiv)  
summary(aov_test)

p_load("agricolae")

hsd_test <- TukeyHSD(aov_test)
summary(hsd_test)
hsd_res <- HSD.test(aov_test, "Treatment", group = TRUE)$groups  
hsd_res


kruskal.test(data_shannon ~ Treatment, data = data_alphadiv)
kruskal.test(data_shannon ~ Strain, data = data_alphadiv)
p_load("FSA")

tmp <- data_alphadiv

tmp$Treatment <- gsub("Modified-microbiome kit", "Modified_microbiome kit",tmp$Treatment)

PT <- dunnTest(data_shannon ~ Treatment, data = tmp, method="bh") # require the FSA package     
PT2 <- PT$res  
p_load("rcompanion")
cldList(comparison = PT2$Comparison, p.value = PT2$P.adj, threshold  = 0.05) # require the rcompanion package   

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")


biocLite("ComplexHeatmap")
library(ComplexHeatmap)

beta <- read.table("/home/colinl/metaG/Git/metaG_EukDepletion/results/illumina/diversity_reports/beta_diversity.txt", row.names=1, col.names =c(";","H2_02_MBX_D" ,"F3_01_MBX_P" ,"H2_03_BT_PBS" ,"F003_03_B_PBS" ,"H2_02_MB_DESS" ,"F003_02_B_PBS" ,"F003_01_MB_PBS" ,"F003_01_BT_DESS" ,"H2_03_MB_PBS" ,"H2_03_BT_DESS" ,"F003_02_BT_PBS" ,"F003_03_MB_PBS" ,"H2_02_BT_DESS" ,"H2_03_MBX_D" ,"H2_03_B_DESS" ,"F3_03_MBX_D" ,"H2_03_MBX_P" ,"H2_01_MB_PBS" ,"H2_01_MBX_D" ,"F003_03_BT_PBS" ,"F003_03_BT_DESS" ,"H2_02_BT_PBS" ,"H2_01_MBX_P" ,"F003_01_BT_PBS" ,"H2_03_MBX_D" ,"H2_02_MB_PBS" ,"H2_01_B_PBS" ,"F003_02_BT_DESS" ,"F003_02_MB_PBS" ,"F3_02_MBX_P" ,"F003_01_B_PBS" ,"H2_02_MBX_P" ,"H2_01_BT_DESS" ,"H2_01_BT_PBS" ,"F3_01_MBX_D" ,"H2_01_MB_DESS" ,"F3_02_MBX_D" ,"H2_02_B_DESS" ,"H2_01_B_DESS" ,"F3_03_MBX_P")  )
beta <- read.table("/home/colinl/metaG/Git/metaG_EukDepletion/results/illumina/kraken2/corals/diversity_reports/beta_diversity.txt", row.names=1, col.names =c(";","H2_02_MBX_D" ,"F3_01_MBX_P" ,"H2_03_BT_PBS" ,"F003_03_B_PBS" ,"H2_02_MB_DESS" ,"F003_02_B_PBS" ,"F003_01_MB_PBS" ,"F003_01_BT_DESS" ,"H2_03_MB_PBS" ,"H2_03_BT_DESS" ,"F003_02_BT_PBS" ,"F003_03_MB_PBS" ,"H2_02_BT_DESS" ,"H2_03_MBX_D" ,"H2_03_B_DESS" ,"F3_03_MBX_D" ,"H2_03_MBX_P" ,"H2_01_MB_PBS" ,"H2_01_MBX_D" ,"F003_03_BT_PBS" ,"F003_03_BT_DESS" ,"H2_02_BT_PBS" ,"H2_01_MBX_P" ,"F003_01_BT_PBS" ,"H2_03_MBX_D" ,"H2_02_MB_PBS" ,"H2_01_B_PBS" ,"F003_02_BT_DESS" ,"F003_02_MB_PBS" ,"F3_02_MBX_P" ,"F003_01_B_PBS" ,"H2_02_MBX_P" ,"H2_01_BT_DESS" ,"H2_01_BT_PBS" ,"F3_01_MBX_D" ,"H2_01_MB_DESS" ,"F3_02_MBX_D" ,"H2_02_B_DESS" ,"H2_01_B_DESS" ,"F3_03_MBX_P")  )
# c("H2_02_MBX_D" ,"F3_01_MBX_P" ,"H2_03_BT_PBS" ,"F003_03_B_PBS" ,"H2_02_MB_DESS" ,"F003_02_B_PBS" ,"F003_01_MB_PBS" ,"F003_01_BT_DESS" ,"H2_03_MB_PBS" ,"H2_03_BT_DESS" ,"F003_02_BT_PBS" ,"F003_03_MB_PBS" ,"H2_02_BT_DESS" ,"H2_03_MBX_D" ,"H2_03_B_DESS" ,"F3_03_MBX_D" ,"H2_03_MBX_P" ,"H2_01_MB_PBS" ,"H2_01_MBX_D" ,"F003_03_BT_PBS" ,"F003_03_BT_DESS" ,"H2_02_BT_PBS" ,"H2_01_MBX_P" ,"F003_01_BT_PBS" ,"H2_03_MBX_D" ,"H2_02_MB_PBS" ,"H2_01_B_PBS" ,"F003_02_BT_DESS" ,"F003_02_MB_PBS" ,"F3_02_MBX_P" ,"F003_01_B_PBS" ,"H2_02_MBX_P" ,"H2_01_BT_DESS" ,"H2_01_BT_PBS" ,"F3_01_MBX_D" ,"H2_01_MB_DESS" ,"F3_02_MBX_D" ,"H2_02_B_DESS" ,"H2_01_B_DESS" ,"F3_03_MBX_P")
beta
length(beta$H2_02_MBX_D)
length(beta)
beta_clean <- beta[2:41, ]
plot(as.matrix(beta_clean))

beta_clean
heatmap(as.matrix(beta_clean))
heatmap(data.matrix(beta_clean))
hclust((beta_clean))


customized_read_tsv <- function(file){
  read_tsv(file, show_col_types = FALSE) %>% # , col_names = c("read_count", "taxid", "species", "genus", "family", "order", "class", "phylum", "superkingdom")
    mutate(filePath = file)
}


sample_c <- read.table("/home/colinl/metaG/Git/metaG_EukDepletion/results/illumina/kraken2/corals/diversity_reports/all_alphadiv.txt", header = TRUE )
sample_a <- read.table("/home/colinl/metaG/Git/metaG_EukDepletion/results/illumina/kraken2/aip/diversity_reports/all_alphadiv.txt", header = TRUE )
sample <- rbind(sample_c, sample_a)

sample <- sample %>% separate(Sample, , into = c("Strain","replicate","Treatment","Buffer"), sep = "_", remove = FALSE)
# sample <- separate(sample, sample, into = c("Strain","replicate","Treatment","Buffer"), sep = "_", remove = TRUE) # nolint: line_length_linter.
sample$Strain <- gsub ("F3", "F003", sample$Strain)
sample$Strain <- gsub ("Ac", "Acropora", sample$Strain)
sample$Strain <- gsub ("Acroporaro", "Acropora", sample$Strain)

sample$Strain <- gsub ("Po", "Pocillopora", sample$Strain)
sample$Strain <- gsub ("Pocilloporar", "Porites", sample$Strain)
sample$Strain <- gsub ("Pocilloporaci", "Pocillopora", sample$Strain)

sample$Strain <- gsub ("Pr", "Porites", sample$Strain)
unique(sample$Strain)

sample$Buffer <- gsub ("D", "DESS", sample$Buffer )
sample$Buffer <- gsub ("DESSESS", "DESS", sample$Buffer )
sample$Buffer <- gsub ("P", "PBS", sample$Buffer )
sample$Buffer <- gsub ("PBSBS", "PBS", sample$Buffer )
unique(sample$Buffer)

sample <- as.data.frame(sample) 
as.data.frame(sample) %>% 
  pivot_longer(
    cols = Shannon_diversity:Fisher_alpha,
    values_to = "diversity"
)


# for (i in c("Shannon_diversity", "Bp_diversity", "Simpson_diversity", "Simpson_RIndex", "Fisher_alpha")) {
#     print(paste("Plotting:", i, sep = " "))
#       print(
sample_2 <- sample %>%
            pivot_longer(
            cols = Shannon_diversity:Fisher_alpha,
            values_to = "diversity" ) 
            

sample_2$Strain <- factor(sample_2$Strain, levels = rev(c("Porites", "Pocillopora", "Acropora", "F003", "H2")))
sample_2$Buffer <- factor(sample_2$Buffer, levels = rev(c("PBS", "DESS")))


            sample_2 %>% 
            ggplot(aes(x = Treatment, y = diversity)) +
        geom_boxplot(aes(x = Treatment, y = diversity, fill = Treatment)) +
        geom_jitter(color="black", size = 0.4, alpha = 0.9) +
        theme_bw(base_size = 16) +
        facet_grid(name ~ Strain+Buffer, scale = "free") + 
          theme(legend.position = "bottom",
            strip.background = element_rect(color = NA, fill = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black")) +
          scale_x_discrete(labels = function(x)
            str_wrap(x, width = 1, whitespace_only = FALSE)) +
        ggtitle("plot diversity")  +
        xlab("") +
        ylab("") 
#         )
# }
names(sample) #"Shannon_diversity" "Bp_diversity"      "Simpson_diversity" "Simpson_RIndex"    "Fisher_alpha"

for (i in unique(sample_2$Strain)) {
        print(paste("Plotting:", i, sep = " "))
        # print(ggarrange(sample_df %>%
        #         filter(Strain == i) %>%
        #         ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
        #         geom_bar(stat = "identity", position="fill") +
        #         scale_fill_manual(values = cols) +
        #         theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
        #         facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
        #         xlab("") +
        #         ylab("") +
        #         theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
        #         strip.placement = "top", # move strip to outside
        #         strip.background = element_blank(), # remove strip background
        #         panel.grid.major = element_blank(),
        #         panel.grid.minor = element_blank(),
        #         axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
        #         ggtitle(paste(i, "contigs illumina", sep = " - ")) +
        #         scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE)), 
        print(        
        sample_2 %>%
          filter(Strain == i) %>%
            ggplot(aes(x = Treatment, y = diversity)) +
        geom_boxplot(aes(x = Treatment, y = diversity, fill = Treatment)) +
        geom_jitter(color="black", size = 0.4, alpha = 0.9) +
        theme_bw(base_size = 16) +
        facet_grid(name ~ Strain+Buffer, scale = "free") + 
          theme(legend.position = "bottom",
            strip.background = element_rect(color = NA, fill = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black")) +
          scale_x_discrete(labels = function(x)
            str_wrap(x, width = 1, whitespace_only = FALSE)) +
        ggtitle("plot diversity")  +
        xlab("") +
        ylab("") 
        )
                # sample_df %>%       
                # filter(Strain == i) %>%
                # # filter(Kingdom != "Other") %>%
                # # ggplot(aes(x = sample_short, y = Total_percentage, fill = Kingdom)) +
                # # geom_bar(stat = "identity", position = "fill") +
                # # scale_fill_manual(values = cols) +
                # # theme_bw(base_size = 10) + # optional argument not to tweak all the sizes independently "base_size = 22"
                # # facet_wrap(Buffer ~ Treatment, nrow = 1, scale = "free_x") +
                # # xlab("") +
                # # ylab("") +
                # # theme(panel.spacing = unit(0.5, "lines"), # reduce panel spacing
                # # strip.placement = "top", # move strip to outside
                # # strip.background = element_blank(), # remove strip background
                # # panel.grid.major = element_blank(),
                # # panel.grid.minor = element_blank(),
                # # axis.text.x = element_text(angle = -0, hjust = 0.5, vjust = 1)) + 
                # # ggtitle(paste(i, "contigs illumina (no other)", sep = " - ")) +
                # # scale_x_discrete(labels = function(x) str_wrap(x, width = 1, whitespace_only = FALSE)), ncol = 1))
        # invisible(capture.output(suppressMessages(ggsave(paste(i, "contigs-illumina.svg", sep = "_"), device = "svg", width = 40, height = 30, units = "cm")))) 
        # invisible(capture.output(suppressMessages(ggsave(paste(i, "contigs-illumina.png", sep = "_"), device = "png", width = 40, height = 30, units = "cm"))))
}

sample_2 %>%
  filter(Strain != "H2", Strain != "F003") %>%
  ggplot(aes(x = Treatment, y = diversity)) +
    geom_boxplot(aes(x = Treatment, y = diversity, fill = Treatment)) +
    geom_jitter(color="black", size = 0.4, alpha = 0.9) +
    theme_bw(base_size = 16) +
    facet_grid(name ~ Strain+Buffer, scale = "free") + 
      theme(legend.position = "bottom",
        strip.background = element_rect(color = NA, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black")) +
      scale_x_discrete(labels = function(x)
        str_wrap(x, width = 1, whitespace_only = FALSE)) +
    ggtitle("boxplot coral 𝛼-diversity")  +
    xlab("") +
    ylab("") 
invisible(capture.output(suppressMessages(ggsave(paste("coral", "alpha-diversity-k2.svg", sep = "_"), device = "svg", width = 40, height = 30, units = "cm")))) 
invisible(capture.output(suppressMessages(ggsave(paste("coral", "alpha-diversity-k2.png", sep = "_"), device = "png", width = 40, height = 30, units = "cm"))))

# Acropora    Pocillopora Porites     F003        H2


sample_2 %>%
  filter(Strain != "Acropora", Strain != "Pocillopora", Strain != "Porites") %>%
  ggplot(aes(x = Treatment, y = diversity)) +
    geom_boxplot(aes(x = Treatment, y = diversity, fill = Treatment)) +
    geom_jitter(color="black", size = 0.4, alpha = 0.9) +
    theme_bw(base_size = 16) +
    facet_grid(name ~ Strain+Buffer, scale = "free") + 
      theme(legend.position = "bottom",
        strip.background = element_rect(color = NA, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black")) +
      scale_x_discrete(labels = function(x)
        str_wrap(x, width = 1, whitespace_only = FALSE)) +
    ggtitle("boxplot aiptasia alpha-diversity")  +
    xlab("") +
    ylab("") 
invisible(capture.output(suppressMessages(ggsave(paste("aiptasia", "alpha-diversity-k2.svg", sep = "_"), device = "svg", width = 40, height = 30, units = "cm")))) 
invisible(capture.output(suppressMessages(ggsave(paste("aiptasia", "alpha-diversity-k2.png", sep = "_"), device = "png", width = 40, height = 30, units = "cm"))))
