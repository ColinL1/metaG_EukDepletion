library(pacman)
library(tidyverse)

# Read the TSV file (adjust file path as needed)
#path 
folder_path_coass <- "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/binning/nf_core_MAG/results/GenomeBinning/QC/"
folder_path <- "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/binning/nf_core_MAG/results_no_coass/GenomeBinning/QC/"
# mags <- read_tsv("/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/binning/nf_core_MAG/results/GenomeBinning/QC/gunc_checkm_summary.tsv")
mags_coass <- read_tsv(paste0(folder_path_coass,"gunc_checkm_summary.tsv"))
mags <- read_tsv(paste0(folder_path,"gunc_checkm_summary.tsv"))


mags_parsed <- mags %>%
    separate(genome, into = c("assembler", "binner", "temp"), sep = "-", remove = FALSE) %>%
    separate(temp, into = c("strain", "replicate", "treatment", "buffer_bin"), sep = "_", remove = FALSE) %>%
    separate(buffer_bin, into = c("buffer", "bin_id"), sep = "\\.", remove = FALSE)


# # Summarize number of bins per strain, treatment, and buffer
summary_mags <- mags_parsed %>%
    filter(checkM.completeness >= 50) %>%
    filter(checkM.contamination <= 10) %>%
    group_by(strain, treatment, buffer) %>%
    summarize(
        bins_count = n(),
        lineage_values = paste(checkM.lineage, collapse = "; "),
        completeness_values = paste(checkM.completeness, collapse = "; "),
        .groups = "drop"
    )

# Parse the "genome" column . # for coassembled 
mags_parsed_coass <- mags_coass %>%
    separate(genome, into = c("assembler", "binner", "group", "temp"), sep = "-", remove = FALSE) %>%
    separate(temp, into = c("strain", "treatment", "buffer_bin"), sep = "_", remove = FALSE) %>%
    separate(buffer_bin, into = c("buffer", "bin_id"), sep = "\\.", remove = FALSE)

# Summarize number of bins per strain, treatment, and buffer
summary_mags_coass <- mags_parsed_coass %>%
  filter(checkM.completeness >= 50) %>%
  filter(checkM.contamination <= 10) %>%
  group_by(strain, treatment, buffer) %>%
  summarize(
    bins_count = n(),
    lineage_values = paste(checkM.lineage, collapse = "; "),
    completeness_values = paste(checkM.completeness, collapse = "; "),
    .groups = "drop"
  )

summary_mags %>%
    arrange(buffer, strain, treatment)

summary_mags_coass %>%
    arrange(buffer, strain, treatment)
# Write the summary to a TSV file
write_tsv(summary_mags, paste0(folder_path,"gunc_checkm_summary_summary.tsv"))
write_tsv(summary_mags_coass, paste0(folder_path_coass,"gunc_checkm_summary_summary.tsv"))

#add origin info column  
summary_mags$origin <- "co-binning"
mags_parsed$origin <- "co-binning"
summary_mags_coass$origin <- "co-assembly"
mags_parsed_coass$origin <- "co-assembly"

# Combine the two summaries
combined_summary <- bind_rows(summary_mags, summary_mags_coass)
combined_mags <- bind_rows(mags_parsed, mags_parsed_coass)

# Write the combined summary to a TSV file
# write_tsv(combined_summary, paste0(folder_path,"gunc_checkm_summary_combined_summary.tsv"))

combined_mags_small <- combined_mags %>%
  group_by(strain, treatment, buffer, bin_id, origin) %>%
  summarize(
    completeness = checkM.completeness,
    contamination = checkM.contamination,
    lineage = checkM.lineage,
    .groups = "drop"
  ) %>%
  ungroup()

#clean up the data
# combined_mags_small$buffer <- ifelse(combined_mags_small$buffer == "D", "DESS", combined_mags_small$buffer)
# combined_mags_small$buffer <- ifelse(combined_mags_small$buffer == "P", "PBS", combined_mags_small$buffer)
# unique(combined_mags_small$buffer)
combined_mags_small$treatment <- ifelse(combined_mags_small$treatment == "B", "Benzonase", combined_mags_small$treatment)
# combined_mags_small$treatment <- ifelse(combined_mags_small$treatment == "Benzo", "Benzonase", combined_mags_small$treatment)
combined_mags_small$treatment <- ifelse(combined_mags_small$treatment == "BT", "Blood & Tissue Kit", combined_mags_small$treatment)
combined_mags_small$treatment <- ifelse(combined_mags_small$treatment == "MB", "Microbiome Kit", combined_mags_small$treatment)
combined_mags_small$treatment <- ifelse(combined_mags_small$treatment == "MBX", "Microbiome Kit + Bead beating", combined_mags_small$treatment)
combined_mags_small$treatment <- ifelse(combined_mags_small$treatment == "MBS", "Microbiome Kit + spinning", combined_mags_small$treatment)
unique(combined_mags_small$treatment)
combined_mags_small$strain <- ifelse(combined_mags_small$strain == "Acro", "Acropora", combined_mags_small$strain)
combined_mags_small$strain <- ifelse(combined_mags_small$strain == "Por", "Porites", combined_mags_small$strain)
combined_mags_small$strain <- ifelse(combined_mags_small$strain == "Poci", "Pocillopora", combined_mags_small$strain)
unique(combined_mags_small$strain)

#add grouping info
combined_mags_small$group <- ifelse(combined_mags_small$strain == "Acropora", "Coral", combined_mags_small$strain)
combined_mags_small$group <- ifelse(combined_mags_small$strain == "Porites", "Coral", combined_mags_small$group)
combined_mags_small$group <- ifelse(combined_mags_small$strain == "Pocillopora", "Coral", combined_mags_small$group)
combined_mags_small$group <- ifelse(combined_mags_small$strain == "H2", "Aiptasia", combined_mags_small$group)
combined_mags_small$group <- ifelse(combined_mags_small$strain == "F003", "Aiptasia", combined_mags_small$group)

# factor the data
combined_mags_small$buffer <- factor(combined_mags_small$buffer, levels = rev(c("DESS", "PBS")))
combined_mags_small$treatment <- factor(combined_mags_small$treatment, levels = c("Benzonase", "Blood & Tissue Kit", "Microbiome Kit", "Microbiome Kit + Bead beating", "Microbiome Kit + spinning"))
combined_mags_small$strain <- factor(combined_mags_small$strain, levels = rev(c("Acropora", "Porites", "Pocillopora", "F003", "H2")))

# make lineage factors alphabetically
#make column lineage_clean
combined_mags_small$lineage_clean <- gsub("[gpscfko]__", "", gsub(" \\(UID[0-9]{4}\\)", "",combined_mags_small$lineage))
combined_mags_small$lineage_clean <- gsub(" \\(UID[0-9]{3}\\)", "", combined_mags_small$lineage_clean)
combined_mags_small$lineage_clean <- gsub(" \\(UID[0-9]{2}\\)", "", combined_mags_small$lineage_clean)
combined_mags_small$lineage_clean <- gsub(" \\(UID[0-9]{1}\\)", "", combined_mags_small$lineage_clean)

# sort combined_mags_small by lineage_clean
# combined_mags_small$lineage_clean <- factor(combined_mags_small$lineage_clean, levels = rev(sort(unique(combined_mags_small$lineage_clean))))
combined_mags_small$lineage <- factor(combined_mags_small$lineage, levels = rev(sort(unique(combined_mags_small$lineage))))

# gsub("c__", "", unique(combined_mags_small$lineage))
#  unique(combined_mags_small$lineage))
# taxonomy <- unique(gsub("[gpscfko]__", "", gsub(" \\(UID[0-9]{4}\\)", "", unique(combined_mags_small$lineage))))
# taxonomy <- unique(gsub(" \\(UID[0-9]{3}\\)", "", taxonomy))
# taxonomy <- unique(gsub(" \\(UID[0-9]{2}\\)", "", taxonomy))
# taxonomy <- unique(gsub(" \\(UID[0-9]{1}\\)", "", taxonomy))

# unique(taxonomy)
# Define custom colors for each taxonomy
taxonomy_colors <- c(
    "Bacteria" = "#1f77b4",
    "Actinomycetales" = "#ff7f0e",
    "Rhizobiales" = "#2ca02c",
    "Rhodospirillales" = "#d62728",
    "Gammaproteobacteria" = "#9467bd",
    "algicola" = "#8c564b",
    "Mollicutes" = "#e377c2",
    "Rhodobacteraceae" = "#7f7f7f",
    "root" = "#bcbd22",
    "Archaea" = "#17becf",
    "Cyanobacteria" = "#aec7e8",
    "Alphaproteobacteria" = "#ffbb78",
    "Flavobacteriaceae" = "#98df8a",
    "Proteobacteria" = "#ff9896",
    "Deltaproteobacteria" = "#c5b0d5",
    "Vibrio" = "#c49c94",
    "Burkholderiaceae" = "#f7b6d2",
    "Bradyrhizobium" = "#c7c7c7",
    "Mycoplasma" = "#dbdb8d",
    "Cytophagales" = "#9edae5",
    "Bacteroidetes" = "#ffbb78"
)

# Apply custom colors to the plot
# Create a function to match taxonomy with colors
get_taxonomy_color <- function(lineage, taxonomy_colors) {
    for (taxon in names(taxonomy_colors)) {
        if (grepl(taxon, lineage, ignore.case = FALSE)) {
            return(taxonomy_colors[[taxon]])
        }
    }
    return(NA)  # Return NA if no match is found
}

# Apply the function to get colors for each lineage
combined_mags_small$lineage_color <- sapply(combined_mags_small$lineage, get_taxonomy_color, taxonomy_colors)
# rm(temp_list)
# rm(taxon_color)
# taxon_color <- c()
# for (i in 1:length(combined_mags_small$lineage)) {
#     temp_list <- paste('"',combined_mags_small$lineage[i],'"', " = ", '"', combined_mags_small$lineage_color[i],'"', sep = "")
#     taxon_color <- c(taxon_color, temp_list)
# }
# Convert taxon_color to a named vector
taxon_color_named <- setNames(combined_mags_small$lineage_color, combined_mags_small$lineage)
p_load("lemon")
p_load("ggforce")

combined_mags_small$grouping <- paste(combined_mags_small$buffer, combined_mags_small$strain, sep = " - ")

unique(combined_mags_small$grouping)

# set completeness value 
completeness_value <- 50

big_plot <- combined_mags_small %>%
        filter(completeness >= completeness_value) %>%
        filter(contamination <= 10) %>%
        ggplot(aes(x = treatment, fill = lineage)) +
        geom_bar(stat = "count") +
        facet_grid( ~ grouping , scales = "free") +#, nrow = 1) +
        theme_minimal(base_size = 12) +
        scale_fill_manual(values = taxon_color_named) 
big_plot + facet_wrap(origin ~ grouping, scales = "free_x", nrow = 1) 

extract_legend <- function(plot) {
    g <- ggplotGrob(plot)
    legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
    return(legend)
}

legend <- extract_legend(big_plot)
mainplot <- list()
for (o in c("co-assembly", "co-binning")) {
    # Filter the data for the current origin
        plot_list <- list()
        x <- 1
        for (i in c("PBS - H2","DESS - H2", "PBS - F003","DESS - F003","DESS - Acropora", "DESS - Pocillopora", "DESS - Porites")){

        data_plot <- combined_mags_small %>%
            filter(completeness >= completeness_value) %>%
            filter(origin == o) %>%
            filter(contamination <= 10) %>%
            filter(grouping == i)

        if(unique(data_plot$group) == "Aiptasia") {
            faceting_var <- " ~ grouping"
            scale_limit <- if (completeness_value == 50) ylim(0, 52) else ylim(0, 30)
            width_bar <- 0.8
        }
        else {
            faceting_var <- " ~ strain"
            scale_limit <- if (completeness_value == 50) ylim(0, 12) else ylim(0, 10)
            width_bar <- 0.5
        }
        p <- data_plot %>%
            ggplot(aes(x = treatment, fill = lineage)) +
            geom_bar(stat = "count", width = width_bar) +
            facet_wrap(faceting_var, scales = "free", nrow = 1) +
            theme_minimal(base_size = 12) +
            labs(
                title = "",
                x = "",
                y = "",
                fill = "Lineage",
            ) +
            scale_fill_manual(values = taxon_color_named) +
            scale_limit +
            scale_x_discrete(labels = ~ str_wrap(.x, width = 10)) +
            theme(strip.text.y = element_blank(),
            strip.text.x = element_text(size = 14, face = ifelse(i == "Aiptasia", "plain", "italic")),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            line = element_line(linewidth = 1),
            axis.line = element_line(colour = "black"),
            axis.ticks.length = unit(0.2 , "cm"))

        plot_list[[x]] <- p
        x <- x + 1
        }
        # httpgd::hgd()

        # p_load("ggpubr")
        fplot <- annotate_figure(ggarrange(
            ggarrange(plot_list[[1]] + theme(legend.position = "none"),
                        plot_list[[2]] + theme(legend.position = "none", axis.text.y = element_blank(), axis.line.y = element_blank()),
                        plot_list[[3]] + theme(legend.position = "none", axis.text.y = element_blank(), axis.line.y = element_blank()),
                        plot_list[[4]] + theme(legend.position = "none", axis.text.y = element_blank(), axis.line.y = element_blank()), nrow = 1),
            ggarrange(plot_list[[5]] + theme(legend.position = "none") ,# + geom_bar(width = 0.02),
                        plot_list[[6]] + theme(legend.position = "none", axis.text.y = element_blank(), axis.line.y = element_blank()),
                        plot_list[[7]] + theme(legend.position = "none", axis.text.y = element_blank(), axis.line.y = element_blank()), legend, nrow = 1), nrow = 2),
            fig.lab = paste0("MAGs summary (", o, ")"), left = "Number of MAGs", bottom = "Treatment")

    mainplot[[o]] <- fplot

}

mainplot
# Save to file 
out_path <- "/home/colinl/metaG/Git/metaG_EukDepletion/plots/"
ggsave(mainplot[[1]], filename = paste0(out_path, "MAGs_summary_",completeness_value,"_co-assembly.pdf"), width = 24, height = 16) # A4 dimensions in inches
ggsave(mainplot[[2]], filename = paste0(out_path, "MAGs_summary_",completeness_value,"_co-binning.pdf"), width = 24, height = 16)

# origin <- c("co-assembly", "co-binning")
# Plot with custom colors
p <- combined_mags_small %>%
    filter(completeness >= completeness_value) %>%
    filter(contamination <= 10) %>%
    # filter(strain != "F003") %>%
    # filter(strain != "H2") %>%
    ggplot(aes(x = treatment, fill = lineage)) +
    geom_bar(stat = "count") +
    facet_wrap(strain ~ buffer + origin, scales = "free_y", nrow = 6, ncol = 3) +
    theme_minimal() +
    labs(
        title = "Number of MAGs",
        x = "",
        y = "Number of MAGs",
        fill = "Lineage"
    ) +
    scale_fill_manual(values = taxon_color_named) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

shift_legend2(p+ coord_flip())
# Custom function to shift legend into empty facet panels
shift_legend2 <- function(p) {
    # check if p is a valid object
    if(!(inherits(p, "gtable"))){
        if(inherits(p, "ggplot")){
        gp <- ggplotGrob(p) # convert to grob
        } else {
        message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
        return(p)
        }
    } else {
        gp <- p
    }

    # check for unfilled facet panels
    facet.panels <- grep("^panel", gp[["layout"]][["name"]])
    empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]), 
                                USE.NAMES = F)
    empty.facet.panels <- facet.panels[empty.facet.panels]

    if(length(empty.facet.panels) == 0){
        message("There are no unfilled facet panels to shift legend into. Returning original plot.")
        return(p)
    }

    # establish name of empty panels
    empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
    names <- empty.facet.panels$name

    # return repositioned legend
    lemon::reposition_legend(p, 'center', panel=names)
}

# Get number of bins per strain, treatment, and buffer with completeness above 50
summary_mags <- combined_mags_small %>%
        filter(completeness >= 90) %>%
        group_by(strain, treatment, buffer, lineage) %>%
        summarize(
                bins_count = n(),
                .groups = "drop"
        )

ggplot(summary_mags, aes(x = treatment, y = bins_count, fill = treatment)) +
    geom_boxplot() +
    facet_wrap(~ buffer + strain, scales = "free_x", nrow = 1) +
    theme_minimal() +
    labs(
        title = "Completeness of MAGs",
        x = "",
        y = "Number of Bins"
    )

# test for significance and add it to the plot
summary_mags <- combined_mags_small %>%
    filter(completeness >= 50) %>%
    # filter(strain != "F003") %>%
    # filter(strain != "H2") %>%
    # filter(strain != "Acropora") %>%
    # filter(strain != "Porites") %>%
    # filter(strain != "Pocillopora") %>%
    group_by(strain, treatment, buffer, lineage) %>%
    summarize(
        bins_count = n(),
        .groups = "drop"
    )

# Perform statistical test for significance of treatment on bins_count
anova_results <- aov(bins_count ~  treatment + buffer, data = summary_mags)
summary(anova_results)

# lineage as factor sorted alphabetically
# summary_mags$lineage <- factor(summary_mags$lineage, levels = rev(sort(unique(summary_mags$lineage))))

# summary_mags %>%
#     ggplot(aes(x = treatment, y = bins_count, fill = lineage)) +
#     geom_bar(stat = "identity") +
#     facet_wrap(~ buffer + strain, scales = "free_x", nrow = 1) +
#     theme_minimal() +
#     labs(
#         title = "Completeness of MAGs",
#         x = "",
#         y = "Number of Bins",
#         fill = "Lineage"
#     ) +
#     scale_fill_manual(values = taxon_color_named) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))



summary_mags %>%
    # filter(completeness >= 90) %>%
    # filter(completeness >= 50) %>%
    # filter(strain != "F003") %>%
    # filter(strain != "H2") %>%
    ggplot(aes( y = bins_count, fill = treatment)) +
    geom_boxplot() + 
    facet_wrap(~ buffer + strain, scales = "free_x", nrow = 1) 


# Plot summary of MAGs as boxplot
combined_mags_small %>%
    # filter(completeness >= 90) %>%
    filter(completeness >= 50) %>%
    # filter(strain != "F003") %>%
    # filter(strain != "H2") %>%
ggplot(aes(x = treatment, fill = lineage)) +
    geom_bar(stat = "count") +
    facet_wrap(~ buffer + strain, scales = "free_x", nrow = 1) +
    theme_minimal() +
    labs(
        title = "Completeness of MAGs",
        x = "",
        y = "Number of Bins"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


###---- check virus summary ----###
path_pattern <- "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/binning/nf_core_MAG/results/VirusIdentification/geNomad/*"


# Read the TSV file 
read_and_combine_virus_summaries <- function(path_pattern) {
  # List all TSV files matching the pattern
  files <- list.files(path = dirname(path_pattern), pattern = "*_virus_summary.tsv", full.names = TRUE, recursive = TRUE)
  
  # Function to read a single TSV file and add the "name" column
  read_virus_summary <- function(file) {
    data <- read_tsv(file)
    data$name <- gsub("_virus_summary\\.tsv$", "", basename(file))
    return(data)
  }
  
  # Read and combine all TSV files
  combined_data <- bind_rows(lapply(files, read_virus_summary))
  
  return(combined_data)
}

# Use the function to read and combine virus summaries
viruses <- read_and_combine_virus_summaries(path_pattern)
# viruses <- read_tsv("/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/binning/nf_core_MAG/results/VirusIdentification/geNomad/group-Acro_B_DESS/Acro_concat_B_D_summary/Acro_concat_B_D_virus_summary.tsv")
head(viruses)

min(viruses$n_hallmarks)
summary_viruses <- viruses %>%
    count(name) %>%
        separate(name, into = c("strain", "ignore", "treatment", "buffer"), sep = "_", remove = FALSE) %>%
        mutate(buffer = ifelse(is.na(buffer), "D", buffer))
summary_viruses



clean_viruses <- viruses[,c("name",  "n_genes", "genetic_code", "virus_score", "n_hallmarks", "taxonomy")] %>%
    separate(name, into = c("strain", "ignore", "treatment", "buffer"), sep = "_", remove = FALSE) %>%
    mutate(buffer = ifelse(is.na(buffer), "D", buffer))


#clean up the data
clean_viruses$buffer <- ifelse(clean_viruses$buffer == "D", "DESS", clean_viruses$buffer)
clean_viruses$buffer <- ifelse(clean_viruses$buffer == "P", "PBS", clean_viruses$buffer)
unique(clean_viruses$buffer)
clean_viruses$treatment <- ifelse(clean_viruses$treatment == "B", "Benzo", clean_viruses$treatment)
clean_viruses$treatment <- ifelse(clean_viruses$treatment == "Benzo", "Benzonase", clean_viruses$treatment)
clean_viruses$treatment <- ifelse(clean_viruses$treatment == "BT", "Blood & Tissue Kit", clean_viruses$treatment)
clean_viruses$treatment <- ifelse(clean_viruses$treatment == "MB", "Microbiome Kit", clean_viruses$treatment)
clean_viruses$treatment <- ifelse(clean_viruses$treatment == "MBX", "Microbiome Kit + Bead beating", clean_viruses$treatment)
clean_viruses$treatment <- ifelse(clean_viruses$treatment == "MBS", "Microbiome Kit + spinning", clean_viruses$treatment)
unique(clean_viruses$treatment)
clean_viruses$strain <- ifelse(clean_viruses$strain == "Acro", "Acropora", clean_viruses$strain)
clean_viruses$strain <- ifelse(clean_viruses$strain == "Por", "Porites", clean_viruses$strain)
clean_viruses$strain <- ifelse(clean_viruses$strain == "Poci", "Pocillopora", clean_viruses$strain)
unique(clean_viruses$strain)

# factor the data
clean_viruses$buffer <- factor(clean_viruses$buffer, levels = c("DESS", "PBS"))
clean_viruses$treatment <- factor(clean_viruses$treatment, levels = c("Benzonase", "Blood & Tissue Kit", "Microbiome Kit", "Microbiome Kit + Bead beating", "Microbiome Kit + spinning"))
clean_viruses$strain <- factor(clean_viruses$strain, levels = c("Acropora", "Porites", "Pocillopora", "F003", "H2"))

# httpgd::hgd()

# Plot summary of viruses as boxplot
ggplot(clean_viruses, aes(x = treatment, fill = treatment)) +
    geom_histogram(stat = "count") +
    facet_wrap(~ buffer + strain, scales = "free") +
    theme_minimal() +
    labs(
        title = "Virus count by GeNomad",
        x = "",
        y = "Number of identified viruses"
    )
