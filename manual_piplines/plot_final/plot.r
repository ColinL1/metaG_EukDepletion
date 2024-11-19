#! /usr/bin/env Rscript

cat("Loading libraries...\n")
invisible(capture.output((library(pacman))))
p_load("tidyr", "tidyfast", "dplyr", "ggplot2", "ggpubr", "tidyverse", "scales", "RColorBrewer", "optparse", "pacman","httpgd", "grid", "stringr") # nolint: line_length_linter, cSpell. install lemon? 

## %% load custom functions
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

# Custom labeller function to italicise species names and leave buffer as is
custom_labeller <- function(labels) {
    if (is.atomic(labels)) {
        # If the input is a vector, process it directly
        labels <- sapply(labels, function(label) {
            parts <- unlist(strsplit(as.character(label), " - ", fixed = TRUE))
            if (parts[1] %in% c("H2", "F003")) {
                label
            } else {
                bquote(italic(.(parts[1])) ~ "- " ~ .(parts[2]))
            }
        })
    } else if (is.data.frame(labels)) {
        # If the input is a data frame, modify the value column
        labels$value <- sapply(labels$value, function(label) {
            parts <- unlist(strsplit(as.character(label), " - ", fixed = TRUE))
            if (parts[1] %in% c("H2", "F003")) {
                label
            } else {
                bquote(italic(.(parts[1])) ~ "- " ~ .(parts[2]))
            }
        })
    }
    return(labels)
}

#custom function to extract legend from ggplot object
extract_legend <- function(plot) {
    g <- ggplotGrob(plot)
    legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
    return(legend)
}

#custom function to read in multiple files and have the file name as a column
read_table_many <- function(file_path, ...) {
    data <- read.table(file_path, ...)
    data$file <- basename(file_path)
    return(data)
}

#start hgd  server
hgd()

# %% Read the file
data_reads <- read.table("/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/mapping/stats_reads_sanitised.txt", header = TRUE, sep = "\t")

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

# extractions methods  B   BT  MB  MBX MBS to Benzonase "Blood & Tissue" "Microbiome" "Microbiome & beadbeating" "Microbiome & spinning"
data_reads$extraction <- gsub("\\bB\\b", "Benzonase", data_reads$extraction)
data_reads$extraction <- gsub("\\bBT\\b", "Blood & Tissue", data_reads$extraction)
data_reads$extraction <- gsub("\\bMB\\b", "Microbiome", data_reads$extraction)
data_reads$extraction <- gsub("\\bMBX\\b", "Microbiome & bead-beating", data_reads$extraction)
data_reads$extraction <- gsub("\\bMBS\\b", "Microbiome & spinning", data_reads$extraction)

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

data$extraction <- factor(data$extraction, levels = (c("Blood & Tissue", "Benzonase", "Microbiome",  "Microbiome & bead-beating", "Microbiome & spinning")))

## make extractions groups to sort sample
data$extraction_groups <- data$extraction
data$extraction_groups <- gsub("Microbiome & bead-beating", "2", data$extraction_groups)
data$extraction_groups <- gsub("Microbiome & spinning", "2", data$extraction_groups)

data$extraction_groups <- gsub("Benzonase", "1", data$extraction_groups)
data$extraction_groups <- gsub("Blood & Tissue", "1", data$extraction_groups)
data$extraction_groups <- gsub("Microbiome", "1", data$extraction_groups)

data$extraction <- factor(data$extraction, levels = rev(c("Blood & Tissue", "Benzonase", "Microbiome",  "Microbiome & bead-beating", "Microbiome & spinning")))
# data$extraction <- factor(data$extraction, levels = rev(c("Blood & Tissue", "Benzonase", "Microbiome", " ",  "Microbiome & bead-beating", "Microbiome & spinning")))

## add number of replicate uniquer to each grouping and treatment to the df.
n_samples <- data %>%
    filter(origin == "reads") %>%
    group_by(grouping, extraction) %>%
    summarise(n_sample = n_distinct(replicate))

data <- left_join(data, n_samples, by = c("grouping", "extraction"))

# %% plot figure 2
spacer <- ggplot() + theme_void()

col_palette_2 <- c("Host" = "#D64933", "Symbiodiniaceae" = "#2A7F62", "Bacteria" = "#3498DB", "Other" = "#414288")# #EAC435 for host? 

SIZE_EXTRA_TEXT <- 6

plot_list <- list()
x <-1
for (i in unique(data$species)) {
    print(paste(x,": ",i, sep = "" ))
    p <- data %>%
        filter(origin == "reads") %>%
        filter(species == i) %>%
        filter(mapping != "total") %>%
        # filter(extraction != " ") %>%  # Exclude the last extraction value for group 2
        # filter(!(extraction_groups == 2 & extraction == " ")) %>%  # Exclude the last extraction value for group 2
        ggplot(aes(x = extraction, y = num_seqs, colour = mapping, fill = mapping)) +
        geom_bar(stat = "identity", position = "fill", width = 0.85) +  # Adjust width based on extraction value
        facet_grid(extraction_groups ~ grouping, scales = "free_y", space = "free", labeller = as_labeller(custom_labeller, label_parsed)) +  # Adjust facet by extraction_groups and grouping
        scale_y_continuous(labels = scales::percent) +  # Convert y-axis to percentage
        scale_fill_manual(values = col_palette_2, breaks = c("Host", "Symbiodiniaceae", "Other", "Bacteria")) + 
        scale_color_manual(values = col_palette_2, breaks = c("Host", "Symbiodiniaceae", "Other", "Bacteria")) +
        labs(x = "Extraction", y = "Percentage", colour = "Mapping", fill = "Mapping") +
        theme(plot.title = element_text(hjust = 0.5)) +
        coord_flip() +
        theme_minimal() +
        theme(
            strip.text.y = element_blank(),
            strip.text.x = element_text(size = 12),
            panel.grid.major = element_blank(),  # Hide major grid lines
            panel.grid.minor = element_blank(),  # Hide minor grid lines
            legend.key.size = unit(1.5, "lines"),  # Increase legend key size
            legend.text = element_text(margin = margin(t = 5, b = 5, l = 5, r = 5)),  # Increase padding in legend text
            panel.spacing.y = unit(1.2, "lines")  # Increase spacing between vertical panels
        ) +
        geom_text(aes(label = ifelse(extraction != " ", paste0("n=", n_sample), "")), angle = -90, y = 1.03, hjust = 0.5, size = (SIZE_EXTRA_TEXT - 2), color = "black")  # Add n_sample at the top of each bar, skip if extraction is " "
    
    plot_list[[x]] <- p
    x <- x + 1
}

#manually add NA to the missing extraction value in each subplot
# F003 - DESS
plot_list[[2]] <- plot_list[[2]] +
    geom_text(data = subset(data, grouping == "F003 - DESS" & extraction_groups == 2), aes(label = "NA", y = 0.5, x = "Microbiome & spinning"), color = "black", angle = 0, size = SIZE_EXTRA_TEXT) +
    geom_text(data = subset(data, grouping == "F003 - DESS" & extraction_groups == 1), aes(label = "NA", y = 0.5, x = "Microbiome"), color = "black", angle = 0, size = SIZE_EXTRA_TEXT) +
    geom_text(data = subset(data, grouping == "F003 - DESS" & extraction_groups == 1), aes(label = "NA", y = 0.5, x = "Benzonase"), color = "black", angle = 0, size = SIZE_EXTRA_TEXT)

#H2 - DESS
plot_list[[3]] <- plot_list[[3]] + 
    geom_text(data = subset(data, grouping == "H2 - DESS" & extraction_groups == 2), aes(label = "NA", y = 0.5, x = "Microbiome & spinning"), color = "black", angle = 0, size = SIZE_EXTRA_TEXT)

figure_2 <- ggarrange(
    annotate_figure(
        ggarrange(
            plot_list[[3]] + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank()),
            plot_list[[2]] + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
                                axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()),
            ncol = 2, nrow = 1, common.legend = TRUE, legend = "none", widths = c(1/2, 1/2.4)),
        fig.lab = "A",
        fig.lab.pos = "top.left",
        fig.lab.size = 14,
        fig.lab.face = "bold"
        ),
    spacer,
    annotate_figure(
        ggarrange(
            plot_list[[1]] + theme(axis.title.x = element_blank()),
            plot_list[[5]] + theme(axis.title.x = element_blank(),
                                axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()),
            plot_list[[4]] + theme(axis.title.x = element_blank(),
                                axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()),
            ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom", widths = c(1/3, 1/4, 1/4)),
        fig.lab = "B",
        fig.lab.pos = "top.left",
        fig.lab.size = 14,
        fig.lab.face = "bold"
        ), 
    ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom", heights = c(0.9,0.1,1))
figure_2

# ggsave(figure_2, filename = "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/plot_final/figure_2.png", width = 25, height = 15, units = "cm", dpi = 300)
# ggsave(figure_2, filename = "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/plot_final/figure_2.pdf", width = 25, height = 15, units = "cm", dpi = 300)

# Plot for figure 2 complete!

# %% plot replicate plots 
# invert order of extractions form before
data$extraction <- factor(data$extraction, levels = (c("Blood & Tissue", "Benzonase", "Microbiome",  "Microbiome & bead-beating", "Microbiome & spinning")))
# Loop through the different plots
plot_list_replicate <- list()
x <- 1
for (i in unique(data$grouping)) {
    # Filter data based on origin and mapping
    plot_data <- data %>%
        filter(origin == "reads") %>%
        filter(mapping != "total") %>%
        filter(extraction != " ") %>% # Exclude the last extraction value for group 2
        filter(grouping == i) 
    
    # Create the plot
    plot <- ggplot(plot_data, aes(x = replicate, y = num_seqs, colour = mapping, fill = mapping)) +
        geom_bar(stat = "identity", position = "fill")  +
        facet_wrap(vars(grouping, extraction), scales = "free_y", nrow = 1, labeller = labeller(grouping = as_labeller(custom_labeller, label_parsed), .default = label_value)) +
        theme_minimal() +
        # theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
        scale_fill_manual(values = col_palette_2, breaks = c("Host", "Symbiodiniaceae", "Other", "Bacteria")) + 
        scale_color_manual(values = col_palette_2, breaks = c("Host", "Symbiodiniaceae", "Other", "Bacteria")) +
        scale_x_discrete(labels = ~ str_wrap(.x, width = 20)) +
        scale_y_continuous(labels = scales::percent) +  # Convert y-axis to percentage
        labs(title = str_wrap(" ", width = 60), x = "Sample number", y = " ") +
        theme(plot.title = element_text(hjust = 0.5), legend.box = "horizontal") +
        guides(colour = guide_legend(title = "Mapping"), fill = guide_legend(title = "Mapping"))
    
    plot_list_replicate[[x]] <- plot + coord_flip()
    x <- x + 1
}

annotate_figure(
    ggarrange(
                        plot_list_replicate[[5]],
                        plot_list_replicate[[4]],
                        plot_list_replicate[[3]],
                        plot_list_replicate[[2]],
                        plot_list_replicate[[1]],
                        plot_list_replicate[[7]],
                        plot_list_replicate[[6]],
                        ncol = 1, nrow = 7, common.legend = TRUE, legend = "bottom", heights = c(1,1,1,1,1,1)),
        top = text_grob("Stacked bar plot of read percentages, using mapping as fill and displaying extractions, species, and buffer as facets.", face = "bold", size = 14))



###---- figure 3 ----###
# %% load data figure 3 # co-assembly plot 
#read co-assembly data ! not filtered by mapping
data_co_assembly <- read.table("/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/CAT/results_download/NR/summary/superkingdom.csv", header = TRUE, sep = ",")

#cleanup data_co_assembly file column removing concat
data_co_assembly$File.Name <- gsub("_concat_", "_", data_co_assembly$File.Name)

#replace NA with Unkclassified
data_co_assembly$Superkingdom <- gsub("<NA>", "Unclassified", data_co_assembly$Superkingdom)
data_co_assembly[is.na(data_co_assembly$Superkingdom),]$Superkingdom <- "Unclassified" 

#split file coluimn at the _ to get sample informaiton
data_co_assembly <- data_co_assembly %>% separate(File.Name, c("species", "extraction", "buffer"), sep = "_", remove = FALSE)

# fix species names
data_co_assembly$species <- gsub("\\bF3\\b", "F003", data_co_assembly$species)
data_co_assembly$species <- gsub("\\bAcro\\b", "Acropora", data_co_assembly$species)
data_co_assembly$species <- gsub("\\bAc\\b", "Acropora", data_co_assembly$species)
data_co_assembly$species <- gsub("\\bPo\\b", "Pocillopora", data_co_assembly$species)
data_co_assembly$species <- gsub("\\bPoci\\b", "Pocillopora", data_co_assembly$species)
data_co_assembly$species <- gsub("\\bPor\\b", "Porites", data_co_assembly$species)
data_co_assembly$species <- gsub("\\bPr\\b", "Porites", data_co_assembly$species)

# fix extraction names
data_co_assembly$extraction <- gsub("\\bB\\b", "Benzonase", data_co_assembly$extraction)
data_co_assembly$extraction <- gsub("\\bBT\\b", "Blood & Tissue", data_co_assembly$extraction)
data_co_assembly$extraction <- gsub("\\bMB\\b", "Microbiome", data_co_assembly$extraction)
data_co_assembly$extraction <- gsub("\\bMBX\\b", "Microbiome & bead-beating", data_co_assembly$extraction)
data_co_assembly$extraction <- gsub("\\bMBS\\b", "Microbiome & spinning", data_co_assembly$extraction)

# fix buffer names
data_co_assembly$buffer <- gsub("\\bP\\b", "PBS", data_co_assembly$buffer)
data_co_assembly$buffer <- gsub("\\bD\\b", "DESS", data_co_assembly$buffer)

#set buffer, methods, species as factors
data_co_assembly$species <- factor(data_co_assembly$species, levels = rev(c("Porites", "Pocillopora", "Acropora", "F003", "H2")))
data_co_assembly$buffer <- factor(data_co_assembly$buffer, levels = rev(sort(unique(data_co_assembly$buffer))))
data_co_assembly$extraction <- factor(data_co_assembly$extraction, levels = rev(sort(unique(data_co_assembly$extraction))))

#create grouping column
data_co_assembly$grouping <- paste(data_co_assembly$species, data_co_assembly$buffer, sep = " - ")

#set grouping as factor
data_co_assembly$grouping <- factor(data_co_assembly$grouping, levels = c("H2 - PBS", "H2 - DESS", "F003 - PBS", "F003 - DESS", "Acropora - DESS", "Porites - DESS", "Pocillopora - DESS"))

data_co_assembly$Superkingdom <- factor(data_co_assembly$Superkingdom, levels = rev(c("Bacteria", "Eukaryota", "Viruses", "Archaea", "no support", "Unclassified")))


# %% plot figure 3
col_palette_3 <- c("Bacteria" = "#5C415D", "Eukaryota" = "#1B998B", "Viruses" = "#D5C67A", "Archaea" ="#04A777" , "no support" = "#D8E4FF", "Unclassified" = "#034748" )# #EAC435 for host? 
breaks <- c("Bacteria", "Eukaryota", "Viruses", "Archaea", "no support", "Unclassified")

#calculate xend for each species buffer to use in geom_segment for line under facet names
data_co_assembly <- data_co_assembly %>%
                    filter(extraction != "Benzoase") %>%
                    filter(extraction != "Blood & Tissue") %>%
                    group_by(species, buffer) %>%
                    mutate(xend = length(unique(extraction)))

# plot for figure 3
plot_list_co_assembly <- list()
x <- 1
for (i in unique(data_co_assembly$species)) {
    print(paste(x, ": ", i, sep = ""))
    p <- data_co_assembly %>%
        filter(extraction != "Benzoase") %>%
        filter(extraction != "Blood & Tissue") %>%
        filter(species == i) %>%
        ggplot(aes(x = extraction, y = Number.of.Contigs, colour = Superkingdom, fill = Superkingdom)) +
        geom_bar(stat = "identity", position = "fill", width = 0.85) +
        facet_grid(buffer ~ species, scales = "free_y") +
        theme_minimal() +
        scale_fill_manual(
        values = col_palette_3,
        breaks = breaks) +
        scale_color_manual(
        values = col_palette_3,
        breaks = breaks) +
        scale_x_discrete(labels = ~ str_wrap(.x, width = 20)) +
        scale_y_continuous(labels = scales::percent) +
        labs(title = str_wrap(" ", width = 60), x = " ", y = " ") +
        theme(
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 14, face = ifelse(i %in% c("F003", "H2"), "plain", "italic")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(1.5, "lines"),
        panel.spacing.y = unit(1.2, "lines"),
        legend.text = element_text(margin = margin(t = 5, b = 5, l = 5, r = 5)),
        plot.title = element_text(hjust = 0.5),
        legend.box = "horizontal",
        strip.background = element_blank(),
        # plot.margin = margin(t = 10, r = 10, b = 10, l = 20)  # Increase margin on x-axis side
        ) +
        coord_flip() +
        guides(colour = guide_legend(title = "Mapping"), fill = guide_legend(title = "Mapping"))

    if (i %in% c("F003", "H2")) {
        p <- p + geom_segment(aes(x = ifelse(grouping == "F003 - DESS", 0.75, 1), xend = ifelse(grouping == "F003 - DESS", 1.25, xend), y = 1.05, yend = 1.05), color = "black", size = 0.85) + theme(strip.text.y = element_text(angle = 0, size = 14))
    }

    plot_list_co_assembly[[x]] <- p
    x <- x + 1
}

# Extract the legend from plot_list_co_assembly[[1]]
legend <- extract_legend(plot_list_co_assembly[[1]] + theme(legend.position = "bottom"))

figure_3 <- ggarrange(
            ggarrange(
                annotate_figure(
                    plot_list_co_assembly[[2]] + theme(axis.title.x = element_blank(), legend.position = "none"),
                    fig.lab = "A",
                    fig.lab.pos = "top.left",
                    fig.lab.size = 14,
                    fig.lab.face = "bold"),
                annotate_figure(
                    plot_list_co_assembly[[4]]+ theme(axis.title.x = element_blank(), legend.position = "none",
                                        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()),
                    fig.lab = "B",
                    fig.lab.pos = "top.left",
                    fig.lab.size = 14,
                    fig.lab.face = "bold"),
                    ncol = 2, nrow = 1, widths = c(1,0.8)),
            ggarrange(
                annotate_figure(
                    plot_list_co_assembly[[5]] + theme(axis.title.x = element_blank(), legend.position = "none"),
                fig.lab = "C",
                    fig.lab.pos = "top.left",
                    fig.lab.size = 14,
                    fig.lab.face = "bold"),
                annotate_figure(
                    plot_list_co_assembly[[3]] + theme(axis.title.x = element_blank(), legend.position = "none",
                                        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()),
                fig.lab = "D",
                    fig.lab.pos = "top.left",
                    fig.lab.size = 14,
                    fig.lab.face = "bold"),
                annotate_figure(
                    plot_list_co_assembly[[1]] + theme(axis.title.x = element_blank(), legend.position = "none",
                                        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()),
                    fig.lab = "E",
                    fig.lab.pos = "top.left",
                    fig.lab.size = 14,
                    fig.lab.face = "bold"),
                ncol = 3, nrow = 1, widths = c(1,0.9,0.9)),
            legend,
            ncol = 1, nrow = 3, heights = c(1,0.9,0.1))
figure_3
# ggsave(figure_3, filename = "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/plot_final/figure_3.png", width = 25, height = 15, units = "cm", dpi = 300)
# ggsave(figure_3, filename = "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/plot_final/figure_3.pdf", width = 25, height = 15, units = "cm", dpi = 300)
# Plot for figure 3 complete!
#TODO: Make figure 3 corals less wide. fix heights proportion

# %% load data figure 4
# non filtered microbial only files
files <- list.files(path = "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/CAT/results_download/NR/summary/bacteria_only_contig_class", pattern = "*.summary.txt", full.names = TRUE)

# Create an empty data frame to store the results
CAT_data_family <- data.frame()

# Loop through each file and read the data into a data frame
for (x in files) {
    # Read the data from the file
    data <- read_table_many(x, header = TRUE, sep = "\t", col.names = c("rank","clade", "number of contigs", "number of ORFs", "number of positions"))    
    # Append the data to the results data frame
    CAT_data_family <- rbind(CAT_data_family, data)
}

# clean up names in file
CAT_data_family$file <- gsub(".summary.txt", "", CAT_data_family$file)
CAT_data_family$file <- gsub("_concat_", "_", CAT_data_family$file)

# filter for rank == "family"
CAT_data_family <- CAT_data_family[CAT_data_family$rank == "family", ]

#keep only clade number.of.contigss columns
CAT_data_family <- CAT_data_family[, c("file", "clade", "number.of.contigs")]

# rename clade to species
colnames(CAT_data_family)[2] <- "Family"

#rename number.of.contigs to number_of_contigs
colnames(CAT_data_family)[3] <- "Number_of_contigs"

#rename file to sample
colnames(CAT_data_family)[1] <- "Sample"

# Calculate the total number of contigs for each Family
family_totals <- CAT_data_family %>%
    group_by(Family) %>%
    summarise(total_contigs = sum(Number_of_contigs)) %>%
    arrange(desc(total_contigs))

# Get the top 20 families
top_families <- family_totals$Family[1:23]

# Filter for the top 10 families and group the rest as "Other"
CAT_data_family_curated <- CAT_data_family %>%
    mutate(Family = ifelse(Family %in% top_families, Family, "Other"))

# rename "no support" in Family to "no match"
CAT_data_family_curated$Family <- gsub("no support", "No match", CAT_data_family_curated$Family)

# rename NA in Family to "ambigous match"
CAT_data_family_curated$Family[is.na(CAT_data_family_curated$Family)] <- "Ambigous match"

# Summarise the number of contigs by sample and Family
CAT_data_family_curated <- CAT_data_family_curated %>%
    group_by(Sample, Family) %>%
    summarise(total_contigs = sum(Number_of_contigs)) %>%
    arrange(Sample, desc(total_contigs))

CAT_data_family_curated <- separate(CAT_data_family_curated, Sample, into = c("Strain", "Treatment", "Buffer"), sep = "_", remove = FALSE)

#rename Treatment 
CAT_data_family_curated$Treatment <- gsub("^MBS$", "Microbiome kit + Spinning", CAT_data_family_curated$Treatment)
CAT_data_family_curated$Treatment <- gsub("^MBX$", "Microbiome kit + Bead beating", CAT_data_family_curated$Treatment)
CAT_data_family_curated$Treatment <- gsub("^MB$", "Microbiome kit", CAT_data_family_curated$Treatment)
CAT_data_family_curated$Treatment <- gsub("^B$", "Benzonase", CAT_data_family_curated$Treatment)
CAT_data_family_curated$Treatment <- gsub("^BT$", "Blood & Tissue", CAT_data_family_curated$Treatment)

CAT_data_family_curated$Strain <- gsub("\\bPor\\b", "Porites", CAT_data_family_curated$Strain)
CAT_data_family_curated$Strain <- gsub("\\bPoci\\b", "Pocillopora", CAT_data_family_curated$Strain)
CAT_data_family_curated$Strain <- gsub("\\bAcro\\b", "Acropora", CAT_data_family_curated$Strain)
CAT_data_family_curated$Strain <- gsub("\\bH2\\b", "H2", CAT_data_family_curated$Strain)
CAT_data_family_curated$Strain <- gsub("\\bF003\\b", "F003", CAT_data_family_curated$Strain)

CAT_data_family_curated$Buffer <- gsub("\\bD\\b", "DESS", CAT_data_family_curated$Buffer)
CAT_data_family_curated$Buffer <- gsub("\\bP\\b", "PBS", CAT_data_family_curated$Buffer)

# Sort Family by number of contigs, then put "Other" and "ambigous match" last
family_order <- CAT_data_family_curated %>%
    group_by(Family) %>%
    summarise(total_contigs = sum(total_contigs)) %>%
    arrange(desc(total_contigs)) %>%
    pull(Family)

# Move "Other" and "ambigous match" to the end
family_order <- c(setdiff(family_order, c("Other", "Ambigous match", "No match")), "Other", "Ambigous match", "No match")

# Reorder the factor levels
CAT_data_family_curated$Family <- factor(CAT_data_family_curated$Family, levels = family_order)

# Reorder factor levels for buffer
CAT_data_family_curated$Buffer <- factor(CAT_data_family_curated$Buffer, levels = c("PBS", "DESS"))

# add plotting category called strain_treatment to plot by strain and treatment
CAT_data_family_curated["grouping"] <- paste(CAT_data_family_curated$Strain, CAT_data_family_curated$Buffer, sep = " - ")

# %% plot figure 4
# color_palette <- c("#CB7A79", "#ED816F", "#FFA7A6", "#FFB94F", "#FFDF41", "#EAC88D", "#D6C3A2", "#A2CB7E", "#C0DE75", "#5E90A1", "#74B2C9", "#80B2A9", "#A0E0E0", "#99AED0", "#9C80B9", "#A372A0", "#B4A0EB", "#D0AAEC", "#D7BCE7", "#D897B1", "#F2A8C5")
color_palette <- c("#CB7A79", "#ED816F", "#FFA7A6", "#FFB94F", "#FFDC7B", "#FFDF41", "#EAC88D", "#D6C3A2", "#BED29D", "#A2CB7E", "#C0DE75", "#5E90A1", "#74B2C9", "#80B2A9", "#A0E0E0", "#ACE6FB", "#99AED0", "#9C80B9", "#A372A0", "#B4A0EB", "#D0AAEC", "#D7BCE7", "#D897B1", "#F2A8C5")
# color_palette <- c("#CB7A79", "#ED816F", "#FFA7A6", "#FFB94F", "#FFDC7B", "#FFDF41", "#EAC88D", "#D6C3A2", "#BED29D", "#A2CB7E", "#c4c5c3", "#5E90A1", "#74B2C9", "#80B2A9", "#A0E0E0", "#ACE6FB", "#99AED0", "#9C80B9", "#A372A0", "#B4A0EB", "#D0AAEC", "#D7BCE7", "#D897B1", "#F2A8C5")
# Use the color palette in the plot


# do plot list for all strains
plot_family_list <- list()
n <- 1
for (i in sort(unique(CAT_data_family_curated$Strain))) {
    print(i)
        p <- CAT_data_family_curated %>%
            filter(Strain == i) %>%
            filter(Family != "No match") %>%
            filter(Treatment != "Benzonase") %>%
            filter(Treatment != "Blood & Tissue") %>%
            filter(grouping != "H2 - DESS") %>%
            filter(grouping != "F003 - DESS") %>%
            ggplot(aes(x = Treatment, y = total_contigs, fill = Family)) +
            geom_bar(stat = "identity", position = "fill", width = 1, col = "black") +
            facet_wrap(~ Strain, scales = "free", nrow = 1) +
            labs(title = "", x = "", y = "") +
            scale_fill_manual(values = color_palette) +
            scale_y_continuous(labels = scales::percent) +
            theme_minimal() +
            theme(
                legend.position = "none",
                # panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                #strip.text = element_text(size = 12)
                strip.text.x = element_text(size = 12, face = ifelse(i %in% c("F003", "H2"), "plain", "italic")),
                axis.text.x = element_text(size = 12, angle = 45, hjust = 0.7)
                ) +
            scale_x_discrete(labels = function(x) str_wrap(x, width = 15))
    
    plot_family_list[[n]] <- p
    n <- n + 1
}

legend <- extract_legend(plot_family_list[[1]] + theme(legend.position = "right", legend.key.size = unit(1.5, "lines"), legend.text = element_text(size = 12), legend.title = element_text(size = 14)))

figure_4 <- ggarrange(
        annotate_figure(
            ggarrange(
                ggarrange(
                    annotate_figure(
                        plot_family_list[[3]],
                        fig.lab = "A",
                        fig.lab.pos = "top.left",
                        fig.lab.size = 14,
                        fig.lab.face = "bold"),
                    annotate_figure(
                        plot_family_list[[2]],
                        fig.lab = "B",
                        fig.lab.pos = "top.left",
                        fig.lab.size = 14,
                        fig.lab.face = "bold"),
                    ncol = 2, nrow = 1, widths = c(0.8, 0.8)),
                ggarrange(
                    annotate_figure(
                        plot_family_list[[1]],
                        fig.lab = "C",
                        fig.lab.pos = "top.left",
                        fig.lab.size = 14,
                        fig.lab.face = "bold"),
                    annotate_figure(
                        plot_family_list[[4]],
                        fig.lab = "D",
                        fig.lab.pos = "top.left",
                        fig.lab.size = 14,
                        fig.lab.face = "bold"),
                    annotate_figure(
                        plot_family_list[[5]],
                        fig.lab = "E",
                        fig.lab.pos = "top.left",
                        fig.lab.size = 14,
                        fig.lab.face = "bold"),
                    ncol = 3, nrow = 1),
            ncol = 1, nrow = 2, heights = c(1,1)),
            top = text_grob("Top 20 bacteria families", face = "bold", size = 14)),
        legend,
        ncol = 2, nrow = 1, widths= c(1.5,0.4))
figure_4


#TODO: Make figure 4 less tall 


# %% load data figure 4 dots (increases the number of family kept) !! not FINAL!!
# non filtered microbial only files
files <- list.files(path = "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/CAT/results_download/NR/summary/bacteria_only_contig_class", pattern = "*.summary.txt", full.names = TRUE)

# Create an empty data frame to store the results
CAT_data_family <- data.frame()

# Loop through each file and read the data into a data frame
for (x in files) {
    # Read the data from the file
    data <- read_table_many(x, header = TRUE, sep = "\t", col.names = c("rank","clade", "number of contigs", "number of ORFs", "number of positions"))    
    # Append the data to the results data frame
    CAT_data_family <- rbind(CAT_data_family, data)
}

# clean up names in file
CAT_data_family$file <- gsub(".summary.txt", "", CAT_data_family$file)
CAT_data_family$file <- gsub("_concat_", "_", CAT_data_family$file)

# filter for rank == "family"
CAT_data_family <- CAT_data_family[CAT_data_family$rank == "family", ]

#keep only clade number.of.contigss columns
CAT_data_family <- CAT_data_family[, c("file", "clade", "number.of.contigs")]

# rename clade to species
colnames(CAT_data_family)[2] <- "Family"

#rename number.of.contigs to number_of_contigs
colnames(CAT_data_family)[3] <- "Number_of_contigs"

#rename file to sample
colnames(CAT_data_family)[1] <- "Sample"

# Calculate the total number of contigs for each Family
family_totals <- CAT_data_family %>%
    group_by(Family) %>%
    summarise(total_contigs = sum(Number_of_contigs)) %>%
    arrange(desc(total_contigs))

# Get the top 40 families
top_families <- family_totals$Family[1:43]

# Filter for the top 10 families and group the rest as "Other"
CAT_data_family_curated <- CAT_data_family %>%
    mutate(Family = ifelse(Family %in% top_families, Family, "Other"))

# rename "no support" in Family to "no match"
CAT_data_family_curated$Family <- gsub("no support", "No match", CAT_data_family_curated$Family)

# rename NA in Family to "ambigous match"
CAT_data_family_curated$Family[is.na(CAT_data_family_curated$Family)] <- "Ambigous match"

# Summarise the number of contigs by sample and Family
CAT_data_family_curated <- CAT_data_family_curated %>%
    group_by(Sample, Family) %>%
    summarise(total_contigs = sum(Number_of_contigs)) %>%
    arrange(Sample, desc(total_contigs))

CAT_data_family_curated <- separate(CAT_data_family_curated, Sample, into = c("Strain", "Treatment", "Buffer"), sep = "_", remove = FALSE)

#rename Treatment 
CAT_data_family_curated$Treatment <- gsub("^MBS$", "Microbiome kit + Spinning", CAT_data_family_curated$Treatment)
CAT_data_family_curated$Treatment <- gsub("^MBX$", "Microbiome kit + Bead beating", CAT_data_family_curated$Treatment)
CAT_data_family_curated$Treatment <- gsub("^MB$", "Microbiome kit", CAT_data_family_curated$Treatment)
CAT_data_family_curated$Treatment <- gsub("^B$", "Benzonase", CAT_data_family_curated$Treatment)
CAT_data_family_curated$Treatment <- gsub("^BT$", "Blood & Tissue", CAT_data_family_curated$Treatment)

CAT_data_family_curated$Strain <- gsub("\\bPor\\b", "Porites", CAT_data_family_curated$Strain)
CAT_data_family_curated$Strain <- gsub("\\bPoci\\b", "Pocillopora", CAT_data_family_curated$Strain)
CAT_data_family_curated$Strain <- gsub("\\bAcro\\b", "Acropora", CAT_data_family_curated$Strain)
CAT_data_family_curated$Strain <- gsub("\\bH2\\b", "H2", CAT_data_family_curated$Strain)
CAT_data_family_curated$Strain <- gsub("\\bF003\\b", "F003", CAT_data_family_curated$Strain)

CAT_data_family_curated$Buffer <- gsub("\\bD\\b", "DESS", CAT_data_family_curated$Buffer)
CAT_data_family_curated$Buffer <- gsub("\\bP\\b", "PBS", CAT_data_family_curated$Buffer)

# Sort Family by number of contigs, then put "Other" and "ambigous match" last
family_order <- CAT_data_family_curated %>%
    group_by(Family) %>%
    summarise(total_contigs = sum(total_contigs)) %>%
    arrange(desc(total_contigs)) %>%
    pull(Family)

# Move "Other" and "ambigous match" to the end
family_order <- c(setdiff(family_order, c("Other", "Ambigous match", "No match")), "Other", "Ambigous match", "No match")

# Reorder the factor levels
CAT_data_family_curated$Family <- factor(CAT_data_family_curated$Family, levels = rev(family_order))

# Reorder factor levels for buffer
CAT_data_family_curated$Buffer <- factor(CAT_data_family_curated$Buffer, levels = c("PBS", "DESS"))

# add plotting category called strain_treatment to plot by strain and treatment
CAT_data_family_curated["grouping"] <- paste(CAT_data_family_curated$Strain, CAT_data_family_curated$Buffer, sep = " - ")

# %% plot figure 4 dot version 
CAT_data_family_curated$Strain <- factor(CAT_data_family_curated$Strain, levels = c("F003", "H2", "Acropora", "Porites", "Pocillopora"))

spacer <- ggplot() + theme_void()

# define a 5 color palette with pastel tones
color_palette_x <- c("#20B2AA", "#BBBE64", "#2A2A72")

plot_family_list_dot <- list()
n <- 1
for (i in sort(unique(CAT_data_family_curated$Strain))) {
    print(i)
        p <- CAT_data_family_curated %>%
        filter(Strain == i) %>%
        filter(Family != "No match") %>%
        filter(Treatment != "Benzonase") %>%
        filter(Treatment != "Blood & Tissue") %>%
        filter(grouping != "H2 - DESS") %>%
        filter(grouping != "F003 - DESS") %>%
        ggplot(aes(x = Treatment, y = Family, color = Treatment)) +
        geom_point(aes(size = log10(total_contigs) +1)) +
        facet_wrap(~ Strain, ncol = 1) +
        labs(title = "", x = "", y = "") +
        scale_color_manual(values = color_palette_x) +
        # scale_y_continuous(labels = scales::percent) +
        theme_minimal() +
        theme(
            legend.position = "none",
            panel.grid.minor = element_blank(),
            strip.text.x = element_text(size = 12, face = ifelse(i %in% c("F003", "H2"), "plain", "italic")),
            axis.text.x = element_text(size = 12, angle = 45, hjust = 0.7)
        ) +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 15))

    plot_family_list_dot[[n]] <- p
    n <- n + 1
}

legend <- extract_legend(plot_family_list_dot[[1]] + theme(legend.position = "bottom", legend.key.size = unit(1.5, "lines"), legend.text = element_text(size = 12), legend.title = element_text(size = 14)))

ggarrange(
    # ggarrange(
        plot_family_list_dot[[2]],
        plot_family_list_dot[[1]],
        # ncol = 2, nrow = 1, widths = c(1, 1)),
    # ggarrange(
        plot_family_list_dot[[3]],
        plot_family_list_dot[[5]],
        plot_family_list_dot[[4]],
        # ncol = 3, nrow = 1),
    spacer, spacer, legend,
    ncol = 5, nrow = 2, heights = c(1.2, 0.1))


# CAT_data_family_curated %>%
#         # filter(Strain == i) %>%
#         filter(Family != "No match") %>%
#         filter(Treatment != "Benzonase") %>%
#         filter(Treatment != "Blood & Tissue") %>%
#         filter(grouping != "H2 - DESS") %>%
#         filter(grouping != "F003 - DESS") %>%
#         ggplot(aes(x = Treatment, y = Family, color = Treatment)) +
#         geom_point(aes(size = (total_contigs))) +
#         facet_wrap(~ Strain, nrow = 1) +
#         labs(title = "", x = "", y = "") +
#         scale_color_manual(values = color_palette_x) +
#         # scale_y_continuous(labels = scales::percent) +
#         theme_minimal() +
#         theme(
#             legend.position = "none",
#             panel.grid.minor = element_blank(),
#             strip.text.x = element_text(size = 12, face = ifelse(grepl("F003|H2", CAT_data_family_curated$Strain), "plain", "italic")),
#             axis.text.x = element_text(size = 12, angle = 45, hjust = 0.7)
#         ) +
#         scale_x_discrete(labels = function(x) str_wrap(x, width = 15))


#%% qPCR -  figure 1
library(tidyr)
library(dplyr)
library(ggplot2)
    library(reshape2)
library(egg)
# install.packages("stringr")          # Install stringr package
library("stringr")                   # Load stringr only to fit labels in two row 

### read the table
# setwd("~/PhD/Meta_G_qPCR/R_paper")

Aip_PBS <-read.table("/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/qPCR/MetaG_qPCR_Aip_PBS.csv" , header = T, sep =",") 

PBSvsDESS <-read.table("/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/qPCR/MetaG_qPCR_PBSvsDESS.csv" , header = T, sep =",") 
PBSvsDESS <- PBSvsDESS %>% 
    filter(Buffer == "DESS")

Corals <-read.table("/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/qPCR/MetaG_qPCR_corals.csv" , header = T, sep =",") 

data <- rbind(Aip_PBS, PBSvsDESS, Corals)
data$Ct <- as.numeric(data$Ct)

data_summary <- dcast(formula = Sample_name~Gene, fun.aggregate = mean,data = data)

data_summary <- data_summary[complete.cases(data_summary), ]

data_summary$dCT_16S <- -1*(data_summary$`16S` - data_summary$`b-actin`)
data_summary <- tidyr::separate(data_summary,Sample_name,into =c("Strain", "replicate", "Buffer", "Treatment"),sep = "_",remove = FALSE,extra = "merge")

#change form abbreviation to full name extractions
data_summary$Treatment <- gsub("MB_", "Microbiome Kit_1", data_summary$Treatment)
data_summary$Treatment <- gsub("SP_", "Sonication + PMA_", data_summary$Treatment)
data_summary$Treatment <- gsub("SB_", "Sonication + Benzonase_", data_summary$Treatment)
data_summary$Treatment <- gsub("BT_", "Blood & Tissue Kit_", data_summary$Treatment)
data_summary$Treatment <- gsub("P_", "PMA_", data_summary$Treatment)
data_summary$Treatment <- gsub("B_", "Benzonase_", data_summary$Treatment)

#fix strain names
data_summary$Strain <- gsub("\\bPori\\b", "Porites", data_summary$Strain)
data_summary$Strain <- gsub("\\bPoci\\b", "Pocillopora", data_summary$Strain)
data_summary$Strain <- gsub("\\bAcro\\b", "Acropora", data_summary$Strain)

data_summary <- tidyr::separate(data_summary, Treatment, into =c("Treatment", "Run"),sep = "_", remove = FALSE, extra = "merge")

#order treatmens in the way we want them on the plot (native first)
data_summary$Treatment <- factor(data_summary$Treatment, levels=c("Blood & Tissue Kit", "Microbiome Kit", "Benzonase", "Sonication + Benzonase", "PMA","Sonication + PMA"))
factor(data_summary$Treatment) #to check whether it worked          "Benzonase", "Blood & Tissue Kit", "Microbiome Kit", "PMA", "Sonication + Benzonase", "Sonication + PMA"

#calculate ddCT mean and standard error per group (box plots in ggplot already calculate mean and sd, for bar plots we have to do it manually)
summary_16S <- data_summary %>% group_by (Strain, Buffer, Run, Treatment) %>% 
  summarise(ratio_mean = mean(dCT_16S), 
            ratio_se = sd(dCT_16S)/sqrt(length(dCT_16S)))
summary_16S$Strain <- factor(summary_16S$Strain, levels=c( "H2", "F003", "Acropora", "Pocillopora", "Porites" ))
# backup_summary_16S <- summary_16S
#add grouping column
summary_16S["grouping"] <- paste(summary_16S$Strain, summary_16S$Buffer, sep = " - ")

# compare differnces between summary_16S and backup_summary_16S
summary_16S %>% anti_join(backup_summary_16S)

#%% plot figure 1
color_palette_1 <- c("Blood & Tissue Kit" = "#FF9000", "Microbiome Kit" = "#D05353", "Benzonase" = "#3A6EA5", "Sonication + Benzonase" = "#88B7B5", "PMA" = "#8FAD88", "Sonication + PMA" = "#5D4A66")
breaks <- c("Blood & Tissue Kit", "Benzonase", "Sonication + Benzonase", "PMA", "Sonication + PMA", "Microbiome Kit")
summary_16S$Treatment <- factor(summary_16S$Treatment, levels = breaks)

qPCR_plot_list <- list()
n <- 1
for (i in unique(summary_16S$grouping)) {
    print(paste("Plot ", n,": ", i, sep = ""))
    p <- summary_16S %>%
            filter(grouping == i) %>%
        ggplot(aes(x = Treatment, y = ratio_mean, fill = Treatment)) +
            geom_col(width=0.5)+
            geom_errorbar(aes(
                ymin = ratio_mean - ratio_se, 
                ymax = ratio_mean + ratio_se),
                width = .2, 
                position=position_dodge(.9), 
                color="black") +
            facet_wrap( ~ grouping, labeller = as_labeller(custom_labeller, label_parsed))+
            geom_hline(yintercept = 0,  linetype="dashed", color = "black")+
            scale_fill_manual(values = color_palette_1) +
            scale_x_discrete(labels = ~ str_wrap(.x, width = 10), breaks = breaks) +
            theme_bw() +
            labs(x= "", y= "- dCT (16S - b-actin)") +
            theme(legend.position="none", line= element_line(size = 1),
                axis.line = element_line(colour = "black"),
                axis.ticks.length = unit(0.2 , "cm"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                strip.background = element_blank(), 
                strip.text.x = element_text(color = "black", size = 12, angle = 0, hjust = 0, vjust = 0.5, face = "plain"),
                panel.spacing = unit(3, "lines")) +
                coord_cartesian (ylim=c(-15,6))

    
    qPCR_plot_list[[n]] <- p
    n <- n + 1
}

legend  <- extract_legend(qPCR_plot_list[[5]] + theme(legend.position = "bottom", legend.key.size = unit(1.5, "lines"), legend.text = element_text(size = 12), legend.title = element_text(size = 14), legend.direction = "horizontal")+ guides(fill = guide_legend(nrow = 1)))

figure_1 <- ggarrange(
        ggarrange(
            annotate_figure(
                qPCR_plot_list[[5]],
                fig.lab = "A",
                fig.lab.pos = "top.left",
                fig.lab.size = 14,
                fig.lab.face = "bold"),
            annotate_figure(
                qPCR_plot_list[[3]] + labs(y = ""),
                fig.lab = "B",
                fig.lab.pos = "top.left",
                fig.lab.size = 14,
                fig.lab.face = "bold"),
            ncol = 2, nrow = 1, widths = c(1, 1)),
        ggarrange(
            annotate_figure(
                qPCR_plot_list[[4]],
                fig.lab = "C",
                fig.lab.pos = "top.left",
                fig.lab.size = 14,
                fig.lab.face = "bold"),
            annotate_figure(
                qPCR_plot_list[[2]] + labs(y = ""),
                fig.lab = "D",
                fig.lab.pos = "top.left",
                fig.lab.size = 14,
                fig.lab.face = "bold"),
            ncol = 2, nrow = 1, widths = c(1, 1)),
        ggarrange(
            annotate_figure(
                qPCR_plot_list[[1]],
                fig.lab = "E",
                fig.lab.pos = "top.left",
                fig.lab.size = 14,
                fig.lab.face = "bold"),
            annotate_figure(
                qPCR_plot_list[[6]] + labs(y = ""),
                fig.lab = "F",
                fig.lab.pos = "top.left",
                fig.lab.size = 14,
                fig.lab.face = "bold"),
            annotate_figure(
                qPCR_plot_list[[7]] + labs(y = ""),
                fig.lab = "G",
                fig.lab.pos = "top.left",
                fig.lab.size = 14,
                fig.lab.face = "bold"),
            ncol = 3, nrow = 1, widths = c(1, 1)),
        legend,
        ncol = 1, nrow = 4, heights = c(0.8,0.8,1,0.1))

    annotate_figure(figure_1, 
    bottom = text_grob("*Ct Values from multiple qPCR run have been combined in plotting when required", color = "black",
            hjust = 1, x = 1, face = "italic", size = 10))


annotate_figure(
    ggarrange(
        ggarrange(
            annotate_figure(
                qPCR_plot_list[[5]],
                fig.lab = "A",
                fig.lab.pos = "top.left",
                fig.lab.size = 14,
                fig.lab.face = "bold"),
            annotate_figure(
                qPCR_plot_list[[3]] + labs(y = ""),
                fig.lab = "B",
                fig.lab.pos = "top.left",
                fig.lab.size = 14,
                fig.lab.face = "bold"),
            ncol = 2, nrow = 1, widths = c(1, 1)),
        ggarrange(
            annotate_figure(
                qPCR_plot_list[[4]],
                fig.lab = "C",
                fig.lab.pos = "top.left",
                fig.lab.size = 14,
                fig.lab.face = "bold"),
            annotate_figure(
                qPCR_plot_list[[2]] + labs(y = ""),
                fig.lab = "D",
                fig.lab.pos = "top.left",
                fig.lab.size = 14,
                fig.lab.face = "bold"),
            annotate_figure(
                qPCR_plot_list[[1]] + labs(y = ""),
                fig.lab = "E",
                fig.lab.pos = "top.left",
                fig.lab.size = 14,
                fig.lab.face = "bold"),
            annotate_figure(
                qPCR_plot_list[[6]] + labs(y = ""),
                fig.lab = "F",
                fig.lab.pos = "top.left",
                fig.lab.size = 14,
                fig.lab.face = "bold"),
            annotate_figure(
                qPCR_plot_list[[7]] + labs(y = ""),
                fig.lab = "G",
                fig.lab.pos = "top.left",
                fig.lab.size = 14,
                fig.lab.face = "bold"),
            ncol = 5, nrow = 1, widths = c(1,1,1,1,1)),
        legend,
        ncol = 1, nrow = 3, heights = c(0.8,0.8,0.1)),
        bottom = text_grob("*Ct Values from multiple qPCR run have been combined in plotting when required", color = "black",
                                   hjust = 1, x = 1, face = "italic", size = 10))

##### old plots 
###plots

annotate_figure(grid.arrange(a,
                             b,
                             c, ncol = 1) ,
                top = text_grob("", color = "black", face = "bold", size = 14),
                left = text_grob("-(Ct 16S - Ct B-actin)",color = "black", rot = 90, size = 16, vjust = 0.6),
                bottom = text_grob("*Ct Values from multiple qPCR run have been combined in plotting when required", color = "black",
                                   hjust = 1, x = 1, face = "italic", size = 10),
                right = "",  )




#%% other stuff to 
# ####### random shit #######
# head(data)


# ## test for signifcance data num_seq per mapping
# # data %>%
# #     filter(origin == "assembly") %>%
# #     filter(mapping != "total") %>%
# #     group_by(mapping) %>%
# #     summarise(mean_num_seqs = mean(num_seqs), sd_num_seqs = sd(num_seqs), n = n()) %>%
# #     mutate(se = sd_num_seqs / sqrt(n)) %>%
# #     mutate(ci = qt(0.975, n - 1) * se) %>%
# #     mutate(lwr = mean_num_seqs - ci, upr = mean_num_seqs + ci) %>%
# #     arrange(mean_num_seqs)


# # data <- rbind(data_assembly, data_reads)

# # head(data)
# # data["grouping"] <- paste(data$species, data$buffer, sep = " - ")
# # data$grouping <- factor(data$grouping, levels = c("H2 - PBS", "H2 - DESS", "F003 - PBS", "F003 - DESS", "Acropora - DESS", "Porites - DESS", "Pocillopora - DESS"))


# library(vegan)

# # Create a community matrix (for simplicity, we'll use num_seqs as OTU counts)
# community_matrix <- as.matrix(data$num_seqs)
# rownames(community_matrix) <- data$sample

# # Calculate Bray-Curtis distance
# dist_matrix <- vegdist(community_matrix, method = "bray")

# # Run PERMANOVA
# adonis_result <- adonis2(dist_matrix ~ extraction + grouping, data = data, permutations = 9999, by = "margin")

# # Print the results
# print(adonis_result)




# # Assuming 'counts' is a vector of read counts
# glm_model <- glm(counts ~ group + other_covariates, 
#                     family = negative.binomial(theta = 1),
#                     data = sample_data)
# summary(glm_model)

# # Load necessary library
# library(MASS)  # for negative binomial GLM

# # Fit a negative binomial GLM
# glm_model <- glm.nb(num_seqs ~ grouping + extraction + origin, data = filter(data, mapping == "bacteria"))

# # View summary of the model
# summary(glm_model)

# # Perform analysis of deviance to assess significance of factors
# anova(glm_model, test = "Chisq")




# # Loop through the different plots
# plot_list <- list()
# for (i in 1:6) {
#     # Filter data based on origin and mapping
#     plot_data <- data %>%
#         filter(origin == ifelse(i <= 3, "assembly", "reads")) %>%
#         filter(mapping != "total")
#     based_on <- ifelse(i <= 3, "assembly", "reads")
#     # Set the y variable and plot type based on the value of i
#     if ((i %% 3 == 1) & (i <= 3)) {
#         y_var <- "sum_len"
#         plot_type <- geom_boxplot()
#         plot_title <- paste("Boxplot of",y_var,"column with extractions and strain as facets, based on", based_on, sep = " ")
#         y_label <- "Sum Length"
#     } else if ((i %% 3 == 1) & (i >= 3)) {
#         y_var <- "num_seqs"
#         plot_type <- geom_boxplot()
#         plot_title <- paste("Boxplot of",y_var,"column with extractions and strain as facets, based on", based_on, sep = " ")
#         y_label <- "num_seqs"
#     } else if (i %% 3 == 2) {
#         y_var <- "sum_len"
#         plot_type <- geom_bar(stat = "identity", position = "fill") 
#         plot_title <- paste("Stacked Barplot of",y_var,"with mapping as fill and extractions and species as facets, based on", based_on, sep = " ")
#         y_label <- "Sum Length"
#     } else {
#         y_var <- "num_seqs"
#         plot_type <- geom_bar(stat = "identity", position = "fill")  
#         plot_title <- paste("Stacked Barplot of",y_var,"with mapping as fill and extractions and species as facets, based on", based_on, sep = " ")
#         y_label <- "Num seqs"
#     }
    
#     # Create the plot
#     plot <- ggplot(plot_data, aes(x = extraction, y = !!sym(y_var), colour = mapping, fill = mapping)) +
#         plot_type +
#         facet_wrap(~ grouping, scale = "free", nrow = 2) +
#         theme_minimal() +
#         # theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
#         scale_x_discrete(labels = ~ str_wrap(.x, width = 20)) +
#         labs(title = str_wrap(plot_title, width = 60), x = "Extraction", y = y_label) +
#         theme(plot.title = element_text(hjust = 0.5)) +
#         guides(colour = guide_legend(title = "Mapping"), fill = guide_legend(title = "Mapping"))

#     if (i %% 3 == 1) {
#         plot_list[[i]] <- plot
#     } else {
#         plot_list[[i]] <- plot #+ coord_flip()
#     }
# }


# plot_list
# # plot_list[[2]]
# shift_legend2(plot_list[[2]])
# # Combine the plots into a single plot using grid.arrange from the gridExtra package
# gridExtra::grid.arrange(grobs = plot_list[1:3], ncol = 1)
# gridExtra::grid.arrange(grobs = plot_list[4:6], ncol = 1)


# # Combine the plots into a single plot using grid.arrange from the gridExtra package
# gridExtra::grid.arrange(grobs = plot_list[1:3], mylegend, ncol = 3, common.legend = TRUE, legend="bottom")
# gridExtra::grid.arrange(grobs = plot_list[4:6], ncol = 3, common.legend = TRUE, legend="bottom")


# Aip_plot <- data %>%
#         filter(origin == "reads") %>% #reads 
#         filter(mapping != "total") %>%
#         filter(species != "Pocillopora") %>%
#         filter(species != "Porites") %>%
#         filter(species != "Acropora") %>%
#             ggplot( aes(x = extraction, y = !!sym(y_var), colour = mapping, fill = mapping)) +
#                 plot_type +
#                 facet_wrap(buffer ~ species, scale = "free_x", ncol = 1) +
#                 theme_minimal() +
#                 # theme(axis.text.y = element_text(angle = 45, hjust = 1)) +
#                 scale_x_discrete(labels = ~ str_wrap(.x, width = 20)) + 
#                 labs(title = str_wrap(plot_title, width = 60), x = "Extraction", y = y_label) +
#                 theme(plot.title = element_text(hjust = 0.5))
# Aip_plot
# # #extract legend
# # #https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
# # g_legend<-function(a.gplot){
# #     tmp <- ggplot_gtable(ggplot_build(a.gplot))
# #     leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
# #     legend <- tmp$grobs[[leg]]
# #     return(legend)}

# # mylegend <- g_legend(plot_list[[1]])

# # gridExtra::grid.arrange(arrangeGrob(plot_list[[1]] + theme(legend.position="none"),
# #                                     plot_list[[2]] + theme(legend.position="none"),
# #                                     plot_list[[3]] + theme(legend.position="none"),
# #                         nrow=1), mylegend , nrow=2, heights=c(10, 1))
    

# # grid_arrange_shared_legend <- function(...) {
# #     plots <- list(...)
# #     g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
# #     legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
# #     lheight <- sum(legend$height)
# #     grid.arrange(
# #         do.call(arrangeGrob, lapply(plots, function(x)
# #             x + theme(legend.position = "none"))),
# #         legend,
# #         ncol = 3,
# #         heights = unit.c(unit(1, "npc") - lheight, lheight))
# # }

# # lemon::grid_arrange_shared_legend(plot_list[[1]],plot_list[[2]],plot_list[[3]], nrow = 1, position = "bottom", plot = TRUE)
# # lemon::grid_arrange_shared_legend(plot_list[[4]],plot_list[[5]],plot_list[[6]], nrow = 1, position = "bottom", plot = TRUE)

# # lemon::grid_arrange_shared_legend(plot_list[[1]],plot_list[[4]], nrow = 1, position = "bottom", plot = TRUE)
# # lemon::grid_arrange_shared_legend(plot_list[[6]],plot_list[[2]], nrow = 1, position = "bottom", plot = TRUE)

# # plot_list[[2]] + coord_flip()
# # shift_legend2(plot_list[[2]])
# # shift_legend2(plot_list[[6]])
# # plot_list

#     # grobs = plot_list[4:6], ncol = 3, common.legend = TRUE, legend="bottom")

# # #test if statemnts
# # for (i in 1:6) {
# #     # Filter data based on origin and mapping
# #     plot_data <- data %>%
# #         filter(origin == ifelse(i <= 3, "assembly", "reads")) %>%
# #         filter(mapping != "total")
# #     based_on <- ifelse(i <= 3, "assembly", "reads")
# #     # Set the y variable and plot type based on the value of i
# #     if ((i %% 3 == 1) & (i <= 3)) {
# #         y_var <- "sum_len"
# #         plot_type <- geom_boxplot()
# #         plot_title <- paste("Boxplot of",y_var,"column with extractions and strain as facets, based on", based_on, sep = " ")
# #         y_label <- "Sum Length"
# #         print(i)
# #         print(y_var)
# #         print(based_on)
# #         print(y_label)
# #     } else if ((i %% 3 == 1) & (i >= 3)) {
# #         y_var <- "num_seqs"
# #         plot_type <- geom_boxplot()
# #         plot_title <- paste("Boxplot of",y_var,"column with extractions and strain as facets, based on", based_on, sep = " ")
# #         y_label <- "num_seqs"
# #         print(i)
# #         print(y_var)
# #         print(based_on)
# #         print(y_label)
# #     } #else if (i %% 3 == 2) {
# #     #     y_var <- "sum_len"
# #     #     plot_type <- geom_bar(stat = "identity", position = "fill") 
# #     #     plot_title <- paste("Stacked Barplot of sum_len with mapping as fill and extractions and species as facets, based on", based_on, sep = " ")
# #     #     y_label <- "Sum Length"
# #     #     print(i)
# #     #     print(y_var)
# #     #     print(based_on)
# #     #     print(y_label)
# #     # } else {
# #     #     y_var <- "num_seqs"
# #     #     plot_type <- geom_bar(stat = "identity", position = "fill")  
# #     #     plot_title <- paste("Stacked Barplot of num_seqs with mapping as fill and extractions and species as facets, based on", based_on, sep = " ")
# #     #     y_label <- "Num seqs"
# #     #     print(i)
# #     #     print(y_var)
# #     #     print(based_on)
# #     #     print(y_label)
# #     # }
# # }

