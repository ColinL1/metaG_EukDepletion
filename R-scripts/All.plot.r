#! /usr/bin/env Rscript

# =============================================================================
# SETUP AND CONFIGURATION
# =============================================================================

cat("Loading libraries...\n")
invisible(capture.output((library(pacman))))
p_load("tidyr", "tidyfast", "dplyr", "ggplot2", "ggpubr", "tidyverse", 
       "scales", "RColorBrewer", "optparse", "pacman", "httpgd", "grid", 
       "stringr", "reshape2", "ggh4x")

rm(list = ls())

# Define base paths
base_path <- "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines"
output_path <- "/home/colinl/Proj/metaG_EukDepletion/plots"

# Define color palettes
col_palette_mapping <- c(
    "Host" = "#D64933", 
    "Symbiodiniaceae" = "#2A7F62", 
    "Bacteria" = "#3498DB", 
    "Other" = "#414288"
)

col_palette_treatment <- c(
    "BT" = "#FFB94F", 
    "BenzoBT" = "#D97474", 
    "MB" = "#9AE0EC", 
    "BBMB" = "#87ACE9", 
    "SMB" = "#9C80BA"
)

col_palette_qpcr <- c(
    "Blood & Tissue Kit" = "#FFB94F", 
    "Microbiome Kit" = "#9AE0EC", 
    "Benzonase + BT" = "#D97474", 
    "Sonication + Benzonase + BT" = "#F2A7C5", 
    "PMA + BT" = "#82A37A", 
    "Sonication + PMA + BT" = "#AFDD84"
)

col_palette_mb_treatments <- c(
    "Microbiome Kit" = "#9AE0EC", 
    "Bead-beating + Microbiome Kit" = "#87ACE9", 
    "Spinning + Microbiome Kit" = "#9C80BA"
)

# =============================================================================
# CUSTOM FUNCTIONS
# =============================================================================

# Data cleaning helper functions
clean_species_names <- function(species_col) {
    species_col %>%
        gsub("\\bF3\\b", "F003", .) %>%
        gsub("\\bAcro\\b|\\bAc\\b", "Acropora", .) %>%
        gsub("\\bPo\\b|\\bPoci\\b", "Pocillopora", .) %>%
        gsub("\\bPor\\b|\\bPr\\b", "Porites", .)
}

clean_extraction_names <- function(extraction_col) {
    extraction_col %>%
        gsub("\\bB\\b", "BenzoBT", .) %>%
        gsub("\\bBT\\b", "BT", .) %>%
        gsub("\\bMB\\b", "MB", .) %>%
        gsub("\\bMBX\\b", "BBMB", .) %>%
        gsub("\\bMBS\\b", "SMB", .)
}

clean_buffer_names <- function(buffer_col) {
    buffer_col %>%
        gsub("\\bP\\b", "PBS", .) %>%
        gsub("\\bD\\b", "DESS", .)
}

clean_mapping_names <- function(mapping_col) {
    mapping_col %>%
        gsub("\\bscleractina\\b|\\baiptasia\\b", "Host", .) %>%
        gsub("\\bsymbiodiniaceae\\b", "Symbiodiniaceae", .) %>%
        gsub("\\bbacteria\\b", "Bacteria", .) %>%
        gsub("\\bother\\b", "Other", .)
}

expand_treatment_names <- function(treatment_col) {
    treatment_col %>%
        gsub("^MBS$", "Spinning + Microbiome Kit", .) %>%
        gsub("^MBX$", "Bead-beating + Microbiome Kit", .) %>%
        gsub("^MB$", "Microbiome Kit", .) %>%
        gsub("^B$", "Benzonase + BT", .) %>%
        gsub("^BT$", "Blood & Tissue Kit", .)
}

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

# =============================================================================
# DATA LOADING
# =============================================================================

data_reads <- read.table(
    file.path(base_path, "reads/mapping/stats_reads_sanitised.txt"), 
    header = TRUE, sep = "\t"
)
# Read re-run of Kaiju annotation with updated database 20250730
data_reads_kaiju <- read.table("/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/202505_kaiju_re-check/split_fq/fq/stats_reads_sanitised.txt", header = TRUE, sep = "\t")

# =============================================================================
# DATA PROCESSING AND CLEANING
# =============================================================================

# Replace values from data_reads with data_reads_kaiju if file matches
data_reads <- data_reads %>%
    left_join(data_reads_kaiju, by = "file", suffix = c("", ".kaiju")) %>%
    mutate(num_seqs = ifelse(is.na(num_seqs.kaiju), num_seqs, num_seqs.kaiju),
           sum_len = ifelse(is.na(sum_len.kaiju), sum_len, sum_len.kaiju),
           min_len = ifelse(is.na(min_len.kaiju), min_len, min_len.kaiju),
           avg_len = ifelse(is.na(avg_len.kaiju), avg_len, avg_len.kaiju),
           max_len = ifelse(is.na(max_len.kaiju), max_len, max_len.kaiju)) %>%
    select(-ends_with(".kaiju"))

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

# Clean up column names using helper functions
data_reads$species <- clean_species_names(data_reads$species)
data_reads$extraction <- clean_extraction_names(data_reads$extraction)
data_reads$buffer <- clean_buffer_names(data_reads$buffer)
data_reads$mapping <- clean_mapping_names(data_reads$mapping)

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

data$extraction <- factor(data$extraction, levels = (c("BT", "BenzoBT", "MB",  "BBMB", "SMB")))

## make extractions groups to sort sample
data$extraction_groups <- data$extraction
data$extraction_groups <- gsub("BBMB", "2", data$extraction_groups)
data$extraction_groups <- gsub("SMB", "2", data$extraction_groups)

data$extraction_groups <- gsub("BenzoBT", "1", data$extraction_groups)
data$extraction_groups <- gsub("BT", "1", data$extraction_groups)
data$extraction_groups <- gsub("MB", "1", data$extraction_groups)

data$extraction <- factor(data$extraction, levels = rev(c("BT", "BenzoBT", "MB",  "BBMB", "SMB")))
# data$extraction <- factor(data$extraction, levels = rev(c("BT", "BenzoBT", "MB", " ",  "BBMB", "SMB")))

## add number of replicate uniquer to each grouping and treatment to the df.
n_samples <- data %>%
    filter(origin == "reads") %>%
    group_by(grouping, extraction) %>%
    summarise(n_sample = n_distinct(replicate))

data <- left_join(data, n_samples, by = c("grouping", "extraction"))



# =============================================================================
# FIGURE 2: BOXPLOTS OF READ MAPPING PERCENTAGES
# =============================================================================

test_data <- data
test_data$mapping <- factor(test_data$mapping, levels = rev(c("total", "Other", "Bacteria", "Symbiodiniaceae", "Host")))
test_data$type <- ifelse(test_data$species %in% c("H2", "F003"), "aiptasia", "corals")
test_data$type <- factor(test_data$type, levels = c("aiptasia", "corals")) 

test_data <- test_data %>%
    select(sample, grouping, type, species, replicate, extraction, buffer, mapping, num_seqs) %>%
    filter(mapping != "total")

# Calculate percentage per sample
test_data <- test_data %>%
    group_by(sample) %>%
    mutate(num_seq_percent = num_seqs / sum(num_seqs, na.rm = TRUE) * 100) %>%
    ungroup()

# Statistical comparisons
for (i in c("Host", "Symbiodiniaceae", "Bacteria", "Archaea", "Other")) {
  for (j in c("aiptasia", "corals")) {
    plot_data <- test_data %>%
      filter(mapping != "total") %>%
      filter(mapping == i) %>%
      filter(type == j)
    print(paste0("Running stats for ", i, " in ", j))
    tryCatch({
      print(compare_means(num_seq_percent ~ extraction + grouping,
                          plot_data,
                          method = "t.test",
                          paired = FALSE))
    }, error = function(e) {
      message("Skipping comparison due to error: ", e$message)
    })
  }
}

test_data$mapping <- factor(test_data$mapping, levels = rev(c("total", "Other", "Archaea", "Bacteria", "Symbiodiniaceae", "Host")))
test_data$type <- ifelse(test_data$species %in% c("H2", "F003"), "aiptasia", "corals")
test_data$type <- factor(test_data$type, levels = c("aiptasia", "corals"))
test_data$extraction <- factor(test_data$extraction, levels = c("BT", "BenzoBT", "MB",  "BBMB", "SMB"))
col_palette_2 <- c("Host" = "#D64933", "Symbiodiniaceae" = "#2A7F62", "Bacteria" = "#3498DB", "Archaea" = "#414288", "Other" = "#414288")# #EAC435 for host? 
my_comparisons <- list(c("BT", "BenzoBT"), c("BT", "MB"), c("BT", "BBMB"), c("BT", "SMB"))

col_palette_2 <- c("Host" = "#D64933", "Symbiodiniaceae" = "#2A7F62", "Bacteria" = "#3498DB", "Other" = "#414288")# #EAC435 for host?
col_palette_2 <- c("BT" = "#FFB94F", "BenzoBT" = "#D97474", "MB" = "#9AE0EC", "BBMB" ="#87ACE9", "SMB" = "#9C80BA")# #EAC435 for host?

plot_list <- list()
SCALE <- "fixed" # "fixed" or "free"
# SCALE <- "fixed"# or "free"
# for (j in c("aiptasia", "corals" )) {
rm(stat_results_save)
for (i in sort(unique(test_data$grouping))) {
  for (j in c("Host", "Symbiodiniaceae", "Bacteria")) {
      print(paste0("Running plot for ", j, " in ", i))
      
      # Reset my_comparisons for each new loop
      my_comparisons <- combn(levels(droplevels(test_data$extraction)), 2, simplify = FALSE)
      # my_comparisons <- list(c("BT", "BenzoBT"), c("BT", "MB"), c("BT", "BBMB"), c("BT", "SMB"))

      plot_data <- test_data %>%
          filter(mapping == j)  %>%
          filter(grouping == i)

      # Filter comparison to avoid issues 
      group_counts <- plot_data %>%
        count(extraction) %>%
        filter(n >= 2)

      # Filter comparisons
      valid_comparisons <- my_comparisons[
        sapply(my_comparisons, function(x) all(x %in% group_counts$extraction))
      ]

      # Perform tests and extract significant comparisons
      stat_results <- compare_means(
        num_seq_percent ~ extraction,
        filter(plot_data, extraction %in% group_counts$extraction),
        method = "t.test",
        comparisons = valid_comparisons
      ) 
      # Add fake mapping for plotting
      stat_results$mapping <- j
      stat_results$grouping <- i
      
      if (exists("stat_results_save")) {
      stat_results_save <- rbind(stat_results_save, stat_results)
      } else {
          stat_results_save <- stat_results
      }
      
      stat_results <- stat_results %>%
      filter(group1 == "BT") %>% # only keep comparisons starting with BT
      filter(p.adj < 0.05) 
      

      stat_results$y_lab_max <- max(round(plot_data$num_seq_percent, 0), na.rm = TRUE) + 3 # Add some space for labels
      stat_results$y <- sapply(stat_results$group2, function(grp) {
        max(plot_data$num_seq_percent[plot_data$extraction == grp], na.rm = TRUE) + 3
      })

      # rname stat_results columns for geom_text
      colnames(stat_results)[which(colnames(stat_results) == "group1")] <- "control"
      colnames(stat_results)[which(colnames(stat_results) == "group2")] <- "extraction"

      #TODO: FIXED SCALE SWITCHING OFF

          plot <- plot_data %>%
              ggplot(aes(x = extraction, y = num_seq_percent, fill = extraction)) +
              geom_boxplot() +
                      facet_grid2(
                      grouping ~ mapping,
                      scales = "free_y",
                      independent = "y",
                      space = "free_x",
                      axes = "all",
                      # switch = "y",
                      labeller = labeller(
                        .rows = as_labeller(
                        custom_labeller,
                        label_parsed
                        ),
                        .cols = function(x) paste(x, "reads")
                      )
                      ) + 
              labs(title = "",
              x = "",
              y = "") +
              scale_fill_manual(values = col_palette_2) +
              theme_minimal(base_size = 24) +
              scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
              scale_y_continuous(labels = function(x) paste0(x, "%")) +
              # geom_text(
              #     data = stat_results,
              #     aes(x = group2, y = y_lab_max, label = p.signif),
              #     size = 6
              # ) + 
              theme(
              legend.position = "none",
              line = element_line(linewidth = 1),
              axis.line = element_line(color = "black"),
              axis.ticks.length = unit(0.2, "cm"),
              axis.ticks.x = element_line(color = "black"),
              axis.ticks.y = element_line(color = "black"),
              axis.text.x = element_text(color = "black", size = 16, angle = 0),
              axis.text.y = element_text(color = "black", size = 16),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              strip.background = element_blank(), 
              strip.text.x = element_text(color = "black", angle = 0, hjust = 0, vjust = 0.5, face = "plain"),
              plot.margin = margin(0.1, 0.1, 0.1, 0.1),           # Reduce plot margins (top, right, bottom, left)
              panel.spacing = unit(0.002, "lines")          # Reduce spacing between facets
      )
      if (SCALE == "fixed" & j == "Host") {
              plot <- plot + geom_text(
                  data = stat_results,
                  aes(x = extraction, y = 85, label = p.signif),
                  size = 6
              ) + scale_y_continuous(labels = function(x) paste0(x, "%") , breaks = seq(0, 80, 20), limits = c(0, 93))
      } else if (SCALE == "fixed" & j == "Symbiodiniaceae") {
        if (i == "Porites - DESS") {
                scale_factor <- 10
                plot <- plot + geom_text(
                    data = stat_results,
                    aes(x = extraction, y = scale_factor, label = p.signif),
                    size = 6
                ) + scale_y_continuous(labels = function(x) paste0(x, "%"), breaks = seq(0, scale_factor, 2), limits = c(0, scale_factor + 1))
        } else {
              plot <- plot + geom_text(
                  data = stat_results,
                  aes(x = extraction, y = 85, label = p.signif),
                  size = 6
              ) + scale_y_continuous(labels = function(x) paste0(x, "%") , breaks = seq(0, 80, 20), limits = c(0, 86))
            }
      } else if (SCALE == "fixed" & j == "Bacteria") {
        if (i == "H2 - PBS" | i == "F003 - PBS") {
              plot <- plot + geom_text(
                  data = stat_results,
                  aes(x = extraction, y = 40, label = p.signif),
                  size = 6
              ) + scale_y_continuous(labels = function(x) paste0(x, "%") , breaks = seq(0, 40, 10), limits = c(0, 41))
        } else if (i == "Acropora - DESS") {
              scale_factor <- 60
              plot <- plot + geom_text(
                  data = stat_results,
                  aes(x = extraction, y = scale_factor, label = p.signif),
                  size = 6
              ) + scale_y_continuous(labels = function(x) paste0(x, "%"), breaks = seq(0, scale_factor, 20), limits = c(0, scale_factor + 1))
      } else if (i == "Porites - DESS") {
              plot <- plot + geom_text(
                  data = stat_results,
                  aes(x = extraction, y = 60, label = p.signif),
                  size = 6
              ) + scale_y_continuous(labels = function(x) paste0(x, "%") , breaks = seq(0, 60, 20), limits = c(0, 61))
      } else if (i == "Pocillopora - DESS") {
              scale_factor <- 10
              plot <- plot + geom_text(
                  data = stat_results,
                  aes(x = extraction, y = scale_factor, label = p.signif),
                  size = 6
              ) + scale_y_continuous(labels = function(x) paste0(x, "%") , breaks = seq(0, scale_factor, 2), limits = c(0, scale_factor + 1))
        }
      } else {
        plot <- plot + geom_text(
          data = stat_results,
          aes(x = extraction, y = y_lab_max, label = p.signif),
          size = 6
          ) 
      }
      # plot
      x <- paste0(j, " - ", i)
      plot_list[[x]] <- plot


      # stat_results_save <- rbind(stat_results_save, stat_results)
      }
    }

figure_2 <- ggpubr::ggarrange(
    plot_list$`Host - H2 - PBS` + theme(axis.text.x = element_blank(), strip.text.y = element_blank()),
    plot_list$`Symbiodiniaceae - H2 - PBS`+ theme(axis.text.x = element_blank(), strip.text.y = element_blank()),
    plot_list$`Bacteria - H2 - PBS` + theme(axis.text.x = element_blank()),
    plot_list$`Host - F003 - PBS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank()),
    plot_list$`Symbiodiniaceae - F003 - PBS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank()),
    plot_list$`Bacteria - F003 - PBS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank()),
    plot_list$`Host - Acropora - DESS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank()),
    plot_list$`Symbiodiniaceae - Acropora - DESS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank()),
    plot_list$`Bacteria - Acropora - DESS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank()),
    plot_list$`Host - Porites - DESS` + theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank()),
    plot_list$`Symbiodiniaceae - Porites - DESS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank()),
    plot_list$`Bacteria - Porites - DESS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank()),
    plot_list$`Host - Pocillopora - DESS` + theme(strip.text.x = element_blank(), strip.text.y = element_blank()),
    plot_list$`Symbiodiniaceae - Pocillopora - DESS` + theme(strip.text.x = element_blank(), strip.text.y = element_blank()),
    plot_list$`Bacteria - Pocillopora - DESS` + theme(strip.text.x = element_blank()),
  ncol = 3, nrow = 5, common.legend = TRUE, legend = "bottom",
  align = "v",
  widths = c(1, 1, 1),
  heights = c(0.9, rep(.8, 3), 1.1)
)

figure_2

ggsave(figure_2, 
       filename = file.path(output_path, "figure_2.png"), 
       width = 16, height = 22, dpi = 600)
ggsave(figure_2, 
       filename = file.path(output_path, "figure_2.svg"), 
       width = 16, height = 22, dpi = 600)



# =============================================================================
# FIGURE 1: qPCR ANALYSIS
# =============================================================================

Aip_PBS <- read.table(
    file.path(base_path, "qPCR/MetaG_qPCR_Aip_PBS.csv"), 
    header = TRUE, sep = ","
) 

PBSvsDESS <-read.table("/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/qPCR/MetaG_qPCR_PBSvsDESS.csv" , header = T, sep =",") 
PBSvsDESS <- PBSvsDESS %>% 
    filter(Buffer == "DESS")

Corals <-read.table("/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/qPCR/MetaG_qPCR_corals.csv" , header = T, sep =",") 

data <- rbind(Aip_PBS, PBSvsDESS, Corals)
data$Ct <- as.numeric(data$Ct)

data_summary <- dcast(
    formula = Sample_name ~ Gene, 
    fun.aggregate = mean, 
    data = data
) %>%
    filter(complete.cases(.))

data_summary$dCT_16S <- -1*(data_summary$`16S` - data_summary$`b-actin`)
data_summary <- tidyr::separate(
    data_summary, 
    Sample_name, 
    into = c("Strain", "replicate", "Buffer", "Treatment"), 
    sep = "_", 
    remove = FALSE, 
    extra = "merge"
)

# Expand treatment abbreviations
data_summary$Treatment <- data_summary$Treatment %>%
    gsub("\\bMB_", "Microbiome Kit_1", .) %>%
    gsub("\\bBT_", "Blood & Tissue Kit_", .) %>%
    gsub("\\bSP_", "Sonication + PMA + BT_", .) %>%
    gsub("\\bSB_", "Sonication + Benzonase + BT_", .) %>%
    gsub("\\bP_", "PMA + BT_", .) %>%
    gsub("\\bB_", "Benzonase + BT_", .)

# Clean species names
data_summary$Strain <- clean_species_names(data_summary$Strain)

data_summary <- tidyr::separate(
    data_summary, 
    Treatment, 
    into = c("Treatment", "Run"), 
    sep = "_", 
    remove = FALSE, 
    extra = "merge"
)

# Set factor levels for treatments
data_summary$Treatment <- factor(
    data_summary$Treatment, 
    levels = c("Blood & Tissue Kit", "Microbiome Kit", "Benzonase + BT", 
               "Sonication + Benzonase + BT", "PMA + BT", "Sonication + PMA + BT")
)

# Calculate summary statistics per group
summary_16S <- data_summary %>%
    group_by(Strain, Buffer, Run, Treatment) %>% 
    summarise(
        ratio_mean = mean(dCT_16S), 
        ratio_se = sd(dCT_16S)/sqrt(length(dCT_16S)),
        .groups = "drop"
    )

summary_16S$Strain <- factor(
    summary_16S$Strain, 
    levels = c("H2", "F003", "Acropora", "Pocillopora", "Porites")
)

# Add grouping column
summary_16S["grouping"] <- paste(summary_16S$Strain, summary_16S$Buffer, sep = " - ")
breaks <- c("Blood & Tissue Kit", "Benzonase + BT", "Sonication + Benzonase + BT", "PMA + BT", "Sonication + PMA + BT", "Microbiome Kit")

# Generate qPCR plots

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
                width = .1, 
                position=position_dodge(.9), 
                color="black") +
            facet_wrap( ~ grouping, labeller = as_labeller(custom_labeller, label_parsed))+
            geom_hline(yintercept = 0,  linetype = "dashed", color = "black")+
            scale_fill_manual(values = col_palette_qpcr) +
            scale_x_discrete(
                labels = labels_abbrev,
                breaks = breaks
            ) +
            theme_bw() +
            labs(x= "", y= "- dCT (16S - b-actin)") +
            theme(legend.position = "none", line = element_line(linewidth = 1),
                axis.line = element_line(colour = "black"),
                axis.ticks.length = unit(0.2 , "cm"),
                axis.text = element_text(size = 14, color = "black"),
                axis.title.y = element_text(size = 16, color = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                strip.background = element_blank(), 
                strip.text.x = element_text(color = "black", size = 14, angle = 0, hjust = 0, vjust = 0.5, face = "plain"),
                panel.spacing = unit(3, "lines")) +
                coord_cartesian(ylim=c(-15,6))
    
    qPCR_plot_list[[n]] <- p
    n <- n + 1
}

legend  <- extract_legend(qPCR_plot_list[[5]] + theme(legend.position = "bottom", legend.key.size = unit(1.5, "lines"), legend.text = element_text(size = 12), legend.title = element_text(size = 14), legend.direction = "horizontal")+ guides(fill = guide_legend(nrow = 1)))
spacer <- ggplot() + theme_void()

figure_1 <- ggarrange(
    spacer, ggarrange(
        ggarrange(
            annotate_figure(
                qPCR_plot_list[[5]],
                fig.lab = "A",
                fig.lab.pos = "top.left",
                fig.lab.size = 14,
                fig.lab.face = "bold"),
            annotate_figure(
                qPCR_plot_list[[4]] + labs(y = ""),
                fig.lab = "B",
                fig.lab.pos = "top.left",
                fig.lab.size = 14,
                fig.lab.face = "bold"),
            ncol = 2, nrow = 1, widths = c(2, 1)),
        ggarrange(
            annotate_figure(
                qPCR_plot_list[[3]],
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
            ncol = 2, nrow = 1, widths = c(2, 1)),
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
            ncol = 3, nrow = 1, widths = c(1, 1, 1)),
        legend,
        ncol = 1, nrow = 4, heights = c(1,1,1,0.1)), 
        ncol = 2, nrow = 1, widths = c(0.1,1.6))

    # annotate_figure(figure_1, 
    # bottom = text_grob("*Ct Values from multiple qPCR run have been combined in plotting when required", color = "black",
    #         hjust = 1, x = 1, face = "italic", size = 10))

figure_1_annotated <- annotate_figure(
        figure_1,
        bottom = text_grob("*Ct Values from multiple qPCR run have been combined in plotting when required", color = "black",
            hjust = 1, x = 1, face = "italic", size = 10),
        left = grid::grobTree(
            grid::textGrob("Increased bacteria proportion", x = 4.5, y = 0.75, rot = 90, gp = grid::gpar(fontsize = 14), just = "center"),
            grid::linesGrob(x = c(6,6), y = c(0.6, 0.9), arrow = grid::arrow(type = "open", length = unit(0.5, "cm")), gp = grid::gpar(lwd = 2)),
            grid::textGrob("Increased host proportion", x = 4.5, y = 0.35, rot = 90, gp = grid::gpar(fontsize = 14), just = "center"),
            grid::linesGrob(x = c(6,6), y = c(0.5, 0.2), arrow = grid::arrow(type = "open", length = unit(0.5, "cm")), gp = grid::gpar(lwd = 2))
        )
    )
figure_1_annotated


# =============================================================================
# SAVE FIGURES
# =============================================================================

ggsave(figure_1_annotated, 
       filename = file.path(output_path, "figure_1.png"), 
       width = 18, height = 22, dpi = 600)
ggsave(figure_1_annotated, 
       filename = file.path(output_path, "figure_1.svg"), 
       width = 18, height = 22, dpi = 600)




# Prepare plotting data
CAT_data_family_curated$Strain <- factor(
    CAT_data_family_curated$Strain, 
    levels = c("F003", "H2", "Acropora", "Porites", "Pocillopora")
)

spacer <- ggplot() + theme_void()

# Define labels for treatments (use previously defined labels_abbrev)

# Generate family plots

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

# Calculate the top 20 families per species
top_families <- CAT_data_family %>%
    separate(Sample, into = c("Strain", "Treatment", "Buffer"), sep = "_", remove = FALSE) %>%
    group_by(Strain, Family) %>%
    summarise(total_contigs = sum(Number_of_contigs), .groups = "drop") %>%
    arrange(Strain, desc(total_contigs)) %>%
    group_by(Strain) %>%
    slice_head(n = 20) %>%
    pull(Family) %>%
    unique()

    # Count unique families for Pocillopora ("Poci") samples
    num_families_poci <- CAT_data_family %>%
        filter(grepl("^Poci|^Pocillopora", Sample, ignore.case = TRUE)) %>%
        pull(Family) %>%
        unique() %>%
        length()
    print(paste("Number of unique families for Pocillopora samples:", num_families_poci))


    # CAT_data_family_curated %>%
    #     filter(grepl("^Pocillopora", Strain, ignore.case = TRUE)) %>%
    #     pull(Family) %>%
    #     unique() %>%
    #     length()
    # print(paste("Number of unique families for Pocillopora samples:", num_families_poci))
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

# Rename Treatment using helper function
CAT_data_family_curated$Treatment <- expand_treatment_names(CAT_data_family_curated$Treatment)

# Clean column names using helper functions
CAT_data_family_curated$Strain <- clean_species_names(CAT_data_family_curated$Strain)
CAT_data_family_curated$Buffer <- clean_buffer_names(CAT_data_family_curated$Buffer)

# Sort Family by number of contigs, then put "Other" and "ambigous match" last
family_order <- CAT_data_family_curated %>%
    group_by(Family) %>%
    summarise(total_contigs = sum(total_contigs)) %>%
    arrange(Family) %>%
    pull(Family)

family_order <- sort(family_order)
# Move "Other" and "ambigous match" to the end
family_order <- c(setdiff(family_order, c("Other", "Ambigous match", "No match")), "Other", "Ambigous match", "No match")

# Reorder the factor levels
CAT_data_family_curated$Family <- factor(CAT_data_family_curated$Family, levels = rev(family_order))

# Reorder factor levels for buffer
CAT_data_family_curated$Buffer <- factor(CAT_data_family_curated$Buffer, levels = c("PBS", "DESS"))

# add plotting category called strain_treatment to plot by strain and treatment
CAT_data_family_curated["grouping"] <- paste(CAT_data_family_curated$Strain, CAT_data_family_curated$Buffer, sep = " - ")

# =============================================================================
# FIGURE 3: CAT FAMILY CLASSIFICATION (DOT PLOTS)
# =============================================================================

# Load non-filtered microbial only files
files <- list.files(
    path = file.path(base_path, "CAT/results_download/NR/summary/bacteria_only_contig_class"),
    pattern = "*.summary.txt", 
    full.names = TRUE
)


plot_family_list_dot <- list()
n <- 1
for (i in sort(unique(CAT_data_family_curated$Strain))) {
    print(i)
        data_to_plot <- CAT_data_family_curated %>%
            filter(Strain == i) %>%
            filter(Family != "No match") %>%
            filter(Treatment != "Benzonase + BT") %>%
            filter(Treatment != "Blood & Tissue Kit") %>%
            filter(grouping != "H2 - DESS") %>%
            filter(grouping != "F003 - DESS")

        p <- data_to_plot %>%
            ggplot(aes(x = Treatment, y = Family, color = Treatment)) +
            geom_point(aes(size = log10(total_contigs) + 1), shape = 16, show.legend = TRUE) +
            facet_wrap(~ Strain, ncol = 1) +
            labs(title = "", x = "", y = "") +
            scale_color_manual(values = color_palette_x) +
            theme_minimal(base_size = 16) +
            theme(
            legend.position = "bottom",
            legend.box = "vertical",
            panel.grid.minor = element_blank(),
            strip.text.x = element_text(size = 12, face = ifelse(i %in% c("F003", "H2"), "plain", "italic")),
            axis.text.x = element_text(
                size = 12,
                angle = 90,
                vjust = 0.5
            ),
            axis.text.y = ggtext::element_markdown(size = 12) # for ggtext if you want color
            ) + scale_x_discrete(
                labels = labels_abbrev
            )

    plot_family_list_dot[[n]] <- p
    n <- n + 1
}

# legend <- extract_legend(plot_family_list_dot[[1]] + theme(legend.position = "right", legend.key.size = unit(1.5, "lines"), legend.title = element_text(size = 16) ) + guides(fill = guide_legend(nrow = 1), color = guide_legend(override.aes = list(size = 20)))) 

figure_3  <- ggarrange(
                # ggarrange(
                    plot_family_list_dot[[2]] + 
                        theme(legend.position = "bottom", 
                                legend.key.size = unit(1.2, "lines"),
                                legend.title = element_text(size = 18),
                                legend.text = element_text(size = 16)) + 
                                guides(fill = guide_legend(nrow = 1), color = guide_legend(override.aes = list(size = 8))),
                    plot_family_list_dot[[1]],
                    # ncol = 2, nrow = 1, widths = c(1, 1)),
                # ggarrange(
                    plot_family_list_dot[[3]],
                    plot_family_list_dot[[5]],
                    plot_family_list_dot[[4]],
                    # ncol = 3, nrow = 1),
                # spacer, spacer, legend, spacer, spacer,
                ncol = 5, nrow = 1, common.legend = TRUE, legend = "bottom") # heights = c(1,0.1)
figure_3

ggsave(figure_3, 
       filename = file.path(output_path, "figure_3.png"), 
       width = 18, height = 14, dpi = 600)
ggsave(figure_3, 
       filename = file.path(output_path, "figure_3.svg"), 
       width = 18, height = 14, dpi = 600)
