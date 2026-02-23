
library(ggplot2)
library(dplyr)
library(ggpubr)
library(stringr)

library(ggh4x)

test_data <- data
test_data$mapping <- factor(test_data$mapping, levels = rev(c("total", "Other", "Bacteria", "Symbiodiniaceae", "Host")))
test_data$type <- ifelse(test_data$species %in% c("H2", "F003"), "aiptasia", "corals")
test_data$type <- factor(test_data$type, levels = c("aiptasia", "corals"))
col_palette_2 <- c("Host" = "#D64933", "Symbiodiniaceae" = "#2A7F62", "Bacteria" = "#3498DB", "Other" = "#414288")# #EAC435 for host? 

test_data <- test_data[,c("sample", "grouping", "type", "species", "replicate", "extraction", "buffer", "mapping", "num_seqs")]

# # Step 1: Pivot wider so each mapping becomes a column
# df_wide <- test_data %>%
#   pivot_wider(names_from = mapping, values_from = num_seqs) # values_fill = 0

# # Step 2: Compute "arch"
# df_with_arch <- df_wide %>%
#   mutate(Archaea = total - (Host + Symbiodiniaceae + Bacteria + Other))

# # Step 3: Pivot back to long format (adding arch as new mapping)
# df_final <- df_with_arch %>%
#   pivot_longer(cols = c(total, Host, Symbiodiniaceae, Bacteria, Archaea, Other),
#                names_to = "mapping",
#                values_to = "num_seqs") %>%
#   arrange(sample, mapping)  # Optional: sort by sample then mapping


# test_data <- df_final %>%
test_data <- test_data %>%
    filter(mapping != "total")


test_data <- test_data %>%
    group_by(sample) %>%
    mutate(num_seq_percent = num_seqs / sum(num_seqs, na.rm = TRUE) * 100) %>%
    ungroup()

for (i in c("Host", "Symbiodiniaceae", "Bacteria", "Archaea", "Other")) {
    # test_data$mapping[test_data$mapping == i] <- i
    for (j in c("aiptasia", "corals")) {
        plot_data <- test_data %>%
            filter(mapping != "total") %>%
            filter(mapping == i) %>%
            filter(type == j)
    print(paste0("Running stats for ", i, " in ", j))
    print(compare_means(num_seq_percent ~ extraction + grouping, plot_data, method = "t.test", paired = FALSE))
    }
}


test_data$mapping <- factor(test_data$mapping, levels = rev(c("total", "Other", "Archaea", "Bacteria", "Symbiodiniaceae", "Host")))
test_data$type <- ifelse(test_data$species %in% c("H2", "F003"), "aiptasia", "corals")
test_data$type <- factor(test_data$type, levels = c("aiptasia", "corals"))
test_data$extraction <- factor(test_data$extraction, levels = c("Blood & Tissue", "Benzonase", "Microbiome",  "Microbiome & bead-beating", "Microbiome & spinning"))
col_palette_2 <- c("Host" = "#D64933", "Symbiodiniaceae" = "#2A7F62", "Bacteria" = "#3498DB", "Archaea" = "#414288", "Other" = "#414288")# #EAC435 for host? 
my_comparisons <- list(c("Blood & Tissue", "Benzonase"), c("Blood & Tissue", "Microbiome"), c("Blood & Tissue", "Microbiome & bead-beating"), c("Blood & Tissue", "Microbiome & spinning"))

col_palette_2 <- c("Host" = "#D64933", "Symbiodiniaceae" = "#2A7F62", "Bacteria" = "#3498DB", "Other" = "#414288")# #EAC435 for host?
col_palette_2 <- c("Blood & Tissue" = "#FFB94F", "Benzonase" = "#D97474", "Microbiome" = "#9AE0EC", "Microbiome & bead-beating" ="#87ACE9", "Microbiome & spinning" = "#9C80BA")# #EAC435 for host?

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
      # my_comparisons <- list(c("Blood & Tissue", "Benzonase"), c("Blood & Tissue", "Microbiome"), c("Blood & Tissue", "Microbiome & bead-beating"), c("Blood & Tissue", "Microbiome & spinning"))

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
      filter(group1 == "Blood & Tissue") %>% # only keep comparisons starting with Blood & Tissue
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
              axis.text.x = element_text(color = "black", size = 16, angle = 90),
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

ggpubr::ggarrange(
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
  ncol = 3, nrow = 5, common.legend = TRUE, legend = "none",
  align = "v",
  widths = c(1, 1, 1),
  heights = c(0.9, rep(.8, 3), 1.1)
)

# if switch is used in facet_grid2, strip.text.y needs to be removed from all but first column
# ggpubr::ggarrange(
#     plot_list$`Host - H2 - PBS` + theme(axis.text.x = element_blank()),
#     plot_list$`Symbiodiniaceae - H2 - PBS`+ theme(axis.text.x = element_blank(), strip.text.y = element_blank()),
#     plot_list$`Bacteria - H2 - PBS` + theme(axis.text.x = element_blank(), strip.text.y = element_blank()),
#     plot_list$`Host - F003 - PBS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank()),
#     plot_list$`Symbiodiniaceae - F003 - PBS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank()),
#     plot_list$`Bacteria - F003 - PBS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank()),
#     plot_list$`Host - Acropora - DESS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank()),
#     plot_list$`Symbiodiniaceae - Acropora - DESS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank()),
#     plot_list$`Bacteria - Acropora - DESS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank()),
#     plot_list$`Host - Porites - DESS` + theme(strip.text.x = element_blank(), axis.text.x = element_blank()),
#     plot_list$`Symbiodiniaceae - Porites - DESS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank()),
#     plot_list$`Bacteria - Porites - DESS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank()),
#     plot_list$`Host - Pocillopora - DESS` + theme(strip.text.x = element_blank()),
#     plot_list$`Symbiodiniaceae - Pocillopora - DESS` + theme(strip.text.x = element_blank(), strip.text.y = element_blank()),
#     plot_list$`Bacteria - Pocillopora - DESS` + theme(strip.text.x = element_blank(), strip.text.y = element_blank()),
#   ncol = 3, nrow = 5, common.legend = TRUE, legend = "none",
#   align = "v",
#   widths = c(1, 1, 1),
#   heights = c(0.9, rep(.8, 3), 1.1)
# )



write.table(stat_results_save, file = "results/mapping/stats_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# ggpubr::ggarrange(
#     plot_list$`Host - H2 - PBS` + theme(axis.text.x = element_blank(), strip.text.y = element_blank()) + ylim(c(0, 80)),
#     plot_list$`Symbiodiniaceae - H2 - PBS`+ theme(axis.text.x = element_blank(), strip.text.y = element_blank())+ ylim(c(0, 80)),
#     plot_list$`Bacteria - H2 - PBS` + theme(axis.text.x = element_blank())+ ylim(c(0, 60)),
#     plot_list$`Host - F003 - PBS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank())+ ylim(c(0, 80)),
#     plot_list$`Symbiodiniaceae - F003 - PBS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank())+ ylim(c(0, 80)),
#     plot_list$`Bacteria - F003 - PBS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank())+ ylim(c(0, 60)),
#     plot_list$`Host - Acropora - DESS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank())+ ylim(c(0, 80)),
#     plot_list$`Symbiodiniaceae - Acropora - DESS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank())+ ylim(c(0, 80)),
#     plot_list$`Bacteria - Acropora - DESS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank())+ ylim(c(0, 60)),
#     plot_list$`Host - Porites - DESS` + theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank())+ ylim(c(0, 80)),
#     plot_list$`Symbiodiniaceae - Porites - DESS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank(), strip.text.y = element_blank())+ ylim(c(0, 80)),
#     plot_list$`Bacteria - Porites - DESS`+ theme(strip.text.x = element_blank(), axis.text.x = element_blank())+ ylim(c(0, 60)),
#     plot_list$`Host - Pocillopora - DESS` + theme(strip.text.x = element_blank(), strip.text.y = element_blank())+ ylim(c(0, 80)),
#     plot_list$`Symbiodiniaceae - Pocillopora - DESS` + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + ylim(c(0, 80)),
#     plot_list$`Bacteria - Pocillopora - DESS` + theme(strip.text.x = element_blank()) + ylim(c(0, 60)),
#   ncol = 3, nrow = 5, common.legend = TRUE, legend = "none",
#   align = "v",
#   widths = c(1, 1, 1),
#   heights = c(0.9, rep(.8, 3), 1.1) 


# plot_list$corals
# plot_list$aiptasia
# ggarrange(plot_list$`Host - H2 - PBS`,
#             plot_list$`Symbiodiniaceae - H2 - PBS`,
#             plot_list$`Bacteria - H2 - PBS`,
#             plot_list$`Host - F003 - PBS`,
#             plot_list$`Symbiodiniaceae - F003 - PBS`,
#             plot_list$`Bacteria - F003 - PBS`,
#           ncol = 3, nrow = 2, common.legend = TRUE, legend = "none")

# ggarrange(plot_list$`Host - Acropora - DESS`,
#             plot_list$`Symbiodiniaceae - Acropora - DESS`,
#             plot_list$`Bacteria - Acropora - DESS`,
#             plot_list$`Host - Porites - DESS`,
#             plot_list$`Symbiodiniaceae - Porites - DESS`,
#             plot_list$`Bacteria - Porites - DESS`,
#             plot_list$`Host - Pocillopora - DESS`,
#             plot_list$`Symbiodiniaceae - Pocillopora - DESS`,
#             plot_list$`Bacteria - Pocillopora - DESS`,
#           ncol = 3, nrow = 3, common.legend = TRUE, legend = "none")

# ggarrange(plot_list$aiptasia + labs(title = "", x = "", y = ""),
#           plot_list$corals +
#               labs(title = "", x = "", y = "") +
#               theme(strip.text.x = element_blank()),
#           heights = c(0.4, 0.7),
#           vjust = -10,
#           ncol = 1, nrow = 2, common.legend = TRUE, legend = "none")

### start of solution to stats missing in faceting. 
library(ggpubr)
library(dplyr)
# Only keep mapping types of interest
plot_data_filtered <- plot_data %>%
  filter(mapping %in% c("Host", "Bacteria", "Symbiodiniaceae")) %>%
  filter(extraction %in% c("Blood & Tissue", "Microbiome"))

# Calculate comparisons per facet
pvals <- compare_means(
  num_seq_percent ~ extraction,
  data = plot_data,
  group.by = c("mapping", "grouping"),
  method = "t.test",
  ref.group = "Blood & Tissue"
)
