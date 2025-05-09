#! /usr/bin/env Rscript
## start from read data of script "/home/colinl/metaG/Git/metaG_EukDepletion/R-scripts/All.plot.r"
library(dplyr)

out_path_test <- "/home/colinl/metaG/Git/metaG_EukDepletion/plots/tests/"

test_data <- data
head(test_data)

# # sort as factor test_data$mapping
# test_data$mapping <- factor(test_data$mapping, levels = c("total", "Other", "Bacteria", "Symbiodiniaceae", "Host"))
# test_data$mapping <- factor(test_data$mapping, levels = c("total", "Other", "Bacteria", "Symbiodiniaceae", "Host"))
test_data$extraction <- factor(test_data$extraction, levels = (c("Blood & Tissue", "Benzonase", "Microbiome",  "Microbiome & bead-beating", "Microbiome & spinning")))

for (i in c("Host", "Symbiodiniaceae", "Bacteria", "Other" )) {
    # test_data$mapping[test_data$mapping == i] <- i
    for (j in c("aiptasia", "corals" )) {
        plot_data <- test_data %>%
            filter(mapping != "total") %>%
            filter(mapping == i) %>%
            filter(type == j)
    print(compare_means(num_seqs ~ extraction , plot_data, method = "wilcox.test", paired = FALSE))
    }
}

my_comparisons <- list(c("Blood & Tissue", "Benzonase"), c("Blood & Tissue", "Microbiome"), c("Blood & Tissue", "Microbiome & bead-beating"), c("Blood & Tissue", "Microbiome & spinning"))

#   stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#   stat_compare_means(label.y = 50)     # Add global p-value

col_palette_2 <- c("Host" = "#D64933", "Symbiodiniaceae" = "#2A7F62", "Bacteria" = "#3498DB", "Other" = "#414288")# #EAC435 for host? 

plot_list <- list()
plot_final_list <- list()
for (i in c("Host", "Symbiodiniaceae", "Bacteria", "Other")) {
    # test_data$mapping[test_data$mapping == i] <- i
    for (j in c("aiptasia", "corals" )) {
        plot_data <- test_data %>%
            filter(mapping != "total") %>%
            filter(mapping == i) %>%
            filter(type == j)

        plot <- plot_data %>%
            ggplot(aes(x = extraction, y = num_seqs, fill = mapping)) +
            geom_boxplot() +
            facet_wrap(mapping ~ species + buffer, nrow = 1) +
            labs(title = "Number of Sequences by Extraction Method and Mapping",
                x = "Extraction Method",
                y = "Number of Sequences") +
            scale_fill_manual(values = col_palette_2) +
            theme_minimal() +
            scale_x_discrete(labels = ~ str_wrap(.x, width = 10)) +
            stat_compare_means(method = "anova", label.y = max(plot_data$num_seqs, na.rm = TRUE) * 1.2)+      # Add global p-values
            stat_compare_means(label = "p.signif", method = "t.test",
                        bracket.size = 0.5, ref.group = "Blood & Tissue")  # Add pairwise comparisons p-value
            # stat_compare_means(method = "t.test", comparisons = my_comparisons)+ 
            # stat_compare_means(label.y = -10 ) # max(plot_data$num_seqs, na.rm = TRUE) * 1.2
        plot_list[[j]] <- plot
    }
        plot_f <- ggarrange(plot_list[[1]],
                            plot_list[[2]],
                            ncol = 1, nrow = 2)

    ggsave(plot_f, filename = paste0(out_path_test, "test_boxplot-stats_", i, ".pdf"), width =  20, height = 11) # A4 dimensions in inches

    plot_final_list[[i]] <- plot_f
}
plot_final_list$Host #best?
plot_final_list$Symbiodiniaceae #best?

# p_load(ggpubr)
# arrange(plot_list, ncol = 2, nrow = 2) +
#     theme(plot.title = element_text(hjust = 0.5))

unique(test_data$mapping)

### test host and symbiont together
for (j in c("aiptasia", "corals")) {
    plot_data <- test_data %>%
        filter(mapping != "total") %>%
        filter(mapping == "Symbiodiniaceae" | mapping == "Host") %>%
        filter(type == j)

    plot <- plot_data %>%
        ggplot(aes(x = extraction, y = num_seqs, fill = mapping)) +
        geom_boxplot() +
        facet_wrap(~ species + buffer, nrow = 1) +
        labs(title = "",
            x = "Extraction Method",
            y = "Number of Sequences") +
        scale_fill_manual(values = col_palette_2) +
        theme_minimal() +
        scale_x_discrete(labels = ~ str_wrap(.x, width = 10)) +
        # stat_compare_means(method = "anova", label.y = max(plot_data$num_seqs, na.rm = TRUE) * 1.2)+      # Add global p-values
        stat_compare_means(label = "p.signif", method = "t.test",
                    bracket.size = 0.5, ref.group = "Blood & Tissue")  # Add pairwise comparisons p-value
        # stat_compare_means(method = "t.test", comparisons = my_comparisons)+ 
        # stat_compare_means(label.y = -10 ) # max(plot_data$num_seqs, na.rm = TRUE) * 1.2
    plot_list[[j]] <- plot
}

        plot_f <- ggarrange(plot_list[[1]],
                            plot_list[[2]],
                            ncol = 1, nrow = 2)

ggsave(plot_f, filename = paste0(out_path_test, "test_boxplot-stats_", "sym-host", ".pdf"), width =  22, height = 11) # A4 dimensions in inches
