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

names(test_data)
sum(check_quick$num_seqs) - 39791238

test_data <- test_data %>%
    filter(mapping != "total")

test_data <- test_data[,c("sample", "grouping", "type", "species", "replicate", "extraction", "buffer", "mapping", "num_seqs")]

test_data <- test_data %>%
    group_by(sample) %>%
    mutate(num_seq_percent = num_seqs / sum(num_seqs, na.rm = TRUE) * 100) %>%
    ungroup()

for (i in c("Host", "Symbiodiniaceae", "Bacteria", "Other")) {
    # test_data$mapping[test_data$mapping == i] <- i
    for (j in c("aiptasia", "corals")) {
        plot_data <- test_data %>%
            filter(mapping != "total") %>%
            filter(mapping == i) %>%
            filter(type == j)
    print(compare_means(num_seq_percent ~ extraction , plot_data, method = "wilcox.test", paired = FALSE))
    }
}

my_comparisons <- list(c("Blood & Tissue", "Benzonase"), c("Blood & Tissue", "Microbiome"), c("Blood & Tissue", "Microbiome & bead-beating"), c("Blood & Tissue", "Microbiome & spinning"))

#   stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#   stat_compare_means(label.y = 50)     # Add global p-value
test_data$type <- ifelse(test_data$species %in% c("H2", "F003"), "aiptasia", "corals")
test_data$type <- factor(test_data$type, levels = c("aiptasia", "corals"))
col_palette_2 <- c("Host" = "#D64933", "Symbiodiniaceae" = "#2A7F62", "Bacteria" = "#3498DB", "Other" = "#414288")# #EAC435 for host? 

plot_list <- list()
plot_final_list <- list()
for (j in c("aiptasia", "corals" )) {
    for (i in c("Host", "Symbiodiniaceae", "Bacteria")) {
    # test_data$mapping[test_data$mapping == i] <- i

        plot_data <- test_data %>%
            filter(mapping != "total") %>%
            filter(mapping != "Other") %>%
            filter(mapping == i) %>%
            filter(type == j)

        # if j == "aiptasia" then filter out dess buffer
        if (j == "aiptasia") {
            plot_data <- plot_data %>%
                filter(buffer != "DESS")
        }
        # # Calculate percentage of num_seqs within each species+buffer group
        # plot_data <- plot_data %>%
        #     group_by(species, buffer) %>%
        #     mutate(num_seqs_pct = num_seqs / sum(num_seqs, na.rm = TRUE) * 100) %>%
        #     ungroup()

      plot <- plot_data %>%
        ggplot(aes(x = extraction, y = num_seq_percent, fill = mapping)) +
        geom_boxplot() +
        facet_wrap(
          ~ grouping + mapping,
          nrow = ifelse(j == "corals", 3, 1),
          ncol = ifelse(j == "corals", 3, 6),
          scales = "free"
        ) +
        labs(
          title = "Percentage of Sequences by Extraction Method and Mapping",
          x = "Extraction Method",
          y = "Percentage of Sequences"
        ) +
        scale_fill_manual(values = col_palette_2) +
        theme_minimal() +
        scale_x_discrete(labels = ~ str_wrap(.x, width = 10)) +
        stat_compare_means(
          method = "anova",
          label.y = max(plot_data$num_seq_percent, na.rm = TRUE) * 1.2
        ) + # Add global p-values
        # stat_compare_means(
        #   label = "p.signif",
        #   method = "t.test",
        #   bracket.size = 0.5,
        #   ref.group = "Blood & Tissue"
        # )  # Add pairwise comparisons p-value
        stat_compare_means(
          method = "t.test",
          comparisons = my_comparisons,
          label = "p.signif",
          hide.ns = TRUE
        )
        plot
            # stat_compare_means(label.y = -10 ) # max(plot_data$num_seqs_pct, na.rm = TRUE) * 1.2
        plot_list[[i]] <- plot
        }
        plot_list
        plot_f <- ggarrange(plot_list[[1]],
                            plot_list[[2]],
                            plot_list[[3]],
                            ncol = 1, nrow = 1)

#     # ggsave(plot_f, filename = paste0(out_path_test, "test_boxplot-stats_", i, ".pdf"), width =  20, height = 11) # A4 dimensions in inches
#     plot_f
      plot_final_list[[j]] <- plot_f
}
plot_final_list$Host #best?
plot_final_list$Symbiodiniaceae #best?

names(plot_final_list)

my_comparisons_2 <- list(
    c("Blood & Tissue", "Benzonase"),
    c("Blood & Tissue", "Microbiome"),
    c("Blood & Tissue", "Microbiome & bead-beating"),
    c("Blood & Tissue", "Microbiome & spinning"),
    c("Blood & Tissue", "Benzonase"),
    c("Blood & Tissue", "Microbiome"),
    c("Blood & Tissue", "Microbiome & bead-beating"),
    c("Blood & Tissue", "Microbiome & spinning"),
    c("Blood & Tissue", "Benzonase"),
    c("Blood & Tissue", "Microbiome"),
    c("Blood & Tissue", "Microbiome & bead-beating"),
    c("Blood & Tissue", "Microbiome & spinning"),
    c("Blood & Tissue", "Benzonase"),
    c("Blood & Tissue", "Microbiome"),
    c("Blood & Tissue", "Microbiome & bead-beating"),
    c("Blood & Tissue", "Microbiome & spinning")
)

### end boxplot code ###


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

        # # Calculate percentage of num_seqs within each species+buffer group
        # plot_data <- plot_data %>%
        #     group_by(species, buffer, mapping, extraction ) %>%
        #     mutate(num_seqs_pct = num_seqs / sum(num_seqs, na.rm = TRUE) * 100) %>%
        #     ungroup()

    plot <- plot_data %>%
        ggplot(aes(x = extraction, y = num_seq_percent, fill = mapping)) +
        geom_boxplot() +
        facet_wrap(~ species + buffer, nrow = 1) +
        labs(title = "",
            x = "Extraction Method",
            y = "Number of Sequences") +
        scale_fill_manual(values = col_palette_2) +
        theme_minimal() +
        scale_x_discrete(labels = ~ str_wrap(.x, width = 10)) +
        # stat_compare_means(method = "anova", label.y = max(plot_data$num_seqs_pct, na.rm = TRUE) * 1.2)+      # Add global p-values
        stat_compare_means(label = "p.signif", method = "t.test",
                    bracket.size = 0.5, ref.group = "Blood & Tissue")  # Add pairwise comparisons p-value
        # stat_compare_means(method = "t.test", comparisons = my_comparisons)+ 
        # stat_compare_means(label.y = -10 ) # max(plot_data$num_seqs_pct, na.rm = TRUE) * 1.2
    plot_list[[j]] <- plot
}

        plot_f <- ggarrange(plot_list[[1]],
                            plot_list[[2]],
                            ncol = 1, nrow = 2)

ggsave(plot_f, filename = paste0(out_path_test, "test_boxplot-stats_", "sym-host", ".pdf"), width =  22, height = 11) # A4 dimensions in inches


# Simple barplot: mean num_seqs per extraction method
library(ggplot2)

unique(test_data$mapping)
test_data$mapping <- factor(test_data$mapping, levels = (c("total", "Other", "Bacteria", "Symbiodiniaceae", "Host")))

    test_data %>%
    filter(sample == "Ac_01_MBS_DESS")

    plot_start_data <- test_data %>%
        filter(mapping != "total") %>%
            group_by(sample) %>%
            reframe(
                mapping = mapping,
                extraction = extraction,
                buffer = buffer,
                species = species,
                total = sum(num_seqs),
                num_seqs = num_seqs)

plot_start_data <- plot_start_data %>%
        group_by(species, extraction, buffer, mapping) %>%
            summarise(total = sum(total), num_seqs = sum(num_seqs), percent = sum(num_seqs)/sum(total) * 100 ) %>%
            ungroup() #%>%
            # mutate(mapping = "total") -> plot_start_data_total
plot_start_data$grouping <- paste(plot_start_data$species, plot_start_data$buffer, sep = " - ")
plot_start_data$grouping <- factor(plot_start_data$grouping, levels = c("H2 - PBS", "H2 - DESS", "F003 - PBS", "F003 - DESS", "Acropora - DESS", "Porites - DESS", "Pocillopora - DESS"))

theme <- theme_minimal(base_size = 14) +
        theme(
            strip.text.y = element_blank(),
            strip.text.x = element_text(size = 12),
            panel.grid.minor = element_blank(),  # Hide minor grid lines
            panel.grid.major = element_blank(),  # Hide minor grid lines
            legend.key.size = unit(2, "lines"),  # Increase legend key size
            panel.spacing.y = unit(1.2, "lines"),  # Increase spacing between vertical panels
            panel.spacing.x = unit(1.2, "lines"),  # Increase spacing between vertical panels
            panel.border = element_blank(),        # Remove full box border
            axis.line = element_line(color = "black") # Only x and y axis lines
        )


        # %>%
        # summarise(average_num_seq = mean(num_seqs))

shift_legend2(
    plot_start_data %>%
        mutate(percent = ifelse(mapping == "Host", -1 * percent, percent)) %>%
        mutate(percent = ifelse(mapping == "Symbiodiniaceae", -1 * percent, percent)) %>%
        # filter(mapping != "total") %>%
        mutate(mapping = factor(mapping, levels = c("Host","Symbiodiniaceae", "Bacteria", "Other"))) %>%
        ggplot(aes(x = extraction, y = percent, fill = mapping)) +
            geom_bar(stat = "identity") +
            facet_wrap(~ grouping, ncol = 4, scale = "free_x",labeller = as_labeller(custom_labeller, label_parsed)) +
            # facet_grid(extraction_groups ~ grouping, scales = "free_y", space = "free", labeller = as_labeller(custom_labeller, label_parsed)) +  # Adjust facet by extraction_groups and grouping
            labs(title = "Percentage of Sequences per Extraction Method",
                x = "Extraction Method",
                y = "Percentage of Sequences",
                fill = "Mapping") +  # Capital letter for legend title
            theme_minimal() +
            theme +
            coord_flip() +
            scale_x_discrete(labels = ~ stringr::str_wrap(.x, width = 10)) +
            scale_fill_manual(values = col_palette_2, breaks = c("Host", "Symbiodiniaceae", "Bacteria", "Other")) +
            theme(legend.position = "right") +
            scale_y_continuous(
                labels = function(x) paste0(ifelse(x < 0, abs(x), x), "%"),
                limits = c(-100, 100),
                breaks = seq(-100, 100, by = 40)
            ) + 
            geom_hline(yintercept = 0, color = "black", linetype = "dashed") 
    )
