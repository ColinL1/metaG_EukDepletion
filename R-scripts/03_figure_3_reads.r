#! /usr/bin/env Rscript
# ==============================================================================
# Figure 3 — share of reads assigned to host / Symbiodiniaceae / bacteria per
# extraction protocol, one row per species x buffer, with the significance
# stars from 01_stats_anova.r.
# Inputs : read stats (00_setup.r) + ../results/anova_posthoc_contrasts.csv
# Output : ../plots/{png,svg,pdf}/Figure_3.*
# ==============================================================================
source("00_setup.r")
p_load(ggh4x, patchwork)

# ---- PARAMS (defaults = paper version) ---------------------------------------
DEFAULTS <- list(base_size = 24, axis_size = 18, strip_size = 18, legend_size = 22)
PARAMS   <- DEFAULTS

# y-axis upper limit and break step per mapping (and grouping, where needed);
# NA = keep the default axis without stars
y_axis <- function(mapping, grouping) {
    switch(mapping,
        Host            = c(82, 20),
        Symbiodiniaceae = if (grouping == "Porites - DESS") c(10, 2) else c(97, 20),
        Bacteria        = switch(grouping,
                                 "H2 - PBS" = , "F003 - PBS" = c(42, 10),
                                 "Acropora - DESS"    = c(62, 20),
                                 "Porites - DESS"     = c(64, 20),
                                 "Pocillopora - DESS" = c(12, 2),
                                 c(NA, NA)))
}

# ---- data ---------------------------------------------------------------------
if (!file.exists(paths$anova_csv)) stop("Run 01_stats_anova.r first: ", paths$anova_csv, " not found", call. = FALSE)
posthoc_results <- read_csv(paths$anova_csv, show_col_types = FALSE)[, -1]

test_data <- load_read_stats() %>%
    filter(mapping != "total") %>%
    group_by(sample) %>%
    mutate(num_seq_percent = num_seqs / sum(num_seqs, na.rm = TRUE) * 100) %>%
    ungroup() %>%
    mutate(grouping   = factor(grouping, levels = grouping_levels),
           extraction = factor(extraction, levels = extraction_levels))

# significance stars as labels above the boxes, coloured by the reference protocol
plot_sig <- function(plot, ymax, breaks = 20, stats, palette = col_protocol) {
    if (nrow(stats |> filter(p.signif != "")) == 0) {
        return(plot + scale_y_continuous(labels = function(x) paste0(x, "%"), breaks = seq(0, ymax, breaks), limits = c(0, ymax)))
    }
    stats_sig <- stats %>%
        filter(p.signif != "") %>%
        group_by(extraction) %>%
        arrange(extraction, base) %>%
        mutate(y_pos = ymax + 3 + (row_number() - 1) * 2) %>%
        ungroup()
    plot +
        geom_label(data = stats_sig,
                   aes(x = extraction, y = (ymax) / 10, label = p.signif, fill = NULL, color = base),
                   size = 6, border.color = "#00ffff00",
                   position = position_stack(vjust = 10, reverse = TRUE), show.legend = FALSE) +
        scale_color_manual(values = palette) +
        scale_y_continuous(labels = function(x) paste0(x, "%"), breaks = seq(0, ymax, breaks), limits = c(0, max(stats_sig$y_pos) + 5))
}

my_theme <- theme_minimal(base_size = PARAMS$base_size) +
    theme(
        legend.position = "none",
        line = element_line(linewidth = 0.8),
        axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size = PARAMS$axis_size),
        axis.text.y = element_text(color = "black", size = PARAMS$axis_size),
        strip.text.x = element_text(color = "black", size = PARAMS$strip_size, hjust = 0, vjust = 0.5),
        strip.text.y = element_text(color = "black", size = PARAMS$strip_size, hjust = 0, vjust = 0.5),
        strip.clip = "off",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(1, 1, 1, 1),
        panel.spacing = unit(2, "pt")
    )

# ---- one boxplot per mapping x grouping ---------------------------------------------
plot_list <- list()
for (i in sort(unique(test_data$grouping))) {
    for (j in c("Host", "Symbiodiniaceae", "Bacteria")) {
        plot_data <- test_data %>% filter(mapping == j, grouping == i)

        stats <- separate(posthoc_results, col = contrast, into = c("base", "extraction"), sep = " - ", remove = FALSE) |>
            filter(grouping == i & mapping == j) |>
            mutate(p.signif = case_when(
                p.value <= 0.0001 ~ "****",
                p.value <= 0.001  ~ "***",
                p.value <= 0.01   ~ "**",
                p.value <= 0.05   ~ "*",
                TRUE ~ ""
            ))

        plot <- plot_data %>%
            ggplot(aes(x = extraction, y = num_seq_percent, fill = extraction)) +
            geom_boxplot(linewidth = 0.5) +
            facet_grid2(grouping ~ mapping, scales = "free_y", independent = "y", space = "free_x", axes = "all",
                        labeller = labeller(
                            .rows = as_labeller(function(x) custom_labeller(x, pbs_fresh = TRUE), label_parsed),
                            .cols = function(x) paste(x, "reads"))) +
            labs(title = "", x = "", y = "") +
            scale_fill_manual(values = col_protocol) +
            theme_minimal(base_size = 24) +
            scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
            scale_y_continuous(labels = function(x) paste0(x, "%")) +
            my_theme + theme(axis.title = element_blank())

        ax <- y_axis(j, i)
        if (!is.na(ax[1])) plot <- plot_sig(plot, ax[1], breaks = ax[2], stats = stats)
        plot_list[[paste0(j, " - ", i)]] <- plot
    }
}

# ---- assemble: rows A-E = H2, F003, Acropora, Pocillopora, Porites -------------
fig_row <- function(plot1, plot2, plot3, lab) {
    p1 <- plot1 + labs(tag = lab) + theme(strip.text.y = element_blank(), axis.text.x = element_text(angle = 0))
    p2 <- plot2 + theme(strip.text.y = element_blank(), axis.text.x = element_text(angle = 0))
    p3 <- plot3 + theme(axis.text.x = element_text(angle = 0))
    (p1 + p2 + p3 + plot_layout(ncol = 3)) &
        theme(plot.tag = element_text(size = PARAMS$axis_size, face = "bold"),
              plot.tag.position = c(0.03, 0.9),
              plot.margin = margin(0, 0, 0, 0))
}
row_for <- function(g, lab) {
    fig_row(plot_list[[paste0("Host - ", g)]], plot_list[[paste0("Symbiodiniaceae - ", g)]],
            plot_list[[paste0("Bacteria - ", g)]], lab)
}

legend <- extract_legend(plot_list$`Host - H2 - PBS` +
    theme(legend.position = "bottom", legend.text = element_text(size = PARAMS$legend_size), legend.title = element_blank()))
sig_caption <- wrap_elements(full = gridtext::richtext_grob(
    "BT <span style='color:#FFB94F'>\\*\\*\\*</span>&nbsp;&nbsp;&nbsp;MB <span style='color:#9AE0EC'>\\*\\*\\*</span>",
    gp = grid::gpar(fontsize = PARAMS$legend_size)))
spacer <- ggplot() + theme_void()

figure <- (row_for("H2 - PBS", "A") / row_for("F003 - PBS", "B") / row_for("Acropora - DESS", "C") /
           row_for("Pocillopora - DESS", "D") / row_for("Porites - DESS", "E") /
           (spacer | spacer | wrap_elements(legend) | spacer | spacer | sig_caption | spacer | spacer)) +
    plot_layout(heights = c(1.1, 1, 1, 1, 1, 0.2))

out <- output_target("Figure_3", PARAMS, DEFAULTS, kind = "main")
save_figure(figure, out$name, out$kind, width = 16, height = 24, bg = c(png = "white", svg = NA, pdf = NA))
