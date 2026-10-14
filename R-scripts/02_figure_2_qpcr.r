#! /usr/bin/env Rscript
# ==============================================================================
# Figure 2 — qPCR: bacterial 16S relative to host b-actin (-dCt) per protocol.
# Output: ../plots/{png,svg,pdf}/Figure_2.*
# ==============================================================================
source("00_setup.r")
p_load(reshape2)

# ---- PARAMS (defaults = paper version) ---------------------------------------
DEFAULTS <- list(ylim = c(-15, 6))
PARAMS   <- DEFAULTS

# ---- data ---------------------------------------------------------------------
qpcr_files <- file.path(paths$qpcr, c("MetaG_qPCR_Aip_PBS.csv", "MetaG_qPCR_PBSvsDESS.csv", "MetaG_qPCR_corals.csv"))
check_inputs(qpcr_files)

Aip_PBS   <- read.table(qpcr_files[1], header = TRUE, sep = ",")
PBSvsDESS <- read.table(qpcr_files[2], header = TRUE, sep = ",") %>% filter(Buffer == "DESS")
Corals    <- read.table(qpcr_files[3], header = TRUE, sep = ",")

data <- rbind(Aip_PBS, PBSvsDESS, Corals)
data$Ct <- as.numeric(data$Ct)

data_summary <- dcast(formula = Sample_name ~ Gene, fun.aggregate = mean, data = data)
data_summary <- data_summary[complete.cases(data_summary), ]
data_summary$dCT_16S <- -1 * (data_summary$`16S` - data_summary$`b-actin`)
data_summary <- tidyr::separate(data_summary, Sample_name, into = c("Strain", "replicate", "Buffer", "Treatment"),
                                sep = "_", remove = FALSE, extra = "merge")

# protocol codes -> full names (the "_<run>" suffix is split off below)
data_summary$Treatment <- gsub("\\bMB_", "Microbiome Kit_1", data_summary$Treatment)
data_summary$Treatment <- gsub("\\bBT_", "Blood & Tissue Kit_", data_summary$Treatment)
data_summary$Treatment <- gsub("\\bSP_", "Sonication + PMA + BT_", data_summary$Treatment)
data_summary$Treatment <- gsub("\\bSB_", "Sonication + Benzonase + BT_", data_summary$Treatment)
data_summary$Treatment <- gsub("\\bP_", "PMA + BT_", data_summary$Treatment)
data_summary$Treatment <- gsub("\\bB_", "Benzonase + BT_", data_summary$Treatment)
data_summary$Strain <- canon_species(data_summary$Strain)
data_summary <- tidyr::separate(data_summary, Treatment, into = c("Treatment", "Run"), sep = "_", remove = FALSE, extra = "merge")

# mean and standard error per group
summary_16S <- data_summary %>%
    group_by(Strain, Buffer, Run, Treatment) %>%
    summarise(ratio_mean = mean(dCT_16S),
              ratio_se   = sd(dCT_16S) / sqrt(length(dCT_16S)), .groups = "drop_last")
summary_16S$Strain <- factor(summary_16S$Strain, levels = c("H2", "F003", "Acropora", "Pocillopora", "Porites"))
summary_16S$grouping <- paste(summary_16S$Strain, summary_16S$Buffer, sep = " - ")

breaks <- c("Blood & Tissue Kit", "Benzonase + BT", "Sonication + Benzonase + BT", "PMA + BT", "Sonication + PMA + BT", "Microbiome Kit")
summary_16S$Treatment <- factor(summary_16S$Treatment, levels = breaks)
labels_abbrev <- c("Blood & Tissue Kit" = "BT", "Benzonase + BT" = "BenzoBT", "Sonication + Benzonase + BT" = "SonBenzoBT",
                   "PMA + BT" = "PMABT", "Sonication + PMA + BT" = "SonPMABT", "Microbiome Kit" = "MB")

# ---- one panel per species x buffer -------------------------------------------
qPCR_plot_list <- list()
for (i in unique(summary_16S$grouping)) {
    qPCR_plot_list[[length(qPCR_plot_list) + 1]] <- summary_16S %>%
        filter(grouping == i) %>%
        ggplot(aes(x = Treatment, y = ratio_mean, fill = Treatment)) +
        geom_col(width = 0.5) +
        geom_errorbar(aes(ymin = ratio_mean - ratio_se, ymax = ratio_mean + ratio_se),
                      width = .1, position = position_dodge(.9), color = "black") +
        facet_wrap(~ grouping, labeller = as_labeller(function(x) custom_labeller(x, aip_raw = TRUE), label_parsed)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        scale_fill_manual(name = "Protocols", values = col_qpcr, breaks = breaks, labels = labels_abbrev) +
        scale_x_discrete(labels = labels_abbrev, breaks = breaks) +
        theme_bw(base_size = 20) +
        labs(x = "", y = expression(paste(Delta, "CT"))) +
        theme(legend.position = "none", line = element_line(linewidth = 1),
              axis.line = element_line(colour = "black"),
              axis.ticks.length = unit(0.2, "cm"),
              axis.text = element_text(size = 15, color = "black"),
              axis.title.y = element_text(size = 18, color = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_text(color = "black", size = 16, angle = 0, hjust = 0, vjust = 0.5, face = "plain"),
              panel.spacing = unit(3, "lines")) +
        coord_cartesian(ylim = PARAMS$ylim)
}

# ---- assemble -------------------------------------------------------------------
tag <- function(p, lab) annotate_figure(p, fig.lab = lab, fig.lab.pos = "top.left", fig.lab.size = 16, fig.lab.face = "bold")

legend <- extract_legend(qPCR_plot_list[[5]] +
    theme(legend.position = "bottom", legend.key.size = unit(1.2, "lines"), legend.text = element_text(size = 14),
          legend.title = element_text(size = 16), legend.direction = "horizontal") +
    guides(fill = guide_legend(nrow = 2)))
spacer <- ggplot() + theme_void()

figure <- ggarrange(
    spacer,
    ggarrange(
        ggarrange(tag(qPCR_plot_list[[5]], "A"), tag(qPCR_plot_list[[4]] + labs(y = ""), "B"),
                  ncol = 2, nrow = 1, widths = c(2, 1)),
        ggarrange(tag(qPCR_plot_list[[3]], "C"), tag(qPCR_plot_list[[2]] + labs(y = ""), "D"),
                  ncol = 2, nrow = 1, widths = c(2, 1)),
        ggarrange(tag(qPCR_plot_list[[1]], "E"), tag(qPCR_plot_list[[6]] + labs(y = ""), "F"),
                  tag(qPCR_plot_list[[7]] + labs(y = ""), "G"),
                  ncol = 3, nrow = 1, widths = c(1, 1, 1)),
        legend,
        ncol = 1, nrow = 4, heights = c(1, 1, 1, 0.2)),
    ncol = 2, nrow = 1, widths = c(0.1, 1.6))

# direction arrows on the left margin
figure <- annotate_figure(
    figure,
    left = grid::grobTree(
        grid::textGrob("Increased bacterial proportion", x = 4.5, y = 0.75, rot = 90, gp = grid::gpar(fontsize = 15), just = "center"),
        grid::linesGrob(x = c(6, 6), y = c(0.6, 0.9), arrow = grid::arrow(type = "open", length = unit(0.5, "cm")), gp = grid::gpar(lwd = 2)),
        grid::textGrob("Increased host proportion", x = 4.5, y = 0.35, rot = 90, gp = grid::gpar(fontsize = 15), just = "center"),
        grid::linesGrob(x = c(6, 6), y = c(0.5, 0.2), arrow = grid::arrow(type = "open", length = unit(0.5, "cm")), gp = grid::gpar(lwd = 2))
    )
)

out <- output_target("Figure_2", PARAMS, DEFAULTS, kind = "main")
save_figure(figure, out$name, out$kind, width = 18, height = 22, scale = 0.7, limitsize = FALSE,
            bg = c(png = "white", svg = NA, pdf = "white"))
