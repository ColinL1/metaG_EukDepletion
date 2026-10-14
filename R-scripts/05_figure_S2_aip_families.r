#! /usr/bin/env Rscript
# ==============================================================================
# Figure S2 — bacterial families detected in Aiptasia (H2 / F003) per protocol
# and buffer; dot size = relative coverage. Families outside each strain's
# top TOP_N are collapsed into "Other".
# Presence calls come from prepare_cat_families() (00_setup.r).
# Output: ../plots/supplementary/Figure_S2.*
# ==============================================================================
source("00_setup.r")
p_load(ggtext)

# ---- PARAMS (defaults = paper version) ---------------------------------------
DEFAULTS <- list(MODE = "persample", RAREFY = TRUE, PRESENCE_FRAC = 0.9, MIN_BREADTH = NULL,
                 COV_METHOD = "trimmed_mean", TOP_N = 35)
PARAMS <- DEFAULTS

# ---- data ---------------------------------------------------------------------
cat_data <- prepare_cat_families(mode = PARAMS$MODE, rarefy = PARAMS$RAREFY, presence_frac = PARAMS$PRESENCE_FRAC,
                                 min_breadth = PARAMS$MIN_BREADTH, cov_method = PARAMS$COV_METHOD, top_n = PARAMS$TOP_N)

labels_abbrev <- setNames(names(protocol_full), protocol_full)
dtp_supp <- cat_data$plot_dat %>%
    filter(Strain %in% c("H2", "F003")) %>%
    mutate(Treatment = factor(Treatment, levels = unname(protocol_full)))

# ---- one panel per strain, faceted by protocol ----------------------------------
plot_aip_dot <- list()
for (i in c("H2", "F003")) {
    d_i <- dtp_supp %>% filter(Strain == i)
    if (nrow(d_i) == 0) next
    plot_aip_dot[[i]] <- ggplot(d_i, aes(x = Buffer, y = Family, color = Treatment)) +
        geom_point(aes(size = log10(rel * 1e6 + 1)), shape = 16, show.legend = TRUE) +
        facet_wrap(~ Treatment, nrow = 1, labeller = as_labeller(labels_abbrev)) +
        labs(title = "", x = "", y = "") +
        scale_color_manual(name = "Protocols", values = col_protocol_full, labels = labels_abbrev) +
        theme_minimal(base_size = 14) +
        theme(legend.position = "bottom", legend.box = "vertical",
              panel.grid.minor = element_blank(),
              strip.text.x = element_text(size = 12, face = "plain"),
              axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5),
              axis.text.y = element_markdown(size = 12))
}

figure <- ggarrange(
    annotate_figure(
        plot_aip_dot[["H2"]] + labs(title = "H2") +
            theme(plot.title = element_text(size = 16), legend.position = "none"),
        fig.lab = "A", fig.lab.pos = "top.left", fig.lab.size = 15, fig.lab.face = "bold"),
    annotate_figure(
        plot_aip_dot[["F003"]] + labs(title = "F003") +
            guides(size  = guide_legend(title.position = "top", nrow = 1),
                   color = guide_legend(title.position = "top", nrow = 1, override.aes = list(size = 8))) +
            theme(plot.title = element_text(size = 16),
                  legend.position = "bottom", legend.box = "horizontal",
                  legend.key.size = unit(1.2, "lines"),
                  legend.title = element_text(size = 10, hjust = 0.5),
                  legend.text = element_text(size = 10)),
        fig.lab = "B", fig.lab.pos = "top.left", fig.lab.size = 15, fig.lab.face = "bold"),
    ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom"
)

out <- output_target("Figure_S2", PARAMS, DEFAULTS, kind = "supp")
save_figure(figure, out$name, out$kind, width = 10, height = 17)
