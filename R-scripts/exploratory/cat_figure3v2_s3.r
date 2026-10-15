#! /usr/bin/env Rscript
# ==============================================================================
# Exploratory CAT figures (not in the paper), from the same family table as
# Figure 4 / Figure S2:
#   extra_02_family_dots  - family dot plot, MB/BBMB/SMB, Aiptasia (Fresh) vs corals
#   extra_03_family_upset - upset plot of families shared across species x buffer, SMB only
# Run from R-scripts/:  Rscript exploratory/cat_figure3v2_s3.r
# Output: ../plots/extra/extra_02_family_dots.png, ../plots/extra/extra_03_family_upset.png
# ==============================================================================
source("00_setup.r")
p_load(ggtext, ComplexUpset)

# ---- PARAMS -------------------------------------------------------------------------
MODE <- "persample"; RAREFY <- TRUE; PRESENCE_FRAC <- 0.9; MIN_BREADTH <- NULL
COV_METHOD <- "trimmed_mean"; TOP_N <- 35
UPSET_TREATMENTS <- c("Spinning + Microbiome Kit")   # protocols pooled into the upset sets

cat_data <- prepare_cat_families(mode = MODE, rarefy = RAREFY, presence_frac = PRESENCE_FRAC,
                                 min_breadth = MIN_BREADTH, cov_method = COV_METHOD, top_n = TOP_N)
fam <- cat_data$fam
plot_dat <- cat_data$plot_dat
labels_abbrev <- setNames(names(protocol_full), protocol_full)

# ---- extra_02: dot plot, species-split (Aiptasia vs coral) ---------------------------
plot_family_list <- list()
for (v in c("Aip", "Coral")) {
    dtp <- plot_dat %>%
        filter(!Treatment %in% c("Benzonase + BT", "Blood & Tissue Kit"),
               !grouping %in% c("H2 - DESS", "F003 - DESS"))
    if (v == "Aip") {
        dtp <- dtp %>% filter(Strain %in% c("H2", "F003"))
        dtp$grouping <- factor(dtp$grouping, levels = c("H2 - Fresh", "H2 - DESS", "F003 - Fresh", "F003 - DESS"))
    } else {
        dtp <- dtp %>% filter(!Strain %in% c("H2", "F003"))
    }
    if (nrow(dtp) == 0) next

    plot_family_list[[v]] <- ggplot(dtp, aes(x = Treatment, y = Family, color = Treatment)) +
        geom_point(aes(size = log10(rel * 1e6 + 1)), shape = 16) +
        facet_wrap(~ grouping, nrow = 1, scale = "free_x",
                   labeller = as_labeller(function(x) custom_labeller(x, style = "stacked"), label_parsed)) +
        labs(title = "", x = "", y = "") +
        scale_color_manual(name = "Protocols", values = col_protocol_full, labels = labels_abbrev) +
        scale_x_discrete(limits = unname(protocol_full[c("MB", "BBMB", "SMB")]), labels = labels_abbrev) +
        theme_minimal(base_size = 16) +
        theme(legend.position = "bottom", legend.box = "vertical",
              panel.grid.minor = element_blank(), strip.background = element_blank(),
              strip.clip = "off",
              strip.text.x = element_text(size = 12, face = ifelse(v == "Aip", "plain", "italic")),
              axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5),
              axis.text.y = element_markdown(size = 12))
}
family_dots <- ggarrange(plot_family_list[["Aip"]], plot_family_list[["Coral"]],
                         common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 1, widths = c(1.5, 2))
save_figure(family_dots, "extra_02_family_dots", "extra", width = 10, height = 12, formats = "png")

# ---- extra_03: upset of families shared across species x buffer -----------------------
sets <- fam %>%
    filter(Treatment %in% UPSET_TREATMENTS) %>%
    distinct(grouping, Family) %>%
    group_by(grouping) %>%
    summarise(Family = list(as.character(Family)), .groups = "drop") %>%
    { setNames(.$Family, .$grouping) }

all_fams <- unique(unlist(sets))
upset_df <- as.data.frame(sapply(sets, function(x) as.integer(all_fams %in% x)))
rownames(upset_df) <- all_fams

group_info <- tibble(grouping = names(sets)) %>%
    mutate(category = if_else(str_detect(grouping, "^(Acropora|Pocillopora|Porites)"), "Corals", "Aiptasia"),
           label_expr = case_when(
               str_detect(grouping, "^Porites")     ~ "italic('P. lobata') ~ plain('- DESS')",
               str_detect(grouping, "^Acropora")    ~ "~italic('A. hemprichii') ~ plain('- DESS')",
               str_detect(grouping, "^Pocillopora") ~ "~italic('P. verrucosa') ~ plain('- DESS')",
               TRUE ~ paste0("plain('", grouping, "')")))

set_colors   <- setNames(c(Corals = "#F48876", Aiptasia = "#D6AB7C")[group_info$category], gsub("-", "_", group_info$grouping, fixed = TRUE))
label_lookup <- setNames(group_info$label_expr, group_info$grouping)
parse_labeller <- function(x) as.expression(lapply(x, function(z) parse(text = z)[[1]]))
set_order <- intersect(rev(c("H2 - Fresh", "F003 - Fresh", "H2 - DESS", "F003 - DESS",
                             "Acropora - DESS", "Pocillopora - DESS", "Porites - DESS")), names(sets))

family_upset <- ComplexUpset::upset(
    upset_df,
    intersect = set_order,
    encode_sets = FALSE,
    name = "Family sharing",
    sort_sets = FALSE,
    labeller = function(x) parse_labeller(unname(label_lookup[x])),
    set_sizes = upset_set_size(geom = geom_bar(aes(fill = group), width = 0.6)) +
        scale_fill_manual(values = set_colors, guide = "none"),
    themes = upset_modify_themes(list(
        "intersections_matrix" = theme(axis.text.y = element_text(size = 11)),
        "overall_sizes" = theme(axis.text.x = element_text(size = 9)),
        "Intersection size" = theme(axis.title.y = element_text(hjust = 0.5, vjust = -30, margin = margin(r = 0, l = 0)),
                                    panel.grid.major = element_blank(), panel.grid.minor = element_blank()))),
    base_annotations = list(
        "Intersection size" = intersection_size(counts = TRUE, mapping = aes(fill = "bars")) +
            scale_fill_manual(values = c(bars = "#8AAFD6"), guide = "none"))
)
save_figure(family_upset, "extra_03_family_upset", "extra", width = 8, height = 5, formats = "png")
