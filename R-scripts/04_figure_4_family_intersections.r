#! /usr/bin/env Rscript
# ==============================================================================
# Figure 4 — bacterial families shared between the Microbiome-Kit protocols
# (MB / BBMB / SMB). One panel per Aiptasia strain (top) and coral species
# (bottom); each panel is a flipped upset: intersection matrix | relative
# coverage boxplots per protocol | number of families coloured by phylum.
# Presence calls come from prepare_cat_families() (00_setup.r).
# Output: ../plots/{png,svg,pdf}/Figure_4.*
# ==============================================================================
source("00_setup.r")
p_load(patchwork, cowplot, RColorBrewer)

# ---- PARAMS (defaults = paper version) ---------------------------------------
DEFAULTS <- list(
    # family table (see prepare_cat_families)
    MODE = "persample", RAREFY = TRUE, PRESENCE_FRAC = 0.9, MIN_BREADTH = NULL,
    COV_METHOD = "trimmed_mean", TOP_N = 35,
    # layout
    TXT           = 1.5,     # global text multiplier
    AIP_DROP_DESS = TRUE,    # TRUE: Aiptasia Fresh (PBS) panels only
    KEEP_ALL_ROWS = TRUE,    # TRUE: every panel shows all intersection rows in the same order
    SHARED_X      = FALSE,   # FALSE: each boxplot panel gets its own x range
    MIN_FAM       = 1,       # hide intersections with fewer families than this
    PHY_PAL       = "Set3",  # RColorBrewer palette for phyla
    MAT_W         = 0.35,    # matrix width per set
    W_ABUND       = 1.5,     # relative width of the boxplot panel
    W_BARS        = 1.5,     # relative width of the bar panel
    BOX_W         = 0.2,     # box thickness in row units
    BAR_W         = 0.7,     # bar thickness in row units
    BAR_PAD       = 0.12,    # headroom right of the longest bar
    LEG_SCALE     = 1.2      # extra legend multiplier (on top of TXT)
)
PARAMS <- DEFAULTS
invisible(list2env(PARAMS, envir = environment()))

# ---- data ---------------------------------------------------------------------
cat_data <- prepare_cat_families(mode = MODE, rarefy = RAREFY, presence_frac = PRESENCE_FRAC,
                                 min_breadth = MIN_BREADTH, cov_method = COV_METHOD, top_n = TOP_N)
fam <- cat_data$fam

corals   <- c("Acropora", "Pocillopora", "Porites")
aips     <- c("H2", "F003")
ab       <- c("MB", "BBMB", "SMB")
trts     <- unname(protocol_full[ab])
trt_cols <- col_protocol[ab]
nudge    <- c(MB = -0.25, BBMB = 0, SMB = 0.25)

coral_title <- list(
    Acropora    = quote(italic("A. hemprichii")),
    Pocillopora = quote(italic("P. verrucosa")),
    Porites     = quote(italic("P. lobata")))

aip_spec <- expand_grid(strain = aips, buffer = c("Fresh", "DESS")) %>%
    mutate(title = paste(strain, "-", buffer))
if (AIP_DROP_DESS) aip_spec <- filter(aip_spec, buffer == "Fresh")
coral_spec <- tibble(strain = corals, buffer = NA_character_)

get_dat <- function(strain, buffer) {
    x <- fam %>% filter(Strain == .env$strain, Treatment %in% trts)
    if (!is.na(buffer)) x <- x %>% filter(Buffer == .env$buffer)
    x %>% mutate(Trt = factor(Treatment, levels = trts, labels = ab),
                 Family = as.character(Family))
}

all_spec <- bind_rows(aip_spec, coral_spec)
all_dat  <- map2_dfr(all_spec$strain, all_spec$buffer, ~ get_dat(.x, .y) %>% mutate(.buf = .y))

# fixed intersection order: all 3, then pairs, then singles
canon_combos <- function(sets) {
    unlist(lapply(rev(seq_along(sets)), function(k) combn(sets, k, paste, collapse = " + ", simplify = TRUE)))
}

# ---- shared phylum lookup + palette ------------------------------------------------
tax <- cat_data$d %>% distinct(family, phylum) %>%
    filter(!is.na(phylum), phylum != "") %>%
    group_by(family) %>% slice(1) %>% ungroup() %>%
    rename(Family = family)

phy_all <- all_dat %>% distinct(Family) %>% left_join(tax, by = "Family") %>%
    mutate(phylum = forcats::fct_lump_n(factor(replace_na(phylum, "Unknown")), 8))

phy_lv    <- levels(phy_all$phylum)
non_other <- setdiff(phy_lv, "Other")
stopifnot("more phyla than colours in this palette" = length(non_other) <= brewer.pal.info[PHY_PAL, "maxcolors"])
pal_cols  <- brewer.pal(max(3, length(non_other)), PHY_PAL)[seq_along(non_other)]
phy_cols  <- c(setNames(pal_cols, non_other), if ("Other" %in% phy_lv) c(Other = "#999999"))

# common abundance range (used only when SHARED_X = TRUE)
ab_all <- all_dat %>%
    group_by(Strain, .buf, Trt) %>% mutate(n_samp = n_distinct(sample)) %>%
    group_by(Strain, .buf, Trt, Family, n_samp) %>%
    summarise(rel_sum = sum(rel), .groups = "drop") %>%
    mutate(v = log10(rel_sum / n_samp * 1e6 + 1))
xlim_shared <- range(ab_all$v, na.rm = TRUE) + c(-0.1, 0.1)

# ---- theme: all text sizes derive from TXT ---------------------------------------
base <- theme_minimal(base_size = 11 * TXT) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3),
          axis.text  = element_text(size = 10 * TXT, color = "black"),
          axis.title = element_text(size = 11 * TXT),
          legend.position = "none")
no_y <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())

# ---- one flipped upset (matrix | boxplots | bars) ---------------------------------
make_panel <- function(strain, buffer, title) {
    dat <- get_dat(strain, buffer)

    n_s <- dat %>% distinct(Trt, sample) %>% count(Trt, name = "n_samp")
    ab_wide <- dat %>%
        group_by(Trt, Family) %>% summarise(rel_sum = sum(rel), .groups = "drop") %>%
        left_join(n_s, by = "Trt") %>%
        mutate(v = log10(rel_sum / n_samp * 1e6 + 1)) %>%
        select(Family, Trt, v) %>%
        pivot_wider(names_from = Trt, values_from = v, names_prefix = "ab_")

    m <- dat %>% distinct(Trt, Family) %>% mutate(v = 1L) %>%
        pivot_wider(names_from = Trt, values_from = v, values_fill = 0L) %>%
        left_join(ab_wide, by = "Family") %>%
        left_join(phy_all %>% filter(Family %in% dat$Family), by = "Family")
    sets <- ab[ab %in% names(m)]
    if (length(sets) < 2) {
        message("skipping ", strain, " ", buffer, ": < 2 of MB/BBMB/SMB present (", nrow(dat), " rows)")
        return(NULL)
    }

    set_n <- colSums(m[, sets])
    m$combo <- apply(m[, sets], 1, function(r) paste(sets[r == 1], collapse = " + "))
    m$deg   <- rowSums(m[, sets])

    keep <- m %>% count(combo) %>% filter(n >= MIN_FAM) %>% pull(combo)
    m <- m %>% filter(combo %in% keep)
    if (nrow(m) == 0) {
        message("skipping ", strain, " ", buffer, ": MIN_FAM removed every intersection")
        return(NULL)
    }

    combos <- canon_combos(sets)
    if (!KEEP_ALL_ROWS) combos <- combos[combos %in% m$combo]
    combo_tab <- tibble(combo = combos) %>%
        mutate(deg  = str_count(combo, fixed("+")) + 1,
               n    = map_int(combo, ~ sum(m$combo == .x)),
               ypos = n() - row_number() + 1)
    K <- nrow(combo_tab)
    m$ypos <- combo_tab$ypos[match(m$combo, combo_tab$combo)]

    yscale <- scale_y_continuous(limits = c(0.5, K + 0.5), breaks = seq_len(K), expand = c(0, 0))

    # matrix
    grid <- expand_grid(ypos = combo_tab$ypos, set = sets) %>%
        left_join(combo_tab %>% select(ypos, combo), by = "ypos") %>%
        mutate(x  = match(set, sets),
               on = str_detect(paste0(" + ", combo, " + "), fixed(paste0(" + ", set, " + "))))
    segs <- grid %>% filter(on) %>% group_by(ypos) %>%
        summarise(xmin = min(x), xmax = max(x), .groups = "drop") %>% filter(xmin < xmax)

    p_mat <- ggplot(grid, aes(x, ypos)) +
        geom_segment(data = segs, aes(x = xmin, xend = xmax, y = ypos, yend = ypos),
                     linewidth = 0.9, inherit.aes = FALSE) +
        geom_point(aes(color = on), size = 2.5 * TXT^0.7) +
        scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "grey85"), guide = "none") +
        scale_x_continuous(breaks = seq_along(sets), position = "top",
                           labels = paste0(sets, "\n(", set_n, ")"),
                           limits = c(0.5, length(sets) + 0.5), expand = c(0, 0)) +
        yscale + labs(x = NULL) + base + no_y +
        theme(axis.text.x = element_text(size = 10 * TXT, color = "black"))

    # boxplots
    long <- m %>%
        select(Family, ypos, deg, starts_with("ab_")) %>%
        pivot_longer(starts_with("ab_"), names_to = "Trt", names_prefix = "ab_", values_to = "v") %>%
        filter(!is.na(v)) %>%
        mutate(Trt = factor(Trt, levels = ab),
               y   = ypos + ifelse(deg > 1, nudge[as.character(Trt)], 0),
               grp = paste(ypos, Trt))

    xr <- if (SHARED_X) {
        xlim_shared
    } else {
        st <- long %>% group_by(grp) %>%
            summarise(lo = boxplot.stats(v)$stats[1], hi = boxplot.stats(v)$stats[5], .groups = "drop")
        r <- c(min(st$lo), max(st$hi))
        r + c(-1, 1) * max(diff(r) * 0.05, 0.05)
    }

    p_abund <- ggplot(long, aes(v, y, group = grp, fill = Trt)) +
        geom_boxplot(orientation = "y", width = BOX_W, position = "identity",
                     outlier.shape = NA, color = "grey25", linewidth = 0.3, alpha = 0.75, na.rm = TRUE) +
        scale_fill_manual(values = trt_cols, name = "Protocol") +
        scale_x_continuous(labels = scales::label_number(accuracy = 0.1, drop0trailing = TRUE),
                           breaks = scales::breaks_pretty(n = 4)) +
        scale_y_continuous(breaks = seq_len(K), expand = c(0, 0)) +
        coord_cartesian(xlim = xr, ylim = c(0.5, K + 0.5), expand = FALSE) +
        labs(x = expression(log[10] * "(rel.cov.×10"^6 * "+1)")) + base + no_y

    # bars
    p_bars <- ggplot(m, aes(y = ypos, fill = phylum)) +
        geom_bar(orientation = "y", width = BAR_W, color = "white", linewidth = 0.2) +
        geom_text(data = filter(combo_tab, n > 0), aes(x = n, y = ypos, label = n),
                  hjust = -0.25, size = 3.3 * TXT, inherit.aes = FALSE) +
        scale_fill_manual(values = phy_cols, name = "Phylum", drop = FALSE) +
        scale_x_continuous(limits = c(0, max(combo_tab$n) * (1 + BAR_PAD)),
                           breaks = scales::breaks_pretty(n = 3), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0.5, K + 0.5), breaks = seq_len(K), expand = c(0, 0)) +
        coord_cartesian(clip = "off") +
        labs(x = "Families") + base + no_y

    row <- (p_mat | p_abund | p_bars) + plot_layout(widths = c(length(sets) * MAT_W, W_ABUND, W_BARS))

    list(K = K,
         plot = wrap_elements(full = patchwork::patchworkGrob(row)) +
             ggtitle(title) +
             theme(plot.title = element_text(size = 14 * TXT, hjust = 0.5)))
}

# ---- build panels ----------------------------------------------------------------
res_aip <- lapply(seq_len(nrow(aip_spec)), function(i) make_panel(aip_spec$strain[i], aip_spec$buffer[i], aip_spec$title[i]))
res_cor <- lapply(corals, function(s) make_panel(s, NA_character_, coral_title[[s]]))
res_aip <- Filter(Negate(is.null), res_aip)
res_cor <- Filter(Negate(is.null), res_cor)
stopifnot("no Aiptasia panels were built" = length(res_aip) > 0)

# ---- one shared legend ---------------------------------------------------------------
LS <- LEG_SCALE * TXT
leg_plot <- ggplot() +
    geom_col(data = data.frame(phylum = factor(phy_lv, levels = phy_lv), x = 1, y = 1), aes(x, y, fill = phylum)) +
    geom_point(data = data.frame(Trt = factor(ab, levels = ab), x = 1, y = 1), aes(x, y, color = Trt), shape = 15, size = 4) +
    scale_color_manual(name = "Protocol", values = trt_cols) +
    scale_fill_manual(name = "Phylum", values = phy_cols, drop = FALSE) +
    guides(color = guide_legend(order = 1, override.aes = list(size = 4 * LS)),
           fill  = guide_legend(order = 2, ncol = 1)) +
    theme(legend.position = "right", legend.box = "horizontal", legend.box.just = "top",
          legend.text     = element_text(size = 10 * LS),
          legend.title    = element_text(size = 11 * LS, face = "bold"),
          legend.key.size = unit(1.1 * LS, "lines"))
legend_el <- wrap_elements(full = cowplot::get_legend(leg_plot))

# ---- assemble: top row = Aiptasia + legend, bottom row = corals -----------------------
top_cells <- c(lapply(res_aip, `[[`, "plot"), list(legend_el))
bot_cells <- lapply(res_cor, `[[`, "plot")
n_cols <- max(length(top_cells), length(bot_cells))
pad_row <- function(cells) {
    while (length(cells) < n_cols) cells <- c(cells, list(plot_spacer()))
    cells
}
h_top <- max(vapply(res_aip, `[[`, numeric(1), "K")) + 4
h_bot <- max(vapply(res_cor, `[[`, numeric(1), "K")) + 4
figure <- wrap_plots(wrap_plots(pad_row(top_cells), nrow = 1), wrap_plots(pad_row(bot_cells), nrow = 1),
                     ncol = 1, heights = c(h_top, h_bot))

out <- output_target("Figure_4", PARAMS, DEFAULTS, kind = "main")
save_figure(figure, out$name, out$kind, width = 22, height = 14)
