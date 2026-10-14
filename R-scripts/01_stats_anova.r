#! /usr/bin/env Rscript
# ==============================================================================
# Read-fraction statistics behind the significance stars of Figure 3.
#
# Per type (aiptasia / corals) and per mapping category: two-way ANOVA
# num_seq_percent ~ extraction * grouping, then emmeans pairwise contrasts of
# extraction within each grouping (Tukey-adjusted), keeping the a-priori
# contrasts (each protocol vs BT, BBMB/SMB vs MB).
#
# Output: ../results/anova_posthoc_contrasts.csv  (read by 03_figure_3_reads.r)
# ==============================================================================
source("00_setup.r")
p_load(emmeans)

reads <- load_read_stats()

test_data <- reads %>%
    filter(mapping != "total") %>%
    group_by(sample) %>%
    mutate(num_seq_percent = num_seqs / sum(num_seqs, na.rm = TRUE) * 100) %>%
    ungroup() %>%
    mutate(type       = factor(ifelse(species %in% c("H2", "F003"), "aiptasia", "corals"), levels = c("aiptasia", "corals")),
           grouping   = factor(grouping, levels = grouping_levels),
           extraction = relevel(factor(extraction, levels = extraction_levels), ref = "BT"))

wanted <- c("BT - BenzoBT", "BT - MB", "BT - BBMB", "BT - SMB",
            "MB - BBMB", "BBMB - MB", "MB - SMB", "SMB - MB")

anova_results <- list()
posthoc_results <- list()
for (i in levels(test_data$type)) {
    sub_data <- test_data %>% filter(type == i)
    for (m in unique(sub_data$mapping)) {
        sub_df <- sub_data %>% filter(mapping == m)
        fit <- aov(num_seq_percent ~ extraction * grouping, data = sub_df)
        anova_results[[i]][[m]] <- summary(fit)
        em <- emmeans(fit, pairwise ~ extraction | grouping, adjust = "tukey")
        posthoc_results[[i]][[m]] <- as.data.frame(em$contrasts) %>%
            filter(contrast %in% wanted) %>%
            mutate(mapping = m, type = i)
    }
}
print(anova_results)

# report MB-vs-(BBMB|SMB) contrasts as "MB - X" (sign of the estimate flipped to match)
results <- bind_rows(lapply(posthoc_results, bind_rows)) %>%
    select(type, mapping, grouping, contrast, estimate, SE, df, t.ratio, p.value) %>%
    mutate(
        contrast = gsub("\\s+", " ", trimws(contrast)),
        swap_mb  = grepl("^(SMB|BBMB)\\s*-\\s*MB$", contrast),
        estimate = ifelse(swap_mb, -estimate, estimate),
        contrast = ifelse(swap_mb, sub("^(SMB|BBMB)\\s*-\\s*MB$", "MB - \\1", contrast), contrast),
        contrast = factor(contrast, levels = c("BT - BenzoBT", "BT - MB", "BT - BBMB", "BT - SMB", "MB - BBMB", "MB - SMB"))
    ) %>%
    arrange(mapping, grouping, contrast) %>%
    select(-swap_mb)

dir.create(paths$results, showWarnings = FALSE, recursive = TRUE)
write.csv(results, paths$anova_csv)
message("wrote ", paths$anova_csv)
