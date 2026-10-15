# R scripts: figures, statistics and tables for the manuscript

Requires R with the `pacman` package; each script loads (and installs if missing) the rest with `pacman::p_load`.

Run `main.nf` and the two [nf-core/mag](../mag_recovery/README.md) runs first. With the default folders no paths need editing (see [Inputs](#inputs)). Then run everything in order:

```bash
Rscript R-scripts/run_all.r       # each script in its own R session
```

or run single scripts from this folder:

```bash
cd R-scripts
Rscript 01_stats_anova.r          # must run before 03_figure_3_reads.r
Rscript 02_figure_2_qpcr.r
Rscript 03_figure_3_reads.r
Rscript 04_figure_4_family_intersections.r
Rscript 05_figure_S2_aip_families.r
Rscript fit_to_letter.r        # optional: after the figure scripts
Rscript tables/table_counts.r
Rscript tables/table_mags.r
```

## Inputs

Every path is defined once in `00_setup.r`. The two input folders are set at its top:

```r
PIPELINE_RESULTS <- file.path("..", "results")        # main.nf --outdir
MAG_RESULTS      <- file.path("..", "mag_recovery")   # folder holding results_cobinning/ and results_coassembly/
```

The defaults (relative to `R-scripts/`) match a fresh run from the repository root: `main.nf` with `--outdir results`, then `run_cobinning.sh` and `run_coassembly.sh` from `mag_recovery/`. If you used another `--outdir` or moved the folders, edit these two lines only.

| Input | Location |
|---|---|
| read stats (Figure 3, stats, tables) | `<PIPELINE_RESULTS>/reads_based/stats/stats_reads.tsv` |
| CAT quantification (Figure 4, S2) | `<PIPELINE_RESULTS>/coassembly_based/quantification/<mode>/contig_table.<mode>.tsv.gz` |
| qPCR (Figure 2) | `../qPCR/MetaG_qPCR_*.csv` (see [qPCR/README.md](../qPCR/README.md)) |
| MAG quality (table_mags) | `<MAG_RESULTS>/results_{cobinning,coassembly}/GenomeBinning/QC/gunc_checkm_summary.tsv` |

Read categories from main.nf: host, Symbiodiniaceae, bacteria and archaea are kept as they are; viruses and non-mapped reads are pooled as "Other".

## Outputs

| Script | Output |
|---|---|
| `01_stats_anova.r` | `../results/anova_posthoc_contrasts.csv` |
| `02_figure_2_qpcr.r` | `../plots/{png,svg,pdf}/Figure_2.*` |
| `03_figure_3_reads.r` | `../plots/{png,svg,pdf}/Figure_3.*` |
| `04_figure_4_family_intersections.r` | `../plots/{png,svg,pdf}/Figure_4.*` |
| `05_figure_S2_aip_families.r` | `../plots/supplementary/Figure_S2.*` |
| `fit_to_letter.r` | `../plots/pdf/letter/Figure_{2,3,4}.pdf`: copies of the main figures scaled onto a US Letter page |
| `tables/table_counts.r` | `../results/tables/table_counts_by_{mapping,species}.docx` |
| `tables/table_mags.r` | `../results/tables/table_mags_{cobinning,full}.docx` |
| `exploratory/*` | `../plots/extra/extra_NN_*`: products that are not in the paper |

Figure 1 (the workflow figure) is not produced by these scripts.

## Alternative versions of a figure

Each figure script starts with a `DEFAULTS` list, which is the version in the paper, followed by `PARAMS <- DEFAULTS`. To make a variant, change a value in `PARAMS`, for example `PARAMS$RAREFY <- FALSE`. Variants are written to `../plots/extra/` with the changed parameters in the file name, so the paper files are never overwritten.

## Folder content

```text
R-scripts
├── README.md
├── 00_setup.r                          # paths, vocabulary, palettes, shared helpers and data loaders
├── run_all.r                           # runs all scripts below in order
├── 01_stats_anova.r
├── 02_figure_2_qpcr.r
├── 03_figure_3_reads.r
├── 04_figure_4_family_intersections.r
├── 05_figure_S2_aip_families.r
├── fit_to_letter.r
├── tables/
│   ├── table_counts.r
│   └── table_mags.r
└── exploratory/                     # not in the paper; outputs go to ../plots/extra/
    ├── saturation_curves.r          # extra_01
    ├── cat_figure3v2_s3.r           # extra_02, extra_03
    └── run_pavian.r
```
