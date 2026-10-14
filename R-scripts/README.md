# R scripts: figures, statistics and tables for the manuscript

All scripts are run from this folder, in the `metaG_EukDepletion` conda env:

```bash
cd R-scripts
conda activate metaG_EukDepletion
Rscript 01_stats_anova.r          # must run before 03_figure_3_reads.r
Rscript 02_figure_2_qpcr.r
Rscript 03_figure_3_reads.r
Rscript 04_figure_4_family_intersections.r
Rscript 05_figure_S2_aip_families.r
Rscript fit_to_letter.r        # after the figure scripts
Rscript tables/table_counts.r
Rscript tables/table_mags.r
```

## Inputs

Every path is defined once in `00_setup.r`.

| Input | Location |
|---|---|
| main.nf results | `PIPELINE_RESULTS` at the top of `00_setup.r` (currently `../results_binac2/results`). Change this one line when the results folder moves. |
| read stats (Figure 3, stats, tables) | `<PIPELINE_RESULTS>/reads_based/stats/stats_reads.tsv` |
| CAT quantification (Figure 4, S2) | `<PIPELINE_RESULTS>/coassembly_based/quantification/<mode>/contig_table.<mode>.tsv.gz` |
| qPCR (Figure 2) | `../qPCR/MetaG_qPCR_*.csv` |
| MAG quality (table_mags) | `../manual_piplines/binning/nf_core_MAG/gunc_checkm_summary_{co,non_co}.tsv` |

Read categories from main.nf: host, Symbiodiniaceae, bacteria and archaea are kept as they are; viruses and non-mapped reads are pooled as "Other".

`RSCRIPTS_INPUT=legacy Rscript <script>` runs a script on the pre-pipeline files under `manual_piplines/` instead. This mode exists only to check that the scripts reproduce the submitted figures.

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
