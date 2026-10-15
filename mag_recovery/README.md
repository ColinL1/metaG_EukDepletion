# MAG recovery (nf-core/mag)

Metagenome-assembled genomes (MAGs) are recovered with [nf-core/mag](https://nf-co.re/mag/3.3.1) **v3.3.1**, run on the output of [`main.nf`](../main.nf). nf-core/mag does not re-assemble. It takes the trimmed reads and the MEGAHIT assemblies that `main.nf` already produced.

There are two runs, which differ only in which assemblies they use:

| Run | Assemblies (from `main.nf`) | Output | Table used by R |
|---|---|---|---|
| **Co-binning** | one per sample, `assembly_based/assembly/` | `results_cobinning/` | `GenomeBinning/QC/gunc_checkm_summary.tsv` (`mags_non_co`) |
| **Co-assembly** | one per group (species + extraction + buffer, replicates pooled), `coassembly_based/assembly/` | `results_coassembly/` | `GenomeBinning/QC/gunc_checkm_summary.tsv` (`mags_co`) |

In both runs, reads from every sample in a group are mapped back to that group's assemblies for binning (`binning_map_mode: group`). A group has the same name as `main.nf`'s co-assembly, e.g. `Ahemp_BBMB_DESS`.

## Settings

Both runs use the same analysis settings ([params/](params/)):

| Step | Setting |
|---|---|
| Read QC | nf-core/mag defaults (fastp, PhiX removal) |
| Binning | MetaBAT2, MaxBin2, CONCOCT; contigs ≥ 1000 bp |
| Bin refinement | DAS Tool; only refined bins are used downstream |
| Bin QC | CheckM + GUNC |
| Taxonomy | GTDB-Tk (release 220, nf-core default) |
| Other | geNomad virus identification, MetaEuk eukaryotic gene prediction, Prokka, QUAST |

The MAG tables (`R-scripts/tables/table_mags.r`) count bins with CheckM completeness ≥ 50 % (medium) or ≥ 90 % (high) and contamination ≤ 10 %.

## Requirements

- Nextflow and Apptainer (or Singularity: use `-profile singularity`)
- A finished `main.nf` run (its `--outdir`)
- Databases, set in both `params/*.json`:
  - `checkm_db`: CheckM reference data ([download](https://data.ace.uq.edu.au/public/CheckM_databases/))
  - `metaeuk_db`: MetaEuk database (we used `MetaEuk_Tara_2020` from the MMseqs2 databases)
  - GTDB-Tk and GUNC databases are downloaded by the pipeline

## Run

From the `mag_recovery/` folder:

```bash
# 1. samplesheets from the main.nf output (writes samplesheets/)
python3 make_samplesheets.py --results ../results

# 2. edit the database paths in params/*.json and the resource limits in nextflow.config

# 3. run nf-core/mag (each run takes several days)
./run_cobinning.sh
./run_coassembly.sh
```

`make_samplesheets.py` options:

| Option | Effect |
|---|---|
| `--run cobinning` / `--run coassembly` | only write the two sheets that run needs (the reads sheet is shared) |
| `--out DIR` | output folder (default `samplesheets/`, which the params files expect) |
| `--no-check` | do not check that the listed files exist (e.g. published symlinks into a work dir on another machine) |

Extra Nextflow options are passed through by the run scripts, e.g. `./run_cobinning.sh -profile singularity`. The scripts set `NXF_SYNTAX_PARSER=v1`, because the nf-core/mag 3.3.1 config does not parse with the strict parser that is the default in recent Nextflow releases.

`nextflow.config` (auto-loaded when launched from this folder) sets the resource limits. `process_resources.config` raises the time and memory limits of the slowest steps (Bowtie2 mapping, CONCOCT, DAS Tool, MetaEuk, geNomad).

## Into the R analysis

`R-scripts/00_setup.r` reads the two `gunc_checkm_summary.tsv` files. After both runs finish, switch `mags_co` / `mags_non_co` to the commented `mag_recovery/results_*` paths below them, then run `Rscript tables/table_mags.r` from `R-scripts/`.
