# metaG_EukDepletion

[![DOI](static/zenodo.19912536.svg)](https://doi.org/10.5281/zenodo.19912536)

Code used for analysis of kingdom taxonomy to assess the development of a cnidarian metagenomic protocol for efficient eukaryotic DNA removal

---

## Flowchart

```mermaid
flowchart LR
 subgraph input["Extraction methods"]
        n7["BT | BenzoBT | MB | BBMB | SMB"]
        A["Input Reads"]
  end
 subgraph subGraph0["reads based"]
        C@{ label: "<div style='font-weight:bold'>Bowtie2</div> Filter Host Reads" }
        D@{ label: "<div style='font-weight:bold'>Bowtie2</div> Filter Symbiodiniaceae Reads" }
        E@{ label: "<div style='font-weight:bold'>Kaiju</div> Classify Reads" }
        E2@{ label: "<div style='font-weight:bold'>seqkit</div> Split: Bacteria / Archaea / Viruses / Non-mapped" }
        S1@{ label: "<div style='font-weight:bold'>seqkit stats</div> Reads per Fraction" }
  end
 subgraph subGraph1["assembly based"]
        F@{ label: "<div style='font-weight:bold'>Megahit</div> Assemble Reads" }
        G@{ label: "<div style='font-weight:bold'>Minimap2</div> Filter Host Sequences" }
        H@{ label: "<div style='font-weight:bold'>Minimap2</div> Filter Symbiodiniaceae Sequences" }
        J@{ label: "<div style='font-weight:bold'>Kaiju</div> Classify Sequences" }
        J2@{ label: "<div style='font-weight:bold'>seqkit</div> Split: Bacteria / Archaea / Viruses / Non-mapped" }
        S2@{ label: "<div style='font-weight:bold'>seqkit stats</div> Contigs per Fraction" }
  end
 subgraph subGraph2["coassembly based"]
        K["Concatenate Reads per Group"]
        L@{ label: "<div style='font-weight:bold'>Megahit</div> Co-assemble Reads" }
        M@{ label: "<div style='font-weight:bold'>CAT</div> Classify Contigs" }
        BF["Filter Bacterial Contigs"]
        T1@{ label: "<div style='font-weight:bold'>MMseqs2 linclust</div> Pool and Dereplicate Contigs per Species" }
        T2@{ label: "<div style='font-weight:bold'>Minimap2</div> Map Reads (Per-sample / Pooled Reference)" }
        T3@{ label: "<div style='font-weight:bold'>CoverM</div> Contig Coverage" }
        T4["Merge CAT Family and Coverage"]
  end
 subgraph mags["MAGs recovery"]
        N["nf-core/mag"]
        O@{ label: "<div style='font-weight:bold'>METABAT2 / MAXBIN2 / CONCOCT</div> Binning modes (multi-sample, co-assembly)" }
  end
 subgraph qPCR["qPCR"]
        Q["qPCR Ct data"]
        R["Outlier detection"]
  end
    n7 --> A
    A ==> B["Trim Reads"]
    B ==> C & F & K
    C --> D
    D --> E
    E --> E2
    B & C & D & E2 --> S1
    F --> G
    G --> H
    H --> J
    J --> J2
    F & G & H & J2 --> S2
    K --> L & T2
    L --> M & BF
    M --> BF
    BF --> T1 & T2 & T4
    T1 --> T2
    T2 --> T3
    T3 --> T4
    T4 --> T5@{ label: "<div style='font-weight:bold'>vegan drarefy</div> Family Presence Call per Species" }
    S1 & S2 & T5 --> P@{ label: "<div style='font-weight:bold'>R tidyr / ggplot2 / emmeans</div> Statistics and Figures" }
    B & F & K & L -.-> N
    N --> O
    O --> P
    Q --> R
    R --> P

    A@{ shape: procs}
    Q@{ shape: procs}
    N@{ shape: paper-tape}
    classDef external stroke-dasharray: 5 5
    class T5,P,N,O,Q,R external
```

Dashed nodes are run outside this Nextflow pipeline (R analysis, MAG recovery — see [mag_recovery/](mag_recovery/README.md), qPCR).

---

## Prerequisites

- Nextflow
- Conda / mamba for environment management
- Kaiju database files:
  - `nodes.dmp`
  - `kaiju_db_refseq.fmi`
- Pre-built reference genome indexes (Bowtie2 for reads-based; Minimap2 for assembly-based):
  - Host - Aiptasiidae or Scleractinia
  - Symbiont - Symbiodiniaceae

> Taxon IDs used in the manuscript: `42823` (Aiptasiidae), `6127` (Acropora), `46730` (Pocillopora), `46719` (Porites), `252141` (Symbiodiniaceae)

---

## Pipeline overview

The pipeline runs **three parallel branches** from trimmed reads:

| Branch | Subworkflow | Steps |
|---|---|---|
| Reads-based | `MAP_CONSE_PE` | Bowtie2 (host) → Bowtie2 (Symbiodiniaceae) → Kaiju → bacteria / archaea / viruses / non-mapped split → seqkit stats |
| Assembly-based | `MAP_CONSE_ASSEMBLY` | MEGAHIT → Minimap2 (host) → Minimap2 (Symbiodiniaceae) → Kaiju → bacteria / archaea / viruses / non-mapped split → seqkit stats |
| Co-assembly / CAT | `CAT_WORKFLOW` | Reads grouped by species+method+buffer → CAT contig classification → read mapping of bacterial contigs (Minimap2, CoverM) → family tables |

### Input file naming assumptions

FASTQ files follow this naming pattern:

```
{species}_{replicate}_{method}_{buffer}_{1,2}.fq.gz
```

Species codes currently recognised by the pipeline:

| Code | Mapped to |
|---|---|
| `Ac`, `Acro`, `Ahemp` | Acropora |
| `Po`, `Por`, `Plob` | Porites |
| `Poci`, `Pr`, `Pver` | Pocillopora |
| `F003`, `F3` | F003 (Aiptasia) |
| `H2` | H2 (Aiptasia) |

> Edit input channel logic in [main.nf](main.nf) to add/remove.
---

## Step-by-step setup

### 1. Clone the repository

```bash
git clone https://github.com/ColinL1/metaG_EukDepletion.git
cd metaG_EukDepletion
```

### 2. Build the Kaiju database - [GitHub](https://github.com/bioinformatics-centre/kaiju)

Pre-built indexes can be downloaded directly from the Kaiju [downloads page](https://bioinformatics-centre.github.io/kaiju/downloads.html):

```bash
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_nr_2024-08-25.tgz
tar xvf kaiju_db_nr_2024-08-25.tgz
```

Or build from scratch with:

```bash
kaiju-makedb -s refseq
```

### 3. Download reference genomes - [NCBI Datasets](https://github.com/ncbi/datasets)

```bash
# Example for Acropora (taxon 6127)
datasets download genome taxon 6127 --reference
unzip ncbi_dataset.zip
cat ncbi_dataset/*/*.fna > 6127_genome.fna && pigz 6127_genome.fna
```

Repeat for each reference taxon (host and symbiont).

### 4. Build reference indexes

**Bowtie2** (reads-based branch):

```bash
bowtie2-build 6127_genome.fna.gz 6127_genome -p number_threads --large-index
```

**Minimap2** (assembly-based branch):

```bash
minimap2 -x asm5 -d 6127_genome.mmi 6127_genome.fna.gz -t number_threads
```

### 5. Configure and run

Copy and edit the example config:

```bash
cp nextflow.config.example nextflow.config
# Edit paths for references, Kaiju DB, and resource limits
```

```bash
nextflow run main.nf -profile conda \
  --input path/to/reads/ \
  --outdir path/to/results/
```

Use `-profile apptainer` instead to run every tool from its container (Apptainer/Singularity) rather than conda.

> **Important:** Adjust CPU and RAM limits in `nextflow.config` and `process_resources.config` to match your machine.

### 6. MAG recovery (nf-core/mag)

MAGs are recovered by running nf-core/mag on the `main.nf` output (trimmed reads, per-sample assemblies and co-assemblies). Settings, samplesheet builder and run scripts are in [mag_recovery/](mag_recovery/README.md).
