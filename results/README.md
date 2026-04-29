# Output folder for the pipeline

All nextflow pipeline outputs are written here. Each subfolder corresponds to a processing stage:

```text
results
├── trimmed_fastq/         # Trimmomatic output — quality-trimmed paired-end reads
├── concatenated_fastq/    # Per-sample concatenated reads (multi-lane merging)
├── mapping/               # Bowtie2 / Minimap2 BAM files and split mapped/unmapped reads
├── assembly/              # MEGAHIT contig assemblies
└── CAT/                   # CAT contig classification output
```

> Outputs are written as symlinks by default to avoid duplicating large files. Set `mode: 'copy'` in the relevant module if hard copies are needed.