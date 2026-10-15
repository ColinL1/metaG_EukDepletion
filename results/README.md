# Pipeline output

Default `--outdir` of `main.nf` (not in git). Outputs are grouped by branch: reads based, assembly based and co-assembly based.

```text
results
├── trimmed_fastq/<species>/       # Trimmomatic output, input to all three branches
├── reads_based/
│   ├── mapping/
│   │   ├── <species>/             # Bowtie2 host BAM, mapped / unmapped reads
│   │   ├── symbiodiniaceae/       # Bowtie2 Symbiodiniaceae BAM, mapped / unmapped reads
│   │   └── bacteria-kaiju/        # Kaiju output, taxon names, summary; bacteria / archaea / viruses / non-mapped reads
│   └── stats/                     # stats_reads.tsv (seqkit stats of every fraction), kaiju_summary_reads.tsv
├── assembly_based/
│   ├── assembly/                  # MEGAHIT assembly per sample
│   ├── mapping/
│   │   ├── <species>/             # Minimap2 host BAM, mapped / unmapped contigs
│   │   ├── symbiodiniaceae/       # Minimap2 Symbiodiniaceae BAM, mapped / unmapped contigs
│   │   └── bacteria-kaiju/        # Kaiju output, taxon names, summary; bacteria / archaea / viruses / non-mapped contigs
│   └── stats/                     # stats_assembly.tsv, kaiju_summary_assembly.tsv
└── coassembly_based/
    ├── concatenated_fastq/        # Reads concatenated per species, method and buffer
    ├── assembly/                  # MEGAHIT co-assembly per group
    ├── CAT/<species>/<group>/     # CAT classification, taxonomy names, bacteria/ subset
    └── quantification/
        ├── persample/             # Reads mapped to the group's own bacterial contigs
        ├── pooled/                # Reads mapped to the species catalogue (catalogue/)
        └── mapping_indexes/       # Minimap2 indexes of both references
```

`persample/` and `pooled/` each hold `family_table`, `contig_table`, `mapping_summary`, `coverm/` and `flagstat/`.

> Reads, BAMs, assemblies and CAT files are published as symlinks into the Nextflow `work/` folder, so they disappear if `work/` is deleted. Tables (stats, Kaiju summaries, quantification, catalogue) are copied. Set `mode: 'copy'` in a module's `publishDir` to copy all files.
