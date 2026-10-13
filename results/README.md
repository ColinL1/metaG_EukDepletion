# Output folder for the pipeline

All nextflow pipeline outputs are written here, grouped by the three branches of the pipeline:

```text
results
├── trimmed_fastq/                 # Trimmomatic output, input to all three branches
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

> Outputs are written as symlinks by default to avoid duplicating large files. Set `mode: 'copy'` in the relevant module if hard copies are needed.