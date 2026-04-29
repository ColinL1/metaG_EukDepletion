# Collection of modules used in the pipeline

## Module descriptions

- [**trimmomatic.nf**](../modules/trimmomatic.nf) - Trims raw paired-end Illumina reads (adapter removal, quality filtering)
- [**concat_fastq.nf**](../modules/concat_fastq.nf) - Concatenates multi-lane FASTQ files into a single R1/R2 pair per sample
- [**bowtie2.nf**](../modules/bowtie2.nf) - Maps PE reads to the host reference genome (Aiptasia or Scleractinia / Symbiodiniaceae); outputs mapped and unmapped reads separately
- [**minimap2.nf**](../modules/minimap2.nf) - Maps assembled contigs to host / Symbiodiniaceae references
- [**samtools_split.nf**](../modules/samtools_split.nf) - Splits a BAM file into mapped and unmapped reads (separate processes for host and Symbiodiniaceae)
- [**assembly_megahit.nf**](../modules/assembly_megahit.nf) - De novo assembly of paired-end reads with MEGAHIT
- [**kaiju.nf**](../modules/kaiju.nf) - Taxonomic classification of reads against the Kaiju reference database (single-end and paired-end processes)
- [**cat.nf**](../modules/cat.nf) - Classifies assembled contigs with CAT against GTDB and adds taxonomy names to the output


## Folder content

```text
modules
├── README.md                   <-- you are here.
├── adapters
│   └── Truseq_V.3_edited.fa    # edited variant of adapter references used in Trimmomatic
├── assembly_megahit.nf
├── bowtie2.nf
├── bracken.nf                  # Functional, not used in pipeline - superseded by Kaiju
├── bwa_index.nf                # Not used in pipeline - superseded by Bowtie2 / Minimap2
├── bwa.nf                      # Not used in pipeline - superseded by Bowtie2 / Minimap2
├── concat_fastq.nf
├── cat.nf
├── exploratory                 # Not actively used. No guarantee of functionality
│   └── templates
│       └── parse_json.py       # Dropped in favour of the seqkit stats approach
├── K2.nf                       # Functional, not used in pipeline - superseded by Kaiju
├── kaiju.nf
├── load_k2_db_shm.nf
├── minimap2_index.nf
├── minimap2.nf
├── samtools_split.nf
└── trimmomatic.nf
```
