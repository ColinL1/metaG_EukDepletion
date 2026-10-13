# Collection of modules used in the pipeline

## Module descriptions

Used by `main.nf`:

| Module | Processes | Used in | Description |
|---|---|---|---|
| [trimmomatic.nf](trimmomatic.nf) | `TRIMMOMATIC` | `main.nf` | Adapter and quality trimming of paired-end reads |
| [bowtie2.nf](bowtie2.nf) | `BOWTIE2_MAP_HOST`, `BOWTIE2_MAP_SYM` | reads based | Maps read pairs to the host, then to the Symbiodiniaceae reference; writes mapped and unmapped pairs |
| [assembly_megahit.nf](assembly_megahit.nf) | `MEGAHIT_PE` | assembly based, coassembly based | De novo assembly of a sample, or co-assembly of a group |
| [minimap2.nf](minimap2.nf) | `MINIMAP2_MAP_HOST`, `MINIMAP2_MAP_SYM`, `MINIMAP2_MAP_CATALOGUE` | assembly based, coassembly based | Maps contigs to the host and Symbiodiniaceae references; maps reads to the bacterial contig references |
| [samtools_split.nf](samtools_split.nf) | `SPLIT_BAM_HOST`, `SPLIT_BAM_SYM` | assembly based | Splits a contig BAM into mapped and unmapped contigs |
| [kaiju.nf](kaiju.nf) | `KAIJU_PE`, `KAIJU_SE`, `KAIJU_ADD_TAXON_NAMES`, `KAIJU_SUMMARISE` | reads based, assembly based | Kaiju classification of reads or contigs, taxon names, and a per-sample summary with the ids of bacteria, archaea, viruses and non-mapped sequences |
| [seqkit.nf](seqkit.nf) | `SEQKIT_SPLIT_KAIJU_PE`, `SEQKIT_SPLIT_KAIJU_SE`, `SEQKIT_STATS`, `SEQKIT_PREFIX_CONTIGS` | all three | Splits reads or contigs by Kaiju group, counts sequences per file, prefixes contig names with the group id |
| [concat_fastq.nf](concat_fastq.nf) | `CONCAT_FASTQ` | coassembly based | Concatenates the trimmed reads of all samples of a group (species, method, buffer) |
| [cat.nf](cat.nf) | `CAT`, `CAT_ADD_NAMES`, `CAT_FILTER_BACTERIA`, `CAT_PREP_FAMILY`, `CAT_MERGE_TABLES` | coassembly based | CAT classification of co-assembled contigs, taxonomy names, bacterial contig subset, family table, merge with coverage |
| [mmseqs2.nf](mmseqs2.nf) | `POOL_CONTIGS`, `MMSEQS_LINCLUST` | coassembly based | Pools the bacterial contigs per species and dereplicates them into a catalogue |
| [minimap2_index.nf](minimap2_index.nf) | `MINIMAP2_INDEX` | coassembly based | Minimap2 index of a bacterial contig reference |
| [coverm.nf](coverm.nf) | `COVERM_CONTIG` | coassembly based | Per-contig read counts and coverage |

Not used by `main.nf`:

| Module | Description |
|---|---|
| [K2.nf](K2.nf), [bracken.nf](bracken.nf), [load_k2_db_shm.nf](load_k2_db_shm.nf) | Kraken2 / Bracken classification, used by `subworkflows/k2b.nf`. Superseded by Kaiju |
| [bwa.nf](bwa.nf) | BWA mapping. Superseded by Bowtie2 / Minimap2 |
| [exploratory](exploratory) | Outdated or undeveloped modules, no guarantees. See its own README |

## Folder content

```text
modules
├── README.md                   <-- you are here.
├── assembly_megahit.nf
├── bowtie2.nf
├── bracken.nf
├── bwa.nf
├── cat.nf
├── concat_fastq.nf
├── coverm.nf
├── exploratory
├── K2.nf
├── kaiju.nf
├── load_k2_db_shm.nf
├── minimap2_index.nf
├── minimap2.nf
├── mmseqs2.nf
├── samtools_split.nf
├── seqkit.nf
└── trimmomatic.nf
```
