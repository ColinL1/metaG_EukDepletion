# Collection of modules used in the pipeline

Folder content:

```text
modules
├── README.md <-- you are here. 
├── adapters
│   └── Truseq_V.3_edited.fa # differnt variant of adapter references used in trimmomatic
├── assembly_megahit.nf 
├── bowtie2.nf
├── bracken.nf # Functional, not used in pipline in favour of Kaiju
├── bwa_index.nf
├── bwa.nf # Functional, not used in pipline in favour of bowtie2 and minimap2
├── cat_fastq.nf
├── cat.nf
├── exploratory # A subset of modules or modules ideas. Not actively used in the pipeline. No guarantee any of this code is functional
│   └── templates # Droped from updated logic. 
│       └── parse_json.py Function to parse fastp json file out for stats. Droped in favour of a simpler "seqkit stats" arpoach
├── K2.nf # Functional, not used in pipline in favour of Kaiju
├── kaiju.nf
├── load_k2_db_shm.nf 
├── minimap2_index.nf
├── minimap2.nf
├── multiqc.nf
├── samtools_split.nf
└── trimmomatic.nf
```