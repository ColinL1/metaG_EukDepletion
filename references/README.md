# Folder with precomputed Minimap2 and Bowtie2 indexes

## Steps

1. Download reference genome with NCBI Datasets

```bash
datasets download genome taxon 6127 --reference
unzip ncbi_dataset.zip
cat ncbi_dataset/*/*.fna > 6127_genome.fna && pigz 6127_genome.fna
```

2. Build index with Bowtie2
  
```bash
bowtie2-build 6127_genome.fna.gz 6127_genome -p number_threads --large-index
```

3. Build index with Minimap2

```bash
minimap2 -x asm5 -d 6127_genome.mmi 6127_genome.fna.gz -t number_threads
```

>**Note**:
The taxids used for manuscript are: 42823 (Aiptasiidae), 6127 (Acropora),  46730 (Pocillopora), 46719 (Porites), 252141 (Symbiodiniaceae)

#### Kaiju databse

References database indexes created from standard `refseq` reference database with the command:

```bash
kaiju-makedb -s refseq
```

Pre-built indexes for the reference database can be downloaded from the official Kaiju [website](https://bioinformatics-centre.github.io/kaiju/downloads.html). 
</details>
