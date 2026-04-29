# folder with precomputed minimap2 and bowtie2 index

## Steps ~
- download reference genome with ncbi dataset

```bash
datasets download genome taxon 6125 --reference
# unzip ncbi_dataset.zip
# cat ncbi_dataset/*/*.fna > 6125_RefSeq_genome.fna && pigz 6125_RefSeq_genome.fna
```

- build index with bowtie2
  
```bash
bowtie2-build 6125_genome.fna.gz 6125_genome -p number_threads --large-index
```

- build index with minimap2
```bash
minimap2 -x asm5 -d 6125_genome.mmi 6125_genome.fna.gz -t number_threads
```