# Folder with precomputed reference indexes

This folder holds the precomputed indexes required by the pipeline: **Bowtie2** indexes (reads-based branch), **Minimap2** indexes (assembly-based branch), the **Kaiju** database (taxonomic classification), and the **CAT** database (contig taxonomic annotation).

> All paths to these files must be set as `params` in `nextflow.config`. See `nextflow.config.example` for the expected parameter names.

---

## Reference genome indexes (Bowtie2 and Minimap2)

Taxon IDs used in the manuscript: `42823` (Aiptasiidae), `6127` (Acropora), `46730` (Pocillopora), `46719` (Porites), `252141` (Symbiodiniaceae)

### 1. Download reference genome

```bash
# Example for Acropora (taxon 6127) — repeat for each taxon
datasets download genome taxon 6127 --reference
unzip ncbi_dataset.zip
cat ncbi_dataset/*/*.fna > 6127_genome.fna && pigz 6127_genome.fna
```

### 2. Build Bowtie2 index

```bash
bowtie2-build 6127_genome.fna.gz 6127_genome -p number_threads --large-index
```

### 3. Build Minimap2 index

```bash
minimap2 -x asm5 -d 6127_genome.mmi 6127_genome.fna.gz -t number_threads
```

---

## Kaiju database

Pre-built indexes can be downloaded from the Kaiju [downloads page](https://bioinformatics-centre.github.io/kaiju/downloads.html):

```bash
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_nr_2024-08-25.tgz
tar xvf kaiju_db_nr_2024-08-25.tgz
```

Or build from scratch with:

```bash
kaiju-makedb -s refseq
```

---

## CAT database (GTDB)

The [CAT pack](https://github.com/MGXlab/CAT_pack) (Contig Annotation Tool) is used for taxonomic classification of assembled contigs. The pipeline uses a GTDB-based database, which must expose two subdirectories: `db/` and `tax/`.

Set `params.cat_db` to the database root directory in `nextflow.config`.

### Option 1 — Download prebuilt GTDB database

```bash
wget tbb.bio.uu.nl/tina/CAT_pack_prepare/20231120_CAT_gtdb.tar.gz   # GTDB release 214
tar -xvzf 20231120_CAT_gtdb.tar.gz
```

### Option 2 — Build from scratch

```bash
# Download and process GTDB source files
CAT_pack download --db gtdb -o path/to/gtdb_data_dir

# Build the CAT pack database
CAT_pack prepare \
  --db_fasta path/to/gtdb_data_dir/gtdb_proteins_aa_reps.faa \
  --names path/to/gtdb_data_dir/names.dmp \
  --nodes path/to/gtdb_data_dir/nodes.dmp \
  --acc2tax path/to/gtdb_data_dir/acc2taxid.txt.gz \
  --db_dir 20231120_CAT_gtdb
```

The resulting directory will contain `db/` and `tax/` subdirectories required by the pipeline.
