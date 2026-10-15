# References and databases

This folder holds the reference indexes and databases the pipeline needs. None of it is in git. Set each path in `nextflow.config` (copied from [nextflow.config.example](../nextflow.config.example)):

| Param | What | Used by |
|---|---|---|
| `ref_bt2_aip`, `ref_bt2_scl`, `ref_bt2_sym` | Bowtie2 indexes of Aiptasiidae, Scleractinia, Symbiodiniaceae | reads based |
| `ref_mmi_aip`, `ref_mmi_scl`, `ref_mmi_sym` | Minimap2 indexes of the same references | assembly based |
| `nodes`, `names`, `kaiju_db` | Kaiju taxonomy and index | reads based, assembly based |
| `cat_db` | CAT database root (with `db/` and `tax/`) | co-assembly based |

---

## Reference genome indexes (Bowtie2 and Minimap2)

Taxon IDs used in the manuscript: `42823` (Aiptasiidae), `6127` (Acropora), `46730` (Pocillopora), `46719` (Porites), `252141` (Symbiodiniaceae)

### 1. Download reference genome

```bash
# Example for Acropora (taxon 6127). Repeat for each taxon
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

The manuscript run used the NCBI NR index. Pre-built indexes can be downloaded from the Kaiju [downloads page](https://bioinformatics-centre.github.io/kaiju/downloads.html):

```bash
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_nr_2024-08-25.tgz
tar xvf kaiju_db_nr_2024-08-25.tgz
```

Or build it from scratch:

```bash
kaiju-makedb -s nr
```

Set `nodes`, `names` and `kaiju_db` (the `.fmi` file) to the extracted files.

---

## CAT database (NR)

The [CAT pack](https://github.com/MGXlab/CAT_pack) (Contig Annotation Tool) is used for taxonomic classification of assembled contigs. The pipeline uses an NCBI NR-based database, which must expose two subdirectories: `db/` and `tax/`.

### Option 1: download a prebuilt NR database

Prebuilt databases are listed on the [CAT pack download page](https://github.com/MGXlab/CAT_pack#downloading-preconstructed-database-files). The pipeline was run with `20241212_CAT_nr_website`.

### Option 2: build from scratch

```bash
# Download and process NR source files
CAT_pack download --db nr -o path/to/nr_data_dir

# Build the CAT pack database (downloaded files carry a date prefix, e.g. 2024-12-12.nr.gz)
CAT_pack prepare \
  --db_fasta path/to/nr_data_dir/nr.gz \
  --names path/to/nr_data_dir/names.dmp \
  --nodes path/to/nr_data_dir/nodes.dmp \
  --acc2tax path/to/nr_data_dir/prot.accession2taxid.FULL.gz \
  --db_dir CAT_nr
```

The resulting directory will contain `db/` and `tax/` subdirectories required by the pipeline.
