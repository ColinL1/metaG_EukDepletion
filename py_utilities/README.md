# Quick utilities in Python for small tasks

## Folder content

```text
py_utilities
├── README.md        <-- you are here.
├── fetch-tax.py     # Queries NCBI via Biopython to retrieve genus/species names for a list of assembly accessions (GCA/GCF). Outputs a TSV.
└── json_parser.py   # Parses fastp JSON report files from a folder and collects read/base counts per kingdom into a single CSV. Superseded by the seqkit stats approach.
```

## Usage

### fetch-tax.py

Fetches genus and species names from NCBI for a list of genome assembly accessions.

```bash
# From a text file (one accession per line)
python fetch-tax.py accessions.txt > output.tsv

# Or edit the ACCESSIONS list directly in the script and run
python fetch-tax.py > output.tsv
```

> **Note:** Requires `biopython`. Set your email in the script (`Entrez.email`) as required by NCBI.

### json_parser.py

Parses a folder of `fastp` JSON report files and compiles read and base counts into a single CSV, with kingdom assignment based on filename patterns.

```bash
python json_parser.py -j path/to/json/folder/ -o output.csv
```

> **Note:** This script was used in a previous version of the pipeline. It has been superseded by the `seqkit stats` approach used in the R scripts.
