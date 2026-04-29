# R scripts used in downstream processing and visualizations

## Prerequisites

Run Seqkit stats on mapped / unmapped reads to extract the stats used for the plots.

```bash
seqkit stats results/mapping/*/*/*_R1*.gz > stats.txt
# awk 'BEGIN{OFS="\t"} {$1=$1; print}' stats.txt > stats.clean.txt # sanitise to have clear tsv format
```

### Folder content

```text
R-scripts
├── README.md       <-- you are here.
├── All.plot.r      # final version of plots included in the manuscript. 
├── exploratory     # Temporary plots. Kept for future reference
│   └── archived    # archived from previous version. Kept for future reference
├── run_pavian.r
├── table_counts.r  # code for nice table generation
└── table_mags.r    # code for nice table generation
```
