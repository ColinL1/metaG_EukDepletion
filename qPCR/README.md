# Input folder for qPCR data

Contains 16S qPCR Ct values used to assess bacterial abundance across samples and treatments. Each CSV follows the same column structure (`DATE`, `Sample_name`, `Strain`, `Buffer`, `Treatment`, `Gene`, `Ct`, ...) and is read directly by the R scripts for outlier detection and downstream plotting.

Example of folder content (according to what the R code expects):

```text
qPCR
├── README.md                   <-- you are here.
├── MetaG_qPCR_Aip_PBS.csv      # Aiptasia samples preserved in PBS
├── MetaG_qPCR_corals.csv       # Coral samples (Acropora, Porites, Pocillopora)
└── MetaG_qPCR_PBSvsDESS.csv    # Aiptasia samples comparing PBS vs DESS preservation buffers
```
