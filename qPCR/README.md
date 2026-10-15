# qPCR input

qPCR Ct values of the bacterial 16S rRNA gene and the host β-actin gene (`b-actin` in the `Gene` column), used to compare the relative abundance of bacteria to host across extraction methods (Figure 2). `R-scripts/02_figure_2_qpcr.r` reads them from this folder, runs the outlier detection and computes, per sample, the signed ΔCt:

```math
\Delta Ct = -1 \times \left( Ct_{\text{16S}} - Ct_{\beta\text{-actin}} \right)_{\text{extraction method}}
```

Higher values mean more bacterial 16S copies relative to host β-actin. The CSV files are not in git.

The raw qPCR data is published on Zenodo ([10.5281/zenodo.23238935](https://doi.org/10.5281/zenodo.23238935)) as one unified xlsx file. The R scripts expect it split into the three CSVs below, each with the same columns (`DATE`, `Sample_name`, `Strain`, `Buffer`, `Treatment`, `Gene`, `Ct`, ...).

```text
qPCR
├── README.md                   <-- you are here.
├── MetaG_qPCR_Aip_PBS.csv      # Aiptasia samples preserved in PBS
├── MetaG_qPCR_corals.csv       # Coral samples (Acropora, Porites, Pocillopora)
└── MetaG_qPCR_PBSvsDESS.csv    # Aiptasia samples comparing PBS and DESS preservation buffers
```
