# Collection of subworkflows used in the pipeline

## Logic

<details>
  <summary> map_conse_pe </summary>

```mermaid
---
config:
  layout: elk
  theme: neutral
  fontFamily: '''Open Sans Variable'', sans-serif'
  themeVariables:
    fontFamily: '''Open Sans Variable'', sans-serif'
---
flowchart LR
 subgraph take["take"]
        v0["reads"]
  end
 subgraph emit["emit"]
        v4["kaiju_report"]
        v10["bowtie2_sym_unmap_fq"]
        v5["bowtie2_host_bam"]
        v8["bowtie2_sym_bam"]
        v6["bowtie2_host_map_fq"]
        v9["bowtie2_sym_map_fq"]
        v7["bowtie2_host_unmap_fq"]
  end
 subgraph MAP_CONSE_PE["MAP_CONSE_PE"]
        take
        v1(["BOWTIE2_MAP_HOST"])
        v2(["BOWTIE2_MAP_SYM"])
        v3(["KAIJU_PE"])
        emit
  end
    v0 --> v1
    v1 --> v2 & v5 & v6 & v7
    v2 --> v3 & v8 & v9 & v10
    v3 --> v4
```
</details>

<details>
  <summary> map_conse_assembly</summary>

```mermaid
---
config:
  layout: elk
  theme: neutral
  fontFamily: '''Open Sans Variable'', sans-serif'
  themeVariables:
    fontFamily: '''Open Sans Variable'', sans-serif'
---
flowchart LR
 subgraph take["take"]
        v0["reads"]
  end
 subgraph emit["emit"]
        v13["bowtie2_sym_unmap_fq"]
        v12["bowtie2_sym_map_fq"]
        v11["bowtie2_sym_bam"]
        v10["bowtie2_host_unmap_fq"]
        v9["bowtie2_host_map_fq"]
        v8["bowtie2_host_bam"]
        v7["kaiju_report"]
  end
 subgraph MAP_CONSE_ASSEMBLY["MAP_CONSE_ASSEMBLY"]
        take
        v1(["MEGAHIT_PE"])
        v2(["MINIMAP2_MAP_HOST"])
        v3(["SPLIT_BAM_HOST"])
        v4(["MINIMAP2_MAP_SYM"])
        v5(["SPLIT_BAM_SYM"])
        v6(["KAIJU_SE"])
        emit
  end
    v0 --> v1
    v1 --> v2
    v2 --> v3 & v8
    v3 --> v4 & v9 & v10
    v4 --> v5 & v11
    v5 --> v6 & v12 & v13
    v6 --> v7

```

</details>

### Folder content

```text
subworkflows
├── README.md               <-- you are here. 
├── CAT.nf
├── exploratory             # A subset of subworkflow or subworkflow ideas. Not actively used in the pipeline. No guarantee any of this code is functional
│   └── archived_old-logic  # As the folder says. Previous snippets of code from the previous version of the pipeline.
├── k2b.nf
├── map_conse_assembly.nf
└── map_conse_pe.nf
```
