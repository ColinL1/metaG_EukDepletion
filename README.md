# metaG_EukDepletion
Code used in for analysis of kingdom taxonomy to assess the development of a cnidarian metagenomic protocol for efficient eukaryotic DNA removal

# MetaG methods paper

Testing the efficiency of different DNA extractions protocol to enrich the bacterial fraction of Coral samples. 

### DSL2 Nextflow pipeline combining mmseqs2, BlastX, and MEGAN to investigate read-based analysis on ONT long reads sequencing.

---

Workflow includes: (Plotting under constructions :construction_worker: ) 
 - Adapters trimming with porechop
 - Appending sample information (from sample name) to sequence identifier
 - Run against a specified reference database using mmseqs2
    - Possibility of plotting taxonomy information when relevant
 - Run Run against a specified reference database using diamond BlastX
    - Prepare file for MEGAN with "daa-meganizer"
    - Add Taxonomy information to text BlastX output and plot directly


![Alt text](diagram_ID_main.png "Workflow diagram")

---

 Coming soon: 
 - Support for tweaking mmseqs2 parameter directly as nextflow params.  &#10003;
 - Direct plots of taxonomy and GO terms form protein ID. 

---
### Additional resources 
<details><summary> (Old scripts no longer updated or maintained) </summary><p>

Various Nextflow scripts (in [DSL1](DSL1/)) to do:
- BlastX on ONT reads 
- run kaiju on ONT reads
- kraken2 on ONT reads.nf
- Assembly illumina paired-end data with megahit (included optional "name cleanup" to run mmseqs2 more easily)
- A pipeline to run mmseqs2 (various databases) on on ONT reads (including adapter trimming and cleanup)

R script to plot and analyse qPCR ratio ([R](qPCR/))

A bash based pipeline to run BlastX on ONT data get results in by the bp (all the scripts can be found [here](blastX/))

##### Older

<s>
 - Kraken2: script [here](kraken2_ONT.nf)

 - Kaiju: 1) standard script [here](kaiju_ONT.nf) 2) custom under construction
 </s>


<s>
Added nextflow pipeline to get blastX results form ONT reads.

Including:
* Adapter removal (porechop)
* BlastX (diamond blastX)
* Custom parsing scripts (partially based on [neavemj/mappd](https://github.com/neavemj/mappd))
* Nextflow script for [kaiju](kraken2_ONT.nf) and [kraken2](kaiju_ONT.nf) operation
* Nextflow script for fasta file manipulation and more in [extra_scripts](extra_scripts/)</s>

more...
</p></details>